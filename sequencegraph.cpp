#include "sequencegraph.h"
#include <iostream>
#include <tuple>
#include <set>
#include "utils.h"
using namespace std;

const vector<pair<int, int>>& SequenceNode::GetOverlaps(char nextBase) const
{
    const vector<pair<int, int>>* result;
    switch (nextBase)
    {
        case 'A': result = &overlapsA; break;
        case 'C': result = &overlapsC; break;
        case 'G': result = &overlapsG; break;
        default: result = &overlapsT; break;
    }
    if (result->size() > 0)
    {
        return *result;
    }
    else
    {
        return overlaps;
    }
}

void SequenceGraph::LoadReads(const string &filename)
{
    FASTA f(filename);
    Sequence s;

    while (f >> s)
    {
        forward.emplace(forward.end(), s);
        reverse.emplace(reverse.end(), s, true);
    }
}

void SequenceGraph::GetNodes(int idRead, vector<SequenceNode> &nodes)
{
    nodes.clear();
    vector<SequenceNode> referenceNodes;
    vector<long long> sequences;
    map<int, int> offsetRelations;
    for (auto &neighbor : overlapGraph.adjacency[idRead])
    {
        int idNeighbor = neighbor.first;

        bool reversed = neighbor.second[0].reversed;
        int meanOffset = 0;
        offsetRelations.clear();

        for (Match &m : neighbor.second)
        {
            // check match
            // Match swappedPositions = m;
            // swap(swappedPositions.pos1, swappedPositions.pos2);
            // string seq = (m.reversed ? reverse : forward)[idNeighbor].GetData().substr(m.pos2, m.length);
            // string seq2 = GetSequence(idRead, m, false).first;
            // if (GetSequence(idRead, m, false).first != seq)
            // {
            //     cout << "megapipkos jemine, nezhoda v zhode" << endl;
            // }
            // if (m.reversed != reversed)
            // {
            //     cout << "megapipkos jemine, pomiesane forwardy a reversy";
            // }
            // sequences from the reference
            GetSequences(idRead, m, false, sequences);
            for (unsigned int i = 0; i < sequences.size(); i++)
            {
                nodes.push_back({sequences[i], m.pos1 + (int)i, true, {}, false, params.nodeLength, 4});
            }
            // referenceNodes.push_back({GetSequence(idRead, m, false).first, m.pos1, true, {}, false});
            //referenceNodes.insert(m.pos1);

            offsetRelations[m.pos2] = m.pos1;
            // meanOffset += m.pos1 - m.pos2;
        }
        // meanOffset /= (int)neighbor.second.size();

        for (auto &neighbor2 : overlapGraph.adjacency[idNeighbor])
        {
            if (neighbor2.first == idRead)
            {
                continue;
            }
            for (Match &m : neighbor2.second)
            {
                // auto seqPos = GetSequence(idNeighbor, m, reversed);
                // nodes.push_back({seqPos.first, seqPos.second + meanOffset, false, {}, reversed});

                int pos = GetSequences(idNeighbor, m, reversed, sequences);
                int refPos = pos;
                auto closestRight = offsetRelations.upper_bound(pos);
                if (closestRight == offsetRelations.end())
                {
                    refPos += offsetRelations.rbegin()->second - offsetRelations.rbegin()->first; 
                }
                else if (closestRight == offsetRelations.begin())
                {
                    refPos += offsetRelations.begin()->second - offsetRelations.begin()->first; 
                }
                else
                {
                    auto closestLeft = closestRight;
                    closestLeft--;
                    refPos += (closestLeft->second + closestRight->second - closestLeft->first - closestRight->first) / 2;
                }

                for (unsigned int i = 0; i < sequences.size(); i++)
                {
                    int nodePos = refPos + (int)i;
                    if (nodePos > -100 && nodePos < forward[idRead].GetData().length() + 100)
                    {
                        nodes.push_back({sequences[i], refPos + (int)i, false, {}, reversed, params.nodeLength});
                    }
                }
            }
        }
    }

    // sort(referenceNodes.begin(), referenceNodes.end(), [](const SequenceNode &left, const SequenceNode &right) { 
    //     return left.expectedPos < right.expectedPos || (!(right.expectedPos < left.expectedPos) && (left.sequence.length() > right.sequence.length())); 
    // });
    // int prevBegin = -1;
    // int prevEnd = -1;
    // for (SequenceNode& n : referenceNodes)
    // {
    //     if ((n.expectedPos != prevBegin) && (n.expectedPos + (int)n.sequence.length() > prevEnd))
    //     {
    //         nodes.push_back(n);
    //         prevBegin = n.expectedPos;
    //         prevEnd = n.expectedPos + n.sequence.length();
    //     }
    // }
    Utils::StartTiming();
    RemoveDuplicates(nodes);
    Utils::VerbalResult("remove duplicates");

    Utils::StartTiming();
    FindOverlaps(nodes);
    Utils::VerbalResult("find overlaps");
}

// pair<string, int> SequenceGraph::GetSequence(int idRead, Match &m, bool reversed)
// {
//     Sequence &s((reversed ? reverse : forward)[idRead]);
//     int length = s.GetData().length();
//     int pos = reversed ? (length - m.pos1 - m.length) : m.pos1;
//     return make_pair(s.GetData().substr(pos, m.length), pos);
// }

int SequenceGraph::GetSequences(int idRead, Match& m, bool reversed, vector<long long>& sequences)
{
    Sequence &s((reversed ? reverse : forward)[idRead]);
    int length = s.GetData().length();
    int pos = reversed ? (length - m.pos1 - m.length) : m.pos1;

    sequences.clear();
    long long sequence = 0;
    unsigned long long mask = ((long long)1 << (2 * params.nodeLength)) - 1;
    for (int i = 0; i < m.length; i++)
    {
        long long base = 0;
        switch (s.GetData()[pos + i])
        {
            case 'C': base = 1; break;
            case 'G': base = 2; break;
            case 'T': base = 3; break;
        }
        sequence = ((sequence << 2) + base) & mask;
        if (i >= params.nodeLength - 1)
        {
            sequences.push_back(sequence);
        }
    }
    return pos;
}

void SequenceGraph::RemoveDuplicates(vector<SequenceNode>& nodes)
{
    sort(nodes.begin(), nodes.end(), [](const SequenceNode &left, const SequenceNode &right) { 
        return left.sequence < right.sequence || (!(right.sequence < left.sequence) && left.expectedPos < right.expectedPos); 
    });

    const int maxExpectedPosDistanceForTheSameNode = 500;
    
    vector<SequenceNode> nodesND; // no duplicates
    for (int i = 0; i < (int)nodes.size(); )
    {
        int j = i;
        bool fromReference = false;
        int posFromReference;
        int sumExpectedPos = 0;
        while (j < (int)nodes.size() && nodes[i].sequence == nodes[j].sequence && (!fromReference || !nodes[j].isFromReference || posFromReference == nodes[j].expectedPos)
            && nodes[i].expectedPos + maxExpectedPosDistanceForTheSameNode >= nodes[j].expectedPos)
        {
            if (nodes[j].isFromReference)
            {
                fromReference = true;
                posFromReference = nodes[j].expectedPos;
            }
            sumExpectedPos += nodes[j].expectedPos;
            j++;
        }
        nodesND.push_back({nodes[i].sequence, fromReference ? posFromReference : (sumExpectedPos / (j - i)), fromReference, {}, false, params.nodeLength});
        i = j;
    }
    nodes.swap(nodesND);
}

void SequenceGraph::FindOverlaps(vector<SequenceNode> &nodes)
{
    vector<vector<tuple<long long, bool, int>>> prefixesSufixes(params.nodeLength);
    for (int nodeId = 0; nodeId < nodes.size(); nodeId++)
    {
        SequenceNode &node(nodes[nodeId]);
        for (int i = params.minKmerOverlap; i < params.nodeLength; i++)
        {
            prefixesSufixes[i].push_back(make_tuple(node.sequence >> ((params.nodeLength - i) * 2), true, nodeId));
            unsigned long long mask = ((long long)1 << (2 * i)) - 1;
            prefixesSufixes[i].push_back(make_tuple(node.sequence & mask, false, nodeId));
        }
    }

    for (int i = params.nodeLength - 1; i >= params.minKmerOverlap; i--)
    {
        sort(prefixesSufixes[i].begin(), prefixesSufixes[i].end());

        auto sufBegin = prefixesSufixes[i].begin();
        while (sufBegin != prefixesSufixes[i].end())
        {
            auto prefBegin = sufBegin;
            while (prefBegin != prefixesSufixes[i].end() && get<0>(*prefBegin) == get<0>(*sufBegin) && get<1>(*prefBegin) == false)
            {
                ++prefBegin;
            }
            auto nextSeq = prefBegin;
            while (nextSeq != prefixesSufixes[i].end() && get<0>(*nextSeq) == get<0>(*sufBegin))
            {
                ++nextSeq;
            }
            int numSuffix = prefBegin - sufBegin;
            int numPrefix = nextSeq - prefBegin;
            for (auto it1 = sufBegin; it1 != prefBegin; ++it1)
            {
                for (auto it2 = prefBegin; it2 != nextSeq; ++it2)
                {
                    int idRead1 = get<2>(*it1);
                    int idRead2 = get<2>(*it2);
                    if (abs(nodes[idRead1].expectedPos - nodes[idRead2].expectedPos) < params.overlappingKmersMaxExpectedPosDistance/*150*/)
                    {
                        nodes[idRead1].overlaps.push_back(make_pair(idRead2, i));
                    }
                }
            }
            sufBegin = nextSeq;
        }
    }

    for (auto& n : nodes)
    {
        for (int i = 0; i < n.overlaps.size() && n.overlaps[i].second == params.nodeLength - 1; i++)
        {
            switch (nodes[n.overlaps[i].first][params.nodeLength - 1])
            {
                case 'A': n.overlapsA.push_back(n.overlaps[i]);
                case 'C': n.overlapsC.push_back(n.overlaps[i]);
                case 'G': n.overlapsG.push_back(n.overlaps[i]);
                default: n.overlapsT.push_back(n.overlaps[i]);
            }
        }
        if (params.maxEdgesFromNode >= -1 && n.overlaps.size() > params.maxEdgesFromNode)
        {
            n.overlaps.resize(params.maxEdgesFromNode);
        }
    }
}