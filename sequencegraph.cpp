#include "sequencegraph.h"
#include <iostream>
#include <tuple>
using namespace std;


void SequenceGraph::LoadReads(const string& filename)
{
    FASTQ f(filename);
    Sequence s;

    while (f >> s)
    {
        forward.emplace(forward.end(), s);
        reverse.emplace(reverse.end(), s, true);
    }
}

void SequenceGraph::GetNodes(int idRead, vector<SequenceNode>& nodes)
{
    nodes.clear();
    for (auto& neighbor : overlapGraph.adjacency[idRead])
    {
        int idNeighbor = neighbor.first;
        if (neighbor.second.size() == 0)
        {
            cout << "sused bez zhody" << endl;
        }
        bool reversed = neighbor.second[0].reversed;
        int meanOffset = 0;
        for (Match& m : neighbor.second)
        {
            // check match
            Match swappedPositions = m;
            swap(swappedPositions.pos1, swappedPositions.pos2);
            string seq = (m.reversed ? reverse : forward)[idNeighbor].GetData().substr(m.pos2, m.length);
            string seq2 = GetSequence(idRead, m, false).first;
            if (GetSequence(idRead, m, false).first != seq)
            {
                cout << "megapipkos jemine, nezhoda v zhode" << endl;
            }
            if (m.reversed != reversed)
            {
                cout << "megapipkos jemine, pomiesane forwardy a reversy";
            }
            // sequences from the reference
            nodes.push_back({GetSequence(idRead, m, false).first, m.pos1, true, {}});
            meanOffset += m.pos1 - m.pos2;
        }
        meanOffset /= neighbor.second.size();

        for (auto& neighbor2 : overlapGraph.adjacency[idNeighbor])
        {
            if (neighbor2.first == idRead)
            {
                continue;
            }
            for (Match& m : neighbor2.second)
            {
                auto seqPos = GetSequence(idNeighbor, m, reversed);
                nodes.push_back({seqPos.first, seqPos.second + meanOffset, true, {}});
            }
        }
    }
    FindOverlaps(nodes);
}

pair<string, int> SequenceGraph::GetSequence(int idRead, Match& m, bool reversed)
{
    Sequence& s((reversed ? reverse : forward)[idRead]);
    int length = s.GetData().length();
    int pos = reversed ? (length - m.pos1 - m.length) : m.pos1;
    return make_pair(s.GetData().substr(pos, m.length), pos);
}

void SequenceGraph::FindOverlaps(vector<SequenceNode>& nodes)
{
    vector<tuple<string, bool, int>> prefixesSufixes;
    for (int nodeId = 0; nodeId < nodes.size(); nodeId++)
    {
        SequenceNode& node(nodes[nodeId]);
        for (int i = 1; i < node.sequence.length() - 1; i++)
        {
            prefixesSufixes.push_back(make_tuple(node.sequence.substr(0, i), true, nodeId));
            prefixesSufixes.push_back(make_tuple(node.sequence.substr(node.sequence.length() - i, i), false, nodeId));
        }
    }
    sort(prefixesSufixes.begin(), prefixesSufixes.end());

    auto sufBegin = prefixesSufixes.begin();
    while (sufBegin != prefixesSufixes.end())
    {
        auto prefBegin = sufBegin;
        while (prefBegin != prefixesSufixes.end() && get<0>(*prefBegin) == get<0>(*sufBegin) && get<1>(*prefBegin) == false)
        {
            ++prefBegin;
        }
        auto nextSeq = prefBegin;
        while (nextSeq != prefixesSufixes.end() && get<0>(*nextSeq) == get<0>(*sufBegin))
        {
            ++nextSeq;
        }
        for (auto it1 = sufBegin; it1 != prefBegin; ++it1)
        {
            for (auto it2 = prefBegin; it2 != nextSeq; ++it2)
            {
                int idRead1 = get<2>(*it1);
                int idRead2 = get<2>(*it2);
                if (idRead1 != idRead2)
                {
                    nodes[idRead1].overlaps.push_back(make_pair(idRead2, get<0>(*it1).length()));
                }
            }
        }
        sufBegin = nextSeq;
    }
}