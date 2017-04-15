#include "reconstruction.h"
#include <queue>
#include <vector>
#include <tuple>
#include <iostream>

string Reconstruction::FindPath(const Sequence& read, const vector<SequenceNode>& nodes)
{
    int endPos = read.GetData().length();

    int startNode = -1;
    for (int i = 0; i < nodes.size(); i++)
    {
        if (nodes[i].isFromReference)
        {
            if (startNode == -1 || nodes[i].expectedPos < nodes[startNode].expectedPos)
            {
                startNode = i;
            }
        }
    }
    if (startNode == -1)
    {
        return "";
    }
    

    // [pos in read][node]
    vector<vector<Position>> positions(read.GetData().length() + 1, vector<Position>(nodes.size(), {-1, -1, 0}));

    struct ReconstructionStep
    {
        int distance;
        int posInRead;
        int node;
        int overlap;
        int previousNode;
    };

    auto compare = [](const ReconstructionStep& left, const ReconstructionStep& right)
    {
        return left.distance > right.distance || (!(right.distance > left.distance) && left.posInRead < right.posInRead);
    };

    // {edit distance, pos in read, node, overlap}
    priority_queue<ReconstructionStep, vector<ReconstructionStep>, decltype(compare)> steps(compare);

    steps.push({0, nodes[startNode].expectedPos + (int)nodes[startNode].sequence.length(), startNode, 0, -1});

    string result;

    long long iterationCounter = 0;

    while (!steps.empty())
    {
        cout << iterationCounter++ << " " << steps.size() << endl;
        auto step = steps.top();
        steps.pop();
        int &distance(step.distance);
        int &posInRead(step.posInRead);
        int &node(step.node);
        int &overlap(step.overlap);
        int &previousNode(step.previousNode);
        
        if (positions[posInRead][node].distance != -1)
        {
            continue;
        }

        positions[posInRead][node] = {distance, previousNode, overlap};

        if (posInRead >= endPos)
        {
            result = RecreatePath(positions, nodes, posInRead, node);
            return result;
        }

        for (auto& nextNode : nodes[node].overlaps)
        {
            // {alignment end pos in read, score}
            vector<pair<int, int>> alignmentScores;
            GetAlignmentScores(read, posInRead, nodes[nextNode.first], nextNode.second, alignmentScores);

            int stepLength = GetStepLength(nodes[nextNode.first], nextNode.second, posInRead);
            for (auto& alignment : alignmentScores)
            {
                int newDistance = distance + stepLength + alignment.second;
                steps.push({newDistance, alignment.first, nextNode.first, nextNode.second, node});
            }
        }
    }
    return "";
}

void Reconstruction::GetAlignmentScores(const Sequence& read, int readPos, const SequenceNode& node, int nodePos, vector<pair<int, int>>& scores)
{
    string readSeq = read.GetData().substr(readPos, min(read.GetData().length() - readPos, node.sequence.length() * 2));
    //string nodeSeq = node.sequence.substr(nodePos, node.sequence.length() - nodePos);
    string nodeSeq(node.sequence.begin() + nodePos, node.sequence.end());

    if (!alignmentMatrix.empty() && alignmentMatrix[0].size() < readSeq.length() + 1)
    {
        for (int i = 0; i < alignmentMatrix.size(); i++)
        {
            alignmentMatrix[i].resize(readSeq.length() + 1);
        }
    }

    if (alignmentMatrix.size() < nodeSeq.length() + 1)
    {
        alignmentMatrix.resize(nodeSeq.length() + 1, vector<int>(readSeq.length() + 1));
    }

    // fill trivial solutions (all inserts and all deletes)
    for (int i = 0; i <= nodeSeq.length(); i++)
    {
        alignmentMatrix[i][0] = i;
    }
    for (int i = 1; i <= readSeq.length(); i++)
    {
        alignmentMatrix[0][i] = i;
    }

    // fill alignment matrix
    for (int i = 1; i <= nodeSeq.length(); i++)
    {
        for (int j = 1; j <= readSeq.length(); j++)
        {
            int insertion = alignmentMatrix[i - 1][j] + scoring.insertion;
            int deletion = alignmentMatrix[i][j - 1] + scoring.deletion;
            int substitution = alignmentMatrix[i - 1][j - 1] + scoring.substitution * (nodeSeq[i - 1] != readSeq[j - 1]);
            alignmentMatrix[i][j] = min({insertion, deletion, substitution});
        }
    }

    scores.clear();
    for (int i = 0; i <= readSeq.length(); i++)
    {
        scores.push_back(make_pair(readPos + i, alignmentMatrix[nodeSeq.length()][i]));
    }
}

string Reconstruction::RecreatePath(const vector<vector<Position>>& positions, const vector<SequenceNode>& nodes, int posInRead, int node)
{
    string result;
    while (node >= 0)
    {
        const Position &pos(positions[posInRead][node]);
        string nextPiece = nodes[node].sequence.substr(pos.overlapWithPrevious, int(nodes[node].sequence.length()) - pos.overlapWithPrevious);
        reverse(nextPiece.begin(), nextPiece.end());
        result.append(nextPiece);
        posInRead = posInRead - nodes[node].sequence.length() + pos.overlapWithPrevious;
        node = pos.previousNode;
    }
    reverse(result.begin(), result.end());
    return result;
}

int Reconstruction::GetStepLength(const SequenceNode& node, int overlap, int posInRead)
{
    return (!node.isFromReference) * scoring.notFromReferencePenalty 
        + scoring.overlapPenalty(overlap, node.sequence.length()) 
        + scoring.misplacementPenalty(posInRead - node.expectedPos);
}