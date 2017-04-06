#include "reconstruction.h"
#include <queue>
#include <vector>
#include <tuple>

void Reconstruction::FindPath(Sequence& read, vector<SequenceNode>& nodes)
{
    int startPos;
    int endPos;

    // [pos in read][node]
    vector<vector<Position>> positions(read.GetData().length(), vector<Position>(nodes.size(), {-1, -1, 0}));

    auto compare = [](const tuple<int, int, int, int>& left, const tuple<int, int, int, int>& right)
    {
        return get<0>(left) > get<0>(right) || (!(get<0>(right) > get<0>(left)) && get<1>(left) < get<1>(right));
    };
    // {edit distance, pos in read, node, overlap}
    priority_queue<tuple<int, int, int, int>, vector<tuple<int, int, int, int>>, decltype(compare)> connections(compare);

    for (int i = 0; i < nodes.size(); i++)
    {
        connections.push(make_tuple(0, startPos, 0, i));
    }

    string result;

    while (!connections.empty())
    {
        auto connection = connections.top();
        connections.pop();
        int &distance(get<0>(connection));
        int &posInRead(get<1>(connection));
        int &node(get<2>(connection));
        int &overlap(get<3>(connection));

        if (posInRead >= endPos)
        {
            result = RecreatePath(positions, nodes, posInRead, overlap);
        }
        
        if (positions[posInRead][node].distance != -1)
        {
            continue;
        }

        positions[posInRead][node] = {distance, node, overlap};

        for (auto& nextNode : nodes[node].overlaps)
        {
            // {alignment end pos in read, score}
            vector<pair<int, int>> alignmentScores;
            GetAlignmentScores(read, posInRead, nodes[nextNode.first], nextNode.second, alignmentScores);

            int stepLength = GetStepLength(nodes[nextNode.first], nextNode.second, posInRead);
            for (auto& alignment : alignmentScores)
            {
                int newDistance = distance + stepLength + alignment.second;
                connections.push(make_tuple(newDistance, alignment.first, nextNode.first, nextNode.second));
            }
        }
    }
}

void Reconstruction::GetAlignmentScores(Sequence& read, int readPos, SequenceNode& node, int nodePos, vector<pair<int, int>>& scores)
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

string Reconstruction::RecreatePath(vector<vector<Position>>& positions, vector<SequenceNode>& nodes, int posInRead, int node)
{
    string result;
    while (positions[posInRead][node].distance > 0)
    {
        Position &pos(positions[posInRead][node]);
        string nextPiece = nodes[node].sequence.substr(pos.overlapWithPrevious, int(nodes[node].sequence.length()) - pos.overlapWithPrevious);
        reverse(nextPiece.begin(), nextPiece.end());
        result.append(nextPiece);
    }
    reverse(result.begin(), result.end());
    return result;
}

int Reconstruction::GetStepLength(SequenceNode& node, int overlap, int posInRead)
{
    return node.isFromReference * scoring.notFromReferencePenalty 
        + scoring.overlapPenalty(overlap, node.sequence.length()) 
        + scoring.misplacementPenalty(posInRead - node.expectedPos);
}