#include "reconstruction.h"
#include <queue>
#include <vector>
#include <tuple>
#include <iostream>
#include <climits>
#include "utils.h"
#include <set>

pair<string, int> Reconstruction::FindPath(const Sequence& read, const vector<SequenceNode>& nodes, int begin, int end)
{
    cout << "node count: " << nodes.size() << endl;

    auto compare = [end](const Position& left, const Position& right)
    {
        return ((left.node == end) < (right.node == end)) || 
        (!((left.node == end) > (right.node == end)) && (left.distance > right.distance || (!(right.distance > left.distance) && left.posInRead < right.posInRead)));
    };

    set<tuple<int, int, int>> visited;
    vector<Position> positions;

    priority_queue<Position, vector<Position>, decltype(compare)> steps(compare);
    SequenceNode beginNode = nodes[begin];
    SequenceNode endNode = nodes[end];

    steps.push({0, beginNode.expectedPos + (int)beginNode.sequence.length(), (int)beginNode.sequence.length(), begin, -1});


    long long iterationCounter = 0;
    int furthestInRead = 0;
    while (!steps.empty())
    {
        iterationCounter++;;
        //cout << iterationCounter++ << " " << steps.size() << endl;
        auto step = steps.top();
        steps.pop();

        if (step.posInRead >= endNode.expectedPos + endNode.sequence.length())
        {
            continue;
        }

        if (visited.find(make_tuple(step.posInRead, step.node, step.posInNode)) != visited.end())
        {
            continue;
        }
        visited.insert(make_tuple(step.posInRead, step.node, step.posInNode));
        int currentPosition = positions.size();
        positions.push_back({step.distance, step.posInRead, step.posInNode, step.node, step.previousPos});

        if (step.node == end && step.posInRead == endNode.expectedPos + step.posInNode)
        {
            //end
            cout << "iterations: " << iterationCounter << endl;
            return make_pair(RecreatePath(positions, nodes, currentPosition), step.distance);
        }

        furthestInRead = max(furthestInRead, step.posInRead);
        if (furthestInRead - 10 > step.posInRead)
        {
            continue;
        }

        if (step.posInNode == nodes[step.node].sequence.length())
        {
            for (auto& overlapNextNode : nodes[step.node].overlaps)
            {
                const SequenceNode& nextNode = nodes[overlapNextNode.first];
                if (!nextNode.isFromReference
                 || (overlapNextNode.first == end && step.posInRead == endNode.expectedPos + overlapNextNode.second))
                {
                    int newDistance = step.distance
                        + scoring.overlapPenalty(overlapNextNode.second, nextNode.sequence.length())
                        + scoring.misplacementPenalty(step.posInRead - nextNode.expectedPos - overlapNextNode.second);
                    steps.push({
                        newDistance, 
                        step.posInRead, 
                        overlapNextNode.second, 
                        overlapNextNode.first, 
                        currentPosition
                    });
                }
            }
        }
        else
        {
            bool match = read.GetData()[step.posInRead] == nodes[step.node].sequence[step.posInNode];
            //match/mismatch
            steps.push({
                step.distance + (!match) * scoring.substitution,
                step.posInRead + 1,
                step.posInNode + 1,
                step.node,
                currentPosition
            });
            if (!match)
            {
                //insertion
                steps.push({
                    step.distance + scoring.insertion,
                    step.posInRead + 1,
                    step.posInNode,
                    step.node,
                    currentPosition
                });
                //deletion
                steps.push({
                    step.distance + scoring.deletion,
                    step.posInRead,
                    step.posInNode + 1,
                    step.node,
                    currentPosition
                });
            }
        }
    }
    cout << "iterations: " << iterationCounter << endl;
    return make_pair("", -1);
}

string Reconstruction::RecreatePath(const vector<Position>& positions, const vector<SequenceNode>& nodes, int currPos)
{
    string result;
    int endNodeBases = positions[currPos].posInNode + 1;
    int prevPosInNode = -1;
    int prevNode = -1;
    while (currPos >= 0)
    {
        const Position& position = positions[currPos];
        if ((position.node != prevNode || position.posInNode != prevPosInNode) 
            && position.posInNode < nodes[position.node].sequence.length())
        {
            result += nodes[position.node].sequence[position.posInNode];
        }
        prevPosInNode = position.posInNode;
        prevNode = position.node;

        currPos = position.previousPos;
    }
    reverse(result.begin(), result.end());
    return result.substr(0, result.length() - endNodeBases);
/*    string result;
    const SequenceNode& endNode = nodes[node];
    int posInEndNode = posInNode;
    while (posInRead != -1)
    {
        const Position& position = positions[posInRead][node][posInNode];
        if ((position.previousNode != node || position.previousPosInNode != posInNode) 
            && position.previousPosInNode < nodes[position.previousNode].sequence.length())
        {
            result += nodes[position.previousNode].sequence[position.previousPosInNode];
        }
        posInRead = position.previousPosInRead;
        posInNode = position.previousPosInNode;
        node = position.previousNode;
    }

    reverse(result.begin(), result.end());
    return result.substr(0, result.length() - posInEndNode);*/
}

string Reconstruction::Reconstruct(const Sequence& read, const vector<SequenceNode>& nodes)
{
    cout << "reconstruct" << endl;
    const int skipsAllowed = 2;

    vector<int> nodesFromReference;
    for (size_t i = 0; i < nodes.size(); i++)
    {
        if (nodes[i].isFromReference)
        {
            nodesFromReference.push_back(i);
        }
    }
    cout << "ndoes from reference: " << nodesFromReference.size() << endl;
    sort(nodesFromReference.begin(), nodesFromReference.end(), [&nodes](int left, int right) {
        return nodes[left].expectedPos < nodes[right].expectedPos;
    });
    cout << "sorted" << endl;

    vector<string> paths(nodesFromReference.size());

    // distance, previous node
    vector<pair<int, int>> distanceToMatch(nodesFromReference.size(), make_pair(INT_MAX, -1));
    distanceToMatch[0].first = 0;

    for (int i = 1; i < nodesFromReference.size(); i++)
    {
        for (int jump = 1; jump <= skipsAllowed && i - jump >= 0; jump++)
        {
            cout << nodesFromReference[i - jump] << " " << nodesFromReference[i] << endl;
            cout << " expected pos: " << nodes[nodesFromReference[i - jump]].expectedPos << " " <<
            nodes[nodesFromReference[i]].expectedPos << " " << nodes[nodesFromReference[i - jump]].sequence << " " << nodes[nodesFromReference[i]].sequence << endl;
            Utils::StartTiming();
            auto path = FindPath(read, nodes, nodesFromReference[i - jump], nodesFromReference[i]);
            Utils::VerbalResult("Find path");
            cout << "from " << i - jump << " to " << i << " distance " << path.second << " " << path.first << endl;
            cout << endl << endl;
            if (path.second != -1 && distanceToMatch[i].first > path.second + distanceToMatch[i - jump].first)
            {
                distanceToMatch[i].first = path.second + distanceToMatch[i - jump].first;
                distanceToMatch[i].second = i - jump;
                path.first.swap(paths[i]);
            }
            
        }
    }
    vector<int> bestPathNodes;
    for (int pos = (int)distanceToMatch.size() - 1; pos != -1; pos = distanceToMatch[pos].second)
    {
        bestPathNodes.push_back(pos);
    }
    reverse(bestPathNodes.begin(), bestPathNodes.end());

    string finalSequence = nodes[nodesFromReference[bestPathNodes[0]]].sequence;
    for (int i = 1; i < bestPathNodes.size(); i++)
    {
        finalSequence += paths[bestPathNodes[i]];
        finalSequence += nodes[nodesFromReference[bestPathNodes[i]]].sequence;
    }
    return finalSequence;
}