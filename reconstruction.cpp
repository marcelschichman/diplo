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
    //cout << "node count: " << nodes.size() << endl;

    auto compare = [end](const Position& left, const Position& right)
    {
        // return ((left.node == end) < (right.node == end)) || 
        // (!((left.node == end) > (right.node == end)) && (left.distance > right.distance || (!(right.distance > left.distance) && left.posInRead < right.posInRead)));
        return left.distance > right.distance || (!(right.distance > left.distance) && left.posInRead < right.posInRead);
    };

    set<tuple<int, int, int>> visited;
    vector<Position> positions;

    priority_queue<Position, vector<Position>, decltype(compare)> steps(compare);
    SequenceNode beginNode = nodes[begin];
    SequenceNode endNode = nodes[end];

    steps.push({0, beginNode.expectedPos + (int)beginNode.length, (int)beginNode.length, begin, -1, pruning.minMatchesWindowSize, (1 << pruning.minMatchesWindowSize) - 1});


    long long iterationCounter = 0;
    int furthestInRead = 0;
    long long pruned = 0;
    while (!steps.empty())
    {
        iterationCounter++;;
        //cout << iterationCounter++ << " " << steps.size() << endl;
        auto step = steps.top();
        steps.pop();
        // if (begin == 85 && end == 5898)
        // {
        //     cout << step.distance << " " << step.posInRead << " " << step.posInNode << " " << step.node << " " << nodes[step.node].GetSequence() << endl;
        // }

        if (step.posInRead >= endNode.expectedPos + endNode.length)
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
            // cout << "pruned: " << pruned << endl;
            // cout << "iterations: " << iterationCounter << endl;
            return make_pair(RecreatePath(positions, nodes, currentPosition), step.distance);
        }

        furthestInRead = max(furthestInRead, step.posInRead);
        if (furthestInRead - pruning.maxDistanceFromFurthest > step.posInRead)
        {
            continue;
        }
        if (step.numMatches < pruning.minMatches)
        {
            pruned++;
            continue;
        }

        if (step.posInNode == nodes[step.node].length)
        {
            for (auto& overlapNextNode : nodes[step.node].GetOverlaps(read.GetData()[step.posInRead]))
            {
                const SequenceNode& nextNode = nodes[overlapNextNode.first];
                // if (!nextNode.isFromReference
                //  || (overlapNextNode.first == end && step.posInRead == endNode.expectedPos + overlapNextNode.second))
                // {
                    int jumpSize = nextNode.length - overlapNextNode.second - 1;
                    int newDistance = step.distance
                        + jumpSize * scoring.jumpLinearPenalty + jumpSize * jumpSize * scoring.jumpQuadraticPenalty
                        + (step.posInRead - nextNode.expectedPos - overlapNextNode.second > scoring.maxKmerMisplacement) * 50;
                    steps.push({
                        newDistance, 
                        step.posInRead, 
                        overlapNextNode.second, 
                        overlapNextNode.first, 
                        currentPosition,
                        step.numMatches,
                        step.matchesBits
                    });
                // }
            }
        }
        else
        {
            bool match = read.GetData()[step.posInRead] == nodes[step.node][step.posInNode];
            int numMatches = step.numMatches - ((step.matchesBits & (1 << (pruning.minMatchesWindowSize - 1))) > 0);
            int matchesBits = (step.matchesBits << 1) & ((1 << pruning.minMatchesWindowSize) - 1);
            //match/mismatch
            steps.push({
                step.distance + (!match) * scoring.substitution,
                step.posInRead + 1,
                step.posInNode + 1,
                step.node,
                currentPosition,
                numMatches + match,
                matchesBits + match
            });
            if (!match)
            {
                //insertion
                steps.push({
                    step.distance + scoring.insertion,
                    step.posInRead + 1,
                    step.posInNode,
                    step.node,
                    currentPosition,
                    numMatches,
                    matchesBits
                });
                //deletion
                steps.push({
                    step.distance + scoring.deletion,
                    step.posInRead,
                    step.posInNode + 1,
                    step.node,
                    currentPosition,
                    numMatches,
                    matchesBits
                });
            }
        }
    }
    // cout << "iterations: " << iterationCounter << endl;
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
            && position.posInNode < nodes[position.node].length)
        {
            result += nodes[position.node][position.posInNode];
        }
        prevPosInNode = position.posInNode;
        prevNode = position.node;

        currPos = position.previousPos;
    }
    reverse(result.begin(), result.end());
    return result.substr(0, max(0, (int)result.length() - endNodeBases));
}

void Reconstruction::Reconstruct(const Sequence& read, const vector<SequenceNode>& nodes, vector<string>& result)
{
    result.clear();
    // cout << "reconstruct" << endl;
    // const int skipsAllowed = 2;
    // const int minSequenceLength = 40;

    vector<int> nodesFromReference;
    for (size_t i = 0; i < nodes.size(); i++)
    {
        if (nodes[i].isFromReference)
        {
            nodesFromReference.push_back(i);
        }
    }
    cerr << "ndoes from reference: " << nodesFromReference.size() << endl;
    sort(nodesFromReference.begin(), nodesFromReference.end(), [&nodes](int left, int right) {
        return nodes[left].expectedPos < nodes[right].expectedPos;
    });
    // cout << "sorted" << endl;

    vector<string> paths(nodesFromReference.size());

    // distance, previous node
    struct PathToMatch
    {
        int distance;
        int previousNode;
        int sequenceLength;
        string lastPathSequence;
        int overlap;
    };
    vector<PathToMatch> pathToMatches(nodesFromReference.size(), {INT_MAX, -1, 0, "", 0});

    for (int i = 0; i < nodesFromReference.size(); i++)
    {
        bool isNewStart = true;
        for (int jump = 1; jump <= pruning.skipsAllowed && i - jump >= 0; jump++)
        {
            if (pathToMatches[i - jump].distance != INT_MAX)
            {
                isNewStart = false;
                break;
            }
        }
        if (isNewStart)
        {
            pathToMatches[i].distance = 0;
            continue;
        }

        for (int jump = 1; jump <= pruning.skipsAllowed && i - jump >= 0; jump++)
        {
            if (pathToMatches[i - jump].distance == INT_MAX)
            {
                continue;
            }
            // /**/cout << nodesFromReference[i - jump] << " " << nodesFromReference[i] << endl;
            // /**/cout << " expected pos: " << nodes[nodesFromReference[i - jump]].expectedPos << " " <<
            // /**/nodes[nodesFromReference[i]].expectedPos << " " << nodes[nodesFromReference[i - jump]].GetSequence() << " " << nodes[nodesFromReference[i]].GetSequence() << endl;


            const SequenceNode& leftNode = nodes[nodesFromReference[i - jump]];
            const SequenceNode& rightNode = nodes[nodesFromReference[i]];
            int leftNodeEnd = leftNode.expectedPos + leftNode.length;
            int rightNodeBegin = rightNode.expectedPos;
            int overlap = 0;
            pair<string, int> path;
            if (leftNodeEnd > rightNodeBegin)
            {
                path = make_pair("", 0);
                overlap = leftNodeEnd - rightNodeBegin;
            }
            else
            {
                // /**/Utils::StartTiming();
                path = FindPath(read, nodes, nodesFromReference[i - jump], nodesFromReference[i]);
                // /**/Utils::VerbalResult("Find path");
            } 
            // /**/cout << "from " << i - jump << " to " << i << " distance " << path.second << " " << path.first << endl;
            // /**/cout << endl << endl;
            if (path.second != -1 && pathToMatches[i].distance > path.second + pathToMatches[i - jump].distance)
            {
                pathToMatches[i].distance = path.second + pathToMatches[i - jump].distance;
                pathToMatches[i].previousNode = i - jump;
                pathToMatches[i].sequenceLength = pathToMatches[i - jump].sequenceLength + path.first.length() + rightNode.length - overlap;
                path.first.swap(pathToMatches[i].lastPathSequence);
                pathToMatches[i].overlap = overlap;
            }
        }
    }
    for (int pos = (int)pathToMatches.size() - 1; pos >= 0; pos--)
    {
        if (pathToMatches[pos].sequenceLength >= pruning.minSequenceLength)
        {
            vector<int> bestPathNodes;
            for (int pos2 = pos; pos2 != -1; pos2 = pathToMatches[pos2].previousNode)
            {
                bestPathNodes.push_back(pos2);
                pos = pos2 - 1;
            }
            reverse(bestPathNodes.begin(), bestPathNodes.end());

            string finalSequence = nodes[nodesFromReference[bestPathNodes[0]]].GetSequence();
            for (int i = 1; i < bestPathNodes.size(); i++)
            {
                int overlap = pathToMatches[bestPathNodes[i]].overlap;
                finalSequence += pathToMatches[bestPathNodes[i]].lastPathSequence;
                const string& nodeSequence = nodes[nodesFromReference[bestPathNodes[i]]].GetSequence();
                finalSequence += overlap ? nodeSequence.substr(overlap) : nodeSequence;
            }
            result.push_back(finalSequence);
        }
    }
    /*vector<int> bestPathNodes;
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
    }*/
    //return finalSequence;

}