#pragma once
#include "matchfinder2.h"

struct SequenceNode
{
    string sequence;
    int expectedPos;
    bool isFromReference;
    
    // idNode, overlap
    vector<pair<int, int>> overlaps;
};

class SequenceGraph
{
public:
    SequenceGraph(OverlapGraph& _overlapGraph)
        : overlapGraph(_overlapGraph)
    {}

    void LoadReads(const string& filename);
    void GetNodes(int idRead, vector<SequenceNode>& nodes);
    void FindOverlaps(vector<SequenceNode>& nodes);
    pair<string, int> GetSequence(int idRead, Match& m, bool reversed);



    OverlapGraph& overlapGraph;
    vector<Sequence> forward;
    vector<Sequence> reverse;
};