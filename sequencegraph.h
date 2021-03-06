#pragma once
#include "matchfinder2.h"

struct SequenceNode
{
    long long sequence;
    int expectedPos;
    bool isFromReference;
    
    // idNode, overlap
    vector<pair<int, int>> overlaps;
    bool reversed;
    int length;
    int info;

    char operator[](unsigned int pos) const
    {
        switch (sequence >> ((length - pos - 1) * 2) & 3)
        {
            case 0: return 'A';
            case 1: return 'C';
            case 2: return 'G';
            case 3: return 'T';
        }
    }
    string GetSequence() const
    {
        string result;
        for (unsigned int i = 0; i < length; i++)
        {
            result.push_back((*this)[i]);
        }
        return result;
    }

    vector<pair<int, int>> overlapsA;
    vector<pair<int, int>> overlapsC;
    vector<pair<int, int>> overlapsG;
    vector<pair<int, int>> overlapsT;
    const vector<pair<int, int>>& GetOverlaps(char nextBase) const;
};

struct SequenceGraphParams
{
    int nodeLength;
    int overlappingKmersMaxExpectedPosDistance;
    int minKmerOverlap;
    int maxEdgesFromNode;
};

class SequenceGraph
{
public:
    SequenceGraph(OverlapGraph& _overlapGraph, const SequenceGraphParams& _params)
        : overlapGraph(_overlapGraph)
        , params(_params)
    {}

    void LoadReads(const string& filename);
    void GetNodes(int idRead, vector<SequenceNode>& nodes);


    vector<Sequence> forward;
    vector<Sequence> reverse;

    void FindOverlaps(vector<SequenceNode>& nodes);

protected:
    void RemoveDuplicates(vector<SequenceNode>& nodes);
    // pair<string, int> GetSequence(int idRead, Match& m, bool reversed);

    int GetSequences(int idRead, Match& m, bool reversed, vector<long long>& sequences);
    


    SequenceGraphParams params;
    OverlapGraph& overlapGraph;
    //int nodeLength;
};