#pragma once
#include "matchfinder.h"
using namespace std;

struct OverlapGraph
{
    OverlapGraph()
        : numReads(0)
    {
    }
    OverlapGraph(int _numReads)
        : numReads(_numReads)
    {
        adjacency.resize(numReads);
    }

    int numReads;
    vector<vector<pair<unsigned short, vector<Match>>>> adjacency;
};

class MatchFinder2
{
public:
    MatchFinder2(int _length)
        : length(_length)
        , numReads(0)
//        , matchCounts(NULL)
    {
        countPos.resize(1 << (2 * length));
    }
    int CreateIndex(const string& fastq);


    void GetCounts(const string& fastq);
    void GetKmers(const string& fastq);
    
    void ProcessMatches(OverlapGraph& graph);
    void ProcessMatches(OverlapGraph& graph, int rangeBegin, int rangeEnd);

    static bool CompareToGroupNicely(const pair<unsigned short, Match>& left, const pair<unsigned short, Match>& right);

    void ExtendMatches(vector<pair<unsigned short, Match>>& oneReadMatches);
    void GetOverlapingReadsWithGoodMatches(vector<pair<unsigned short, Match>>& oneReadMatches, vector<pair<unsigned short, vector<Match>>>& neighbors);

    //void GetMatchCounts();

    void Clear();

    int length;
    long long numReads;
    //unsigned char *matchCounts;
    vector<int> countPos;
    vector<pair<unsigned short, short>> kmers;
};