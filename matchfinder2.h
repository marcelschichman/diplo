#pragma once
#include "matchfinder.h"
#include <sstream>
#include <unordered_map>
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

struct MatchingParams
{
    int kmerLength;

    int numProcessedReadsPerIteration;
    int validKmerMaxOccurrences;
    
    int matchesDiagonalWindowSize;
    int overlappingReadsMinScore;
    float overlappingReadsScoreRatio;

    string ToString()
    {
        stringstream s;
        s << kmerLength << "_";
        s << numProcessedReadsPerIteration << "_";
        s << validKmerMaxOccurrences << "_";
        s << matchesDiagonalWindowSize << "_";
        s << overlappingReadsMinScore << "_";
        s << overlappingReadsScoreRatio;
        return s.str();
    }
};

class MatchFinder2
{
public:
    MatchFinder2(const MatchingParams& _params)
        : params(_params)
        , numReads(0)
    {
        //countPos.resize(1 << (2 * params.kmerLength));
    }

    int CreateIndex(const string& fastq);

    void ProcessMatches(OverlapGraph& graph);

    void Clear();

protected:
    void GetCounts(const string& fastq);
    void GetKmers(const string& fastq);
    
    void ProcessMatches(OverlapGraph& graph, int rangeBegin, int rangeEnd);

    static bool CompareToGroupNicely(const pair<unsigned short, Match>& left, const pair<unsigned short, Match>& right);

    void ExtendMatches(vector<pair<unsigned short, Match>>& oneReadMatches);


    struct Window
    {
        int read;
        bool reversed;
        vector<pair<unsigned short, Match>>::iterator windowEnd, windowStart;
        int score;
    };

    void GetOverlapingReadsWithGoodMatches(int idRead, vector<pair<unsigned short, Match>>& oneReadMatches, vector<pair<unsigned short, vector<Match>>>& neighbors);
    int GetOverlapLength(Window& w, int idReadRef);

    //void GetMatchCounts();


    MatchingParams params;

    long long numReads;
    //unsigned char *matchCounts;
    //vector<int> countPos;
    unordered_map<long long, int> countPos;
    vector<pair<unsigned short, short>> kmers;
    vector<int> readLengths;
};