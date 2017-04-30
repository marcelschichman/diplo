#include "matchfinder2.h"
#include <iostream>
using namespace std;

int  MatchFinder2::CreateIndex(const string &fastq)
{
    readLengths.clear();
    GetCounts(fastq);

    // counts to positions
    int sum = 0;
    int temp;
    for (auto& x : countPos)
    {
        temp = x;
        x = sum;
        sum += temp;
    }

    kmers.resize(sum);
    GetKmers(fastq);

    // return positions
    for (int i = (int)countPos.size() - 1; i > 0; i--)
    {
        countPos[i] = countPos[i - 1];
    }
    countPos[0] = 0;
    return numReads;
//    cout << "done" << endl;
//    cin >> sum;
//    GetMatchCounts();
}

void MatchFinder2::GetCounts(const string &fastq)
{
    FASTQ reader(fastq);

    Sequence seq;
    while (reader >> seq)
    {
        readLengths.push_back(seq.GetData().length());

        numReads++;
        unsigned int seed = 0;
        unsigned int reversedSeed = 0;
        unsigned int mask = ((long long)1 << (2 * params.kmerLength)) - 1;
        char *data = seq.ToDalignFromat();
        int genomeLength = (int)seq.GetData().length();

        for (int i = 0; i < genomeLength; i++)
        {
            seed = ((seed << 2) + data[i]) & mask;
            reversedSeed = ((reversedSeed << 2) + 3 - data[genomeLength - i - 1]) & mask;

            if (i >= params.kmerLength - 1)
            {
                countPos[seed]++;
                countPos[reversedSeed]++;
            }
        }
    }
}



void MatchFinder2::GetKmers(const string& fastq)
{
    FASTQ reader(fastq);

    Sequence seq;
    int id = 0;
    while (reader >> seq)
    {
        unsigned int seed = 0;
        unsigned int reversedSeed = 0;
        unsigned int mask = ((long long)1 << (2 * params.kmerLength)) - 1;
        char *data = seq.ToDalignFromat();
        int genomeLength = (int)seq.GetData().length();

        for (int i = 0; i < genomeLength; i++)
        {
            seed = ((seed << 2) + data[i]) & mask;
            reversedSeed = ((reversedSeed << 2) + 3 - data[genomeLength - i - 1]) & mask;

            if (i >= params.kmerLength - 1)
            {
                if (countPos[seed] >= kmers.size())
                {
                    cout << "megapipkos" << endl;
                }
                kmers[countPos[seed]++] = {(unsigned short)id, i - params.kmerLength + 1};
                if (countPos[reversedSeed] >= kmers.size())
                {
                    cout << "megapipkos" << endl;
                }
                kmers[countPos[reversedSeed]++] = {(unsigned short)id, -(i - params.kmerLength + 1) - 1};
            }
        }
        id++;
    }
}
/*
void MatchFinder2::GetMatchCounts()
{
    if (matchCounts != NULL)
    {
        delete matchCounts;
    }
    matchCounts = new unsigned char[numReads * numReads];

    long long numOverflowed = 0;
    for (int i = 0; i < (int)countPos.size(); i++)
    {
        int beginRegion = countPos[i];
        int endRegion = (i < (int)countPos.size() - 1) ? countPos[i + 1] : (int)kmers.size();
        int numPositions = endRegion - beginRegion;
        if (numPositions > 20)
        {
            continue;
        }
        for (int k = beginRegion; k < endRegion; k++)
        {
            for (int l = beginRegion; l < endRegion; l++)
            {
                if (kmers[k].first != kmers[l].first)
                {
                    if (matchCounts[kmers[k].first * numReads + kmers[l].first] == 255)
                    {
                        numOverflowed++;
                    }
                    matchCounts[kmers[k].first * numReads + kmers[l].first]++;
                }
            }
        }        
    }
    cout << "num overflowed: " << numOverflowed << endl;

    long long matchCount = 0;
    for (int i = 0; i < numReads; i++)
    {
        for (int j = 0; j < numReads; j++)
        {
            matchCount += matchCounts[i * numReads + j];
        }
    }
    cout << "total matches: " << matchCount << endl;
    cin >> matchCount;
}*/

void MatchFinder2::ProcessMatches(OverlapGraph& graph)
{
    //const int readsPerIteration = 10000;
    OverlapGraph tempGraph(graph.numReads);
    for (int rangeBegin = 0; rangeBegin < numReads; rangeBegin += params.numProcessedReadsPerIteration)
    {
        int rangeEnd = min(rangeBegin + params.numProcessedReadsPerIteration, (int)numReads);
        cerr << "region: " << rangeBegin << " -> " << rangeEnd << endl;
        ProcessMatches(tempGraph, rangeBegin, rangeEnd);
    }
    graph = tempGraph;
}

void MatchFinder2::ProcessMatches(OverlapGraph& graph, int rangeBegin, int rangeEnd)
{
    //const int validKmerMaxOccurences = 100;

    vector<vector<pair<unsigned short, Match>>> matches (rangeEnd - rangeBegin);

    for (int i = 0; i < (int)countPos.size(); i++)
    {
        int beginKmerRegion = countPos[i];
        int endKmerRegion = (i < (int)countPos.size() - 1) ? countPos[i + 1] : (int)kmers.size();
        int numPositions = endKmerRegion - beginKmerRegion;
        if (numPositions > params.validKmerMaxOccurrences)
        {
            continue;
        }
        for (int k = beginKmerRegion; k < endKmerRegion; k++)
        {
            if (kmers[k].first >= rangeBegin && kmers[k].first < rangeEnd && kmers[k].second >= 0)
            {
                for (int l = beginKmerRegion; l < endKmerRegion; l++)
                {
                    if (kmers[k].first != kmers[l].first)
                    {
                        if (kmers[l].second >= 0)
                        {
                            matches[kmers[k].first - rangeBegin].push_back({kmers[l].first, Match(kmers[k].second, kmers[l].second, false, params.kmerLength)});
                        }
                        else
                        {
                            matches[kmers[k].first - rangeBegin].push_back({kmers[l].first, Match(kmers[k].second, -kmers[l].second - 1, true, params.kmerLength)});
                        }
                    }
                }
            }
        }        
    }

    int readId = rangeBegin;
    for (auto& oneReadMatches : matches)
    {
        sort(oneReadMatches.begin(), oneReadMatches.end(), CompareToGroupNicely);

        ExtendMatches(oneReadMatches);

        GetOverlapingReadsWithGoodMatches(readId, oneReadMatches, graph.adjacency[readId]);

        readId++;
    }
    matches = decltype(matches)();
}

bool MatchFinder2::CompareToGroupNicely(const pair<unsigned short, Match>& left, const pair<unsigned short, Match>& right)
{
    // order by matched read, orientation, diagonal, position within reference read
    int leftDiagonal = left.second.pos2 - left.second.pos1;
    int rightDiagonal = right.second.pos2 - right.second.pos1;
    return 
    left.first != right.first ? 
        (left.first < right.first) 
        :
        (
            left.second.reversed != right.second.reversed ?
                left.second.reversed == false
                :
                (
                    leftDiagonal != rightDiagonal ?
                        leftDiagonal < rightDiagonal
                        :
                        left.second.pos1 < right.second.pos1
                )
        );
}

void MatchFinder2::ExtendMatches(vector<pair<unsigned short, Match>>& oneReadMatches)
{
    if (oneReadMatches.empty())
    {
        return;
    }
    vector<pair<unsigned short, Match>> extendedMatches;
    pair<unsigned short, Match> lastMatch = oneReadMatches[0];
    for (int i = 1; i < (int)oneReadMatches.size(); i++)
    {
        auto& readMatch(oneReadMatches[i]);
        auto& prevMatch(oneReadMatches[i - 1]);

        if (readMatch.first == prevMatch.first &&
            readMatch.second.reversed == prevMatch.second.reversed &&
            readMatch.second.pos1 == prevMatch.second.pos1 + 1 &&
            readMatch.second.pos2 == prevMatch.second.pos2 + 1 )
        {
            lastMatch.second.length++;
        }
        else
        {
            extendedMatches.push_back(lastMatch);
            lastMatch = readMatch;
        }
    }
    extendedMatches.push_back(lastMatch);
    oneReadMatches.swap(extendedMatches);
}


void MatchFinder2::GetOverlapingReadsWithGoodMatches(int idRead, vector<pair<unsigned short, Match>>& oneReadMatches, vector<pair<unsigned short, vector<Match>>>& neighbors)
{
    // const int windowSize = 50;
    // const int overlappingReadsMinScore = 100;
    // const float overlappingReadsScoreRatio = 0.01;

    vector<Window> bestWindowsPerRead;
    
    int currentWindowScore = 0;
    for (auto it1 = oneReadMatches.begin(), it2 = it1; it1 != oneReadMatches.end(); it1++)
    {   
        currentWindowScore += it1->second.length;
                            // continue while:
        while (it1 != it2   // the window has more than one match
            && (it1->first != it2->first // window has matches from different reads
                || it1->second.reversed != it2->second.reversed // has matches of different orientation
                || (it1->second.pos2 - it1->second.pos1) > (it2->second.pos2 - it2->second.pos1) + params.matchesDiagonalWindowSize // match diagonals are too far apart
                )
            )
        {
            currentWindowScore -= it2->second.length;
            it2++;
        }

        if (bestWindowsPerRead.empty() || bestWindowsPerRead.back().read != it1->first)
        {
            bestWindowsPerRead.push_back({it1->first, it1->second.reversed, it1, it2, currentWindowScore});
        }
        else if (bestWindowsPerRead.back().score < currentWindowScore)
        {
            bestWindowsPerRead.back().reversed = it1->second.reversed;
            bestWindowsPerRead.back().score = currentWindowScore;
            bestWindowsPerRead.back().windowEnd = it1;
            bestWindowsPerRead.back().windowStart = it2;
        }
    }

    neighbors.clear();
    for (Window& w : bestWindowsPerRead)
    {
        if (w.score > params.overlappingReadsMinScore + int(GetOverlapLength(w, idRead) * params.overlappingReadsScoreRatio))
        {
            vector<Match> goodMatchesWithNeighbor;
            for (auto it = w.windowStart; it <= w.windowEnd; it++)
            {
                goodMatchesWithNeighbor.push_back(it->second);
            }
            neighbors.push_back({w.windowEnd->first, goodMatchesWithNeighbor});
        }
    }
}

int MatchFinder2::GetOverlapLength(Window& w, int idReadRef)
{
    int offset = (w.windowStart->second.pos1 - w.windowStart->second.pos2 + w.windowEnd->second.pos1 - w.windowEnd->second.pos2) / 2;
    int overlapBegin = max(0, offset);
    int overlapEnd = min(readLengths[idReadRef], readLengths[w.windowEnd->first] + offset);
    return max(0, overlapEnd - overlapBegin);
}

void MatchFinder2::Clear()
{
    vector<int>().swap(countPos);
    decltype(kmers)().swap(kmers);
}