#include "matchfinder2.h"
#include <iostream>
using namespace std;

void MatchFinder2::CreateIndex(const string &fastq)
{
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
    cout << "done" << endl;
    cin >> sum;
    GetMatchCounts();
}

void MatchFinder2::GetCounts(const string &fastq)
{
    FASTQ reader(fastq);

    Sequence seq;
    while (reader >> seq)
    {
        numReads++;
        unsigned int seed = 0;
        unsigned int reversedSeed = 0;
        unsigned int mask = ((long long)1 << (2 * length)) - 1;
        char *data = seq.ToDalignFromat();
        int genomeLength = (int)seq.GetData().length();

        for (int i = 0; i < genomeLength; i++)
        {
            seed = ((seed << 2) + data[i]) & mask;
            reversedSeed = ((reversedSeed << 2) + 3 - data[genomeLength - i - 1]) & mask;

            if (i >= length - 1)
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
        unsigned int mask = ((long long)1 << (2 * length)) - 1;
        char *data = seq.ToDalignFromat();
        int genomeLength = (int)seq.GetData().length();

        for (int i = 0; i < genomeLength; i++)
        {
            seed = ((seed << 2) + data[i]) & mask;
            reversedSeed = ((reversedSeed << 2) + 3 - data[genomeLength - i - 1]) & mask;

            if (i >= length - 1)
            {
                if (countPos[seed] >= kmers.size())
                {
                    cout << "megapipkos" << endl;
                }
                kmers[countPos[seed]++] = {(unsigned short)id, i - length + 1};
                if (countPos[reversedSeed] >= kmers.size())
                {
                    cout << "megapipkos" << endl;
                }
                kmers[countPos[reversedSeed]++] = {(unsigned short)id, -(i - length + 1) - 1};
            }
        }
        id++;
    }
}

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
}

void MatchFinder2::ProcessMatches(OverlapGraph& graph)
{
    const int readsPerIteration = 5000;
    
    for (int rangeBegin = 0; rangeBegin < numReads; rangeBegin += readsPerIteration)
    {
        int rangeEnd = min(rangeBegin + readsPerIteration, (int)numReads);
        ProcessMatches(graph, rangeBegin, rangeEnd);
    }
}

void MatchFinder2::ProcessMatches(OverlapGraph& graph, int rangeBegin, int rangeEnd)
{
    const int validKmerMaxOccurences = 100;

    vector<vector<pair<short, Match>>> matches (rangeEnd - rangeBegin);

    for (int i = 0; i < (int)countPos.size(); i++)
    {
        int beginKmerRegion = countPos[i];
        int endKmerRegion = (i < (int)countPos.size() - 1) ? countPos[i + 1] : (int)kmers.size();
        int numPositions = endKmerRegion - beginKmerRegion;
        if (numPositions > validKmerMaxOccurences)
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
                            matches[kmers[k].first - rangeBegin].push_back({kmers[l].first, Match(kmers[k].second, kmers[l].second, false, length)});
                        }
                        else
                        {
                            matches[kmers[k].first - rangeBegin].push_back({kmers[l].first, Match(kmers[k].second, -kmers[l].second - 1, true, length)});
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

        GetOverlapingReadsWithGoodMatches(oneReadMatches, graph.adjacency[readId++]);
    }
}

bool MatchFinder2::CompareToGroupNicely(const pair<short, Match>& left, const pair<short, Match>& right)
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

void MatchFinder2::ExtendMatches(vector<pair<short, Match>>& oneReadMatches)
{
    //...
}


void MatchFinder2::GetOverlapingReadsWithGoodMatches(vector<pair<short, Match>>& oneReadMatches, vector<pair<short, vector<Match>>>& neighbors)
{

}