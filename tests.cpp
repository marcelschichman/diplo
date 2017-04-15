#include "tests.h"
#include <vector>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include "sequencegraph.h"
#include "reconstruction.h"
using namespace Tests;
using namespace std;

void Tests::AllocationOverhead()
{
    int* ptr;
    for (int i = 0; i < 150000000; i++)
    {
        ptr = new int;
    }
}

void Tests::KmerCountDistribution(const string& filename, int size)
{
    FASTQ reader(filename);

    Sequence seq;
    vector<unsigned int> counts(1 << (2 * size), 0);
    while (reader >> seq)
    {
        unsigned int seed = 0;
        unsigned int reversedSeed = 0;
        unsigned int mask = ((long long)1 << (2 * size)) - 1;
        char *data = seq.ToDalignFromat();
        int genomeLength = (int)seq.GetData().length();

        for (int i = 0; i < genomeLength; i++) 
        {
            seed = ((seed << 2) + data[i]) & mask;
            reversedSeed = ((reversedSeed << 2) + 3 - data[genomeLength - i - 1]) & mask;
            
            if (i >= size - 1) 
            {
                counts[seed]++;
                counts[reversedSeed]++;
            }
        }
    }
    map<unsigned int, unsigned int> countCounts;
    for (unsigned int count : counts)
    {
        if (countCounts.find(count) == countCounts.end())
        {
            countCounts[count] = 0;
        }
        countCounts[count]++;
    }
    ofstream o("identical_kmer_count_distribution.txt");
    for (auto& count : countCounts)
    {
        o << count.first << "\t" << count.second << endl;
    }
    o.close();
}

long long Tests::MemoryConsumptionWithStructure(long long (*memorySize)(long long))
{
    ifstream is("identical_kmer_count_distribution.txt");
    long long sum = 0;
    long long count, countCount;
    while (is >> count >> countCount)
    {
        sum += countCount * memorySize(count);
    }
    return sum;
}

void Tests::CompareDataStructures()
{
    cout << "test sum: " << MemoryConsumptionWithStructure([] (long long count) -> long long { return count; }) << endl;
}

void Tests::MatchFinderSimpleTest()
{
    MatchFinder mf(7);
    Sequence s1("AAATTTCTTTAAAAA");
    Sequence s2("CCCCCCCCCTTTCTTTCCCCCCCCCAAAGAAACCCCCCCC");

    
    mf.AddSequence(1, s1);
    mf.AddSequence(2, s2);
    map<int, vector<Match>> matches;
    mf.GetMatches(1, matches);

    int matchesCount = 0;
    for (auto matchesPair : matches)
    {
        matchesCount += matchesPair.second.size();
        cout << "matches with sequence " << matchesPair.first << endl;
        for (auto match : matchesPair.second)
        {
            cout << match.pos1 << " " << match.pos2 << " " << (match.reversed ? "reversed" : "") << endl;
        }
    }
}

void Tests::TestHashmapSize(const string& filename, int size)
{
    FASTQ reader(filename);

    Sequence seq;
    vector<unsigned int> counts(1 << (2 * size), 0);
    while (reader >> seq)
    {
        unsigned int seed = 0;
        unsigned int reversedSeed = 0;
        unsigned int mask = ((long long)1 << (2 * size)) - 1;
        char *data = seq.ToDalignFromat();
        int genomeLength = (int)seq.GetData().length();

        for (int i = 0; i < genomeLength; i++) 
        {
            seed = ((seed << 2) + data[i]) & mask;
            reversedSeed = ((reversedSeed << 2) + 3 - data[genomeLength - i - 1]) & mask;
            
            if (i >= size - 1) 
            {
                counts[seed]++;
                counts[reversedSeed]++;
            }
        }
    }

    unsigned int maxCount = 0;
    for (auto& x : counts)
    {
        maxCount = max(x, maxCount);
    }
    cout << "max count: " << maxCount << endl;
    cin >> maxCount;
}

void Tests::FindOverlaps()
{
    OverlapGraph og(5);
    SequenceGraph sg(og);

    vector<SequenceNode> nodes;
    string first = "ACTGCATGCGCTCGAGC";
    string second = first.substr(first.length() - 5, 5) + "CCGCACC";
    nodes.push_back({first, 0, false, {}});
    nodes.push_back({second, 0, false, {}});
    nodes.push_back({first.substr(first.length() - 10, 10) + "CGCGCGCGCGCGCGACT", 0, false, {}});
    nodes.push_back({first + "CTG", 0, false, {}});
    sg.FindOverlaps(nodes);
    for (auto& overlap : nodes[0].overlaps)
    {
        cout << overlap.first << "\t" << overlap.second << endl;
    }
    for (auto& overlap : nodes[2].overlaps)
    {
        cout << overlap.first << "\t" << overlap.second << endl;
    }
    int x;
    cin >> x;
}

void Tests::GetAlignmentScores()
{
    Scoring s;
    s.insertion = 1;
    s.deletion = 1;
    s.substitution = 1;
    
    Sequence read("ACTGCATGCGCTCGAGC");
    SequenceNode node1 = {"ATGACGC", 0, false, {}};

    Reconstruction r(s);
    vector<pair<int, int>> scores;
    r.GetAlignmentScores(read, 5, node1, 0, scores);
    for (auto& score : scores)
    {
        cout << score.first << ": " << score.second << endl;
    }
}

void Tests::FindPath()
{
    Sequence read("ACTGCATGCGCTCGAGCTGTAGCT");
    cout << read.GetData() << endl;


    SequenceNode node1 = {"ACTGCATGCG", 0, true, {{1, 3}}};
    SequenceNode node2 = {"GCGCTCGAGC", 10, true, {{2, 5}}};
    SequenceNode node3 = {"CGAGCTGTAGCT", 17, true, {}};


    Scoring s;
    s.insertion = 1;
    s.deletion = 1;
    s.substitution = 1;
    s.notFromReferencePenalty = 10;
    s.misplacementPenalty = [](int distance) { return 0; };
    s.overlapPenalty = [](int overlap, int length) { return 0; };
    Reconstruction r(s);

    cout << r.FindPath(read, {node1, node2, node3}) << endl;
}