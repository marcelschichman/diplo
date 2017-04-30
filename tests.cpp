#include "tests.h"
#include <vector>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include "sequencegraph.h"
#include "reconstruction.h"
#include "utils.h"
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
    SequenceGraphParams params;
    params.nodeLength = 13;
    params.overlappingKmersMaxExpectedPosDistance = 150;
    params.minKmerOverlap = 2;

    OverlapGraph og(13);
    SequenceGraph sg(og, params);

    vector<SequenceNode> nodes;
    string first = "ACTGCATGCGCTC";
    string second = first.substr(first.length() - 5, 5) + "CCGCACCG";
    nodes.push_back({Utils::SeqToLongLong(first), 0, false, {}});
    nodes.push_back({Utils::SeqToLongLong(second), 0, false, {}});
    nodes.push_back({Utils::SeqToLongLong(first.substr(first.length() - 10, 10) + "CGC"), 0, false, {}});
    //nodes.push_back({first + "CTG", 0, false, {}});
    sg.FindOverlaps(nodes);
    for (auto& overlap : nodes[0].overlaps)
    {
        cout << overlap.first << "\t" << overlap.second << endl;
    }
    cout << endl;
    for (auto& overlap : nodes[2].overlaps)
    {
        cout << overlap.first << "\t" << overlap.second << endl;
    }
    int x;
    cin >> x;
}



void Tests::FindPath()
{
/*
    Sequence read("ACTGCATGCGCTGTAGCT");
    cout << read.GetData() << endl;


    SequenceNode node1 = {"ACTGCATGC", 0, true, {{1, 2}}};
    SequenceNode node2 = {"GCGCTCGAGC", 10, false, {{2, 2}}};
    SequenceNode node3 = {"GCTGTAGCT", 15, true, {}};
*/
    //Sequence read("ACTGCATGCGCTGTAGCTATCGCTAGCTAGCCGCGCGTATTATATCCTTGAGCT");
    // string original = "ACTGCATGCGCTGTAGCTATCGCTAGCTAGCCGCGCGTATTATATCCTTGAGCT";
    // //                     ]      
    // Sequence read("ACTGCATGCGCTGAGCTATCGCGCTAGCCGCGCGATTATTCCTTGAGCT");

    // SequenceNode node0 = {"ACTGCATGC", 0, true, {{1, 4}}};
    // SequenceNode node1 = {"ATGCGCTGT", 0, false, {{2, 4}}};
    // SequenceNode node2 = {"CTGTAGCTAT", 0, false, {{3, 3}}};
    // SequenceNode node3 = {"TATCGCTAGCTA", 0, false, {{4, 3}, {9, 2}}};
    // SequenceNode node4 = {"CTAGCCGCGC", 0, true, {{5, 5}}};
    // SequenceNode node5 = {"CGCGCGTATT", 0, false, {{6, 6}}};
    // SequenceNode node6 = {"CGTATTATAT", 0, false, {{7, 5}}};
    // SequenceNode node7 = {"TATATCCTTGAG", 0, false, {{8, 5}}};
    // SequenceNode node8 = {"TTGAGCT", 0, true, {}};

    // SequenceNode node9 = {"TACGCG", 0, false, {{5, 4}}};
    
    // node8.expectedPos = read.GetData().length() - node8.sequence.length();
    // node4.expectedPos = read.GetData().find(node4.sequence);

    // Scoring s;
    // s.insertion = 1;
    // s.deletion = 1;
    // s.substitution = 1;
    // s.misplacementPenalty = [](int distance) { return 0; };
    // s.overlapPenalty = [](int overlap, int length) { return 0; };
    // Reconstruction r(s);
    // vector<SequenceNode> nodes = {node0, node1, node2, node3, node4, node5, node6, node7, node8, node9};

    // vector<string> result;
    // r.Reconstruct(read, nodes, result);
    // Utils::ResultToOStream(result, cout);
    // if (result[0] == original)
    // {
    //     cout << "sukces" << endl;
    // }
    // int x = 5;
    OverlapGraph og(3);
    og.adjacency[0].push_back({1, vector<Match>()});
    og.adjacency[0][0].second.push_back(Match(0, 0, false, 13));
    og.adjacency[0][0].second.push_back(Match(25, 26, false, 11));
    og.adjacency[0][0].second.push_back(Match(44, 44, false, 10));

    og.adjacency[1].push_back({2, vector<Match>()});
    og.adjacency[1][0].second.push_back(Match(5, 5, false, 13));
    og.adjacency[1][0].second.push_back(Match(10, 10, false, 13));
    og.adjacency[1][0].second.push_back(Match(18, 18, false, 13));
    og.adjacency[1][0].second.push_back(Match(20, 20, false, 13));
    og.adjacency[1][0].second.push_back(Match(23, 23, false, 13));
    og.adjacency[1][0].second.push_back(Match(28, 28, false, 13));
    og.adjacency[1][0].second.push_back(Match(31, 31, false, 13));
    og.adjacency[1][0].second.push_back(Match(36, 36, false, 13));
    og.adjacency[1][0].second.push_back(Match(40, 40, false, 13));

    int idRead = 0;

    SequenceGraphParams params;
    params.nodeLength = 10;
    params.overlappingKmersMaxExpectedPosDistance = 150;
    params.minKmerOverlap = 2;

    SequenceGraph seqGraph(og, params);
    seqGraph.LoadReads("reconstruct_test_reads.fasta");
    vector<SequenceNode> nodes;
    Utils::StartTiming();
    seqGraph.GetNodes(idRead, nodes);
    Utils::VerbalResult("get nodes took ");

    for (int i = 0; i < nodes.size(); i++)
    {
        for (auto& overlap : nodes[i].overlaps)
        {
            cout << nodes[i].GetSequence() << "\t" << nodes[overlap.first].GetSequence() << "\t" << overlap.first << "\t" << overlap.second << endl;
        }
        cout << endl;
    }

    Scoring s;
    s.insertion = 1;
    s.deletion = 1;
    s.substitution = 1;
    s.misplacementPenalty = [](int distance) { return 0; };
    s.overlapPenalty = [](int overlap, int length) { return 0; };

    Pruning p; 
    p.maxDistanceFromFurthest = 10;
    p.minMatches = 10;
    p.minMatchesWindowSize = 30;
    p.skipsAllowed = 2;
    p.minSequenceLength = 40;

    Reconstruction r(s, p);

    vector<string> result;
    r.Reconstruct(seqGraph.forward[idRead], nodes, result);
    Utils::ResultToOStream(result, cout);
    if (result[0] == seqGraph.forward[1].GetData())
    {
        cout << "sukces" << endl;
    }

    string x; getline(cin, x);
}