#include <iostream>
#include <vector>
#include <tuple>
#include "sequence.h"
#include "matchfinder2.h"
#include <list>
#include <forward_list>
#include "tests.h"
#include <unordered_set>
#include <set>
#include "sequencegraph.h"
#include "reconstruction.h"
#include "filestorage.h"
#include "utils.h"
using namespace std;

pair<string, string> GetMatchingSequences(Sequence& seq1, Sequence& seq2, Match& m)
{
    string subseq1 = seq1.GetData().substr(m.pos1, m.length);
    string subseq2 = Sequence(seq2, m.reversed).GetData().substr(m.pos2, m.length);
    return make_pair(subseq1, subseq2);
}

int main()
{
    //Tests::FindPath();
    //return 0;

    string filename = "/home/marcel/programming/data/PacBio_10kb_CLR.fastq";
    //string filename = "/home/marcel/programming/data/test1.fastq";
    //string filename = "/home/marcel/programming/data/overlap.fastq";

    OverlapGraph graph;
    if (!FileStorage::Load("overlapgraph.txt", graph))
    {
        MatchFinder2 mf(13);
        int numReads = mf.CreateIndex(filename);

        graph = OverlapGraph(numReads);
        mf.ProcessMatches(graph);
        mf.Clear();
        FileStorage::Save("overlapgraph.txt", graph);
    }

    /*for (auto& x : graph.adjacency[35146])
    {
        int score = 0;
        cout << x.first << endl;
        for (auto& m : x.second)
        {
            cout << m.pos1 << "\t" << m.pos2 << "\t" << m.length << "\t" << (m.pos1 - m.pos2) << endl;
            score += m.length;
        }
        cout << "score: " << score << endl;
        cout << endl;
    }*/
/*
    unsigned long long sumOverlapping = 0;
    unsigned long long totalMatches = 0;
    int readId = 0;
    for (auto& x :graph.adjacency)
    {
        sumOverlapping += x.size();
        for (auto& y : x)
        {
            totalMatches += y.second.size();
            for (auto& z : y.second)
            {
                cout << readId << "\t" << y.first << "\t" << z.pos1 << "\t" << z.pos2 << endl;
            }
        }
        readId++;
    }
    cout << "average num overlapping: " << (double(sumOverlapping) / graph.adjacency.size()) << endl;
    cout << "total matches: " << totalMatches << endl;

    cin >> x;*/

    int idRead = 1001;

    SequenceGraph seqGraph(graph);
    seqGraph.LoadReads(filename);
    vector<SequenceNode> nodes;
    Utils::StartTiming();
    seqGraph.GetNodes(idRead, nodes);
    Utils::VerbalResult("get nodes took ");

    cout << "reference length: " << seqGraph.forward[idRead].GetData().length() << endl;
    cout << "num overlapping reads: " << graph.adjacency[idRead].size() << endl;
    Utils::NodesToFasta(nodes);
    //return 0;

    Scoring s;
    s.insertion = 1;
    s.deletion = 1;
    s.substitution = 1;
    s.misplacementPenalty = [](int distance) { return 0; };
    s.overlapPenalty = [](int overlap, int length) { return 0; };

    cout << seqGraph.forward[idRead].GetData().length() << endl;

    Reconstruction r(s);

    cout << r.Reconstruct(seqGraph.forward[idRead], nodes) << endl;
    int x = 5;
}
//509 3144