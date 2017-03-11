#include <iostream>
#include <vector>
#include <tuple>
#include "sequence.h"
#include "matchfinder2.h"
#include <list>
#include <forward_list>
#include "tests.h"
using namespace std;

pair<string, string> GetMatchingSequences(Sequence& seq1, Sequence& seq2, Match& m)
{
    string subseq1 = seq1.GetData().substr(m.pos1, m.length);
    string subseq2 = Sequence(seq2, m.reversed).GetData().substr(m.pos2, m.length);
    return make_pair(subseq1, subseq2);
}

int main()
{
    string filename = "/home/marcel/programming/data/PacBio_10kb_CLR.fastq";
    //string filename = "/home/marcel/programming/data/test1.fastq";
    //string filename = "/home/marcel/programming/data/overlap.fastq";

    MatchFinder2 mf(14);
    int numReads = mf.CreateIndex(filename);

    OverlapGraph graph(numReads);
    mf.ProcessMatches(graph);
    mf.Clear();

    for (auto& x : graph.adjacency[2])
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
    }
    int x;
    cin >> x;
}