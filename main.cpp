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

    MatchFinder2 mf(13);
    int numReads = mf.CreateIndex(filename);

    OverlapGraph graph(numReads);
    mf.ProcessMatches(graph);
}