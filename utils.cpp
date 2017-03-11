#include "utils.h"

using namespace std;


pair<string, string> Utils::GetMatchingSequences(Sequence& seq1, Sequence& seq2, Match& m)
{
    string subseq1 = seq1.GetData().substr(m.pos1, m.length);
    string subseq2 = Sequence(seq2, m.reversed).GetData().substr(m.pos2, m.length);
    return make_pair(subseq1, subseq2);
}