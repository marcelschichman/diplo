#include <vector>
#include <tuple>
#include "sequence.h"
#include <map>
using namespace std;

struct Match {
    Match(int _pos1, int _pos2, bool _reversed)
        : pos1(_pos1)
        , pos2(_pos2)
        , reversed(_reversed)
        {}
    int pos1;
    int pos2;
    bool reversed;
    int kmer;
    operator pair<int, int>() {
        return make_pair(pos1, pos2);
    }
};

class MatchFinder
{
public:
    MatchFinder(int _length)
        : length(_length)
        {
            kmers.resize(1 << (2 * length));
        }
    void AddSequence(int id, Sequence & seq);
    void Process();
    void GetMatches(int id1, int id2, vector<Match> & matches);
    void GetMatches(int id, map<int, vector<Match>> & matches);

protected:
    int length;

    // tuples are <sequence id, position within sequence, direction>
    vector<vector<pair<short, short>>*> kmers;

    map<pair<short, short>, vector<Match>> pairs; 
    //unordered_set<long long, vector<tuple<int, int, bool>>> kmers;
    void AddKmer(unsigned int kmer, short id, short pos);
};