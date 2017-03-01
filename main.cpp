#include <iostream>
#include <vector>
#include <tuple>
#include "sequence.h"
#include "matchfinder.h"
#include <list>
#include <forward_list>
using namespace std;
int main()
{
    MatchFinder mf(13);
/*    Sequence s1("AAATTTCTTTAAAAA");
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

    return 0;*/

    string filename = "/home/marcel/programming/data/PacBio_10kb_CLR.fastq";

    cout << sizeof(forward_list<int>) << endl;


    FASTQ fqReader(filename);
    int idSequence = 0;
    long long sumLengths = 0;
    while (fqReader)
    {
        Sequence s;
        fqReader >> s;
        mf.AddSequence(idSequence++, s);
        sumLengths += s.GetData().length();
        if (idSequence % 1000 == 0)
        {
            cout << "added: " << idSequence << endl;
        }

        if (idSequence > 20000)
            break;
    }
    cout << "num sequences: " << idSequence << endl;
    cout << "sum lengths: " << sumLengths << endl;
        
    map<int, vector<Match>> matches;
    mf.GetMatches(50, matches);

    int matchesCount = 0;

    for (auto matchesPair : matches)
    {
        if (matchesPair.second.size() < 10)
        {
            continue;
        }
        matchesCount += matchesPair.second.size();

        cout << "matches with sequence " << matchesPair.first << endl;
        for (auto match : matchesPair.second)
        {
            cout << match.pos1 << " " << match.pos2 << " " << (match.reversed ? "reversed" : "") << " " << match.length << endl;
        }
    }
    cout << "done" << endl;
}