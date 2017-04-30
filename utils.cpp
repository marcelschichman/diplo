#include "utils.h"
#include <iostream>
#include <fstream>

using namespace std;

vector<chrono::time_point<chrono::system_clock>> tpts;

pair<string, string> Utils::GetMatchingSequences(Sequence& seq1, Sequence& seq2, Match& m)
{
    string subseq1 = seq1.GetData().substr(m.pos1, m.length);
    string subseq2 = Sequence(seq2, m.reversed).GetData().substr(m.pos2, m.length);
    return make_pair(subseq1, subseq2);
}

void Utils::StartTiming() {
    tpts.push_back(chrono::system_clock::now());
}

double Utils::GetTimerResult() {
    double result = ((chrono::duration<double>)(chrono::system_clock::now() - tpts.back())).count();
    tpts.pop_back();
    return result;
}

double Utils::VerbalResult(const string& name) {
    double result = GetTimerResult();
    cout << name << " took " << result << "s" << endl;
    return result;
}

void Utils::NodesToFasta(vector<SequenceNode>& nodes)
{
    ofstream of("dobremiesta.fasta");
    
    for (SequenceNode& n : nodes)
    {
        of << ">MD_" << (int)n.isFromReference << "_" << (int)n.reversed << "_" << n.info << "_" << (int)n.expectedPos << endl;
        of << n.GetSequence() << endl;
    }
}

void Utils::ResultToOStream(vector<string>& result, ostream& stream)
{
    for (int i = 0; i < result.size(); i++)
    {
        stream << ">result" << i << endl;
        for (int j = 0; j < result[i].length(); j += 80)
        {
            stream << result[i].substr(j, min(80, (int)result[i].length() - j)) << endl;
        }
    }    
}

long long Utils::SeqToLongLong(const string& seq)
{
    long long result = 0;
    for (char x : seq)
    {
        result <<= 2;
        switch (x)
        {
            case 'C': result += 1; break;
            case 'G': result += 2; break;
            case 'T': result += 3; break;
        }
    }
    return result;
}

void Utils::ExportFASTQ(vector<string>& reads, const string& filename)
{
    ofstream ofs(filename);
    for (int i = 0; i < reads.size(); i++)
    {
        ofs << "@read" << i << endl;
        ofs << reads[i] << endl;
        ofs << "+" << endl;
        ofs << endl;
    }
}

void Utils::ExportReadsWithOverlaps(OverlapGraph& og, int from, int to)
{
    ofstream ofs("my_overlaps.txt");
    for (int i = from; i < to; i++)
    {
        ofs << i;
        for (auto& n : og.adjacency[i])
        {
            ofs << " " << n.first;
        }
        ofs << endl;
    }
}