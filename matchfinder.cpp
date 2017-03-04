#include "matchfinder.h"
#include <algorithm>


void MatchFinder::AddSequence(int id, Sequence & seq)
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
            AddKmer(seed, (short)id, i - length + 1);
            AddKmer(reversedSeed, (short)id, -(i - length + 1) - 1);
        }
    }
}

void MatchFinder::AddKmer(unsigned int kmer, short id, short pos)
{
/*    auto& pVec(kmers[kmer]);
    if (pVec == NULL)
    {
        pVec = new vector<pair<short, short>>();
    }
    pVec->push_back(make_pair(id, pos));*/
    kmers[kmer].push_front(make_pair(id, pos));
}

void MatchFinder::Process()
{
    /*for (unsigned int kmer = 0; kmer < kmers.size(); kmer++)
    {
        auto& kmerRef(kmers[kmer]);
        if (kmerRef != NULL)
        {
            //for (auto  tu pokracuj)
        }
        
    }*/
}

void MatchFinder::GetMatches(int id, map<int, vector<Match>> & matches)
{
    matches.clear();
    for (unsigned int kmer = 0; kmer < kmers.size(); kmer++)
    {
        auto& kmerRef(kmers[kmer]);
        if (!kmerRef.empty())
        {
            int refPos = 0;
            bool found = false;
            for (auto pos : kmerRef)
            {
                if (pos.first == id)
                {
                    refPos = pos.second;
                    found = true;
                    break;
                }
            }
            if (found && refPos >= 0)
            {
                for (auto pos : kmerRef)
                {
                    if (pos.first != id)
                    {
                        if (pos.second >= 0)
                        {
                            matches[pos.first].push_back(Match(refPos, pos.second, false, length));
                        }
                        else
                        {
                            matches[pos.first].push_back(Match(refPos, -pos.second - 1, true, length));
                        }
                    }
                }
            }
        }        
    }
    for (auto& it : matches)
    {
        //ExtendMatches(it.second);
    }
}

void MatchFinder::ExtendMatches(vector<Match>& matches)
{
    sort(matches.begin(), matches.end(), [](const Match& left, const Match& right) -> bool 
    {
        int leftDiagonal = left.pos2 - left.pos1;
        int rightDiagonal = right.pos2 - right.pos1;
        return left.reversed == right.reversed ? ((leftDiagonal == rightDiagonal) ? (left.pos1 < right.pos1) : (leftDiagonal < rightDiagonal))
            : left.reversed == true;
    });
    vector<Match> newMatches;
    for (auto it = matches.begin(); it != matches.end(); )
    {
        int currentLength = it->length;
        auto following = it + 1;
        auto prev = it;
        while (following != matches.end() && following->pos1 == prev->pos1 + 1 && following->pos2 == prev->pos2 + 1 && prev->reversed == following->reversed)
        {
            currentLength++;
            prev = following;
            following++;
        }
        newMatches.push_back(*it);
        newMatches.back().length = currentLength;
        it = following;
    }
    matches.swap(newMatches);
    //newMatches.clear();
}

bool MatchFinder::AreOnTheSameDiagonal(const Match& m1, const Match& m2)
{
    const int minDistance = 50;
    const double diagonalDistanceRatio = 0.05;
    int distance = min(abs(m1.pos1 - m2.pos1), abs(m1.pos2 - m2.pos2));
    int diagonalDistance = abs(m1.pos1 - m2.pos1 - (m1.pos2 - m2.pos2));

    return distance >= minDistance && diagonalDistance <= distance * diagonalDistanceRatio;
}

void MatchFinder::GetGoodMatches(vector<Match>& matches, vector<Match>& goodMatches)
{
    goodMatches.clear();
    // todo: use set
    for (Match& m : matches)
    {
        for (Match& m2 : matches)
        {
            if (AreOnTheSameDiagonal(m, m2))
            {
                goodMatches.push_back(m);
                break;
            }
        }
    }
}