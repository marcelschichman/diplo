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
    auto& pVec(kmers[kmer]);
    if (pVec == NULL)
    {
        pVec = new vector<pair<short, short>>();
    }
    pVec->push_back(make_pair(id, pos));
}

void MatchFinder::Process()
{
    for (unsigned int kmer = 0; kmer < kmers.size(); kmer++)
    {
        auto& kmerRef(kmers[kmer]);
        if (kmerRef != NULL)
        {
            //for (auto /* tu pokracuj*/)
        }
        
    }
}

void MatchFinder::GetMatches(int id, map<int, vector<Match>> & matches)
{
    matches.clear();
    for (unsigned int kmer = 0; kmer < kmers.size(); kmer++)
    {
        auto& kmerRef(kmers[kmer]);
        if (kmerRef != NULL)
        {
            int refPos = 0;
            bool found = false;
            for (auto pos : *kmerRef)
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
                for (auto pos : *kmerRef)
                {
                    if (pos.first != id)
                    {
                        if (pos.second >= 0)
                        {
                            matches[pos.first].push_back(Match(refPos, pos.second, false));
                        }
                        else
                        {
                            matches[pos.first].push_back(Match(refPos, -pos.second - 1, true));
                        }
                    }
                }
            }
        }        
    }
}