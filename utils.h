#pragma once
#include "matchfinder.h"

namespace Utils
{
    pair<string, string> GetMatchingSequences(Sequence& seq1, Sequence& seq2, Match& m);
}