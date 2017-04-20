#pragma once
#include "matchfinder.h"
#include <chrono>
#include "sequencegraph.h"

namespace Utils
{
    pair<string, string> GetMatchingSequences(Sequence& seq1, Sequence& seq2, Match& m);

    void NodesToFasta(vector<SequenceNode>& nodes);

    void StartTiming();
    double GetTimerResult();
    double VerbalResult(const string& name);

    void ResultToOStream(vector<string>& result, ostream& stream);
}