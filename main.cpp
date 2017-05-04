#include <iostream>
#include <vector>
#include <tuple>
#include "sequence.h"
#include "matchfinder2.h"
#include <list>
#include <forward_list>
#include "tests.h"
#include <unordered_set>
#include <set>
#include "sequencegraph.h"
#include "reconstruction.h"
#include "filestorage.h"
#include "utils.h"
#include <tclap/CmdLine.h>
using namespace std;

bool ParseParams(int argc, char** argv, 
    MatchingParams& matchingParams,
    Scoring& scoring,
    Pruning& pruning,
    SequenceGraphParams& sequenceGraphParams,
    string& filename,
    int& numReadsToFix);

int main(int argc, char** argv)
{
    MatchingParams matchingParams;
    Scoring scoring;
    Pruning pruning; 
    SequenceGraphParams sequenceGraphParams;
    
    string filename = "/home/marcel/programming/diplo/tests/test_dataset.fasta";
    int numReadsToFix = 1;

    if (!ParseParams(argc, argv, matchingParams, scoring, pruning, sequenceGraphParams, filename, numReadsToFix))
    {
        return 1;
    }

    //string filename = "/home/marcel/programming/data/PacBio_10kb_CLR.fastq";
    //string filename = "relevant_reads.fastq";
    //string filename = "/home/marcel/programming/data/test1.fastq";
    //string filename = "/home/marcel/programming/data/overlap.fastq";
    //string filename = "/home/marcel/programming/data/pacbio.100X.fastq";

    OverlapGraph graph;
    if (!FileStorage::Load("overlapgraph.txt", graph))
    {
        MatchFinder2 mf(matchingParams);
        int numReads = mf.CreateIndex(filename);

        graph = OverlapGraph(numReads);
        mf.ProcessMatches(graph);
        mf.Clear();
        FileStorage::Save("overlapgraph.txt", graph);
    }

    Utils::ExportReadsWithOverlaps(graph, 0, graph.adjacency.size());
    //return 0;

    SequenceGraph seqGraph(graph, sequenceGraphParams);
    seqGraph.LoadReads(filename);

    int idRead = 11;

    if (numReadsToFix == -1 || numReadsToFix > seqGraph.forward.size())
    {
        numReadsToFix = seqGraph.forward.size();
    }
    for (int idRead = 0; idRead < numReadsToFix; idRead++)
    {
        
        cerr << "repairing read " << idRead << endl;

        vector<SequenceNode> nodes;
        Utils::StartTiming();
        seqGraph.GetNodes(idRead, nodes);
        Utils::VerbalResult("get nodes took ");

        cerr << "reference length: " << seqGraph.forward[idRead].GetData().length() << endl;
        cerr << "num overlapping reads: " << graph.adjacency[idRead].size() << endl;
        Utils::NodesToFasta(nodes);
        //return 0;

        // scoring.insertion = 4;
        // scoring.deletion = 5;
        // scoring.substitution = 5;
        // scoring.misplacementPenalty = [](int distance) { return 0; };
        // scoring.overlapPenalty = [](int overlap, int length) { return (length - overlap - 1) * 4 ; };

        // pruning.maxDistanceFromFurthest = 30;
        // pruning.minMatches = 11;
        // pruning.minMatchesWindowSize = 20;
        // pruning.skipsAllowed = 2;
        // pruning.minSequenceLength = 40;

        //cout << seqGraph.forward[idRead].GetData().length() << endl;


        Reconstruction r(scoring, pruning);
        vector<string> result;

        Utils::StartTiming();
        r.Reconstruct(seqGraph.forward[idRead], nodes, result);
        Utils::VerbalResult("reconstruct");

        Utils::ResultToOStream(result, cout, seqGraph.forward[idRead].GetId());
        // ofstream of1("results/opravene1.fasta");
        // Utils::ResultToOStream(result, of1);
        int x = 5;
    }
    //cin >> x;
}
//509 3144

bool ParseParams(int argc, char** argv, 
    MatchingParams& matchingParams,
    Scoring& scoring,
    Pruning& pruning,
    SequenceGraphParams& sequenceGraphParams,
    string& filename,
    int& numReadsToFix)
{

    TCLAP::CmdLine cmd(
        "This is a read correction tool for long reads with high error rate. "
        , ' ', "0.9");

    // Matching params
    TCLAP::ValueArg<int> kmerLengthArg("", "kmerLength", "Correction seed length. (13)", false, 13, "int");
    TCLAP::ValueArg<int> readsPerIterationArg("", "readsPerIteration", "Number of reads proccessed in a single run. This affects memory consumption. (10000)", false, 10000, "int");
    TCLAP::ValueArg<int> kmerMaxOccurrencesArg("", "kmerMaxOccurrences", "Maximal number of occurrences of single kmer. (100)", false, 100, "int");
    TCLAP::ValueArg<int> diagonalWindowSizeArg("", "diagonalWindowSize", "Matches diagonal interval size. (50)", false, 50, "int");
    TCLAP::ValueArg<int> overlapMinScoreArg("", "overlapMinScore", "Minimal match length sum of overlapping reads. (60)", false, 60, "int");
    TCLAP::ValueArg<float> overlapScoreRatioArg("", "overlapScoreRatio", "Match length sum must be bigger when the overlap is long. This sets the ratio. (0.005)", false, 0.005, "float");
    cmd.add(overlapMinScoreArg);
    cmd.add(diagonalWindowSizeArg);
    cmd.add(kmerMaxOccurrencesArg);
    cmd.add(readsPerIterationArg);
    cmd.add(kmerLengthArg);
    cmd.add(overlapScoreRatioArg);

    // SequenceGraph params
    TCLAP::ValueArg<int> maxKmerDistanceArg("", "maxKmerDistance", "Maximum expected distance between kmers to create edge in sequence graph. (80)", false, 80, "int");
    TCLAP::ValueArg<int> minKmerOverlapArg("", "minKmerOverlap", "Minimal kmer overlap to create edge in sequence graph. (2)", false, 2, "int");
    TCLAP::ValueArg<int> maxEdgesFromNode("", "maxEdgesFromNode", "When non negative, limits the number of edges from node in sequence graph. (-1)", false, 40, "int");
    cmd.add(maxEdgesFromNode);
    cmd.add(maxKmerDistanceArg);
    cmd.add(minKmerOverlapArg);

    // Scoring
    TCLAP::ValueArg<int> insertionPenaltyArg("", "insertionPenalty", "SCORING: penalty for insertion. (4)", false, 4, "int");
    TCLAP::ValueArg<int> deletionPenaltyArg("", "deletionPenalty", "SCORING: penalty for deletion. (5)", false, 5, "int");
    TCLAP::ValueArg<int> substitutionPenaltyArg("", "substitutionPenalty", "SCORING: penalty for substitution. (5)", false, 5, "int");
    TCLAP::ValueArg<int> maxKmerMisplacementArg("", "maxKmerMisplacement", "SCORING: if kmer's distance from its expected position exceeds this value, the path is penalized. (100)", false, 100, "int");
    TCLAP::ValueArg<float> jumpSizeLinearPenaltyArg("", "jumpSizeLinearPenalty", "SCORING: Jump size linear penalty coefficient. (0.25)", false, 0, "float");
    TCLAP::ValueArg<float> jumpSizeQuadraticPenaltyArg("", "jumpSizeQuadraticPenalty", "SCORING: Jump size quadratic penalty coefficient. (0)", false, 0.2, "float");
    cmd.add(insertionPenaltyArg);
    cmd.add(deletionPenaltyArg);
    cmd.add(substitutionPenaltyArg);
    cmd.add(maxKmerMisplacementArg);
    cmd.add(jumpSizeLinearPenaltyArg);
    cmd.add(jumpSizeQuadraticPenaltyArg);
    

    // Pruning
    TCLAP::ValueArg<int> maxDistFromFurthestArg("", "maxDistFromFurthest", "PRUNING: Maximal distance from path end yo furthest reaching path with equal edit distance. (15)", false, 15, "int");
    TCLAP::ValueArg<int> minMatchesArg("", "minMatches", "PRUNING: Minimal number of matches in last 'minMatchesWindowSize' bases. (11)", false, 11, "int");
    TCLAP::ValueArg<int> minMatchesWindowSizeArg("", "minMatchesWindowSize", "PRUNING: Window size for 'minMatches'. (20)", false, 20, "int");
    TCLAP::ValueArg<int> skipsAllowedArg("", "skipsAllowed", "PRUNING: Maximal number of seeds from reference to be skipped. (3)", false, 3, "int");
    TCLAP::ValueArg<int> minSequenceLengthArg("", "minSequenceLength", "Minimal length of exported sequence. (40)", false, 40, "int");
    cmd.add(maxDistFromFurthestArg);
    cmd.add(minMatchesArg);
    cmd.add(minMatchesWindowSizeArg);
    cmd.add(skipsAllowedArg);
    cmd.add(minSequenceLengthArg);

    TCLAP::ValueArg<int> numReadsToFixArg("", "numReadsToFix", "Number of reads to fix. -1 means all reads. (-1)", false, -1, "int");
    TCLAP::UnlabeledValueArg<string> filenameArg("filename", "Input file name in FASTQ format.", false, "asdf", "filename");
    cmd.add(numReadsToFixArg);
    cmd.add(filenameArg);
    

    try
    {
        cmd.parse(argc, argv);

        // acquiring param values
        matchingParams.kmerLength = kmerLengthArg.getValue();
        matchingParams.numProcessedReadsPerIteration = readsPerIterationArg.getValue();
        matchingParams.validKmerMaxOccurrences = kmerMaxOccurrencesArg.getValue();
        matchingParams.matchesDiagonalWindowSize = diagonalWindowSizeArg.getValue();
        matchingParams.overlappingReadsMinScore = overlapMinScoreArg.getValue();
        matchingParams.overlappingReadsScoreRatio = overlapScoreRatioArg.getValue();

        sequenceGraphParams.nodeLength = kmerLengthArg.getValue();
        sequenceGraphParams.overlappingKmersMaxExpectedPosDistance = maxKmerDistanceArg.getValue();
        sequenceGraphParams.minKmerOverlap = minKmerOverlapArg.getValue();
        sequenceGraphParams.maxEdgesFromNode = maxEdgesFromNode.getValue();

        scoring.insertion = insertionPenaltyArg.getValue();
        scoring.deletion = deletionPenaltyArg.getValue();
        scoring.substitution = substitutionPenaltyArg.getValue();
        scoring.maxKmerMisplacement = maxKmerMisplacementArg.getValue();
        //scoring.misplacementPenalty = [maxKmerDisplacement](int distance) { return (distance > maxKmerDisplacement) * 50; };
        scoring.jumpLinearPenalty = jumpSizeLinearPenaltyArg.getValue();
        scoring.jumpQuadraticPenalty = jumpSizeQuadraticPenaltyArg.getValue();
        // scoring.overlapPenalty = [linearPenalty, quadraticPenalty](int overlap, int length) {
        //     int jumpSize = length - overlap - 1;
        //     return int(linearPenalty * jumpSize + quadraticPenalty * jumpSize * jumpSize); 
        // };

        pruning.maxDistanceFromFurthest = maxDistFromFurthestArg.getValue();
        pruning.minMatches = minMatchesArg.getValue();
        pruning.minMatchesWindowSize = minMatchesWindowSizeArg.getValue();
        pruning.skipsAllowed = skipsAllowedArg.getValue();
        pruning.minSequenceLength = minSequenceLengthArg.getValue();

        filename = filenameArg.getValue();
        numReadsToFix = numReadsToFixArg.getValue();

        return true;
    }

    catch (TCLAP::ArgException& e)
	{
        cout << "ERROR: " << e.error() << " " << e.argId() << endl;
        return false;
    }
}