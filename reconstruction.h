#include "sequencegraph.h"
using namespace std;

struct Scoring
{
    int insertion;
    int deletion;
    int substitution;

    int (*misplacementPenalty)(int);
    int (*overlapPenalty)(int, int); // overlap, length
    
    int maxKmerMisplacement;
    float jumpLinearPenalty;
    float jumpQuadraticPenalty;

};

struct Pruning
{
    // distance from furthest in read
    int maxDistanceFromFurthest;

    // min matches in last n added bases
    int minMatches;
    int minMatchesWindowSize;

    // skips over good places from reference
    int skipsAllowed;

    int minSequenceLength;

    int maxGoodPlaceDistanceToFindPath;
};

class Reconstruction
{
public:
    Reconstruction(const Scoring& _scoring, const Pruning& _pruning)
        : scoring(_scoring)
        , pruning(_pruning)
    {

    }

    struct Position
    {
        int distance;
        int posInRead;
        int posInNode;
        int node;
        int previousPos;

        int numMatches;
        int matchesBits;
    };

    void Reconstruct(const Sequence& read, const vector<SequenceNode>& nodes, vector<string>& result);

    pair<string, int> FindPath(const Sequence& read, const vector<SequenceNode>& nodes, int begin, int end);

    string RecreatePath(const vector<Position>& positions, const vector<SequenceNode>& nodes, int currPos);

    int GetStepLength(const SequenceNode& node, int overlap, int posInRead);

    

    Scoring scoring;
    Pruning pruning;
    vector<vector<int>> alignmentMatrix;

};