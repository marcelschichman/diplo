#include "sequencegraph.h"
using namespace std;

struct Scoring
{
    int insertion;
    int deletion;
    int substitution;

    int (*misplacementPenalty)(int);
    int (*overlapPenalty)(int, int); // overlap, length
    
};

class Reconstruction
{
public:
    Reconstruction(const Scoring& _scoring)
        : scoring(_scoring)
    {

    }

    struct Position
    {
        int distance;
        int previousNode;
        int overlapWithPrevious;
    };

    string FindPath(const Sequence& read, const vector<SequenceNode>& nodes);

    void GetAlignmentScores(const Sequence& read, int readPos, const SequenceNode& node, int nodePos, vector<pair<int, int>>& scores);


    string RecreatePath(const vector<vector<Position>>& positions, const vector<SequenceNode>& nodes, int posInRead, int node);

    int GetStepLength(const SequenceNode& node, int overlap, int posInRead);

    Scoring scoring;
    vector<vector<int>> alignmentMatrix;

};