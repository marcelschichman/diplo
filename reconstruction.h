#include "sequencegraph.h"
using namespace std;

struct Scoring
{
    int insertion;
    int deletion;
    int substitution;

    int (*misplacementPenalty)(int);
    int (*overlapPenalty)(int, int); // overlap, length
    
    int notFromReferencePenalty;
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

    void FindPath(Sequence& read, vector<SequenceNode>& nodes);

    void GetAlignmentScores(Sequence& read, int readPos, SequenceNode& node, int nodePos, vector<pair<int, int>>& scores);


    string RecreatePath(vector<vector<Position>>& positions, vector<SequenceNode>& nodes, int posInRead, int node);

    int GetStepLength(SequenceNode& node, int overlap, int posInRead);

    Scoring scoring;
    vector<vector<int>> alignmentMatrix;

};