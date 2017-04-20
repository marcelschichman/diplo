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
        int posInRead;
        int posInNode;
        int node;
        int previousPos;
    };

    void Reconstruct(const Sequence& read, const vector<SequenceNode>& nodes, vector<string>& result);

    pair<string, int> FindPath(const Sequence& read, const vector<SequenceNode>& nodes, int begin, int end);

    string RecreatePath(const vector<Position>& positions, const vector<SequenceNode>& nodes, int currPos);

    int GetStepLength(const SequenceNode& node, int overlap, int posInRead);

    

    Scoring scoring;
    vector<vector<int>> alignmentMatrix;

};