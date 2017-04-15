#include "filestorage.h"


void FileStorage::Load(istream& file, OverlapGraph& graph)
{
    file >> graph.numReads;
    graph.adjacency.clear();
    graph.adjacency.resize(graph.numReads);

    for (int i = 0; i < graph.numReads; i++)
    {
        int numNeighbors;
        file >> numNeighbors;
        graph.adjacency[i].resize(numNeighbors);
        for (int j = 0; j < numNeighbors; j++)
        {
            unsigned short neighbor;
            int numMatches;
            file >> neighbor >> numMatches;
            graph.adjacency[i][j].first = neighbor;
            graph.adjacency[i][j].second.resize(numMatches);
            for (int k = 0; k < numMatches; k++)
            {
                Load(file, graph.adjacency[i][j].second[k]);
            }
        }
    }
}

void FileStorage::Load(istream& file, Match& match)
{
    file >> match.pos1 >> match.pos2 >> match.reversed >> match.length;
}

void FileStorage::Save(ostream& file, OverlapGraph& graph)
{
    file << graph.numReads << endl;
    for (int i = 0; i < graph.numReads; i++)
    {
        file << graph.adjacency[i].size() << endl;

        for (int j = 0; j < graph.adjacency[i].size(); j++)
        {
            file << graph.adjacency[i][j].first << " " << graph.adjacency[i][j].second.size() << endl;
            for (int k = 0; k < graph.adjacency[i][j].second.size(); k++)
            {
                Save(file, graph.adjacency[i][j].second[k]);
            }
        }
    }
}

void FileStorage::Save(ostream& file, Match& match)
{
    file << match.pos1 << " " << match.pos2 << " " << match.reversed << " " << match.length << endl;
}