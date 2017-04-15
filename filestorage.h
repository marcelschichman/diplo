#include "matchfinder2.h"
#include <fstream>

class FileStorage
{
public:
    static bool Load(const string& filename, OverlapGraph& graph)
    {
        ifstream ifs(filename);
        if (ifs.is_open())
        {
            Load(ifs, graph);
            return true;
        }
        return false;
    }
    static void Save(const string& filename, OverlapGraph& graph)
    {
        ofstream ofs(filename);
        Save(ofs, graph);
    }

    static void Load(istream& file, OverlapGraph& graph);
    static void Load(istream& file, Match& match);

    static void Save(ostream& file, OverlapGraph& graph);
    static void Save(ostream& file, Match& match);
};