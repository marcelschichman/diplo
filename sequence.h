#pragma once
#include <fstream>
#include <algorithm>

class Sequence {
public:
    Sequence();
    Sequence(const std::string& _data, std::string _id = std::string());
    Sequence(const Sequence& orig, bool remapReverse = false);
    virtual ~Sequence();

    char* ToDalignFromat();
    const std::string &GetData() const {
        return data;
    }
    const std::string &GetId() const {
        return id;
    }

private:
    char* dalignFromat;
    std::string data;
    std::string id;
};

class FASTA {
public:
    FASTA(const std::string& _filename);

    FASTA& operator>>(Sequence& seq);

    operator bool() const {
        return isOk;
    };
private:
    std::ifstream is;
    std::string filename;
    bool isOk;
};

class FASTQ {
public:
    FASTQ(const std::string& filename);

    FASTQ& operator>>(Sequence& seq);

    operator bool() const {
        return isOk;
    };
private:
    std::ifstream is;
    bool isOk;
};

enum Direction {
    EDIR_FORWARD,
    EDIR_BACKWARD
};