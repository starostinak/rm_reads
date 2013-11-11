#ifndef SEQ_H
#define SEQ_H

#include <string>
#include <fstream>

enum ReadType{
    ok,
    adapter = 1,
    n,
    polyG,
    polyC,
    length
};

void init_type_names(int length, int polyG);
const std::string & get_type_name (ReadType type);

class Seq {
public:
    Seq() {}

    bool read_seq(std::ifstream & fin)
    {
        std::string tmp;
        std::getline(fin, id);
        if (!fin) {
            return false;
        }
        std::getline(fin, seq);
        std::getline(fin, tmp);
        std::getline(fin, qual);
        return true;
    }

    void write_seq(std::ofstream & fout)
    {
        fout << id << std::endl;
        fout << seq << std::endl;
        fout << '+' << id.substr(1) << std::endl;
        fout << qual << std::endl;
    }

    void update_id(ReadType type)
    {
        id.insert(1, get_type_name(type) + "__");
    }

    std::string id;
    std::string seq;
    std::string qual;
};

#endif // SEQ_H
