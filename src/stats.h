#ifndef STATS_H
#define STATS_H

#include <map>
#include <fstream>
#include <string>

#include "seq.h"

class Stats
{
public:
    Stats(std::string const & filename) : filename(filename), complete(0), pe(0), se(0) {}

    void update(ReadType type, bool paired = false);

    friend std::ostream & operator << (std::ostream & out, const Stats & stats);

    std::string filename;
    std::map <ReadType, unsigned int> reads;
    unsigned int complete;
    unsigned int pe;
    unsigned int se;
};

std::ostream & operator << (std::ostream & out, const Stats & stats);

#endif // STATS_H
