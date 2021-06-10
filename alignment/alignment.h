#ifndef ALIGNMENT_H_
#define ALIGNMENT_H_

#include <iostream>
#include <string>
#include <getopt.h>
#include <vector>
#include <time.h>
#include <algorithm>

namespace alignment {

enum AlignmentType {
    kLocal,
    kGlobal,
    kSemiGlobal
};

enum Direction {
    kLeft,
    kUp,
    kDiagonal,
    kNone
};

struct Cell {
public:
    std::int32_t val_;
    Direction direction_;
    Cell() : val_(0), direction_(kNone) {}
};

class Alignment{
private:
    const char* query;
    const char* target;
    unsigned int query_len;
    unsigned int target_len;
    std::vector<std::vector<Cell>>& matrix;
    int match, mismatch, gap;

private:
    void generateCell(int row, int col);

public:
    CellGenerator(
        const char* query, unsigned int query_len,
        const char* target, unsigned int target_len,
        std::vector<std::vector<Cell>>& matrix,
        int match, int mismatch, int gap
    );

    void generateAllCells(AlignmentType type);
};

int Align(
    const char* query, unsigned int query_len,
    const char* target, unsigned int target_len,
    AlignmentType type,
    int match,
    int mismatch,
    int gap);
}

#endif
