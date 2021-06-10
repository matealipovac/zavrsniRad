#include "alignment.h"

namespace alignment
{

    Alignment::Alignment(
        const char *query, unsigned int query_len,
        const char *target, unsigned int target_len,
        std::vector<std::vector<Cell>> &matrix,
        int match, int mismatch, int gap) : query(query),
                                            query_len(query_len),
                                            target(target),
                                            target_len(target_len),
                                            matrix(matrix), match(match),
                                            mismatch(mismatch),
                                            gap(gap)
    {
    }

    void Alignment::generateCell(int i, int j)
    {

        int diagonalValue = matrix[i - 1][j - 1].val_;
        diagonalValue += (query[i - 1] == target[j - 1]) ? match : mismatch;
        int topValue = matrix[i - 1][j].val_ + gap;
        int leftValue = matrix[i][j - 1].val_ + gap;

        int maxValue;
        Direction dir;
        if (diagonalValue >= leftValue && diagonal >= topValue)
        {
            dir = kDiagonal;
            maxValue = diagonalValue;
        }
        else if (leftValue >= diagonalValue && leftValue >= topValue)
        {
            dir = kLeft;
            maxValue = leftValue;
        }
        else if (topValue >= leftValue && topValue >= diagonalValue)
        {
            dir = kUp;
            maxValue = topValue;
        }
        matrix[i][j].val_ = maxValue;
        matrix[i][j].direction_ = dir;
    }

    void Alignment::generateAllCells(AlignmentType type)
    {
        for (int i = 1; i < matrix.size(); i++)
        {
            for (int j = 1; j < matrix[0].size(); j++)
            {
                generateCell(i, j);
                if (type == kLocal && matrix[i][j].val_ <= 0)
                {
                    matrix[i][j].val_ = 0;
                    matrix[i][j].direction_ = kNone;
                }
            }
        }
    }

    int Align(
        const char *query, unsigned int query_len,
        const char *target, unsigned int target_len,
        AlignmentType type,
        int match,
        int mismatch,
        int gap)
    {
        std::vector<std::vector<Cell>> matrix = std::vector<std::vector<Cell>>(query_len + 1, std::vector<Cell>(target_len + 1));

        //initialize matrix
        int num_of_rows = matrix.size();
        int num_of_cols = matrix[0].size();
        if (type == kLocal || type == kSemiGlobal)
            gap = 0;
        Direction top_row_direction = kNone;
        Direction first_col_direction = kNone;
        if (type == kGlobal)
        {
            top_row_direction = kLeft;
            first_col_direction = kUp;
        }

        for (int i = 1; i < num_of_rows; i++)
        {
            matrix[i][0].val_ = i * gap;
            matrix[i][0].direction_ = first_col_direction;
        }
        for (int i = 1; i < num_of_cols; i++)
        {
            matrix[0][i].val_ = i * gap;
            matrix[0][i].direction_ = top_row_direction;
        }
        int rowLen = matrix.size();
        int colLen = matrix[0].size();

        CellGenerator generator = CellGenerator (query, query_len, target, target_len, matrix, match, mismatch, gap);
        generator.generateAllCells(type);

        int result;
        int maxRow, maxCol;
        switch (type)
        {
        case kLocal:
            result = matrix[0][0].val_;
            maxRow = 0;
            maxCol = 0;
            for (int i = 0; i < rowLen; i++)
            {
                for (int j = 0; j < colLen; j++)
                {
                    if (matrix[i][j].val_ >= result)
                    {
                        result = matrix[i][j].val_;
                        maxRow = i;
                        maxCol = j;
                    }
                }
            }
            break;

        case kGlobal:
            maxRow = query_len;
            maxCol = target_len;
            result = matrix[query_len][target_len].val_;
            break;

        case kSemiGlobal:
            result = matrix[0][target_len].val_;
            maxRow = 0;
            maxCol = target_len;
            for (int i = 0; i < rowLen; i++)
            {
                if (matrix[i][target_len].val_ >= result)
                {
                    result = matrix[i][target_len].val_;
                    maxRow = i;
                    maxCol = target_len;
                }
            }
            for (int i = 0; i < colLen; i++)
            {
                if (matrix[query_len][i].val_ >= result)
                {
                    result = matrix[query_len][i].val_;
                    maxRow = query_len;
                    maxCol = i;
                }
            }
            break;

        default:
            break;
        }
        return result;
    }

}
