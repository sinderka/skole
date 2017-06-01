#ifndef TOOLS_H
#define TOOLS_H

#include <string>

namespace Tools {
    // allocate
    double* initArray(long length);
    double* randomArray(long length, int seed);

    // print
    void print(auto str);
    void print(std::string str);
    void print(char * str);
    void print(double* mat, long n_row, long n_col);

}

#endif
