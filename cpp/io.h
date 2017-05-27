#ifndef IO_H
#define OI_H

#include <string>

// output functions
void print(auto str);
void print(std::string str);
void print(char const* str);
void printArray(double* mat, long n, long m);

// input functions
std::string getStringFromUser();
double getDoubleFromUser();
long getLongFromUser();
int getIntFromUser();

#endif
