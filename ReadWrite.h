#ifndef READWRITE_H_INCLUDED
#define READWRITE_H_INCLUDED

#include <iostream>
#include <vector>
using namespace std;

class ReadWrite {

public:
    //void readCSV();
    void writeCSV(vector<vector<double > > & vec, const char* filename);


private:

};

#endif // READWRITE_H_INCLUDED
