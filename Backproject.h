#ifndef BACKPROJECT_H_INCLUDED
#define BACKPROJECT_H_INCLUDED

#include <vector>
#include "Siddon.h"
using namespace std;

class Backproject {

public:
    Backproject();
    Backproject(vector<vector<double> > &data);
    Backproject(vector<double> & data);
    void radon(vector<double> & data, vector<vector<double> > & bpAngle, vector<vector<double> > & chordLength, Siddon & s);
    void rotate(vector<double>&);
    vector<vector<double> > getYValues() { return yVal; };
    void setAngle(double & angle);

private:
    vector<vector<double> >* dataTmp;
    vector<vector<double> > yVal;
    double angle;
};

#endif // BACKPROJECT_H_INCLUDED
