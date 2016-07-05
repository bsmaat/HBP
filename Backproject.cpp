#include "Backproject.h"
#include <iostream>
#include "Path.h"
#include "Siddon.h"
#include <math.h>
using namespace std;

Backproject::Backproject() {
}

Backproject::Backproject(vector<vector<double> > &data) {
    //std::cout << data[0][0] << std::endl;
    dataTmp = &data;
    (*dataTmp)[0][0] = 1000;
    //std::cout << (*dataTmp)[0][0] << std::endl;

}

void Backproject::radon(vector<double> & data, vector<vector<double> > & bpAngle, vector<vector<double> > & chordLength, Siddon & s) {

    // first compute the path
    Path * p = new Path(data);
    double phi = M_PI/4;
    //p->rotate(phi);

    vector<vector<double> > & path = p->getPath();
    vector<vector<double> > indices;
    indices = s.getIntersect(path[0], path[1]);
    for (vector<vector<double> >::iterator i = ++path.begin(); i != --path.end(); ++i) {
        yVal = (s.getIntersect(*i, *(i+1)));
        // check if it's the same index and if it is add them together
        if ((indices.back()[0] == yVal[0][0]) && (indices.back()[1] == yVal[0][1])) {
            indices.back()[2] = indices.back()[2] + yVal[0][2];
            indices.insert(indices.end(), yVal.begin()+1, yVal.end());
        } else
            indices.insert(indices.end(), yVal.begin(), yVal.end());
    }

    for (vector<vector<double> > ::const_iterator i = indices.begin(); i != indices.end(); ++i) {
        bpAngle[(*i)[0]-1][(*i)[1]-1] = (*i)[2] * data[4] + bpAngle[(*i)[0]-1][(*i)[1]-1];
        chordLength[(*i)[0]-1][(*i)[1]-1] = chordLength[(*i)[0]-1][(*i)[1]-1] + (*i)[2];
    }

    delete p;

}

Backproject::Backproject(vector<double> & data) {




}
