#include "Siddon.h"
#include <iostream>
#include <algorithm>
#include <math.h>

Siddon::Siddon() {
	d = 1;
	N = 5;
}

Siddon::Siddon(int N, double d) {
	this->N = N;
	this->d = d;
	this->Xp0 = -230;
	this->Yp0 = -230;
    for (int i=0; i<N; i++) {

        Xp.push_back(Xp0 + i*d);
        Yp.push_back(Yp0 + i*d);
    }

}

vector<vector<double> > Siddon::getIntersect(vector<double> & P0, vector<double> & P1) {
    double& X1 = P0[0];
    double& X2 = P1[0];
    double& Y1 = P0[1];
    double& Y2 = P1[1];

    vector<double> ax, ay;
    if ( (X2-X1) != 0 ) {
        //std::cout << X1 << " : " << X2 << std::endl;
        ax.resize(Xp.size());
        std::transform(Xp.begin(), Xp.end(), ax.begin(), bind2nd(std::plus<double>(), -X1));
        std::transform(ax.begin(), ax.end(), ax.begin(), bind1st(std::multiplies<double>(), 1.0/(X2-X1)));
    }

    if ( (Y2-Y1) != 0 ) {
        //std::cout << Y1 << " : " << Y2 << std::endl;

        ay.resize(Yp.size());
        std::transform(Yp.begin(), Yp.end(), ay.begin(), bind2nd(std::plus<double>(), -Y1));
        std::transform(ay.begin(), ay.end(), ay.begin(), bind1st(std::multiplies<double>(), 1.0/(Y2-Y1)));
    }
/*
    if (ax.empty())
        std::cout << "ax empty" << std::endl;
    if (ay.empty())
        std::cout << "ay empty" << std::endl;
*/


    // calculate amin and amax
    double amin, amax;
    vector<double> compare;
    if (ax.empty() == 0 && ay.empty() == 0) {

        compare.push_back(0);
        if (ax[0] < ax.back())
            compare.push_back(ax[0]);
        else
            compare.push_back(ax.back());
        if (ay[0] < ay.back())
            compare.push_back(ay[0]);
        else
            compare.push_back(ay.back());

        amin = *std::max_element(compare.begin(), compare.end());

        compare.clear();
        compare.push_back(1);
        if (ax[0] > ax.back())
            compare.push_back(ax[0]);
        else
            compare.push_back(ax.back());
        if (ay[0] > ay.back())
            compare.push_back(ay[0]);
        else
            compare.push_back(ay.back());

        amax = *std::min_element(compare.begin(), compare.end());

    } else if (ax.empty()) {

        compare.push_back(0);
        if (ay[0] < ay.back()) {
            compare.push_back(ay[0]);
        }
        else {
            compare.push_back(ay.back());
        }

        amin = *std::max_element(compare.begin(), compare.end());

        compare.clear();
        compare.push_back(1);
        if (ay[0] > ay.back())
            compare.push_back(ay[0]);
        else
            compare.push_back(ay.back());

        amax = *std::min_element(compare.begin(), compare.end());

    } else if (ay.empty()) {
        compare.push_back(0);
        if (ax[0] < ax.back())
            compare.push_back(ax[0]);
        else
            compare.push_back(ax.back());

        amin = *std::max_element(compare.begin(), compare.end());

        compare.clear();
        compare.push_back(1);
        if (ax[0] > ax.back())
            compare.push_back(ax[0]);
        else
            compare.push_back(ax.back());

        amax = *std::min_element(compare.begin(), compare.end());
    }


    // calculate imin, imax, ax, jmin, jmax, ay
    double imin, imax;
    if (!ax.empty()) {
        //sort(ax.begin(), ax.end());
        ax.erase(unique(ax.begin(), ax.end()), ax.end());

        if ((X2 - X1) >= 0) {
            imin = N - (Xp.back() - amin * (X2-X1) - X1) / d;
            imax = 1 + (X1 + amax * (X2 - X1) - Xp[0])/d;
            imin = ceil(imin);
            imax = floor(imax);
            //std::cout << imin << ", " << imax << ", " << amin << ", " << amax << std::endl;
            if (imin<=imax)
                ax = std::vector<double>(ax.begin() + (int)imin-1, ax.begin()+(int)imax);
            else
                ax.clear();

        } else if ((X2 - X1) < 0) {
            imin = N - (Xp.back() - amax * (X2 - X1) - X1) / d;
            imax = 1 + (X1 + amin * (X2 - X1) - Xp[0]) / d;
            imin = ceil(imin);
            imax = floor(imax);
            if (imin<=imax)
                ax = std::vector<double>(ax.begin() + (int)imin-1, ax.begin() + (int)imax);
            else
                ax.clear();
            //std::reverse(ax.begin(), ax.end());
        }
    }

    double jmin, jmax;
    if(!ay.empty()) {
        //sort(ay.begin(), ay.end());
        ay.erase(unique(ay.begin(), ay.end()), ay.end());
        if ((Y2 - Y1) >= 0) {
            jmin = N - (Yp.back() - amin * (Y2-Y1) - Y1) / d;
            jmax = 1 + (Y1 + amax * (Y2 - Y1) - Yp[0])/d;
            jmin = ceil(jmin);
            jmax = floor(jmax);
            if (jmin<=jmax)
                ay = std::vector<double>(ay.begin() + (int)jmin-1, ay.begin()+(int)jmax);
            else
                ay.clear();
        } else if ((Y2 - Y1) < 0) {
            jmin = N - (Yp.back() - amax * (Y2 - Y1) - Y1) / d;
            jmax = 1 + (Y1 + amin * (Y2 - Y1) - Yp[0]) / d;
            jmin = ceil(jmin);
            jmax = floor(jmax);
            if (jmin<=jmax)
                ay = std::vector<double>(ay.begin() + (int)jmin-1, ay.begin() + (int)jmax);
            else
                ay.clear();
        }
    }


    vector<double> a;
    a.push_back(amin);
    a.push_back(amax);
    a.insert(a.end(), ax.begin(), ax.end());
    a.insert(a.end(), ay.begin(), ay.end());
    sort(a.begin(), a.end());
    a.erase(std::unique(a.begin(), a.end()), a.end());

    /* debugging
    std::cout << "imin: " << imin << ", imax: " << imax <<
        ", amin: " << amin << ", amax: " << amax << std::endl;
    std::cout << "jmin: " << jmin << ", jmax: " << jmax << std::endl;

    std::cout << "ax" << std::endl;
    writeVector(ax);
    std::cout << "ay" << std::endl;
    writeVector(ay);
    std::cout << "a" << std::endl;
    writeVector(a);
    */

    double amid;
    int i, j;
    double pathLength;
    vector<vector<double> > indices;
    std::vector<double> v;
    v.resize(3);
    /*
    cout << "ax: ";
    writeVector(ax);
    cout << "ay: ";
    writeVector(ay);
    cout << "a: ";
    writeVector(a);
    */
    for(int m = 0; m<a.size()-1; m++) {
        amid = (a[m+1] + a[m])/2.0;
        // i = x position, j = y position
        i = floor( 1 + (X1 + amid*(X2-X1) - Xp[0])/d );
        j = floor( 1 + (Y1 + amid*(Y2-Y1) - Yp[0])/d );
        if ((i>0 && j > 0) && (i <= N-1 && j <= N-1)) {
            pathLength = (a[m+1] - a[m])* sqrt( pow((X2-X1),2) + pow((Y2-Y1),2) );
            v[0] = j;
            v[1] = i;
            v[2] = pathLength;
            indices.push_back(v);
        }
    }

/*
    for (int m = 0; m<indices.size(); m++) {
        writeVector(indices[m]);
    }
*/

    return indices;

}

void Siddon::writeVector(vector<double> vec) {
    for (vector<double>::const_iterator i = vec.begin(); i != vec.end(); ++i)
        std::cout << *i << ' ';
        std::cout << std::endl;
}

