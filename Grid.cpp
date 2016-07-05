#include "Grid.h"

Grid::Grid(int N, double d) {
	this->N = N;
	this->d = d;
	this->Xp0 = -2;
	this->Yp0 = -2;
    for (int i=0; i<N; i++) {
        std::cout << "yo" << std::endl;
        Xp.push_back(Xp0 + i*d);
        Yp.push_back(Yp0 + i*d);
    }
}
