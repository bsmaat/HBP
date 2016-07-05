#ifndef GRID_H_INCLUDED
#define GRID_H_INCLUDED

class Grid {
public:
    Grid(int N, double d);

private:
	int N;  // number of planes
	double d; // spacing between planes
	vector<double> Xp; // position of x planes
	vector<double> Yp; // position of y planes
	double Xp0;
	double Yp0;
}



#endif // GRID_H_INCLUDED
