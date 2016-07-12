#include "Path.h"
#include <iostream>
#include <algorithm>
#include <math.h>

Path::Path() {
    //CubicSpline(1, 2, 3, 4, 5, 6);
}

Path::Path(double & entryPos, double & entryAngle, double & exitPos, double & exitAngle) {
    this->entryPos = entryPos;
    this->entryAngle = entryAngle;
    this->exitPos = exitPos;
    this->exitAngle = exitAngle;
}

Path::Path(vector<double> & row, double & angle) {
    this->entryPos = row[0];
    this->entryAngle = row[1];
    this->exitPos = row[2];
    this->exitAngle = row[3];

    // does this path intersect with ellipse?
    vector<double> A, B;
    //double angle = M_PI/4;


    A.push_back(-230);
    A.push_back(row[0]);
    B.push_back(+230);
    B.push_back(row[2]);

/*
    A.push_back(-230 * cos(angle) - row[0] * sin(angle));
    A.push_back(-230 * sin(angle) + row[0] * cos(angle));

    B.push_back(230 * cos(angle) - row[2] * sin(angle));
    B.push_back(230 * sin(angle) + row[2] * cos(angle));
*/


    //cout << A[0] << ", " << A[1] << " : " << B[0] << ", " << B[1] << endl;

    vector<vector<double> > ellipseCoords =  ellipseIntersect(A, row[1], B, row[3], angle); // will be empty if it doesn't intersect
    //vector<vector<double> > ellipseCoords =  ellipseIntersect(A, entryAngleNew, B, exitAngleNew, angTmp); // will be empty if it doesn't intersect


    vector<double> point;
    point.push_back(-230 * cos(angle) - row[0] * sin(angle));
    point.push_back(-230 * sin(angle) + row[0] * cos(angle));
    //point.push_back(-230);
    //point.push_back(row[0]);
    path.push_back(point);

    // the -angle is kind of a guess
    //entrydydx = sin(row[1]-angle)/cos(row[1]-angle);
    //exitdydx = sin(row[3]-angle)/cos(row[3]-angle);

    entrydydx = sin(row[1])/cos(row[1]);
    exitdydx = sin(row[3])/cos(row[3]);

    // if path would intersect with ellipse, then find CubicSpline stuff
    if (!ellipseCoords.empty()) {
        Vector4d spline = CubicSpline(ellipseCoords[0][0], ellipseCoords[1][0], ellipseCoords[0][1], ellipseCoords[1][1], entrydydx, exitdydx);
        vector<double> u1;

        int numOfSteps = 5;
        //double distanceInEllipse = sqrt(pow(ellipseCoords[1][0] - ellipseCoords[0][0],2) + pow(ellipseCoords[1][1] - ellipseCoords[0][1],2));
        double stepSize = (ellipseCoords[1][0]-ellipseCoords[0][0])/numOfSteps;
        //double stepSize = distanceInEllipse/numOfSteps;
        for (int i = 0; i <= numOfSteps; i++) {
            u1.push_back(ellipseCoords[0][0] + stepSize*i);

        }

        vector<vector<double> > cubicPath = CubicSplinePath(u1, spline, angle);

        path.insert(path.end(), cubicPath.begin(), cubicPath.end() );
    }

    point.clear();
    point.push_back(+230 * cos(angle) - row[2] * sin(angle));
    point.push_back(+230 * sin(angle) + row[2] * cos(angle));
    path.push_back(point);


}

Vector4d Path::CubicSpline() {
    Matrix4d A;
    Vector4d y;
    double x0 = -230;
    double x1 = 230;
    double y0 = entryPos;
    double y1 = exitPos;
    double dy0 = sin(entryAngle)/cos(entryAngle);
    double dy1 = sin(exitAngle)/cos(exitAngle);


    for (int i=0; i<3; i++) {
        A(0, i) = pow(x0, 3-i);
        A(1, i) = pow(x1, 3-i);
        A(2, i) = (3-i)*pow(x0,2-i);
        A(3, i) = (3-i)*pow(x1,2-i);
    }
    A(0, 3) = 1;
    A(1, 3) = 1;
    A(2, 3) = 0;
    A(3, 3) = 0;

    y(0) = y0;
    y(1) = y1;
    y(2) = dy0;
    y(3) = dy1;


    return y;
}

// calculate cubic spline coefficients
Vector4d Path::CubicSpline(double & x0, double & x1, double & y0, double & y1, double & dy0, double & dy1) {
    Matrix4d A;
    Vector4d x, y;

    for (int i=0; i<3; i++) {
        A(0, i) = pow(x0, 3-i);
        A(1, i) = pow(x1, 3-i);
        A(2, i) = (3-i)*pow(x0,2-i);
        A(3, i) = (3-i)*pow(x1,2-i);
    }
    A(0, 3) = 1;
    A(1, 3) = 1;
    A(2, 3) = 0;
    A(3, 3) = 0;


    y(0) = y0;
    y(1) = y1;
    y(2) = dy0;
    y(3) = dy1;

    x = A.inverse() * y;

    //vector<double> coefficients(y.data(), y.data() + y.size());
    //return coefficients;

    return x;
}

vector<vector<double> > Path::CubicSplinePath(vector<double> & x, Vector4d & coeff, double & phi) {
    vector<vector<double> > path;
    vector<double> coord;
    Vector4d argument;
    double y, xTmp, yTmp;
    for(unsigned int i = 0; i<x.size(); i++) {
        argument = Vector4d(pow(x[i],3), pow(x[i], 2), x[i], 1);
        yTmp = coeff.dot(argument);
        xTmp = x[i];
        x[i] = xTmp * cos(phi) - yTmp * sin(phi);
        y = xTmp * sin(phi) + yTmp * cos(phi);
        coord.push_back(x[i]);
        coord.push_back(y);
        path.push_back(coord);
        coord.clear();
    }

    return path;
}

vector<vector<double> > Path::mlp(vector<double> & x, Vector2d & y0, Vector2d & y1, vector<double> & coeff, double & phi) {
    vector<vector<double> > path;
    vector<double> coord;
    vector<double> a(6); // coefficient of polynomial for 1/beta^2 p^2
    double X0 = 361; // radiation length, mm
    double E0 = 1.36; // MeV/mm


    auto fst = [&X0, &E0, &a](double & u1, double & u0) -> double {
        double out = (pow(E0,2) * pow(1 + 0.038*log((u1 - u0)/X0),2)/X0) *
        ( (1/1.0) * (pow(u1,2)*a[0]   ) * (u1 - u0)
        + (1/2.0) * (pow(u1,2)*a[1] - 2*u1*a[0]) * (pow(u1,2) - pow(u0,2))
        + (1/3.0) * (pow(u1,2)*a[2] - 2*u1*a[1] + a[0]) * (pow(u1,3) - pow(u0,3))
        + (1/4.0) * (pow(u1,2)*a[3] - 2*u1*a[2] + a[1]) * (pow(u1,4) - pow(u0,4))
        + (1/5.0) * (pow(u1,2)*a[4] - 2*u1*a[3] + a[2]) * (pow(u1,5) - pow(u0,5))
        + (1/6.0) * (pow(u1,2)*a[5] - 2*u1*a[4] + a[3]) * (pow(u1,6) - pow(u0,6))
        + (1/7.0) * ( - 2 * u1 * a[5] + a[4]) * (pow(u1,7) - pow(u0,7))
        + (1/8.0) * (+ a[5]) * (pow(u1,8) - pow(u0,8)) ) ;
        return out;
    };

    auto fsth = [&X0, &E0, &a](double & u1, double & u0) -> double {
        double out = (pow(E0,2) * pow(1 + 0.038*log((u1-u0)/X0),2)/X0 ) *
        ( a[0] * (u1-u0) + a[1] * (pow(u1,2) - pow(u0,2))/2.0 + a[2] * (pow(u1,3) - pow(u0,3))/3.0
        + a[3] * (pow(u1,4) - pow(u0,4))/4.0 + a[4] * (pow(u1,5) - pow(u0,5))/5.0
        + a[5] * (pow(u1,6) - pow(u0,6))/6.0 );
        return out;
    };

    auto fstth = [&X0, &E0, &a](double & u1, double & u0) -> double {
        double out = (pow(E0,2) * pow(1 + 0.038*log((u1 - u0)/X0),2)/X0) *
        ( (1/1.0) * u1 * a[0] * (u1 - u0)
        + (1/2.0) * (u1 * a[1] - a[0]) * (pow(u1,2) - pow(u0,2))
        + (1/3.0) * (u1 * a[2] - a[1]) * (pow(u1,3) - pow(u0,3))
        + (1/4.0) * (u1 * a[3] - a[2]) * (pow(u1,4) - pow(u0,4))
        + (1/5.0) * (u1 * a[4] - a[3]) * (pow(u1,5) - pow(u0,5))
        + (1/6.0) * (u1 * a[5] - a[4]) * (pow(u1,6) - pow(u0,6))
        + (1/7.0) * (- a[5]) * (pow(u1,7) - pow(u0,7)) );
        return out;
    };

    vector<Vector2d> ymlp(x.size());
    ymlp[0] = y0;
    ymlp.back() = y1;
    double st1, sth1, stth1, st2, sth2, stth2;
    double yTmp;
    double xTmp;
    Matrix2d SIG1, SIG2, R0, R1;
    for(unsigned int i=1; i < x.size()-1; i++) {
        st1 = fst(x[i], x[0]);
        sth1 = fsth(x[i], x[0]);
        stth1 = fstth(x[i], x[0]);

        st2 = fst(x.back(), x[i]);
        sth2 = fsth(x.back(), x[i]);
        stth2 = fstth(x.back(), x[i]);

        SIG1 << st1, stth1,
                stth1, sth1;

        SIG2 << st2, stth2,
                stth2, sth2;

        R0 << 1, x[i]-x[0],
              0, 1;

        R1 << 1, x.back() - x[i],
              0, 1;

        Matrix2d firstTerm;
        Vector2d secondTerm, thirdTerm;
        firstTerm = SIG1.inverse() + R1.transpose() * SIG2.inverse() * R1;
        secondTerm = (SIG1.inverse() * R0) * y0;
        thirdTerm = R1.transpose() * (SIG2.inverse() * y1);

        ymlp[i] = firstTerm.inverse() * (secondTerm + thirdTerm);

        xTmp = x[i];
        yTmp = ymlp[i](0);
        x[i] = xTmp * cos(phi) - yTmp * sin(phi);
        ymlp[i](0) = xTmp * sin(phi) + yTmp * cos(phi);
        coord.push_back(x[i]);
        coord.push_back(ymlp[i](0));
        path.push_back(coord);
        coord.clear();
    }

    return path;


}

// check where the intersection with an ellipse occurs
// return a vector with intersection1, intersection2, and some kind of bool
vector<vector<double> > Path::ellipseIntersect(vector<double> & A, double & entryAngle, vector<double> & B, double & exitAngle, double & phi) {
    vector<vector<double> > intersect;
    vector<double> p0;
    vector<double> p1;
    double a = 70, b = 80;
    vector<double> dirVecA, dirVecB;
    dirVecA.push_back(cos(entryAngle));
    dirVecA.push_back(sin(entryAngle));
    dirVecB.push_back(cos(exitAngle));
    dirVecB.push_back(sin(exitAngle));
    double cosphia2, cosphib2, sinphia2, sinphib2;
    cosphia2 = pow(cos(phi)/a,2);
    cosphib2 = pow(cos(phi)/b,2);
    sinphia2 = pow(sin(phi)/a,2);
    sinphib2 = pow(sin(phi)/b,2);

    double a1, a2, a3, x1, y1, x2, y2;

    a1 = pow(dirVecA[0],2) * ( cosphia2 + sinphib2 ) + pow(dirVecA[1],2) * (sinphia2 + cosphib2)
        - 2 * cos(phi) * sin(phi) * (1.0/pow(a, 2) - 1.0/pow(b,2)) * dirVecA[0] * dirVecA[1];

    a2 = 2 * A[0] * dirVecA[0] * (cosphia2 + sinphib2) + 2 * A[1] * dirVecA[1] * (sinphia2 + cosphib2)
        - 2 * cos(phi) * sin(phi) * (1.0/pow(a,2) - 1.0/pow(b,2)) * (A[0] * dirVecA[1] + dirVecA[0] * A[1]);

    a3 = pow(A[0],2) * ( cosphia2 + sinphib2) + pow(A[1],2) * (sinphia2 + cosphib2) - 2 * cos(phi) * sin(phi)
        * (1.0/pow(a,2) - 1.0/pow(b,2)) * A[0] * A[1] - 1;

    vector<complex<double> > rts;
    rts = quadraticSolver(a1, a2, a3);
    //check that all roots are real
    bool allReal;
    allReal = std::all_of(rts.begin(), rts.end(), [](complex<double> i) { return (i.imag() == 0);} );

    bool allMoreThanZero;
    allMoreThanZero = std::all_of(rts.begin(), rts.end(), [](complex<double> i) { return (i.real() > 0);});

    double d;
    if (allReal && allMoreThanZero) {
        //std::max_element(rts.begin(), rts.end(), real_less); <- don't know why this had a problem!
        d = (*std::min_element(rts.begin(), rts.end(), [](complex<double> const & lhs, complex<double> const & rhs) {return lhs.real() < rhs.real();})).real();
        x1 = A[0] + dirVecA[0] * d;
        y1 = A[1] + dirVecA[1] * d;
        p0.push_back(x1);
        p0.push_back(y1);
    } else {
        return intersect;
    }

    // computations for output line

    a1 = pow(dirVecB[0],2) * ( cosphia2 + sinphib2 ) + pow(dirVecB[1],2) * (sinphia2 + cosphib2)
        - 2 * cos(phi) * sin(phi) * (1.0/pow(a, 2) - 1.0/pow(b,2)) * dirVecB[0] * dirVecB[1];

    a2 = 2 * B[0] * dirVecB[0] * (cosphia2 + sinphib2) + 2 * B[1] * dirVecB[1] * (sinphia2 + cosphib2)
        - 2 * cos(phi) * sin(phi) * (1.0/pow(a,2) - 1.0/pow(b,2)) * (B[0] * dirVecB[1] + dirVecB[0] * B[1]);

    a3 = pow(B[0],2) * ( cosphia2 + sinphib2) + pow(B[1],2) * (sinphia2 + cosphib2) - 2 * cos(phi) * sin(phi)
        * (1.0/pow(a,2) - 1.0/pow(b,2)) * B[0] * B[1] - 1;

    rts.clear();
    rts = quadraticSolver(a1, a2, a3);

    //check that all roots are real
    allReal = std::all_of(rts.begin(), rts.end(), [](complex<double> i) { return (i.imag() == 0);} );

    double allLessThanZero = std::all_of(rts.begin(), rts.end(), [](complex<double> i) { return (i.real() < 0);});

    if (allReal && allLessThanZero) {
        //std::max_element(rts.begin(), rts.end(), real_less); <- don't know why this had a problem!
        d = (*std::max_element(rts.begin(), rts.end(), [](complex<double> const & lhs, complex<double> const & rhs) {return lhs.real() < rhs.real();})).real();
        x2 = B[0] + dirVecB[0] * d;
        y2 = B[1] + dirVecB[1] * d;
        p1.push_back(x2);
        p1.push_back(y2);
    } else {
        return intersect;
    }

    if (x2 <= x1) {
        return intersect;
    }
    intersect.push_back(p0);
    intersect.push_back(p1);
    return intersect;

}

// could be practice to write as a template?
vector<complex<double> > Path::quadraticSolver(double & a, double & b, double & c) {

    vector<complex<double> > roots;
    complex<double> sol1(0,0), sol2(0,0);
    double discriminant = pow(b,2) - 4 * a * c;
    if (discriminant >= 0) {
        sol1.real( (-b + sqrt(discriminant))/(2*a) );
        sol2.real( (-b - sqrt(discriminant))/(2*a) );
    } else {
        sol1.real(-b/(2*a));
        sol2.real(-b/(2*a));
        sol1.imag( (sqrt(-1 * discriminant))/(2*a) );
        sol2.imag( -1 * (sqrt(-1 * discriminant))/(2*a) );
    }

    roots.push_back(sol1);
    roots.push_back(sol2);
    return roots;

}

// unneeded at the moment...
bool Path::real_less(std::complex<double> const & lhs, std::complex<double> const & rhs) {
    return lhs.real() < rhs.real();
}

void Path::testFn(double & a) {

}

void Path::rotate(double & phi) {
    double xTmp, yTmp;
    for(vector<vector<double> >::iterator i = path.begin(); i != path.end(); ++i) {
        xTmp = (*i)[0];
        yTmp = (*i)[1];
        //(*i)[0] = xTmp * cos(phi) - yTmp * sin(phi);
        //(*i)[1] = xTmp * sin(phi) + yTmp * cos(phi);

        (*i)[0] = xTmp * 0.7071 - yTmp * 0.7071;
        (*i)[1] = xTmp * 0.7071 + yTmp * 0.7071;
    }
}
