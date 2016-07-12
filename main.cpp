#include <iostream>

#include <fstream>

#include <stdexcept>
#include "Siddon.h"
#include "Path.h"
#include "CSVRow.h"
#include "Backproject.h"
#include "ReadWrite.h"
#include <thread>
#include <mutex>
#include <numeric>
#include <iomanip>

std::mutex m;
std::mutex mut;
int counter = 0;
int N = 361; // numer of planes...
double d = 0.5; // spacing between planes

vector<MatrixXd> bpAngle(180, MatrixXd::Zero(N-1, N-1));
vector<MatrixXd> chordLength(180, MatrixXd::Zero(N-1, N-1));

//MatrixXd bpAngleBPF;
//MatrixXd chordLengthBPF;


template <class T1, class T2>
T1 lexical_cast(const T2& t2) {
    std::stringstream s;
    s << t2;
    T1 t1;
    if (s >> t1 && s.eof()) {
        return t1;
    } else {
        throw std::runtime_error("bad conversion");
        return T1();
        exit(1);
    }

}

std::istream& operator>>(std::istream& str, CSVRow& data) {
    data.readNextRow(str);
    return str;
}

//find mean and standard deviation of vector
vector<double> stats(vector<double> & v) {
    double sum = std::accumulate(v.begin(), v.end(), 0.0);
    double mean = sum / v.size();


    std::vector<double> diff(v.size());
    std::transform(v.begin(), v.end(), diff.begin(), [mean](double x) { return (x - mean);});
    double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
    double stdev = std::sqrt(sq_sum/v.size());


    vector<double> vec(2);
    vec[0] = mean;
    vec[1] = stdev;
    return vec;
}

// perform sigma filtering on proton data
void filterData(vector<vector<double> > & data) {

    vector<double> relAngle(data.size()), relPos(data.size());

    for (int i =0; i < data.size(); i++) {
        relAngle[i] = data[i][3] - data[i][1];
        relPos[i] = data[i][2] - data[i][0];
    }

    vector<double> angleStats = stats(relAngle);
    vector<double> posStats = stats(relPos);

    // filter elements
    for (int i = relAngle.size() - 1; i >= 0; i--) {
        if ((abs(relAngle[i] - angleStats[0]) > 3 * angleStats[1]) || (abs(relPos[i] - posStats[0]) > 3 * posStats[1])) {
            data.erase(data.begin() + i, data.begin() + i + 1);
        }
    }

}

void backprojectAngle(int id) {

    std::ostringstream filename;
    filename << "/home/billy/geant/EdgeModelx2/proj/proj" << id << ".csv";
    std::string str = filename.str();
    //std::ifstream file("/home/billy/geant/EdgeModelx2/proj/proj45.csv");
	std::ifstream file(str);
	CSVRow row;
    vector<vector<double> > detected;
    vector<double> tmp;
    tmp.resize(5);
    double entryAngle, exitAngle;
	while(file >> row)
	{
        //vec.push_back(row);
        entryAngle = atan(lexical_cast<double>(row[4])/lexical_cast<double>(row[6]));
        exitAngle = atan(lexical_cast<double>(row[12])/lexical_cast<double>(row[14]));

        tmp[0] = lexical_cast<double>(row[1])*1000; // lateral entry position (entryPos) (mm)
        tmp[1] = entryAngle;
        tmp[2] = lexical_cast<double>(row[9])*1000; // lateral exit position (exitPos) (mm)
        tmp[3] = exitAngle;
        tmp[4] = lexical_cast<double>(row[16]); // projection
        detected.push_back(tmp);

	}

	//filter data
	cout << detected.size() << endl;
	filterData(detected);
    cout << detected.size() << endl;

    //vector<vector<double> > bpAngle(460, vector<double>(460, 0));
    //vector<vector<double> > chordLength(460, vector<double>(460,0));

    vector<MatrixXd> bpAngle(180, MatrixXd::Zero(460, 460));
    vector<MatrixXd> chordLength(180, MatrixXd::Zero(460, 460));

    ReadWrite* r = new ReadWrite();
    Siddon* s = new Siddon(461, 1);

    double angle = id*M_PI/180.0;
    for (int i = 0; i < detected.size(); i++) {
        Backproject* bp = new Backproject();
        bp->setAngle(angle);
        bp->radon(detected[i], bpAngle, chordLength, *s, mut);
        delete bp;
    }

/*
    for (int i = 0; i < bpAngle.size(); i++) {
        for (int j = 0; j < bpAngle.size(); j++) {
            if (chordLength[i][j] != 0) {
                bpAngle(i,j) = bpAngle(i,j) / chordLength(i,j);
                //bpAngle[i][j] = bpAngle[i][j] / chordLength[i][j];
                //bpAngle[i][j]= chordLength[i][j];
            }
        }
    }

*/

    for (int k = 0; k < bpAngle.size(); k++) {
    for (int i = 0; i < bpAngle[k].rows(); i++) {
        for (int j = 0; j < bpAngle[k].rows(); j++) {
            if (chordLength[k](i,j) != 0) {
                bpAngle[k](i,j) = bpAngle[k](i,j) / chordLength[k](i,j);
                //bpAngle[i][j] = bpAngle[i][j] / chordLength[i][j];
                //bpAngle[i][j]= chordLength[i][j];
            }

            //fill the air parts
            if (bpAngle[k](i,j) == 0.0) {
                bpAngle[k](i,j) = 0.4930;
            }
        }
    }
    }

    //cout << "Val: " << setprecision(15) << bpAngle[213][162] << endl;
    // fill the air parts...
    /*
    for (vector<vector<double> >::iterator i = bpAngle.begin(); i != bpAngle.end(); ++i) {
        std::replace(i->begin(), i->end(), 0.0, 0.4930);
    }
    */


    filename.str("");
    filename << "output" << id << ".csv";
    str = filename.str();
    const char* fileStr = str.c_str();
    //r->writeCSV(bpAngle, fileStr);
    m.lock();
        counter++;
        std::cout << counter << std::endl;
        std::cout << fileStr << std::endl;
    m.unlock();
}

void backprojectAngle2(int id) {

    std::ostringstream filename;
    filename << "/home/billy/geant/EdgeModelx2/proj/proj" << id << ".csv";
    std::string str = filename.str();
    //std::ifstream file("/home/billy/geant/EdgeModelx2/proj/proj45.csv");
	std::ifstream file(str);
	CSVRow row;
    vector<vector<double> > detected;
    vector<double> tmp;
    tmp.resize(5);
    double entryAngle, exitAngle;
	while(file >> row)
	{
        //vec.push_back(row);
        entryAngle = atan(lexical_cast<double>(row[4])/lexical_cast<double>(row[6]));
        exitAngle = atan(lexical_cast<double>(row[12])/lexical_cast<double>(row[14]));

        tmp[0] = lexical_cast<double>(row[1])*1000; // lateral entry position (entryPos) (mm)
        tmp[1] = entryAngle;
        tmp[2] = lexical_cast<double>(row[9])*1000; // lateral exit position (exitPos) (mm)
        tmp[3] = exitAngle;
        tmp[4] = lexical_cast<double>(row[16]); // projection
        detected.push_back(tmp);

	}

	//filter data
	cout << detected.size() << endl;
	filterData(detected);
    cout << detected.size() << endl;


    ReadWrite* r = new ReadWrite();
    Siddon* s = new Siddon(461, 1);

    double angle = id*M_PI/180.0;
    for (int i = 0; i < 1; i++) {//detected.size(); i++) {
        Backproject* bp = new Backproject();
        bp->setAngle(angle);
        bp->radon(detected[i], bpAngle, chordLength, *s, mut);
        delete bp;
    }



    for (int k = 0; k < bpAngle.size(); k++) {
        for (int i = 0; i < bpAngle[k].rows(); i++) {
            for (int j = 0; j < bpAngle[k].rows(); j++) {
                if (chordLength[k](i,j) != 0) {
                    bpAngle[k](i,j) = bpAngle[k](i,j) / chordLength[k](i,j);
                }

                //fill the air parts
                if (bpAngle[k](i,j) == 0.0) {
                    bpAngle[k](i,j) = 0.4930;
                }
            }
        }
    }



    filename.str("");
    filename << "output" << id << ".csv";
    str = filename.str();
    const char* fileStr = str.c_str();
    //r->writeCSV(bpAngle, fileStr);
    m.lock();
        counter++;
        std::cout << counter << std::endl;
        std::cout << fileStr << std::endl;
    m.unlock();
}

void backprojectAngleBPF(int id) {

    std::ostringstream filename;
    filename << "/home/billy/geant/SlitModelx2/proj/proj" << id << ".csv";
    std::string str = filename.str();
	std::ifstream file(str);
	CSVRow row;
    vector<vector<double> > detected;
    vector<double> tmp;
    tmp.resize(5);
    double entryAngle, exitAngle;
	while(file >> row)
	{
        //vec.push_back(row);
        entryAngle = atan(lexical_cast<double>(row[4])/lexical_cast<double>(row[6]));
        exitAngle = atan(lexical_cast<double>(row[12])/lexical_cast<double>(row[14]));

        tmp[0] = lexical_cast<double>(row[1])*1000; // lateral entry position (entryPos) (mm)
        tmp[1] = entryAngle;
        tmp[2] = lexical_cast<double>(row[9])*1000; // lateral exit position (exitPos) (mm)
        tmp[3] = exitAngle;
        tmp[4] = lexical_cast<double>(row[16]); // projection
        detected.push_back(tmp);

	}

	//filter data
	cout << "Before filtering: " << detected.size() << ", ";
	filterData(detected);
    cout << "After filtering: " << detected.size() << endl;


    ReadWrite* r = new ReadWrite();
    Siddon* s = new Siddon(461, 1);

    double angle = id*M_PI/180.0;
    MatrixXd bpAngleBPF;
    MatrixXd chordLengthBPF;
    bpAngleBPF = MatrixXd::Zero(460,460);
    chordLengthBPF = MatrixXd::Zero(460,460);

    for (int i = 0; i < detected.size(); i++) {
        Backproject* bp = new Backproject();
        bp->setAngle(angle);
        bp->radonBPF(detected[i], bpAngleBPF, chordLengthBPF, *s);

        delete bp;
    }


    for (int i = 0; i < bpAngleBPF.rows(); i++) {
        for (int j = 0; j < bpAngleBPF.rows(); j++) {
            if (chordLengthBPF(i,j) != 0) {
                bpAngleBPF(i,j) = bpAngleBPF(i,j) / chordLengthBPF(i,j);
            }
            //fill the air parts
            if (bpAngleBPF(i,j) == 0.0) {
                bpAngleBPF(i,j) = 0.4930;
            }
        }
    }



    filename.str("");
    filename << "BPFFiles/bpt" << id << ".csv";
    str = filename.str();
    const char* fileStr = str.c_str();
    r->writeCSV(bpAngleBPF, fileStr);
    m.lock();
        counter++;
        std::cout << "Counter: " << counter << std::endl;
    m.unlock();
}

void backprojectAngleHBP(int id) {

    std::ostringstream filename;
    filename << "/home/billy/geant/SlitModelx2/proj/proj" << id << ".csv";
    std::string str = filename.str();
    //std::ifstream file("/home/billy/geant/EdgeModelx2/proj/proj45.csv");
	std::ifstream file(str);
	CSVRow row;
    vector<vector<double> > detected;
    vector<double> tmp;
    tmp.resize(5);
    double entryAngle, exitAngle;
	while(file >> row)
	{
        //vec.push_back(row);
        entryAngle = atan(lexical_cast<double>(row[4])/lexical_cast<double>(row[6]));
        exitAngle = atan(lexical_cast<double>(row[12])/lexical_cast<double>(row[14]));

        tmp[0] = lexical_cast<double>(row[1])*1000; // lateral entry position (entryPos) (mm)
        tmp[1] = entryAngle;
        tmp[2] = lexical_cast<double>(row[9])*1000; // lateral exit position (exitPos) (mm)
        tmp[3] = exitAngle;
        tmp[4] = lexical_cast<double>(row[16]); // projection
        detected.push_back(tmp);

	}

	//filter data
	cout << detected.size() << endl;
	filterData(detected);
    cout << detected.size() << endl;


    ReadWrite* r = new ReadWrite();
    Siddon* s = new Siddon(N, d);

    double angle = id*M_PI/180.0;
    for (int i = 0; i < detected.size(); i++) {
        Backproject* bp = new Backproject();
        bp->setAngle(angle);
        bp->radonHBP(detected[i], bpAngle, chordLength, *s, mut);
        delete bp;
    }

    //filename.str("");
    //filename << "output" << id << ".csv";
    //str = filename.str();
    //const char* fileStr = str.c_str();
    //r->writeCSV(bpAngle, fileStr);
    m.lock();
        counter++;
        std::cout << counter << std::endl;
        //std::cout << fileStr << std::endl;
    m.unlock();

}

//weighting the backprojections...
tuple<MatrixXd, MatrixXd> manipHBP2(vector<MatrixXd> & bpAngle) {

    ReadWrite r;
    MatrixXd bx = MatrixXd::Zero(N-1,N-1);
    MatrixXd by = MatrixXd::Zero(N-1,N-1);

    for (int k = 0; k < bpAngle.size(); k++) {
        for (int i = 0; i < bpAngle[k].rows(); i++) {
            for (int j = 0; j < bpAngle[k].rows(); j++) {
                if (chordLength[k](i,j) != 0) {
                    bpAngle[k](i,j) = bpAngle[k](i,j) / chordLength[k](i,j);
                }
                //fill the air parts
                if (bpAngle[k](i,j) == 0.0) {
                    bpAngle[k](i,j) = 0.4930;
                }
            }
        }
        std::ostringstream filename;
        filename.str("");
        filename << "b" << k << ".csv";
        string str = filename.str();
        const char* fileStr = str.c_str();
        r.writeCSV(bpAngle[k], fileStr);

    }


    for(int i = 0; i < bpAngle.size(); i++) {
        bx = bx + bpAngle[i] * -sin(i*M_PI/180);
        by = by + bpAngle[i] * cos(i * M_PI/180);
    }

    bx = bx * M_PI/180;
    by = by * M_PI/180;


    return make_tuple(bx, by);

}

// don't think we need this method...
tuple<MatrixXd, MatrixXd > manipHBP() {
    for (int k = 0; k < bpAngle.size(); k++) {
        for (int i = 0; i < bpAngle[k].rows(); i++) {
            for (int j = 0; j < bpAngle[k].rows(); j++) {
                if (chordLength[k](i,j) != 0) {
                    bpAngle[k](i,j) = bpAngle[k](i,j) / chordLength[k](i,j);
                }

                //fill the air parts
                if (bpAngle[k](i,j) == 0.0) {
                    bpAngle[k](i,j) = 0.4930;
                }
            }
        }
    }

    MatrixXd bpAngleCos = MatrixXd::Zero(460,460);
    MatrixXd bpAngleSin = MatrixXd::Zero(460,460);

    for (int i = 0; i < bpAngle.size(); i++) {
        bpAngleCos = bpAngleCos + bpAngle[i] * cos(i*M_PI/180.0);
        bpAngleSin = bpAngleSin + bpAngle[i] * sin(i*M_PI/180.0);
    }

    return make_tuple(bpAngleCos, bpAngleSin);
}

void testFn(int id) {
    std::ostringstream filename;
    vector<vector<double> > v(460, vector<double>(460,0));
    filename.str("");
    filename << "test" << id << ".csv";
    const char* fileStr = filename.str().c_str();
    cout << id << " : " << fileStr << endl;
    //ReadWrite* r = new ReadWrite();
    ReadWrite r;
    r.writeCSV(v, fileStr);

}

int main() {
	//Siddon* s = new Siddon(5, 1.0);


// MULTITHREADING CODE, 10 THREADS!


    int numOfThreads = 10;
    std::thread t[numOfThreads];

    string input = "";
    int option;

    while(true) {
        cout << "Enter 1 for HBP, 2 for BPF:\n>";
        getline(cin, input);

        stringstream myStream(input);
        if (myStream >> option)
            break;
        cout << "Invalid number, try again\n";
    }

    if (option == 1) {
        cout << "Weighted projections saved as bx.csv (sin theta) and by.csv (cos theta)\n";


        for (int j = 0; j < 18; j++) {
            for (int i = 0; i < numOfThreads; i++) {
                t[i] = std::thread (backprojectAngleHBP, i + numOfThreads*j);
            }

            for (int i = 0; i< numOfThreads; i++) {
                t[i].join();
            }
        }

        MatrixXd bx, by;
        tie(bx, by) = manipHBP2(bpAngle);

        ReadWrite r;

        r.writeCSV(bx,"bx.csv");
        r.writeCSV(by,"by.csv");
        return 0;
    }

    if (option == 2) {
        cout << "Projections saved as bp<angle>.csv\n";

        for (int j = 0; j < 18; j++) {
            for (int i = 0; i < numOfThreads; i++) {
                t[i] = std::thread (backprojectAngleBPF, i + numOfThreads*j);
            }

            for (int i = 0; i< numOfThreads; i++) {
                t[i].join();
            }
        }


        //backprojectAngleBPF(0);
        return 0;
    }

    if (option == 3) {
        backprojectAngleBPF(0);
    }

}


