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
int counter = 0;

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

    vector<vector<double> > bpAngle(460, vector<double>(460, 0));
    vector<vector<double> > chordLength(460, vector<double>(460,0));

    ReadWrite* r = new ReadWrite();
    Siddon* s = new Siddon(461, 1);

    double angle = id*M_PI/180.0;
    for (int i = 0; i < detected.size(); i++) {
        Backproject* bp = new Backproject();
        bp->setAngle(angle);
        bp->radon(detected[i], bpAngle, chordLength, *s);
        delete bp;
    }

    for (int i = 0; i < bpAngle.size(); i++) {
        for (int j = 0; j < bpAngle.size(); j++) {
            if (chordLength[i][j] != 0) {
                bpAngle[i][j] = bpAngle[i][j] / chordLength[i][j];
                //bpAngle[i][j]= chordLength[i][j];
            }
        }
    }

    //cout << "Val: " << setprecision(15) << bpAngle[213][162] << endl;
    // fill the air parts...
    for (vector<vector<double> >::iterator i = bpAngle.begin(); i != bpAngle.end(); ++i) {
        std::replace(i->begin(), i->end(), 0.0, 0.4930);
    }


    filename.str("");
    filename << "output" << id << ".csv";
    str = filename.str();
    const char* fileStr = str.c_str();
    r->writeCSV(bpAngle, fileStr);
    m.lock();
        counter++;
        std::cout << counter << std::endl;
        std::cout << fileStr << std::endl;
    m.unlock();
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


/* MULTITHREADING CODE, 10 THREADS!
    int numOfThreads = 10;
    std::thread t[numOfThreads];


    for (int j = 0; j < 18; j++) {
        for (int i = 0; i < numOfThreads; i++) {
            //cout << "i: " << i << endl;
            t[i] = std::thread (backprojectAngle, i + numOfThreads*j);
        }

        for (int i = 0; i< numOfThreads; i++) {
            t[i].join();
        }
    }

*/

    backprojectAngle(0);
	return 0;

}


