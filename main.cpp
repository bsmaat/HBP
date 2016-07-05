#include <iostream>

#include <fstream>

#include <stdexcept>
#include "Siddon.h"
#include "Path.h"
#include "CSVRow.h"
#include "Backproject.h"
#include "ReadWrite.h"

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

int main() {
	//Siddon* s = new Siddon(5, 1.0);

    std::ifstream file("/home/billy/geant/EdgeModelx2/proj/proj45.csv");
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

    vector<vector<double> > bpAngle(460, vector<double>(460, 0));
    vector<vector<double> > chordLength(460, vector<double>(460,0));

    ReadWrite* r = new ReadWrite();
    Siddon* s = new Siddon(461, 1);

    for (int i = 0; i < detected.size(); i++) {
        Backproject* bp = new Backproject();
        bp->radon(detected[i], bpAngle, chordLength, *s);
        delete bp;
    }

    //vector<vector<double> > tt;
    //std::transform(bpAngle.begin(), bpAngle.end(), chordLength.begin(), tt.begin(), std::divides<double>());


    for (int i = 0; i < bpAngle.size(); i++) {
        for (int j = 0; j < bpAngle.size(); j++) {
            if (chordLength[i][j] != 0)
                bpAngle[i][j] = bpAngle[i][j] / chordLength[i][j];
        }
    }

    // fill the air parts...
    for (vector<vector<double> >::iterator i = bpAngle.begin(); i != bpAngle.end(); ++i) {
        std::replace(i->begin(), i->end(), 0.0, 0.4950);
    }


    r->writeCSV(bpAngle, "output.csv");


	return 0;

}


