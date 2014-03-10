/*
 * common.h
 *
 *  Created on: Aug 3, 2013
 *      Author: mnandan
 */

#ifndef FILE_INT_H_
#define FILE_INT_H_

#include <vector>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include "common.h"

using namespace std;

class fHand_T {
	ifstream inpFH;
	UINT lineNum;
	bool *delims1, *delims2;
	TrainDat *TrD;
	vector <FeatType> Ftemp;
	UINT D;
	UINT R;
	UINT totVects;
	UINT maxNNZ;
public:
	fHand_T(char *fName, TrainDat *TrD) {
		inpFH.open(fName);
		if(! (inpFH.is_open() && inpFH.good())) {
			cerr << "Cannot open file " << fName << endl;
			exit(1);
		}
		lineNum = 0;

		delims1 = new bool[128];
		for (int i = 0; i < 128; i++)
			delims1[i] = false;
		delims1[' '] = true;
		delims1['\t'] = true;

		delims2 = new bool[128];
		for (int i = 0; i < 128; i++)
			delims2[i] = false;
		delims2[' '] = true;
		delims2['\t'] = true;
		delims2[':'] = true;
		this->TrD = TrD;
		R = TrD->getR();
		D = TrD->getD();
		totVects = TrD->getN();
		maxNNZ = 0;
	}

	~fHand_T() {
		inpFH.close();
		delete [] delims1;
		delete [] delims2;
		//vector<FeatType>().swap(Ftemp);
	}

	UINT getLineNum() {
		return lineNum;
	}

	bool procLine();
	bool readFile();
};
#endif /* FILE_INT_H_ */
