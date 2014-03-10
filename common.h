/*
 * common.h
 *
 *  Created on: Aug 5, 2013
 *      Author: mnandan
 */

#ifndef COMMON_H_
#define COMMON_H_

#include <vector>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <algorithm>

using namespace std;

typedef unsigned int UINT;
typedef char LABEL_T;
#define M_VAL (UINT)2
#define INF HUGE_VAL

typedef struct {
	double dist;
	UINT ind;
} DistDat;

struct FeatType {
	UINT fNum; // feature number
	double fVal;
};

struct DataVect {
	UINT numFeats;
	UINT index; // index of vector in input file
	double nrm; // norm of vector
	FeatType * F; // features
};

class TrainDat {
	DataVect * X;
	UINT D; // largest index of features among all vectors
	UINT R;
	ofstream hOUT, wOUT;
	char aeFile[1024];
	UINT vectNum;
	UINT totNumVects;
    double ** W;	//factor W
	double ** H;	//factor H

public:
	const DataVect * XC;
    const double * const * WC;	//factor W
	const double *const * HC;	//factor H

	TrainDat(UINT D, UINT N, UINT R, char *hFname, char *wFname) {
		this->D = D;
		totNumVects = N;
		this->R = R;
		X = new DataVect[totNumVects];
		H = new double*[totNumVects];
		for (UINT ind = 0; ind < totNumVects; ind++) {
			X[ind].F = NULL;
			X[ind].numFeats = 0;
			H[ind] = new double [R];
			for (UINT ind2 = 0; ind2 < R; ind2++)
				H[ind][ind2] = 0;
		}
		XC = X;
		HC = H;
		vectNum = 0;

		W = new double*[R];
		for (UINT ind = 0; ind < R; ind++) {
			W[ind] = new double[D];
			for (UINT ind2 = 0; ind2 < D; ind2++)
				W[ind][ind2] = 0;
		}
		WC = W;

		hOUT.open(hFname);
		wOUT.open(wFname);
	}

	~TrainDat() {
		hOUT.close();
		wOUT.close();
		for (UINT ind = 0; ind < totNumVects; ind++) {
			delete [] H[ind];
			if (X[ind].F != NULL)
				delete[] X[ind].F;
		}
		delete [] H;
		for (UINT ind = 0; ind < R; ind++)
			delete [] W[ind];
		delete [] W;
	}

	void addVector(vector<FeatType> &F, UINT numFeats) {
		X[vectNum].numFeats = numFeats;
		X[vectNum].index = vectNum;
		X[vectNum].nrm = 0;
		if (numFeats > 0) {
			X[vectNum].F = new FeatType[numFeats];
			for (UINT i = 0; i < numFeats; i++) {
				X[vectNum].F[i] = F[i];
				X[vectNum].nrm += (F[i].fVal) * (F[i].fVal);
			}
		} else {
			X[vectNum].F = NULL;
		}
		vectNum++;
	}

	void copyXW() {	//initialize W matrix with the first R vectors in X
		for(UINT i = 0; i<R; i++) {
			FeatType * F = X[i].F;
			for(UINT j = 0; j < X[i].numFeats; j++)
				W[i][F[j].fNum] = F[j].fVal;
		}
		UINT i = 0;
		while(i < totNumVects) {	//Put X in original order
		    if(i == X[i].index)
		    	i++;
		    else
		    	swap(X[i], X[X[i].index]);
		}
	}

	void swapX(UINT i, UINT j) {
		swap(X[i], X[j]);
	}

	void putX(UINT si, UINT p, DataVect *X2, DistDat *distVals) {
		for (UINT k = 0; k < p; k++) //mass swap
			X2[k] = X[distVals[k].ind];
		for (UINT k = 0; k < p; k++) {
			UINT j = k + si;
			X[j] = X2[k];
		}
	}

	void putH(UINT i, double * lambda) {
		for (UINT k = 0; k < R; k++)
			H[i][k] = lambda[k];
	}

	void putW(UINT i, double * lambda) {
		for (UINT k = 0; k < D; k++)
			W[i][k] = lambda[k];
	}

	void putWval(UINT i, UINT f, double val) {
		W[i][f] = val;
	}

	void putHval(UINT i, UINT f, double val) {
		H[i][f] = val;
	}

	void writeWH() {
		for(UINT ind = 0; ind < totNumVects; ind++) {
			if (hOUT.is_open() && hOUT.good()) {
				hOUT << '1';
				for (UINT k = 0; k < R; k++)
					if(H[ind][k] != 0)
						hOUT <<" "<< k + 1<<':'<<setprecision(8) << H[ind][k];
				hOUT<<endl;
			} else {
				cerr << "Cannot write matrix factor H to file\n";
				exit(1);
			}
		}
		for(UINT ind = 0; ind < R; ind++) {
			if (wOUT.is_open() && wOUT.good()) {
				wOUT << '1';
				for (UINT k = 0; k < D; k++)
					if(W[ind][k] != 0)
						wOUT <<" "<< k + 1<<':'<<setprecision(8) << W[ind][k];
				wOUT<<endl;
			} else {
				cerr << "Cannot write matrix factor W to file\n";
				exit(1);
			}
		}
	}

	double calcNorm() {
		double *hwVect = new double[D];
		double frobNorm = 0;
		for(UINT i = 0;i < totNumVects;i++) {
			for(UINT j = 0;j < D;j++)
				hwVect[j] = 0;
			for(UINT j = 0;j < R;j++) {
				if(H[i][j] != 0 )
					for(UINT k = 0;k < D;k++)
						hwVect[k] += H[i][j]*W[j][k];
			}
			FeatType * F = X[i].F;
			UINT prevFnum = 0;
			for(UINT fID = 0;fID < X[i].numFeats;fID++) {
				for(;prevFnum < F[fID].fNum; prevFnum++)
					frobNorm += hwVect[prevFnum] * hwVect[prevFnum];
				double val = F[fID].fVal - hwVect[F[fID].fNum];
				frobNorm += val*val;
				prevFnum++;
			}
		}
		delete [] hwVect;
		return frobNorm;
	}

	UINT getD() {
		return D;
	}

	UINT getN() {
		return totNumVects;
	}

	UINT getR() {
		return R;
	}

	void clearW() {
		for (UINT ind = 0; ind < R; ind++)
			for (UINT ind2 = 0; ind2 < D; ind2++)
				W[ind][ind2] = 0;
	}

	void clearH() {
		for (UINT ind = 0; ind < totNumVects; ind++)
			for (UINT ind2 = 0; ind2 < R; ind2++)
				H[ind][ind2] = 0;
	}

};

#endif /* COMMON_H_ */
