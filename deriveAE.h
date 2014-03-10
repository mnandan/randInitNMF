/*
 * deriveAE.h
 *
 *  Created on: Aug 6, 2013
 *      Author: mnandan
 */

#ifndef DERIVEAE_H_
#define DERIVEAE_H_

#include "common.h"

#define TAU 1e-12
#define LAMBDA_MAX 1.0
#define LAMBDA_MIN 0

class DeriveAE {
	UINT R, D, N;
	TrainDat *trDat;
	double getDotProductXW(UINT i, UINT j) {
		double dotProduct = 0;
		FeatType *F = X[i].F;
		for(UINT ind = 0; ind < X[i].numFeats; ind++) {
			double wfVal = W[j][F[ind].fNum];
			if(wfVal != 0)
				dotProduct += wfVal*F[ind].fVal;
		}
		return dotProduct;
	}

	double getDotProductW(UINT i, UINT j) {
		double dotProduct = 0;
		for(UINT fInd = 0; fInd < D; fInd++)
			if (W[i][fInd] != 0 && W[j][fInd] != 0)
				dotProduct += W[i][fInd]*W[j][fInd];
		return dotProduct;
	}

	void updateCacheW() {
		for (UINT ind1 = 0; ind1 < R; ind1++) {
			rpCache[ind1][ind1] = getDotProductW(ind1, ind1);
			for (UINT ind2 = ind1; ind2 > 0; ind2--) {
				UINT ind2u = ind2 - 1;
				rpCache[ind2u][ind1] = getDotProductW(ind1, ind2u);
				rpCache[ind1][ind2u] = rpCache[ind2u][ind1];
			}
		}
	}

	double deriveHi(UINT hInd);
protected:
	double ** rpCache, **GM;
	const double * const* W, * const*H;
	UINT *index;
	double *G, *xTz;
	const DataVect* X;
public:
	DeriveAE(TrainDat *trDat) {
		this->trDat = trDat;
		R = trDat->getR();
		N = trDat->getN();
		D = trDat->getD();
		W = trDat->WC;
		H = trDat->HC;
		X = trDat->XC;

		G = new double[R];
		xTz = new double[R];
		index = new UINT[R];
		rpCache = new double*[R]; // size RxR to store hTh
		GM = new double*[R];	//gradient matrix RxD
		for (UINT ind = 0; ind < R; ind++) {
			index[ind] = ind;
			rpCache[ind] = new double[R];
			GM[ind] = new double[D];
		}

	}

	~DeriveAE() {
		for (UINT ind = 0; ind < R; ind++) {
			delete[] rpCache[ind];
			delete[] GM[ind];
		}
		delete [] rpCache;
		delete [] GM;
		delete [] index;
		delete [] G;
		delete [] xTz;
	}


	double deriveW();
	double deriveH();
};


#endif /* DERIVEAE_H_ */
