/*
 * deriveAE.cpp
 *
 *  Created on: Aug 6, 2013
 *      Author: mnandan
 */

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include "deriveAE.h"

using namespace std;

#define REG_PARAM (double) 0.0
double DeriveAE::deriveW( ) {
	double *lambda = new double[D];
	for (UINT ind = 0; ind < R; ind++) {
		for (UINT ind2 = 0; ind2 < R; ind2++)
			rpCache[ind][ind2] = 0;
		for (UINT ind2 = 0; ind2 < D; ind2++)
			GM[ind][ind2] = 0;
	}
	FeatType * F;
	UINT xFnum;
	double tempH;
	for (UINT ind = 0; ind < N; ind++) {
		for (UINT ind2 = 0; ind2 < R; ind2++) {
			double hVal = H[ind][ind2];
			rpCache[ind2][ind2] += hVal*hVal;
			for (UINT ind3 = ind2; ind3 > 0; ind3--) {
				UINT ind3u = ind3 - 1;
				tempH = hVal*H[ind][ind3u];
				rpCache[ind2][ind3u] += tempH;
				rpCache[ind3u][ind2] += tempH;
			}
		}
		F = trDat->XC[ind].F;
		xFnum = trDat->XC[ind].numFeats;
		for (UINT fInd = 0; fInd < xFnum; fInd++)
			for (UINT ind2 = 0; ind2 < R; ind2++) {
				GM[ind2][F[fInd].fNum] -= H[ind][ind2]*F[fInd].fVal;
			}
	}

	for (UINT ind = 0; ind < R; ind++)
		for (UINT fInd = 0; fInd < D; fInd++)
			for (UINT ind2 = 0; ind2 < R; ind2++)
				GM[ind][fInd] += rpCache[ind][ind2]*W[ind2][fInd];

	// optimization step
	double gPod, lambdaNew, delta, gradSum, prevGradSum = INF;
	UINT iter = 0;
	while (iter++ < 1000) {
		for (UINT i = 0; i < R; i++) {
			UINT j = i + rand() % (R - i);
			swap(index[i], index[j]);
		}
		gradSum = 0;
		for (UINT s = 0; s < R; s++) {
			UINT vI = index[s]; //vector index
			if(rpCache[vI][vI] != 0) {
				for (UINT fInd = 0; fInd < D; fInd++) {
					gPod = GM[vI][fInd]/(REG_PARAM + rpCache[vI][vI]);
					lambdaNew = max(W[vI][fInd] - gPod, 0.0); //constraint
					delta = lambdaNew - W[vI][fInd];
					if(delta != 0.0)
						trDat->putWval(vI,fInd,lambdaNew);		//W[vI][fInd] = lambdaNew;
					lambda[fInd] = delta;
				}
				for (UINT gI = 0; gI < R; gI++){ //better cpu cache usage
					for (UINT fInd = 0; fInd < D; fInd++)
						if(lambda[fInd] != 0)
							GM[gI][fInd] += lambda[fInd] *rpCache[gI][vI];
				}

			}
		}
		for (UINT s = 0; s < R; s++)
			for (UINT fInd = 0; fInd < D; fInd++)
				gradSum += GM[s][fInd]*GM[s][fInd];
		gradSum = sqrt(gradSum);
		if(fabs(prevGradSum - gradSum) < 0.001) {
			//cout << "deltaSum = "<< deltaSum <<endl;
			break;
		}
		prevGradSum = gradSum;
	}

	delete [] lambda;

	return gradSum;
}

double DeriveAE::deriveH() {
	updateCacheW();	//update rpCache

// derive lambda for other vectors
	double frobNormSq = 0;
	for (UINT ind = 0; ind < N; ind++) {
		for (UINT rpInd2 = 0; rpInd2 < R; rpInd2++)
			xTz[rpInd2] = getDotProductXW(ind, rpInd2);
		frobNormSq += deriveHi(ind);
	}
	return sqrt(frobNormSq);
}

double DeriveAE::deriveHi(UINT hInd) {
	// initialize
	double PG, Gmax, Gmin;
	double GmaxOld = INF, GminOld = -INF;
	UINT activeSize = R;
	for (UINT i = 0; i < R; i++) {
		G[i] = -xTz[i];
		for (UINT k = 0; k < R; k++)
			G[i] += rpCache[i][k] * H[hInd][k];
	}
	// optimization step
	int iter = 0;
	while (iter++ < 1000) {
		for (UINT i = 0; i < activeSize; i++) {
			UINT j = i + rand() % (activeSize - i);
			swap(index[i], index[j]);
		}
		Gmax = -INF;
		Gmin = INF;
		for (UINT s = 0; s < activeSize; s++) {
			UINT vI = index[s]; //relative vector index
			PG = 0;
			if (H[hInd][vI] == 0) {
				if (G[vI] > GmaxOld) {
					activeSize--;
					swap(index[s], index[activeSize]);
					s--;
					continue;
				} else if (G[vI] < 0)
					PG = G[vI];
			} else
				PG = G[vI];

			if (fabs(PG) > TAU) {
				double lambdaNew = max(H[hInd][vI] - PG / (REG_PARAM + rpCache[vI][vI]), 0.0); //constraint
				double delta = lambdaNew - H[hInd][vI];
				if(delta != 0) {
					trDat->putHval(hInd,vI,lambdaNew);
					for (UINT k = 0; k < R; k++)
						G[k] += rpCache[vI][k] * delta;
				}
			}
			Gmax = max(Gmax, PG);
			Gmin = min(Gmin, PG);
		}
		if (Gmax - Gmin <= 0.1) {
			if (activeSize == R)
				break;
			else {
				activeSize = R;
				//cout << "*";
				GmaxOld = INF;
				GminOld = -INF;
				continue;
			}
		}
		GmaxOld = Gmax;
		GminOld = Gmin;
		if (GmaxOld <= 0)
			GmaxOld = INF;
		if (GminOld >= 0)
			GminOld = -INF;
	}
	double normSum = X[hInd].nrm;
	for(UINT i = 0; i < R; i++) {
		if(H[hInd][i] != 0) {
			normSum -= 2*H[hInd][i]*xTz[i];
			double subSum = 0;
			for(UINT j = 0; j < R; j++)
				subSum += rpCache[i][j]*H[hInd][j];
			normSum += subSum*H[hInd][i];
		}
	}
	return normSum;
}
