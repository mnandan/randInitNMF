#include "fileInt.h"

using namespace std;

bool fHand_T::readFile() {
	while (lineNum < totVects && procLine());

	if(lineNum == totVects)
		return true;
	else
		return false;
}

bool fHand_T:: procLine() {
	string line;
	if (inpFH.good() && getline(inpFH, line)) {
		++lineNum;
		UINT lineSz = (UINT)line.length();
		char * lineC = (char *)line.c_str();
		char *p = lineC;

		UINT fNum;
		double fVal;
		FeatType tempF;
		char *endPtr;
		UINT featInd = 0;

		UINT cInd = 0;
		for(;delims1[(int)lineC[cInd]] && cInd < lineSz; ++cInd);		//remove initial spaces
		for(;!delims1[(int)lineC[cInd]] && cInd < lineSz; ++cInd);		//find label characters
		++cInd;

		while (cInd < lineSz) {
			for(;delims1[(int)lineC[cInd]] && cInd < lineSz; ++cInd);		//remove initial spaces
			p = lineC + cInd;	//start of sub-string
			for(;!delims2[(int)lineC[cInd]] && cInd < lineSz; ++cInd);		//find fNum characters
			lineC[cInd] = 0;	//end of sub-string
			fNum = (UINT) strtoll(p, &endPtr, 10);
			fNum -= 1;		// LIBSVM format features are numbered from 1, need it to be from 0
			if (endPtr == p || *endPtr != '\0') {
				cerr << "Error in input file read\n";
				exit(1);
			}
			++cInd;

			for(;delims1[(int)lineC[cInd]] && cInd < lineSz; ++cInd);		//remove initial spaces
			p = lineC + cInd;
			for(;!delims1[(int)lineC[cInd]] && cInd < lineSz; ++cInd);		//find fVal characters
			lineC[cInd] = 0;
			fVal = strtod(p, &endPtr);
			if (endPtr == p || (*endPtr != '\0' && !isspace(*endPtr))) {
				cerr << "Error in input file read\n";
				exit(1);
			}
			if(featInd < maxNNZ) {
				Ftemp[featInd].fNum = fNum;
				Ftemp[featInd].fVal = fVal;
			} else {
				tempF.fNum = fNum;
				tempF.fVal = fVal;
				Ftemp.push_back(tempF);
				++maxNNZ;
			}
			++featInd;
			++cInd;
		}
		TrD->addVector(Ftemp, featInd);
	}

	if(inpFH.good())
		return true;
	else
		return false;
}
