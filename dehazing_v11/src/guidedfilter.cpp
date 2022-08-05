/*
	The original guided filter is just converted version of He et al.'s matlab code
	The detailed description of the code and the algorithm is presented
	in "http://research.microsoft.com/en-us/um/people/kahe/eccv10/index.html"

	The fast guided filter is approximated filter to reduce the computational 
	complexity of filtering. We only apply the filter to sampled points. The 
	similar works are presented in the paper, J.-Y. Kim, L.-S. Kim, S.-H. Hwang, 
	"An advanced contrast enhancement using partially overlapped subblock
	histogram equalization," IEEE Trans. Circuits and Syst. Video Technol. 
	11 (4) (2001) 475-484. Also, at each window, the gussian weight is applied.

	The guided filter with shiftable window is the enhanced algorithm of original
	guided image filtering by shifting the filter window. 

	Last updated: 2013-02-06
	Author: Jin-Hwan, Kim.
 */
#include "dehazing.h"
#include <assert.h>

/*
	Function: CalcAcoeff
	Description: calculate the coefficent "a" of guided filter (intrinsic function for guide filtering)
	Parameters: 
		pfSigma - Sigma + eps * eye(3) --> see the paper and original matlab code
		pfCov - Cov of image
		nIdx - location of image(index).
	Return:
		pfA1, pfA2, pfA3 - coefficient of "a" at each color channel
 */
void dehazing::CalcAcoeff(double* pfSigma, double* pfCov, double* pfA1, double* pfA2, double* pfA3, int nIdx)
{
	double fDet;
	double fOneOverDeterminant;
	double pfInvSig[9];

	int nIdx9 = nIdx*9;

	// a_k = (sum_i(I_i*p_i-mu_k*p_k)/(abs(omega)*(sigma_k^2+epsilon))
	fDet = ( (pfSigma[nIdx9]*(pfSigma[nIdx9+4] * pfSigma[nIdx9+8] - pfSigma[nIdx9+5] * pfSigma[nIdx9+7]))
		-	( pfSigma[nIdx9+1]*(pfSigma[nIdx9+3] * pfSigma[nIdx9+8] - pfSigma[nIdx9+5] * pfSigma[nIdx9+6]))
		+	( pfSigma[nIdx9+2]*(pfSigma[nIdx9+3] * pfSigma[nIdx9+7] - pfSigma[nIdx9+4] * pfSigma[nIdx9+6])) );
	fOneOverDeterminant = 1.0f/fDet;
    assert(!isnan(fOneOverDeterminant));

	pfInvSig [0] =  (( pfSigma[nIdx9+4]*pfSigma[nIdx9+8] ) - (pfSigma[nIdx9+5]*pfSigma[nIdx9+7]))*fOneOverDeterminant;
	pfInvSig [1] = -(( pfSigma[nIdx9+1]*pfSigma[nIdx9+8] ) - (pfSigma[nIdx9+2]*pfSigma[nIdx9+7]))*fOneOverDeterminant;
	pfInvSig [2] =  (( pfSigma[nIdx9+1]*pfSigma[nIdx9+5] ) - (pfSigma[nIdx9+2]*pfSigma[nIdx9+4]))*fOneOverDeterminant;
	pfInvSig [3] = -(( pfSigma[nIdx9+3]*pfSigma[nIdx9+8] ) - (pfSigma[nIdx9+5]*pfSigma[nIdx9+6]))*fOneOverDeterminant;
	pfInvSig [4] =  (( pfSigma[nIdx9]*pfSigma[nIdx9+8] ) - (pfSigma[nIdx9+2]*pfSigma[nIdx9+6]))*fOneOverDeterminant;
	pfInvSig [5] = -(( pfSigma[nIdx9]*pfSigma[nIdx9+5] ) - (pfSigma[nIdx9+2]*pfSigma[nIdx9+3]))*fOneOverDeterminant;
	pfInvSig [6] =  (( pfSigma[nIdx9+3]*pfSigma[nIdx9+7] ) - (pfSigma[nIdx9+4]*pfSigma[nIdx9+6]))*fOneOverDeterminant;
	pfInvSig [7] = -(( pfSigma[nIdx9]*pfSigma[nIdx9+7] ) - (pfSigma[nIdx9+1]*pfSigma[nIdx9+6]))*fOneOverDeterminant;
	pfInvSig [8] =  (( pfSigma[nIdx9]*pfSigma[nIdx9+4] ) - (pfSigma[nIdx9+1]*pfSigma[nIdx9+3]))*fOneOverDeterminant;

	int nIdx3 = nIdx*3;

	pfA1[nIdx] = pfCov[nIdx3]*pfInvSig[0]+pfCov[nIdx3+1]*pfInvSig[3]+pfCov[nIdx3+2]*pfInvSig[6];
	pfA2[nIdx] = pfCov[nIdx3]*pfInvSig[1]+pfCov[nIdx3+1]*pfInvSig[4]+pfCov[nIdx3+2]*pfInvSig[7];
	pfA3[nIdx] = pfCov[nIdx3]*pfInvSig[2]+pfCov[nIdx3+1]*pfInvSig[5]+pfCov[nIdx3+2]*pfInvSig[8];
}

/*
	Function: BoxFilter
	Description: cummulative function for calculating the integral image (It may apply other arraies.)
	Parameters:
		pfInArray - input array
		nR - radius of filter window
		nWid - width of array
		nHei - height of array
	Return:
		fOutArray - output array (integrated array)
 */
void dehazing::BoxFilter(double* pfInArray, int nR, int nWid, int nHei, double*& fOutArray)
{
    checkdouble(pfInArray, nWid * nHei);
   // for (int i = 0; i < nWid * nHei; i++) {
  //      assert(!isnan(pfInArray[i]));
  //  }

	double* pfArrayCum = newdouble(nWid*nHei);

	//cumulative sum over Y axis
	for ( int nX = 0; nX < nWid; nX++ )
		pfArrayCum[nX] = pfInArray[nX];

	for ( int nIdx = nWid; nIdx <nWid*nHei; nIdx++ )
		pfArrayCum[nIdx] = pfArrayCum[nIdx-nWid] + pfInArray[nIdx];

	//difference over Y axis
	for ( int nIdx = 0; nIdx < nWid*(nR+1); nIdx++)
		fOutArray[nIdx] = pfArrayCum[nIdx+nR*nWid];

    checkdouble(fOutArray, nWid * nHei);

	for ( int nIdx = (nR+1)*nWid; nIdx < (nHei-nR)*nWid; nIdx++ )
		fOutArray[nIdx] = pfArrayCum[nIdx+nR*nWid] - pfArrayCum[nIdx-nR*nWid-nWid];

    checkdouble(fOutArray, nWid * nHei);

	for ( int nY = nHei-nR; nY < nHei; nY++ )
		for ( int nX = 0; nX < nWid; nX++ )
			fOutArray[nY*nWid+nX] = pfArrayCum[(nHei-1)*nWid+nX] - pfArrayCum[(nY-nR-1)*nWid+nX];

    checkdouble(fOutArray, nWid * nHei);
	
	//cumulative sum over X axis
	for ( int nIdx = 0; nIdx < nHei*nWid; nIdx += nWid )
		pfArrayCum[nIdx] = fOutArray[nIdx];

	for ( int nY = 0; nY < nHei*nWid; nY+=nWid )
		for ( int nX = 1; nX < nWid; nX++ )
			pfArrayCum[nY+nX] = pfArrayCum[nY+nX-1]+fOutArray[nY+nX];

	//difference over Y axis
	for ( int nY = 0; nY < nHei*nWid; nY+=nWid )
		for ( int nX = 0; nX < nR+1; nX++ )
			fOutArray[nY+nX] = pfArrayCum[nY+nX+nR];

    checkdouble(fOutArray, nWid * nHei);

	for ( int nY = 0; nY < nHei*nWid; nY+=nWid )
		for ( int nX = nR+1; nX < nWid-nR; nX++ )
			fOutArray[nY+nX] = pfArrayCum[nY+nX+nR] - pfArrayCum[nY+nX-nR-1];

    checkdouble(fOutArray, nWid * nHei);

	for ( int nY = 0; nY < nHei*nWid; nY+=nWid )
		for ( int nX = nWid-nR; nX < nWid; nX++ )
			fOutArray[nY+nX] = pfArrayCum[nY+nWid-1] - pfArrayCum[nY+nX-nR-1];

    checkdouble(fOutArray, nWid * nHei);

	delete []pfArrayCum;
}

/*
	Function: BoxFilter (for 3D array)
	Description: cummulative function for calculating the integral image (It may apply other arraies.)
	Parameters:
		pfInArray1 - input array D1
		pfInArray2 - input array D2
		pfInArray3 - input array D3
		nR - radius of filter window
		nWid - width of array
		nHei - height of array
	Return:
		fOutArray1 - output array D1(integrated array)
		fOutArray1 - output array D2(integrated array)
		fOutArray1 - output array D3(integrated array)
 */
void dehazing::BoxFilter(double* pfInArray1, double* pfInArray2, double* pfInArray3, int nR, int nWid, int nHei, double*& pfOutArray1, double*& pfOutArray2, double*& pfOutArray3)
{
	double* pfArrayCum1 = newdouble(nWid*nHei);
	double* pfArrayCum2 = newdouble(nWid*nHei);
	double* pfArrayCum3 = newdouble(nWid*nHei);

	//cumulative sum over Y axis
	for ( int nX = 0; nX < nWid; nX++ )
	{
		pfArrayCum1[nX] = pfInArray1[nX];
		pfArrayCum2[nX] = pfInArray2[nX];
		pfArrayCum3[nX] = pfInArray3[nX];
	}

	for ( int nIdx = nWid; nIdx < nHei*nWid; nIdx++ )
	{
		pfArrayCum1[nIdx] = pfArrayCum1[nIdx-nWid]+pfInArray1[nIdx];
		pfArrayCum2[nIdx] = pfArrayCum2[nIdx-nWid]+pfInArray2[nIdx];
		pfArrayCum3[nIdx] = pfArrayCum3[nIdx-nWid]+pfInArray3[nIdx];
	}

	//difference over Y axis
	for ( int nIdx = 0; nIdx < (nR+1)*nWid; nIdx++ )
	{
		pfOutArray1[nIdx] = pfArrayCum1[nIdx+nR*nWid];
		pfOutArray2[nIdx] = pfArrayCum2[nIdx+nR*nWid];
		pfOutArray3[nIdx] = pfArrayCum3[nIdx+nR*nWid];
	}

	for ( int nIdx = (nR+1)*nWid; nIdx < (nHei-nR)*nWid; nIdx++ )
	{
		pfOutArray1[nIdx] = pfArrayCum1[nIdx+nR*nWid] - pfArrayCum1[nIdx-(nR+1)*nWid];
		pfOutArray2[nIdx] = pfArrayCum2[nIdx+nR*nWid] - pfArrayCum2[nIdx-(nR+1)*nWid];
		pfOutArray3[nIdx] = pfArrayCum3[nIdx+nR*nWid] - pfArrayCum3[nIdx-(nR+1)*nWid];
	}

	for ( int nY = nHei-nR; nY < nHei; nY++ )
	{
		for ( int nX = 0; nX < nWid; nX++ )
		{
			pfOutArray1[nY*nWid+nX] = pfArrayCum1[(nHei-1)*nWid+nX] - pfArrayCum1[(nY-nR-1)*nWid+nX];
			pfOutArray2[nY*nWid+nX] = pfArrayCum2[(nHei-1)*nWid+nX] - pfArrayCum2[(nY-nR-1)*nWid+nX];
			pfOutArray3[nY*nWid+nX] = pfArrayCum3[(nHei-1)*nWid+nX] - pfArrayCum3[(nY-nR-1)*nWid+nX];
		}
	}

	//cumulative sum over X axis
	for ( int nY = 0; nY < nHei*nWid; nY+=nWid )
	{
		pfArrayCum1[nY] = pfOutArray1[nY];
		pfArrayCum2[nY] = pfOutArray2[nY];
		pfArrayCum3[nY] = pfOutArray3[nY];
	}

	for ( int nY = 0; nY < nHei*nWid; nY+=nWid )
	{
		for ( int nX = 1; nX < nWid; nX++ )
		{
			pfArrayCum1[nY+nX] = pfArrayCum1[nY+nX-1]+pfOutArray1[nY+nX];
			pfArrayCum2[nY+nX] = pfArrayCum2[nY+nX-1]+pfOutArray2[nY+nX];
			pfArrayCum3[nY+nX] = pfArrayCum3[nY+nX-1]+pfOutArray3[nY+nX];
		}
	}

	//difference over Y axis
	for ( int nY = 0; nY < nHei*nWid; nY+=nWid )
	{
		for ( int nX = 0; nX < nR+1; nX++ )
		{
			pfOutArray1[nY+nX] = pfArrayCum1[nY+nX+nR];
			pfOutArray2[nY+nX] = pfArrayCum2[nY+nX+nR];
			pfOutArray3[nY+nX] = pfArrayCum3[nY+nX+nR];
		}
	}

	for ( int nY = 0; nY < nHei*nWid; nY+=nWid )
	{
		for ( int nX = nR+1; nX < nWid-nR; nX++ )
		{	
			pfOutArray1[nY+nX] = pfArrayCum1[nY+nX+nR] - pfArrayCum1[nY+nX-nR-1];
			pfOutArray2[nY+nX] = pfArrayCum2[nY+nX+nR] - pfArrayCum2[nY+nX-nR-1];
			pfOutArray3[nY+nX] = pfArrayCum3[nY+nX+nR] - pfArrayCum3[nY+nX-nR-1];
		}
	}

	for ( int nY = 0; nY < nHei*nWid; nY+=nWid )
	{
		for ( int nX = nWid-nR; nX < nWid; nX++ )
		{
			pfOutArray1[nY+nX] = pfArrayCum1[nY+nWid-1] - pfArrayCum1[nY+nX-nR-1];
			pfOutArray2[nY+nX] = pfArrayCum2[nY+nWid-1] - pfArrayCum2[nY+nX-nR-1];
			pfOutArray3[nY+nX] = pfArrayCum3[nY+nWid-1] - pfArrayCum3[nY+nX-nR-1];
		}
	}

	delete []pfArrayCum1;
	delete []pfArrayCum2;
	delete []pfArrayCum3;
}

/*
	Function: GuidedFilterY
	Description: the original guided filter for gray image. This function is used for image dehazing.
		The video dehazing algorithm uses appoximated filter for fast refinement. 
	Parameter:
		nW - width of array
		nH - height of array
		fEps - epsilon
	(member variable)
		m_pfTransmission - initial transmission (block_based)
		m_pnYImg - guidance image (Y image)
	Return:
		m_pfTransmissionR - filtered transmission
 */
void dehazing::GuidedFilterY(int nW, int nH, double fEps)
{
	double* pfImageY		= newdouble(nW*nH);
	double* pfInitN		= newdouble(nW*nH);
	double* pfInitMeanIp = newdouble(nW*nH);
	double* pfMeanP		= newdouble(nW*nH);
	double* pfA			= newdouble(nW*nH);
	double* pfB			= newdouble(nW*nH);
	double* pfOutA		= newdouble(nW*nH);
	double* pfOutB		= newdouble(nW*nH);
	double* pfN			= newdouble(nW*nH);
	double* pfMeanI		= newdouble(nW*nH);
	double* pfMeanIp		= newdouble(nW*nH);
	double* pfInitvarI	= newdouble(nW*nH);
	double* pfvarI		= newdouble(nW*nH);
	double* pfCovIp		= newdouble(nW*nH);

	int nIdx;

	// Converting to doubleing point array
	for(nIdx = 0 ; nIdx < nW*nH ; nIdx++ )
	{
		pfImageY[nIdx] = (double)m_pnYImg[nIdx];
	}
	//////////////////////////////////////////////////////////////////////////

	// Make an integral image
	for (nIdx = 0; nIdx < nW*nH; nIdx++ )
	{
		pfInitN[nIdx] = 1.0f;
		pfInitMeanIp[nIdx] = pfImageY[nIdx]*m_pfTransmission[nIdx];
	}

	BoxFilter(pfInitN, m_nGBlockSize, nW, nH, pfN);
	BoxFilter(m_pfTransmission, m_nGBlockSize, nW, nH, pfMeanP);
	
	BoxFilter(pfImageY, m_nGBlockSize, nW, nH, pfMeanI);

	BoxFilter(pfInitMeanIp, m_nGBlockSize, nW, nH, pfMeanIp);

	// Variance of (I, pfTrans) in each local patch
	for (nIdx = 0; nIdx < nW*nH; nIdx++ )
	{
		pfInitvarI[nIdx] = pfImageY[nIdx]*pfImageY[nIdx];
	}
	BoxFilter(pfInitvarI, m_nGBlockSize, nW, nH, pfvarI);

	for (nIdx = 0; nIdx < nW*nH; nIdx++ )
	{
		pfMeanI[nIdx] = pfMeanI[nIdx]/pfN[nIdx];
		pfMeanP[nIdx] = pfMeanP[nIdx]/pfN[nIdx];
		pfMeanIp[nIdx] = pfMeanIp[nIdx]/pfN[nIdx];
		pfCovIp[nIdx] = pfMeanIp[nIdx] - pfMeanI[nIdx] * pfMeanP[nIdx];
		pfvarI[nIdx] = pfvarI[nIdx]/pfN[nIdx] - pfMeanI[nIdx]*pfMeanI[nIdx];
	}

	// Calculate coefficient a and coefficient b
	for(int nIdx =0; nIdx<nW*nH; nIdx++)
	{
		// coefficient a
		pfA[nIdx] = (pfCovIp[nIdx])/(pfvarI[nIdx]+fEps);
		// coefficient b
		pfB[nIdx] = pfMeanP[nIdx] - pfA[nIdx]*pfMeanI[nIdx];
	}

	// Transmission refinement at each pixel

	BoxFilter(pfA,m_nGBlockSize,nW,nH,pfOutA);
	BoxFilter(pfB,m_nGBlockSize,nW,nH,pfOutB);

	for ( int nIdx = 0; nIdx < nW*nH; nIdx++ )
	{
        assert(pfN[nIdx] != 0);
		m_pfTransmissionR[nIdx] = ( pfOutA[nIdx]*pfImageY[nIdx] +  pfOutB[nIdx] ) / pfN[nIdx];
        assert(!isnan(m_pfTransmissionR[nIdx]));
	}

	delete []pfInitN;
	delete []pfInitMeanIp;
	delete []pfMeanP;
	delete []pfN;
	delete []pfMeanI;
	delete []pfMeanIp;
	delete []pfCovIp;
	delete []pfInitvarI;
	delete []pfvarI;
	delete []pfA;
	delete []pfB;
	delete []pfOutA;
	delete []pfOutB;
	//////////////////////////////////////////////////////////////////////////

	delete []pfImageY;
}

double* newdouble(int size) {
    double* retv = new double[size];
    memset(retv, 0, size * sizeof(double));
    return retv;
}
int* newint(int size) {
    int* retv = new int[size];
    memset(retv, 0, size * sizeof(int));
    return retv;
}
void checkdouble(double* arr, int size) {
    for (int i = 0; i < size; i++) {
        assert(!isnan(arr[i]));
    }
}

/*
	Function: GuidedFilter
	Description: the original guided filter for rgb color image. This function is used for image dehazing.
		The video dehazing algorithm uses appoximated filter for fast refinement. 
	Parameter:
		nW - width of array
		nH - height of array
		fEps - epsilon
	(member variable)
		m_pfTransmission - initial transmission (block_based)
		m_pnYImg - guidance image (Y image)
	Return:
		m_pfTransmissionR - filtered transmission
 */

void dehazing::GuidedFilter(int nW, int nH, double fEps)
{
	double* pfImageR = newdouble(nW*nH);
	double* pfImageG = newdouble(nW*nH);
	double* pfImageB = newdouble(nW*nH);

	double* pfInitN = newdouble(nW*nH);
	double* pfInitMeanIpR = newdouble(nW*nH);
	double* pfInitMeanIpG = newdouble(nW*nH);
	double* pfInitMeanIpB = newdouble(nW*nH);
	double* pfMeanP = newdouble(nW*nH);

	double* pfN = newdouble(nW*nH);
	double* pfMeanIr = newdouble(nW*nH);
	double* pfMeanIg = newdouble(nW*nH);
	double* pfMeanIb = newdouble(nW*nH);
	double* pfMeanIpR = newdouble(nW*nH);
	double* pfMeanIpG = newdouble(nW*nH);
	double* pfMeanIpB = newdouble(nW*nH);
	double* pfCovIpR = newdouble(nW*nH);
	double* pfCovIpG = newdouble(nW*nH);
	double* pfCovIpB = newdouble(nW*nH);

	double* pfCovEntire = newdouble(nW*nH*3);

	double* pfInitVarIrr = newdouble(nW*nH);
	double* pfInitVarIrg = newdouble(nW*nH);
	double* pfInitVarIrb = newdouble(nW*nH);
	double* pfInitVarIgg = newdouble(nW*nH);
	double* pfInitVarIgb = newdouble(nW*nH);
	double* pfInitVarIbb = newdouble(nW*nH);

	double* pfVarIrr = newdouble(nW*nH);
	double* pfVarIrg = newdouble(nW*nH);
	double* pfVarIrb = newdouble(nW*nH);
	double* pfVarIgg = newdouble(nW*nH);
	double* pfVarIgb = newdouble(nW*nH);
	double* pfVarIbb = newdouble(nW*nH);

	double* pfA1 = newdouble(nW*nH);
	double* pfA2 = newdouble(nW*nH);
	double* pfA3 = newdouble(nW*nH);
	double* pfOutA1 = newdouble(nW*nH);
	double* pfOutA2 = newdouble(nW*nH);
	double* pfOutA3 = newdouble(nW*nH);

	double* pfSigmaEntire = newdouble(nW*nH*9);

	double* pfB = newdouble(nW*nH);
	double* pfOutB = newdouble(nW*nH);

	int nIdx;
	// Converting to double point
	for(nIdx = 0 ; nIdx < nW*nH ; nIdx++ )
	{
		pfImageR[nIdx] = (double)m_pnRImg[nIdx];
		pfImageG[nIdx] = (double)m_pnGImg[nIdx];
		pfImageB[nIdx] = (double)m_pnBImg[nIdx];
	}
	//////////////////////////////////////////////////////////////////////////

	// Make an integral image
	for (nIdx = 0; nIdx < nW*nH; nIdx++ )
	{
		pfInitN[nIdx] = (double)1;

		pfInitMeanIpR[nIdx] = pfImageR[nIdx]*m_pfTransmission[nIdx];
		pfInitMeanIpG[nIdx] = pfImageG[nIdx]*m_pfTransmission[nIdx];
		pfInitMeanIpB[nIdx] = pfImageB[nIdx]*m_pfTransmission[nIdx];
	}

	BoxFilter(pfInitN, m_nGBlockSize, nW, nH, pfN);
	BoxFilter(m_pfTransmission, m_nGBlockSize, nW, nH, pfMeanP);

	BoxFilter(pfImageR, pfImageG, pfImageB, m_nGBlockSize, nW, nH, pfMeanIr, pfMeanIg, pfMeanIb);
		
	BoxFilter(pfInitMeanIpR, pfInitMeanIpG, pfInitMeanIpB, m_nGBlockSize, nW, nH, pfMeanIpR, pfMeanIpG, pfMeanIpB);

	//Covariance of (I, pfTrans) in each local patch
	
	for (nIdx = 0; nIdx < nW*nH; nIdx++ )
	{
		pfMeanIr[nIdx] = pfMeanIr[nIdx]/pfN[nIdx];
		pfMeanIg[nIdx] = pfMeanIg[nIdx]/pfN[nIdx];
		pfMeanIb[nIdx] = pfMeanIb[nIdx]/pfN[nIdx];

		pfMeanP[nIdx] = pfMeanP[nIdx]/pfN[nIdx];

		pfMeanIpR[nIdx] = pfMeanIpR[nIdx]/pfN[nIdx];
		pfMeanIpG[nIdx] = pfMeanIpG[nIdx]/pfN[nIdx];
		pfMeanIpB[nIdx] = pfMeanIpB[nIdx]/pfN[nIdx];

		pfCovIpR[nIdx] = pfMeanIpR[nIdx] - pfMeanIr[nIdx]*pfMeanP[nIdx];
		pfCovIpG[nIdx] = pfMeanIpG[nIdx] - pfMeanIg[nIdx]*pfMeanP[nIdx];
		pfCovIpB[nIdx] = pfMeanIpB[nIdx] - pfMeanIb[nIdx]*pfMeanP[nIdx];

		pfCovEntire[nIdx*3] = pfCovIpR[nIdx];
		pfCovEntire[nIdx*3+1] = pfCovIpG[nIdx];
		pfCovEntire[nIdx*3+2] = pfCovIpB[nIdx];

		pfInitVarIrr[nIdx] = pfImageR[nIdx]*pfImageR[nIdx];
		pfInitVarIrg[nIdx] = pfImageR[nIdx]*pfImageG[nIdx];
		pfInitVarIrb[nIdx] = pfImageR[nIdx]*pfImageB[nIdx];
		pfInitVarIgg[nIdx] = pfImageG[nIdx]*pfImageG[nIdx];
		pfInitVarIgb[nIdx] = pfImageG[nIdx]*pfImageB[nIdx];
		pfInitVarIbb[nIdx] = pfImageB[nIdx]*pfImageB[nIdx];
	}

	// Variance of I in each local patch: the matrix Sigma.
	// 		    rr, rg, rb
	// pfSigma  rg, gg, gb
	//	 	    rb, gb, bb

	BoxFilter(pfInitVarIrr, pfInitVarIrg, pfInitVarIrb, m_nGBlockSize, nW, nH, pfVarIrr, pfVarIrg, pfVarIrb);
	BoxFilter(pfInitVarIgg, pfInitVarIgb, pfInitVarIbb, m_nGBlockSize, nW, nH, pfVarIgg, pfVarIgb, pfVarIbb);

	for (nIdx = 0; nIdx < nW*nH; nIdx++ )
	{
		pfVarIrr[nIdx] = pfVarIrr[nIdx]/pfN[nIdx] - pfMeanIr[nIdx]*pfMeanIr[nIdx];
		pfVarIrg[nIdx] = pfVarIrg[nIdx]/pfN[nIdx] - pfMeanIr[nIdx]*pfMeanIg[nIdx];
		pfVarIrb[nIdx] = pfVarIrb[nIdx]/pfN[nIdx] - pfMeanIr[nIdx]*pfMeanIb[nIdx];
		pfVarIgg[nIdx] = pfVarIgg[nIdx]/pfN[nIdx] - pfMeanIg[nIdx]*pfMeanIg[nIdx];
		pfVarIgb[nIdx] = pfVarIgb[nIdx]/pfN[nIdx] - pfMeanIg[nIdx]*pfMeanIb[nIdx];
		pfVarIbb[nIdx] = pfVarIbb[nIdx]/pfN[nIdx] - pfMeanIb[nIdx]*pfMeanIb[nIdx];

		pfSigmaEntire[nIdx*9+0] = pfVarIrr[nIdx]+fEps*2.0f;
		pfSigmaEntire[nIdx*9+1] = pfVarIrg[nIdx];
		pfSigmaEntire[nIdx*9+2] = pfVarIrb[nIdx];
		pfSigmaEntire[nIdx*9+3] = pfVarIrg[nIdx];
		pfSigmaEntire[nIdx*9+4] = pfVarIgg[nIdx]+fEps*2.0f;
		pfSigmaEntire[nIdx*9+5] = pfVarIgb[nIdx];
		pfSigmaEntire[nIdx*9+6] = pfVarIrb[nIdx];
		pfSigmaEntire[nIdx*9+7] = pfVarIgb[nIdx];
		pfSigmaEntire[nIdx*9+8] = pfVarIbb[nIdx]+fEps*2.0f;
	}
	// Coefficient a와 coefficient b계산
	// Coefficienta
	for(nIdx =0; nIdx<nW*nH; nIdx++)
	{
		CalcAcoeff(pfSigmaEntire, pfCovEntire, pfA1, pfA2, pfA3, nIdx);
	}
    checkdouble(pfA1, nW * nH);

    checkdouble(pfMeanP, nW * nH);
    checkdouble(pfA1, nW * nH);
    checkdouble(pfA2, nW * nH);
    checkdouble(pfA3, nW * nH);
    checkdouble(pfMeanIr, nW * nH);
    checkdouble(pfMeanIg, nW * nH);
    checkdouble(pfMeanIb, nW * nH);

	// Coefficient b
	for (nIdx = 0; nIdx < nW*nH; nIdx++ )
	{
		pfB[nIdx] = pfMeanP[nIdx] - pfA1[nIdx]*pfMeanIr[nIdx] - pfA2[nIdx]*pfMeanIg[nIdx] - pfA3[nIdx]*pfMeanIb[nIdx];
	}
	
	// Transmission refinement at each pixel
	BoxFilter(pfA1,pfA2,pfA3,m_nGBlockSize,nW,nH,pfOutA1,pfOutA2,pfOutA3);

    for (int i = 0; i < nW * nH; i++) {
        assert(!isnan(pfB[i]));
    }
	BoxFilter(pfB,m_nGBlockSize,nW,nH,pfOutB);
    for (int i = 0; i < nW * nH; i++) {
        assert(!isnan(pfOutB[i]));
    }

	for (nIdx = 0; nIdx < nW*nH; nIdx++ )
	{
        assert(pfN[nIdx] != 0);

        double a1 = pfOutA1[nIdx] * pfImageR[nIdx];
        double a2 = pfOutA2[nIdx] * pfImageG[nIdx];
        double a3 = pfOutA3[nIdx] * pfImageB[nIdx];
        double a4 = pfOutB[nIdx];

        double xx = a1 + a2 + a3 + a4;
        double yy = pfN[nIdx];
        m_pfTransmissionR[nIdx] = xx / yy;
        assert(!isnan(m_pfTransmissionR[nIdx]));
    }

	delete []pfInitN;
	delete []pfInitMeanIpR;
	delete []pfInitMeanIpG;
	delete []pfInitMeanIpB;
	delete []pfMeanP;

	delete []pfN;
	delete []pfMeanIr;
	delete []pfMeanIg;
	delete []pfMeanIb;
	delete []pfMeanIpR;
	delete []pfMeanIpG;
	delete []pfMeanIpB;
	delete []pfCovIpR;
	delete []pfCovIpG;
	delete []pfCovIpB;
	delete []pfCovEntire;
	delete []pfInitVarIrr;
	delete []pfInitVarIrg;
	delete []pfInitVarIrb;
	delete []pfInitVarIgg;
	delete []pfInitVarIgb;
	delete []pfInitVarIbb;
	delete []pfVarIrr;
	delete []pfVarIrg;
	delete []pfVarIrb;
	delete []pfVarIgg;
	delete []pfVarIgb;
	delete []pfVarIbb;
	delete []pfA1;
	delete []pfA2;
	delete []pfA3;
	delete []pfOutA1;
	delete []pfOutA2;
	delete []pfOutA3;
	delete []pfSigmaEntire;
	delete []pfB;
	delete []pfOutB;
	//////////////////////////////////////////////////////////////////////////

	delete []pfImageR;
	delete []pfImageG;
	delete []pfImageB;
}

/*
	Function: FastGuidedFilter (downsampled image)
	Description: Transmission refinement based on guided fitlering, but we approximate the filtering 
		using partial window. In addtion, we analyze the guided filter using projection method.
		SSE (SIMD) is applied.
		* Limitation - In this code, we do not consider the border line of filtering, hence we will update the code.
					   The user may regulate the window size that the image size divides into the window size.
	Parameters: 
		nW - width of down-sampled image
		nH - height of down-sampled image	
	(hidden)
		m_pnSmallYImg - down-sampled image
		m_pnSmallTrans - down-sampled initial transmission
	Return:
		m_pfSmallTransR - down-sampled refined transmission
 */
void dehazing::FastGuidedFilterS()
{
	int nW = 320;
	int nH = 240;
	int nStep = m_nStepSize;				
	int nWstep = m_nGBlockSize/nStep;		// Step size of x-axis
	int nHstep = m_nGBlockSize/nStep;		// Step size of y-axis

	int nYb, nXb, nYa, nXa;
	int nIdxA, nIdxB, nI;
	double fDenom, fNormV1, fMeanI;		
	double fAk, fBk;						

	// SIMD is applied
	__m128d sseMean, sseY, ssePk_p, sseDenom, sseInteg, sseNormPk,	
		sseAk, sseBk, sseP, sseNormV1, sseGauss, sseBkV1, sseQ;

	// Initialization
	for ( nYb = 0 ; nYb < nH; nYb++ )
	{
		for ( nXb = 0; nXb < nW; nXb++ )
		{
			m_pfSmallY[nYb*nW+nXb]=(double)m_pnSmallYImg[nYb*nW+nXb];
			m_pfSmallInteg[nYb*nW+nXb]=0;
			m_pfSmallDenom[nYb*nW+nXb]=0;
		}
	}
	// Transmission refinement is applied to sampling pixels
	for ( nYa = 0 ; nYa < (nH-m_nGBlockSize)/nHstep+1 ; nYa++ )
	{
		for ( nXa = 0 ; nXa < (nW-m_nGBlockSize)/nWstep+1 ; nXa++ )
		{
			nIdxA=nYa*nHstep*nW+nXa*nWstep;

			// vector normalize -> fNormV1
			fNormV1 = 1.0f/m_nGBlockSize;
			fMeanI = 0.0f;

			// fNormV1 --> create perpendicular vector(m_pfSmallPk_p)
			// m_pfSmallPk_p == m_pfSmallY - fMeanI;

			sseMean = _mm_setzero_pd();
			sseY = _mm_setzero_pd();

			// Mean value
			for ( nYb = 0 ; nYb < m_nGBlockSize; nYb++ )
			{
				for ( nXb = 0; nXb < m_nGBlockSize; nXb+=4 )
				{
					sseY = _mm_loadu_pd((m_pfSmallY+(nIdxA+nYb*nW+nXb)));
					sseMean = _mm_add_pd( sseMean, sseY );
				}
			}

			for ( nI = 0 ; nI < 4; nI++ )
			{
				fMeanI += *((double*)(&sseMean)+nI);
			}

			fMeanI = fMeanI/(m_nGBlockSize*m_nGBlockSize);

			sseMean = _mm_set1_pd(fMeanI);
			sseY = _mm_setzero_pd();

			ssePk_p;
			
			// m_pfSmallY - fMeanI
			for ( nYb = 0 ; nYb < m_nGBlockSize; nYb++ )
			{
				for ( nXb = 0; nXb < m_nGBlockSize; nXb+=4 )
				{
					sseY = _mm_loadu_pd(m_pfSmallY+(nIdxA+nYb*nW+nXb));
					ssePk_p = _mm_sub_pd( sseY, sseMean );

					for ( nI = 0 ; nI < 4; nI++)
					{
						m_pfSmallPk_p[nYb*m_nGBlockSize+nXb+nI] = *((double*)(&ssePk_p)+nI);
					}
				}
			}

			// p vector normalize -> m_pfSmallNormPk
			// m_pfSmallPk_p^2
			fDenom=0.0f;

			sseDenom =_mm_setzero_pd();

			for ( nIdxB = 0 ; nIdxB < m_nGBlockSize*m_nGBlockSize; nIdxB+=4 )
			{
				ssePk_p = _mm_loadu_pd(m_pfSmallPk_p+nIdxB);
				sseDenom = _mm_add_pd( sseDenom, _mm_mul_pd( ssePk_p, ssePk_p ));
			}

			for ( nI = 0 ; nI < 4; nI++ )
			{
				fDenom += *((double*)(&sseDenom)+nI);
			}

			sseDenom = _mm_set1_pd(sqrt(fDenom));

			sseNormPk;

			fAk = 0.0f;
			fBk = 0.0f;

			sseAk = _mm_setzero_pd();
			sseBk = _mm_setzero_pd();
			sseP;
			
			
			// a = q'norm_p
			// b = q'norm_v1
			// Exception handling (denominator == 0)
			if(fDenom == 0)
			{
				// a = 0;
				fAk = 0;

				sseNormV1 = _mm_set1_pd(fNormV1);
				for ( nYb = 0 ; nYb < m_nGBlockSize; nYb++ )
				{
					for ( nXb = 0; nXb < m_nGBlockSize; nXb+=4 )
					{
						sseP = _mm_loadu_pd(m_pfSmallTrans+(nIdxA+nYb*nW+nXb));
						sseBk = _mm_add_pd( _mm_mul_pd( sseP , sseNormV1 ), sseBk );
					}
				}

				// b = q'norm_v1
				for ( nI = 0 ; nI < 4; nI++ )
					fBk += *((double*)(&sseBk)+nI);
			}
			else
			{
				// m_pfSmallNormPk
				for ( nIdxB = 0 ; nIdxB < m_nGBlockSize*m_nGBlockSize; nIdxB+=4 )
				{
					ssePk_p = _mm_loadu_pd(m_pfSmallPk_p+nIdxB);
					sseNormPk = _mm_div_pd(ssePk_p, sseDenom);

					for ( nI = 0 ; nI < 4; nI++ )
					{
						m_pfSmallNormPk[nIdxB+nI] = *((double*)(&sseNormPk)+nI);
					}
				}

				// a = q'norm_p
				// b = q'norm_v1
				sseNormV1 = _mm_set1_pd(fNormV1);

				for ( nYb = 0 ; nYb < m_nGBlockSize; nYb++ )
				{
					for ( nXb = 0; nXb < m_nGBlockSize; nXb+=4 )
					{
						sseP = _mm_loadu_pd(m_pfSmallTrans+(nIdxA+nYb*nW+nXb));
						sseNormPk = _mm_loadu_pd(m_pfSmallNormPk+(nYb*m_nGBlockSize+nXb));

						sseAk = _mm_add_pd( _mm_mul_pd( sseP , sseNormPk ), sseAk );
						sseBk = _mm_add_pd( _mm_mul_pd( sseP , sseNormV1 ), sseBk );
					}
				}

				for ( nI = 0 ; nI < 4; nI++ )
				{
					fAk += *((double*)(&sseAk)+nI);
					fBk += *((double*)(&sseBk)+nI);
				}
			}

			sseGauss;
			sseInteg;
			sseBkV1 = _mm_set1_pd(fBk*fNormV1);
			sseAk = _mm_set1_pd(fAk);
			// Gaussian weighting
			// Weighted_denom

			for ( nYb = 0 ; nYb < m_nGBlockSize; nYb++ )
			{
				for ( nXb = 0; nXb < m_nGBlockSize; nXb+=4 )
				{ 
					sseNormPk = _mm_loadu_pd(m_pfSmallNormPk+(nYb*m_nGBlockSize+nXb));
					sseGauss = _mm_loadu_pd(m_pfGuidedLUT+(nYb*m_nGBlockSize+nXb));

					sseInteg = _mm_mul_pd(_mm_add_pd( _mm_mul_pd(sseNormPk, sseAk), sseBkV1 ), sseGauss);
					// m_pfSmallInteg
					// m_pfSmallDenom
					for ( nI = 0 ; nI < 4; nI++ )
					{
						m_pfSmallInteg[nIdxA+nYb*nW+nXb+nI] += *((double*)(&sseInteg)+nI);
						m_pfSmallDenom[nIdxA+nYb*nW+nXb+nI] += *((double*)(&sseGauss)+nI);
					}
				}
			}
		}
	}

	// m_pfSmallTransR = m_pfSmallInteg / m_pfSmallDenom
	for( nIdxA = 0 ; nIdxA < nW*nH; nIdxA+=4 )
	{
		sseInteg = _mm_loadu_pd(m_pfSmallInteg+nIdxA);
		sseDenom = _mm_loadu_pd(m_pfSmallDenom+nIdxA);

		sseQ = _mm_div_pd(sseInteg, sseDenom);

		for( nI = 0; nI < 4; nI++)
			m_pfSmallTransR[nIdxA+nI] = *((double*)(&sseQ)+nI);
	}
}

/*
	Function: FastGuidedFilter (original sized image)
	Description: Transmission refinement based on guided fitlering, but we approximate the filtering 
		using partial window. In addtion, we analyze the guided filter using projection method.
		SSE (SIMD) is applied.
		* Limitation - In this code, we do not consider the border line of filtering, hence we will update the code.
					   The user may regulate the window size that the image size divides into the window size.
	Parameters: 
		nW - width of down-sampled image
		nH - height of down-sampled image	
	(hidden)
		m_pnYImg - down-sampled image
		m_pnTransmission - down-sampled initial transmission
	Return:
		m_pfTransmissionR - down-sampled refined transmission
 */
void dehazing::FastGuidedFilter()
{
	int nStep = m_nStepSize;					// step size

	int nEPBlockSizeL = m_nGBlockSize;			
	int nHeiL = m_nHei;							
	int nWidL = m_nWid;							

	double* pfYL = m_pfY;						
	double* pfIntegL = m_pfInteg;
	double* pfDenomL = m_pfDenom;
	int* pnYImgL = m_pnYImg;
	double* pfPk_pL = m_pfPk_p;
	double* pfNormPkL = m_pfNormPk;
	double* pfGuidedLUTL = m_pfGuidedLUT;
	double* pfTransmissionL = m_pfTransmission;

	int nWstep = nEPBlockSizeL/nStep;			
	int nHstep = nEPBlockSizeL/nStep;			

	double fDenom, fNormV1, fMeanI;				
	double fAk = 0.0f;							// fAk : a in "a*I+b", fBk : b in "a*I+b"
	double fBk = 0.0f;
	int nIdxA, nYa, nXa, nYb, nXb, nI, nIdxB;	//  indexing
	
	// SIMD is applied
	__m128d sseMean, sseY, ssePk_p, sseDenom, sseInteg, sseNormPk, 
		sseAk, sseBk, sseP, sseNormV1, sseGauss, sseBkV1,	sseQ;
	
	// Initialization
	for ( nYb = 0 ; nYb < nHeiL; nYb++ )
	{
		for ( nXb = 0; nXb < nWidL; nXb++ )
		{
			pfYL[nYb*nWidL+nXb]=(double)pnYImgL[nYb*nWidL+nXb];
			pfIntegL[nYb*nWidL+nXb]=0;
			pfDenomL[nYb*nWidL+nXb]=0;
		}
	}
	// sampling points
	for ( nYa = 0 ; nYa < (nHeiL-nEPBlockSizeL)/nHstep+1 ; nYa++ )
	{
		for ( nXa = 0 ; nXa < (nWidL-nEPBlockSizeL)/nWstep+1 ; nXa++ )
		{
			nIdxA=nYa*nHstep*nWidL+nXa*nWstep;
			// 1 vector normalize -> fNormV1
			fNormV1 = 1.0f/nEPBlockSizeL;
			fMeanI = 0.0f;

			// Create p vector(m_pfSmallPk_p) which is perpendicular to fNormV1
			// p vector == m_pfSmallY - fMeanI
			sseMean = _mm_setzero_pd();
			sseY = _mm_setzero_pd();

			// Mean value
			for ( nYb = 0 ; nYb < nEPBlockSizeL; nYb++ )
			{
				for ( nXb = 0; nXb < nEPBlockSizeL; nXb+=4 )
				{
					sseY = _mm_loadu_pd((pfYL+(nIdxA+nYb*nWidL+nXb)));
					sseMean = _mm_add_pd( sseMean, sseY );
				}
			}

			for ( nI = 0 ; nI < 4; nI++ )
			{
				fMeanI += *((double*)(&sseMean)+nI);
			}

			fMeanI = fMeanI/(nEPBlockSizeL*nEPBlockSizeL);
			
			sseMean = _mm_set1_pd(fMeanI);
			sseY = _mm_setzero_pd();

			// m_pfSmallY - fMeanI
			for ( nYb = 0 ; nYb < nEPBlockSizeL; nYb++ )
			{
				for ( nXb = 0; nXb < nEPBlockSizeL; nXb+=4 )
				{
					sseY = _mm_loadu_pd(pfYL+(nIdxA+nYb*nWidL+nXb));
					ssePk_p = _mm_sub_pd( sseY, sseMean );

					for ( nI = 0 ; nI < 4; nI++)
					{
						pfPk_pL[nYb*nEPBlockSizeL+nXb+nI] = *((double*)(&ssePk_p)+nI);
					}
				}
			}
			// p vector normalize -> m_pfSmallNormPk
			// m_pfSmallPk_p^2
			fDenom=0.0f;

			sseDenom =_mm_setzero_pd();

			for ( nIdxB = 0 ; nIdxB < nEPBlockSizeL*nEPBlockSizeL; nIdxB+=4 )
			{
				ssePk_p = _mm_loadu_pd(pfPk_pL+nIdxB);
				sseDenom = _mm_add_pd( sseDenom, _mm_mul_pd( ssePk_p, ssePk_p ));
			}

			for ( nI = 0 ; nI < 4; nI++ )
			{
				fDenom += *((double*)(&sseDenom)+nI);
			}

			sseDenom = _mm_set1_pd(sqrt(fDenom));
			
			sseAk = _mm_setzero_pd();
			sseBk = _mm_setzero_pd();
			// a = q'norm_p
			// b = q'norm_v1
			// Exception handling (denominator == 0)
			// a = 0;
			fAk = 0.0f;
			fBk = 0.0f;

			if(fDenom == 0)
			{
				sseNormV1 = _mm_set1_pd(fNormV1);
				for ( nYb = 0 ; nYb < nEPBlockSizeL; nYb++ )
				{
					for ( nXb = 0; nXb < nEPBlockSizeL; nXb+=4 )
					{
						sseP = _mm_loadu_pd(pfTransmissionL+(nIdxA+nYb*nWidL+nXb));
						sseBk = _mm_add_pd( _mm_mul_pd( sseP , sseNormV1 ), sseBk );
					}
				}
				// b = q'norm_v1
				for ( nI = 0 ; nI < 4; nI++ )
					fBk += *((double*)(&sseBk)+nI);
			}
			else
			{
				// m_pfSmallNormPk
				for ( nIdxB = 0 ; nIdxB < nEPBlockSizeL*nEPBlockSizeL; nIdxB+=4 )
				{
					ssePk_p = _mm_loadu_pd(pfPk_pL+nIdxB);
					sseNormPk = _mm_div_pd(ssePk_p, sseDenom);

					for ( int nI = 0 ; nI < 4; nI++ )
					{
						pfNormPkL[nIdxB+nI] = *((double*)(&sseNormPk)+nI);
					}
				}

				sseNormV1 = _mm_set1_pd(fNormV1);
				// a = q'norm_p
				// b = q'norm_v1
				for ( nYb = 0 ; nYb < nEPBlockSizeL; nYb++ )
				{
					for ( nXb = 0; nXb < nEPBlockSizeL; nXb+=4 )
					{
						sseP = _mm_loadu_pd(pfTransmissionL+(nIdxA+nYb*nWidL+nXb));
						sseNormPk = _mm_loadu_pd(pfNormPkL+(nYb*nEPBlockSizeL+nXb));

						sseAk = _mm_add_pd( _mm_mul_pd( sseP , sseNormPk ), sseAk );
						sseBk = _mm_add_pd( _mm_mul_pd( sseP , sseNormV1 ), sseBk );
					}
				}

				for ( nI = 0 ; nI < 4; nI++ )
				{
					fAk += *((double*)(&sseAk)+nI);
					fBk += *((double*)(&sseBk)+nI);
				}
			}
					
			sseBkV1 = _mm_set1_pd(fBk*fNormV1);
			sseAk = _mm_set1_pd(fAk);
			// Gaussian weighting
			// Weighted_denom
			for ( nYb = 0 ; nYb < nEPBlockSizeL; nYb++ )
			{
				for ( nXb = 0; nXb < nEPBlockSizeL; nXb+=4 )
				{ 
					sseNormPk = _mm_loadu_pd(pfNormPkL+(nYb*nEPBlockSizeL+nXb));
					sseGauss = _mm_loadu_pd(pfGuidedLUTL+(nYb*nEPBlockSizeL+nXb));

					sseInteg = _mm_mul_pd(_mm_add_pd( _mm_mul_pd(sseNormPk, sseAk), sseBkV1 ), sseGauss);
					// m_pfSmallInteg
					// m_pfSmallDenom
					for ( nI = 0 ; nI < 4; nI++ )
					{
						pfIntegL[nIdxA+nYb*nWidL+nXb+nI] += *((double*)(&sseInteg)+nI);
						pfDenomL[nIdxA+nYb*nWidL+nXb+nI] += *((double*)(&sseGauss)+nI);
					}
				}
			}
		}
	}

	// m_pfSmallTransR = m_pfSmallInteg / m_pfSmallDenom
	for( int nIdx = 0 ; nIdx < nWidL*nHeiL; nIdx+=4 )
	{
		sseInteg = _mm_loadu_pd(pfIntegL+nIdx);
		sseDenom = _mm_loadu_pd(pfDenomL+nIdx);

		sseQ = _mm_div_pd(sseInteg, sseDenom);

		for( int nI = 0; nI < 4; nI++)
			m_pfTransmissionR[nIdx+nI] = *((double*)(&sseQ)+nI);
	}
}

/*
	Function: GuidedFilterShiftableWindow
	Description: It is a modified guided filter to reduce blur artifacts of original 
		one. The algorithm shift the window which has the minimum error.
		The other scheme is the same with original one.
	Parameters: 
		fEps - the epsilon (the same with original guided filter).
	Return:
		m_pfSmallTransR - refined transmission.		
 */
void dehazing::GuidedFilterShiftableWindow(double fEps)
{
	int nY, nX;
	int nH, nW;
	int nYmin, nXmin, nYmax, nXmax;
	int nYminOpt, nXminOpt, nYmaxOpt, nXmaxOpt;

	double fMeanI, fMeanR, fMeanG, fMeanB;
	double fActualBlock;
	double fMeant;
	double ftrans;
	double fa1, fa2, fa3, fb;
	double fCovRp, fCovGp, fCovBp;
	double afSigma[9], afInv[9];
	double fVRR, fVRG, fVRB, fVGB, fVBB, fVGG;
	double fMRt, fMGt, fMBt;
	double fR, fG, fB;
	double fDet;
	int x1, x2, y1, y2;
	int nItersX, nItersY;
	double fVariance, fOptVariance;

	double *pfImageY		= newdouble(m_nWid*m_nHei);
	double *pfSqImageY	= newdouble(m_nWid*m_nHei);
	double *pfSumTrans	= newdouble(m_nWid*m_nHei);
	
	double *pfIntImageY	= newdouble(m_nWid*m_nHei);
	double *pfSqIntImageY= newdouble(m_nWid*m_nHei);
	double *pfintN		= newdouble(m_nWid*m_nHei);
	double *pfones		= newdouble(m_nWid*m_nHei);

	double *pfWeight		= newdouble(m_nWid*m_nHei);
	double *pfWeiSum		= newdouble(m_nWid*m_nHei);
		
	double *pfVarImage	= newdouble(m_nWid*m_nHei);

	int	*pnN			= newint (m_nWid*m_nHei);
	int *pnVarN			= newint (m_nWid*m_nHei);
		
	int nMax = 0;
	double fMaxvar=0;
	double grey;
	bool bIterationFlag;

	double fMB, fMG, fMR, fMBS, fMGS, fMRS;

	int	nSizes = 31;
	int nNewX, nNewY, nOptNewX, nOptNewY;
	double inverseweight;

	bIterationFlag = TRUE;

	// initialize
	for(nX=0; nX<m_nWid*m_nHei; nX++)
	{
		pfSumTrans[nX] = 0;
		pfintN[nX] = 0;
		pfones[nX] = 1.0f;
		pnVarN[nX] = 0;
		pfWeiSum[nX] = 0;

		pfImageY[nX] = (double)m_pnYImg[nX];
		pfSqImageY[nX] = (double)m_pnYImg[nX]*(double)m_pnYImg[nX];
	}

	BoxFilter(pfImageY, 15, m_nWid, m_nHei, pfIntImageY);
	BoxFilter(pfSqImageY, 15, m_nWid, m_nHei, pfSqIntImageY);
	BoxFilter(pfones, 15, m_nWid, m_nHei, pfintN);

	for(nY=0; nY<m_nHei; nY++)
	{
		for(nX=0; nX<m_nWid; nX++)
		{
			pfVarImage[nX+nY*m_nWid] = pfSqIntImageY[nX+nY*m_nWid]/pfintN[nX+nY*m_nWid] - pfIntImageY[nX+nY*m_nWid]/pfintN[nX+nY*m_nWid]*pfIntImageY[nX+nY*m_nWid]/pfintN[nX+nY*m_nWid];
		}
	}

	int totalsum = 0;

iterative:
	
	for(nY=0; nY<m_nHei; nY++)
	{
		for(nX=0; nX<m_nWid; nX++)
		{
			//find minimum
			fOptVariance = 99999999999999;
			
			for(nItersY=-nSizes/2; nItersY<nSizes/2+1; nItersY+=2)
			{
				for(nItersX=-nSizes/2; nItersX<nSizes/2+1; nItersX+=2)
				{
					nYmin = __max(0, (nY+nItersY) - nSizes/2);
					nYmax = __min(m_nHei, (nY+nItersY) + nSizes/2+1);
					nXmin = __max(0, (nX+nItersX) - nSizes/2);
					nXmax = __min(m_nWid, (nX+nItersX) + (nSizes)/2+1);
				
					nNewX = __min(__max(nX+nItersX, 0), m_nWid-1);
					nNewY = __min(__max(nY+nItersY, 0), m_nHei-1);
				
					inverseweight = 1;//exp(-((double)(nItersX*nItersX)+(double)(nItersY*nItersY))/70.0f);

					if(fOptVariance > pfVarImage[nNewX+nNewY*m_nWid]*inverseweight)
					{
						fOptVariance = pfVarImage[nNewX+nNewY*m_nWid]*inverseweight;
						nYminOpt = nYmin;
						nXminOpt = nXmin;
						nYmaxOpt = nYmax;
						nXmaxOpt = nXmax;
						nOptNewX = nNewX;
						nOptNewY = nNewY;
					}
				}
			}

			fActualBlock = (double)(nYmaxOpt-nYminOpt)*(nXmaxOpt-nXminOpt);
			if(bIterationFlag)
			{
				pnVarN[nOptNewX+nOptNewY*m_nWid]++; 
				totalsum++;
			}

			fVRR=0; fVRG=0; fVRB=0; fVGB=0; fVBB=0; fVGG=0;
			fMeanR=0; fMeanG=0; fMeanB=0; fMeant=0;
			fMRt=0; fMGt=0; fMBt=0;

			// guided filtering at optimal block
			for(nH=nYminOpt; nH<nYmaxOpt; nH++)
			{
				for(nW=nXminOpt; nW<nXmaxOpt; nW++)
				{
					fB = ((double)(m_pnBImg[nW+m_nWid*nH]))/255.0f;
					fG = ((double)(m_pnGImg[nW+m_nWid*nH]))/255.0f;
					fR = ((double)(m_pnRImg[nW+m_nWid*nH]))/255.0f;
					
					fMeanB += fB;
					fMeanG += fG;
					fMeanR += fR;
					ftrans = CLIP_TRS(m_pfTransmission[nW+m_nWid*nH]);

					fMeant += ftrans;

					fMBt += fB*ftrans;
					fMGt += fG*ftrans;
					fMRt += fR*ftrans;

					fVRR += fR*fR;
					fVRG += fR*fG;
					fVRB += fR*fB;
					fVGG += fG*fG;
					fVGB += fG*fB;
					fVBB += fB*fB;
				}
			}

			fMeanB /= fActualBlock;
			fMeanG /= fActualBlock;
			fMeanR /= fActualBlock;
			fMeant /= fActualBlock;
			
			fMRt /= fActualBlock;
			fMGt /= fActualBlock;
			fMBt /= fActualBlock;

			fCovRp = fMRt - fMeanR*fMeant;
			fCovGp = fMGt - fMeanG*fMeant;
			fCovBp = fMBt - fMeanB*fMeant;

			afSigma[0] = fVRR/fActualBlock - fMeanR*fMeanR + fEps*2.0f;
			afSigma[1] = fVRG/fActualBlock - fMeanR*fMeanG;
			afSigma[2] = fVRB/fActualBlock - fMeanR*fMeanB;
			afSigma[3] = fVRG/fActualBlock - fMeanR*fMeanG;
			afSigma[4] = fVGG/fActualBlock - fMeanG*fMeanG + fEps*2.0f;
			afSigma[5] = fVGB/fActualBlock - fMeanG*fMeanB;
			afSigma[6] = fVRB/fActualBlock - fMeanR*fMeanB;
			afSigma[7] = fVGB/fActualBlock - fMeanG*fMeanB;
			afSigma[8] = fVBB/fActualBlock - fMeanB*fMeanB + fEps*2.0f;
			
			//3x3 inverse
			fDet = 1.0f/( (afSigma[0]*(afSigma[4] * afSigma[8] - afSigma[5] * afSigma[7]))
				-	(	afSigma[1]*(afSigma[3] * afSigma[8] - afSigma[5] * afSigma[6]))
				+	(	afSigma[2]*(afSigma[3] * afSigma[7] - afSigma[4] * afSigma[6])) );
			for(nH=0; nH<3; nH++)
			{
				for(nW=0; nW<3; nW++)
				{
					x1 = nH == 0? 1: 0;
					x2 = nH == 2? 1: 2;
					y1 = nW == 0? 1: 0;
					y2 = nW == 2? 1: 2;

					afInv[nH*3+nW] = (( afSigma[y1*3+x1]  *  afSigma[y2*3+x2] )
						-  ( afSigma[y1*3+x2]  *  afSigma[y2*3+x1] )) * fDet;
					if(1==((nW+nH)%2))
						afInv[nH*3+nW] = - afInv[nH*3+nW];
				}
			}
			//
			fa1 = fCovRp*afInv[0] + fCovGp*afInv[3] + fCovBp*afInv[6]; 
			fa2 = fCovRp*afInv[1] + fCovGp*afInv[4] + fCovBp*afInv[7]; 
			fa3 = fCovRp*afInv[2] + fCovGp*afInv[5] + fCovBp*afInv[8]; 

			fb = fMeant - fa1*fMeanR - fa2*fMeanG - fa3*fMeanB;
			if(bIterationFlag)
			{
				for(nH=nYminOpt; nH<nYmaxOpt; nH++)
				{
					for(nW=nXminOpt; nW<nXmaxOpt; nW++)
					{
						pnN[nH*m_nWid+nW]++;
					}
				}
			}
			else
			{
				for(nH=nYminOpt; nH<nYmaxOpt; nH++)
				{
					for(nW=nXminOpt; nW<nXmaxOpt; nW++)
					{
						pfSumTrans[nH*m_nWid+nW] = pfSumTrans[nH*m_nWid+nW] + 
							(fa1*(double)(m_pnRImg[nW+m_nWid*nH])/255.0f
							+ fa2*(double)(m_pnGImg[nW+m_nWid*nH])/255.0f	
							+ fa3*(double)(m_pnBImg[nW+m_nWid*nH])/255.0f + fb)*pfWeight[nOptNewY*m_nWid+nOptNewX];
						//pnN[nH*nWid+nW]++;
						pfWeiSum[nH*m_nWid+nW] += pfWeight[nOptNewY*m_nWid+nOptNewX];
					}
				}
			}
		}
	}
	
	nMax = 0;
	double fmax=0;
	for(nY=0; nY<m_nHei; nY++)
	{
		for(nX=0; nX<m_nWid; nX++)
		{
			if(bIterationFlag)
			{
				if(nMax < pnVarN[nY*m_nWid+nX])
					nMax = pnVarN[nY*m_nWid+nX];
				pfWeight[nX+nY*m_nWid] = (double)pnVarN[nX+nY*m_nWid]/(double)pnN[nX+nY*m_nWid];
			}
			else
			{
                assert(pfWeiSum[nY * m_nWid + nX] != 0);
				m_pfTransmissionR[nY*m_nWid+nX] = pfSumTrans[nY*m_nWid+nX]/pfWeiSum[nY*m_nWid+nX];
                assert(!isnan(m_pfTransmissionR[nY * m_nWid + nX]));
			}
		}
	}

	if(bIterationFlag)
	{
		bIterationFlag = FALSE;
		goto iterative;
	}

	delete [] pfImageY;
	delete [] pfSqImageY;
	delete [] pfSumTrans;
	
	delete [] pfIntImageY;
	delete [] pfSqIntImageY;
	delete [] pfintN;
	delete [] pfones;

	delete [] pfWeight;
	delete [] pfWeiSum;
		
	delete [] pfVarImage;

	delete [] pnN;
	delete [] pnVarN;
}