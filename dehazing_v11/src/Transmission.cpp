/*
	The transmission estimation algorithm is in here. 
	The detailed description of the algorithm is presented
	in "http://mcl.korea.ac.kr/projects/dehazing". See also 
	J.-H. Kim, W.-D. Jang, Y. Park, D.-H. Lee, J.-Y. Sim, C.-S. Kim, "Temporally
	coherent real-time video dehazing," in Proc. IEEE ICIP, 2012.
	

	The transmission is estimated by maximizing the contrast of restored
	image while minimizing the information loss, which denotes the trucated
	pixel value. 

	The transmission estimation is block based approach, and thus its performance
	depends on the size of block. 

	Last updated: 2013-02-07
	Author: Jin-Hwan, Kim.
 */
#include "dehazing.h"

/*
	Function: TransmissionEstimation
	Description: Estiamte the transmission in the frame(Color)
				 Specified size.
	Parameters:
		nFrame - frame no.
		nWid - frame width
		nHei - frame height.
	Return:
		m_pfTransmission
 */
void dehazing::TransmissionEstimationColor(int *pnImageR, int *pnImageG, int *pnImageB, double *pfTransmission,int *pnImageRP, int *pnImageGP, int *pnImageBP, double *pfTransmissionP,int nFrame, int nWid, int nHei)
{
	//int nX, nY, nXstep, nYstep;
	//double fTrans;
	
	if(m_bPreviousFlag == true&&nFrame>0)
	{
		for(int nY=0; nY<nHei; nY+=m_nTBlockSize)
		{
			for(int nX=0; nX<nWid; nX+=m_nTBlockSize)
			{
                double fTrans = NFTrsEstimationPColor(pnImageR, pnImageG, pnImageB, pnImageRP, pnImageGP, pnImageBP, pfTransmissionP, __max(nX, 0), __max(nY, 0), nWid, nHei);
				for(int nYstep=nY; nYstep<nY+m_nTBlockSize; nYstep++)
				{
					for(int nXstep=nX; nXstep<nX+m_nTBlockSize; nXstep++)
					{
						pfTransmission[nYstep*nWid+nXstep] = fTrans;
					}
				}
			}
		}
	}
	else
	{
		for(int nY=0; nY<nHei; nY+=m_nTBlockSize)
		{
			for(int nX=0; nX<nWid; nX+=m_nTBlockSize)
			{
                double fTrans = NFTrsEstimationColor(pnImageR, pnImageG, pnImageB, __max(nX, 0), __max(nY, 0), nWid, nHei);

				for(int nYstep=nY; nYstep<nY+m_nTBlockSize && nYstep < nHei; nYstep++)
				{
					for(int nXstep=nX; nXstep<nX+m_nTBlockSize && nXstep < nWid; nXstep++)
					{
                        int index = nYstep * nWid + nXstep;
                        int imax = nWid * nHei;
                        if (index < imax) {
                            pfTransmission[index] = fTrans;
                        }
                        else {
                            assert(false);
                        }
					}
				}
			}
		}
	}
}

/*
	Function: TransmissionEstimation
	Description: Estiamte the transmission in the frame
				 Specified size.
	Parameters:
		nFrame - frame no.
		nWid - frame width
		nHei - frame height.
	Return:
		m_pfTransmission
 */
void dehazing::TransmissionEstimation(int *pnImageY, double *pfTransmission, int *pnImageYP, double *pfTransmissionP, int nFrame, int nWid, int nHei)
{
	int nX, nY, nXstep, nYstep;
	double fTrans;

	if(m_bPreviousFlag == true&&nFrame>0)
	{
		for(nY=0; nY<nHei; nY+=m_nTBlockSize)
		{
			for(nX=0; nX<nWid; nX+=m_nTBlockSize)
			{
				fTrans = NFTrsEstimationP(pnImageY, pnImageYP, pfTransmissionP, __max(nX, 0), __max(nY, 0), nWid, nHei);
				for(nYstep=nY; nYstep<nY+m_nTBlockSize; nYstep++)
				{
					for(nXstep=nX; nXstep<nX+m_nTBlockSize; nXstep++)
					{
						pfTransmission[nYstep*nWid+nXstep] = fTrans;
					}
				}
			}
		}
	}
	else
	{
		for(nY=0; nY<nHei; nY+=m_nTBlockSize)
		{
			for(nX=0; nX<nWid; nX+=m_nTBlockSize)
			{
				fTrans = NFTrsEstimation(pnImageY, __max(nX, 0), __max(nY, 0), nWid, nHei);

				for(nYstep=nY; nYstep<nY+m_nTBlockSize; nYstep++)
				{
					for(nXstep=nX; nXstep<nX+m_nTBlockSize; nXstep++)
					{
						pfTransmission[nYstep*nWid+nXstep] = fTrans;
					}
				}
			}
		}
	}
}

/*
	Function: NFTrsEstimation
	Description: Estiamte the transmission in the block. 
		The algorithm use exhaustive searching method and its step size
		is sampled to 0.1

	Parameters:
		nStartx - top left point of a block
		nStarty - top left point of a block
		nWid - frame width
		nHei - frame height.
	Return:
		fOptTrs
 */
double dehazing::NFTrsEstimation(int *pnImageY, int nStartX, int nStartY, int nWid, int nHei)
{
	int nCounter;	
	int nX, nY;		
	int nEndX;
	int nEndY;

	int nOut;						// Restored image
	int nSquaredOut;				// Squared value of restored image
	int nSumofOuts;					// Sum of restored image
	int nSumofSquaredOuts;			// Sum of squared restored image
	double fTrans, fOptTrs;			// Transmission and optimal value
	int nTrans;						// Integer transformation 
	int nSumofSLoss;				// Sum of loss info
	double fCost, fMinCost, fMean;	 
	int nNumberofPixels, nLossCount;

	nEndX = __min(nStartX+m_nTBlockSize, nWid); // End point of the block
	nEndY = __min(nStartY+m_nTBlockSize, nHei); // End point of the block

	nNumberofPixels = (nEndY-nStartY)*(nEndX-nStartX);	

	fTrans = 0.3f;	// Init trans is started from 0.3
	nTrans = 427;	// Convert transmission to integer 

	for(nCounter=0; nCounter<7; nCounter++)
	{
		nSumofSLoss = 0;
		nLossCount = 0;
		nSumofSquaredOuts = 0;
		nSumofOuts = 0;
		for(nY=nStartY; nY<nEndY; nY++)
		{
			for(nX=nStartX; nX<nEndX; nX++)
			{
				nOut = ((pnImageY[nY*nWid+nX] - m_nAirlight)*nTrans + 128*m_nAirlight)>>7; // (I-A)/t + A --> ((I-A)*k*128 + A*128)/128
				nSquaredOut = nOut * nOut;

				if(nOut>255)
				{
					nSumofSLoss += (nOut - 255)*(nOut - 255);
					nLossCount++;
				}
				else if(nOut < 0)
				{
					nSumofSLoss += nSquaredOut;
					nLossCount++;
				}
				nSumofSquaredOuts += nSquaredOut;
				nSumofOuts += nOut;
			}
		}
		fMean = (double)(nSumofOuts)/(double)(nNumberofPixels);  
		fCost = m_fLambda1 * (double)nSumofSLoss/(double)(nNumberofPixels) 
			- ((double)nSumofSquaredOuts/(double)nNumberofPixels - fMean*fMean); 
		
		if(nCounter==0 || fMinCost > fCost)
		{
			fMinCost = fCost;
			fOptTrs = fTrans;
		}

		fTrans += 0.1f;
		nTrans = (int)(1.0f/fTrans*128.0f);
	}
	return fOptTrs; 
}

/*
	Function: NFTrsEstimation
	Description: Estiamte the transmission in the block. (COLOR)
		The algorithm use exhaustive searching method and its step size
		is sampled to 0.1

	Parameters:
		nStartx - top left point of a block
		nStarty - top left point of a block
		nWid - frame width
		nHei - frame height.
	Return:
		fOptTrs
 */
double dehazing::NFTrsEstimationColor(int *pnImageR, int *pnImageG, int *pnImageB, int nStartX, int nStartY, int nWid, int nHei)
{
	int nCounter;	
	int nX, nY;		
	int nEndX;
	int nEndY;

	int nOutR, nOutG, nOutB;		
	int nSquaredOut;				
	int nSumofOuts;					
	int nSumofSquaredOuts;			
	double fTrans, fOptTrs;			
	int nTrans;						
	int nSumofSLoss;				
	double fCost, fMinCost, fMean;	
	int nNumberofPixels, nLossCount;

	nEndX = __min(nStartX+m_nTBlockSize, nWid); 
	nEndY = __min(nStartY+m_nTBlockSize, nHei); 

	nNumberofPixels = (nEndY-nStartY)*(nEndX-nStartX) * 3;	

	fTrans = 0.3f;	
	nTrans = 427;

	for(nCounter=0; nCounter<7; nCounter++)
	{
		nSumofSLoss = 0;
		nLossCount = 0;
		nSumofSquaredOuts = 0;
		nSumofOuts = 0;
	
		for(nY=nStartY; nY<nEndY; nY++)
		{
			for(nX=nStartX; nX<nEndX; nX++)
			{
				
				nOutB = ((pnImageB[nY*nWid+nX] - m_anAirlight[0])*nTrans + 128*m_anAirlight[0])>>7; // (I-A)/t + A --> ((I-A)*k*128 + A*128)/128
				nOutG = ((pnImageG[nY*nWid+nX] - m_anAirlight[1])*nTrans + 128*m_anAirlight[1])>>7;
				nOutR = ((pnImageR[nY*nWid+nX] - m_anAirlight[2])*nTrans + 128*m_anAirlight[2])>>7;		

				if(nOutR>255)
				{
					nSumofSLoss += (nOutR - 255)*(nOutR - 255);
					nLossCount++;
				}
				else if(nOutR < 0)
				{
					nSumofSLoss += nOutR * nOutR;
					nLossCount++;
				}
				if(nOutG>255)
				{
					nSumofSLoss += (nOutG - 255)*(nOutG - 255);
					nLossCount++;
				}
				else if(nOutG < 0)
				{
					nSumofSLoss += nOutG * nOutG;
					nLossCount++;
				}
				if(nOutB>255)
				{
					nSumofSLoss += (nOutB - 255)*(nOutB - 255);
					nLossCount++;
				}
				else if(nOutB < 0)
				{
					nSumofSLoss += nOutB * nOutB;
					nLossCount++;
				}
				nSumofSquaredOuts += nOutB * nOutB + nOutR * nOutR + nOutG * nOutG;;
				nSumofOuts += nOutR + nOutG + nOutB;
			}
		}
		fMean = (double)(nSumofOuts)/(double)(nNumberofPixels);  
		fCost = m_fLambda1 * (double)nSumofSLoss/(double)(nNumberofPixels) 
			- ((double)nSumofSquaredOuts/(double)nNumberofPixels - fMean*fMean); 

		if(nCounter==0 || fMinCost > fCost)
		{
			fMinCost = fCost;
			fOptTrs = fTrans;
		}

		fTrans += 0.1f;
		nTrans = (int)(1.0f/fTrans*128.0f);
	}
	return fOptTrs; 
}

/*
	Function: NFTrsEstimationP
	Description: Estiamte the transmission in the block. 
		The algorithm use exhaustive searching method and its step size
		is sampled to 0.1.
		The previous frame information is used to estimate transmission.

	Parameters:
		nStartx - top left point of a block
		nStarty - top left point of a block
		nWid - frame width
		nHei - frame height.
	Return:
		fOptTrs
 */
double dehazing::NFTrsEstimationP(int *pnImageY, int *pnImageYP, double *pfTransmissionP, int nStartX, int nStartY, int nWid, int nHei)
{
	int nCounter;	// for find out transmission 0.1~1.0, 10 iteration 
	int nX, nY;		// variable for index
	int nEndX;
	int nEndY;

	double fMean;

	int nOut;								
	double fPreTrs;							
	int nSquaredOut;						
	int nSumofOuts;							
	int nSumofSquaredOuts;					
	int nTrans;								
	int nSumofSLoss;						
	double fCost, fMinCost, fTrans, fOptTrs;	
	int nNumberofPixels;					

	nEndX = __min(nStartX+m_nTBlockSize, nWid); 
	nEndY = __min(nStartY+m_nTBlockSize, nHei); 

	nNumberofPixels = (nEndY-nStartY)*(nEndX-nStartX);	

	fTrans = 0.3f;	// initial transmission value
	nTrans = 427;	
	fPreTrs = 0;

	double fNewKSum = 0;						// Sum of new kappa which is multiplied the weight
	double fNewK;							// New kappa
	double fWi;								// Weight 
	double fPreJ;							// evade 0 division
	double fWsum = 0;						// Sum of weight
	int nIdx = 0;	
	int nLossCount;
	
	for(nY=nStartY; nY<nEndY; nY++)
	{
		for(nX=nStartX; nX<nEndX; nX++)
		{
			fPreJ = (double)(pnImageYP[nY*nWid+nX]-m_nAirlight);
			if(fPreJ != 0){
				fWi = m_pfExpLUT[abs(pnImageY[nY*nWid+nX]-pnImageYP[nY*nWid+nX])];
				fWsum += fWi;	
				fNewKSum += fWi*(double)(pnImageY[nY*nWid+nX]-m_nAirlight)/fPreJ;
			}
		}
	}
	fNewK = fNewKSum/fWsum;			// Compute new kappa
	fPreTrs = pfTransmissionP[nStartY*nWid+nStartX]*fNewK;	// Update the previous transmission using new kappa

	for(nCounter=0; nCounter<7; nCounter++)
	{
		nSumofSLoss = 0;
		nLossCount = 0;
		nSumofSquaredOuts = 0;
		nSumofOuts = 0;

		for(nY=nStartY; nY<nEndY; nY++)
		{
			for(nX=nStartX; nX<nEndX; nX++)
			{
				nOut = ((pnImageY[nY*nWid+nX] - m_nAirlight)*nTrans + 128*m_nAirlight)>>7; // (I-A)/t + A --> ((I-A)*k*128 + A*128)/128
				nSquaredOut = nOut * nOut;	
				
				if(nOut>255)
				{
					nSumofSLoss += (nOut - 255)*(nOut - 255);
					nLossCount++;
				}
				else if(nOut < 0)
				{
					nSumofSLoss += nSquaredOut;
					nLossCount++;
				}
				nSumofSquaredOuts += nSquaredOut;
				nSumofOuts += nOut;
			}
		}
		fMean = (double)(nSumofOuts)/(double)(nNumberofPixels);
		fCost = m_fLambda1 * (double)nSumofSLoss/(double)(nNumberofPixels) // information loss cost
			- ((double)nSumofSquaredOuts/(double)nNumberofPixels - fMean*fMean)	// contrast cost
			+ m_fLambda2/fPreTrs/fPreTrs*fWsum/(double)nNumberofPixels*((fPreTrs-fTrans)*(fPreTrs-fTrans)*255.0f*255.0f);// fPreTrs/fPreTrs*fWsum/(double)nNumberofPixels

		if(nCounter==0 || fMinCost > fCost)
		{
			fMinCost = fCost;
			fOptTrs = fTrans;
		}

		fTrans += 0.1f;
		nTrans = (int)(1/fTrans*128.0f);
	}
	return fOptTrs;
}
/*
	Function: NFTrsEstimationP(COLOR)
	Description: Estiamte the transmission in the block. 
		The algorithm use exhaustive searching method and its step size
		is sampled to 0.1.
		The previous frame information is used to estimate transmission.

	Parameters:
		nStartx - top left point of a block
		nStarty - top left point of a block
		nWid - frame width
		nHei - frame height.
	Return:
		fOptTrs
 */
double dehazing::NFTrsEstimationPColor(int *pnImageR, int *pnImageG, int *pnImageB, int *pnImageRP, int *pnImageGP, int *pnImageBP, double *pfTransmissionP, int nStartX, int nStartY, int nWid, int nHei)
{
	int nCounter;	
	int nX, nY;		
	int nEndX;
	int nEndY;

	double fMean;

	int nOutR, nOutG, nOutB;				
	double fPreTrs;							
	int nSquaredOut;						
	int nSumofOuts;							
	int nSumofSquaredOuts;					
	int nTrans;								
	int nSumofSLoss;						
	double fCost, fMinCost, fTrans, fOptTrs;	
	int nNumberofPixels;					

	nEndX = __min(nStartX+m_nTBlockSize, nWid); 
	nEndY = __min(nStartY+m_nTBlockSize, nHei); 

	nNumberofPixels = (nEndY-nStartY)*(nEndX-nStartX) * 3;	

	fTrans = 0.3f;	
	nTrans = 427;	
	fPreTrs = 0;

	double fNewKSum = 0;						
	double fNewK;							
	double fWiR, fWiG, fWiB;					
	double fPreJR, fPreJG, fPreJB;			
	double fWsum = 0;						
	int nIdx = 0;	
	int nLossCount;

	for(nY=nStartY; nY<nEndY; nY++)
	{
		for(nX=nStartX; nX<nEndX; nX++)
		{
			fPreJB = (double)(pnImageBP[nY*nWid+nX]-m_anAirlight[0]);
			fPreJG = (double)(pnImageGP[nY*nWid+nX]-m_anAirlight[1]);
			fPreJR = (double)(pnImageRP[nY*nWid+nX]-m_anAirlight[2]);	
			if(fPreJB != 0){
				fWiB = m_pfExpLUT[abs(pnImageB[nY*nWid+nX]-pnImageBP[nY*nWid+nX])]; 
				fWsum += fWiB;	
				fNewKSum += fWiB*(double)(pnImageB[nY*nWid+nX]-m_anAirlight[0])/fPreJB;
			}
			if(fPreJG != 0){
				fWiG = m_pfExpLUT[abs(pnImageG[nY*nWid+nX]-pnImageGP[nY*nWid+nX])];
				fWsum += fWiG;
				fNewKSum += fWiG*(double)(pnImageG[nY*nWid+nX]-m_anAirlight[1])/fPreJG;
			}
			if(fPreJR != 0){
				fWiR = m_pfExpLUT[abs(pnImageR[nY*nWid+nX]-pnImageRP[nY*nWid+nX])];
				fWsum += fWiR;
				fNewKSum += fWiR*(double)(pnImageR[nY*nWid+nX]-m_anAirlight[2])/fPreJR;
			}
		}
	}
	fNewK = fNewKSum/fWsum;			
	fPreTrs = pfTransmissionP[nStartY*nWid+nStartX]*fNewK;	


	for(nCounter=0; nCounter<7; nCounter++)
	{

		nSumofSLoss = 0;
		nLossCount = 0;
		nSumofSquaredOuts = 0;
		nSumofOuts = 0;
		
		for(nY=nStartY; nY<nEndY; nY++)
		{
			for(nX=nStartX; nX<nEndX; nX++)
			{
				
				nOutB = ((pnImageB[nY*nWid+nX] - m_anAirlight[0])*nTrans + 128*m_anAirlight[0])>>7;
				nOutG = ((pnImageG[nY*nWid+nX] - m_anAirlight[1])*nTrans + 128*m_anAirlight[1])>>7;
				nOutR = ((pnImageR[nY*nWid+nX] - m_anAirlight[2])*nTrans + 128*m_anAirlight[2])>>7;		// (I-A)/t + A --> ((I-A)*k*128 + A*128)/128

				if(nOutR>255)
				{
					nSumofSLoss += (nOutR - 255)*(nOutR - 255);
					nLossCount++;
				}
				else if(nOutR < 0)
				{
					nSumofSLoss += nOutR * nOutR;
					nLossCount++;
				}
				if(nOutG>255)
				{
					nSumofSLoss += (nOutG - 255)*(nOutG - 255);
					nLossCount++;
				}
				else if(nOutG < 0)
				{
					nSumofSLoss += nOutG * nOutG;
					nLossCount++;
				}
				if(nOutB>255)
				{
					nSumofSLoss += (nOutB - 255)*(nOutB - 255);
					nLossCount++;
				}
				else if(nOutB < 0)
				{
					nSumofSLoss += nOutB * nOutB;
					nLossCount++;
				}
				nSumofSquaredOuts += nOutB * nOutB + nOutR * nOutR + nOutG * nOutG;
				nSumofOuts += nOutR + nOutG + nOutB;
			}
		}
		fMean = (double)(nSumofOuts)/(double)(nNumberofPixels);
		fCost = m_fLambda1 * (double)nSumofSLoss/(double)(nNumberofPixels)
			- ((double)nSumofSquaredOuts/(double)nNumberofPixels - fMean*fMean)	
			+ m_fLambda2/fPreTrs/fPreTrs*fWsum/(double)nNumberofPixels*((fPreTrs-fTrans)*(fPreTrs-fTrans)*255.0f*255.0f);

		if(nCounter==0 || fMinCost > fCost)
		{
			fMinCost = fCost;
			fOptTrs = fTrans;
		}

		fTrans += 0.1f;
		nTrans = (int)(1/fTrans*128.0f);
	}
	return fOptTrs;
}
