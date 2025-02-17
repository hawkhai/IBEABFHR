﻿/*
	This source file contains dehazing construnctor, destructor, and 
	dehazing functions. 

	The detailed description of the algorithm is presented
	in "http://mcl.korea.ac.kr/projects/dehazing". See also 
	J.-H. Kim, W.-D. Jang, Y. Park, D.-H. Lee, J.-Y. Sim, C.-S. Kim, "Temporally
	coherent real-time video dehazing," in Proc. IEEE ICIP, 2012.

	Last updated: 2013-02-06
	Author: Jin-Hwan, Kim.
 */
#include <assert.h>
#include <iostream>
#include "dehazing.h"
#include "opencv2/highgui/highgui.hpp"

void cvMean_StdDev(const CvArr* image, double* imean, double* isdv, const CvArr* mask) {
    CvScalar mean = cvScalar(imean[0], imean[1], imean[2], imean[3]);
    CvScalar std_dev = cvScalar(isdv[0], isdv[1], isdv[2], isdv[3]);
    cvAvgSdv(image, &mean, &std_dev, mask);
}
IplImage* cvLoadImage(const char* filename,
    int 	iscolor
) {
    assert(iscolor == 1);
    std::string fpath(filename);
    cv::Mat mat = cv::imread(fpath, cv::ImreadModes::IMREAD_COLOR);
    if (mat.cols % 4 != 0) { // 使其内存连续。
        cv::Size size(mat.cols + 4 - (mat.cols % 4), mat.rows);
        cv::Mat temp;
        cv::resize(mat, temp, size);
        mat = temp;
    }

    auto temp = mat.at<cv::Vec3b>(mat.rows - 1, mat.cols - 1);
    printf("%s %d,%d,%d --input \n", filename, temp[0], temp[1], temp[2]);

    cv::Mat* matp = new cv::Mat(mat);
    return new IplImage(cvIplImage(*matp));
}
int cvSaveImage(const char* filename,
    const CvArr* image,
    const int* params
) {
    assert(params == nullptr);
    cv::Mat mat = cv::cvarrToMat(image);

    auto temp = mat.at<cv::Vec3b>(mat.rows-1, mat.cols - 1);
    printf("%s %d,%d,%d \n", filename, temp[0], temp[1], temp[2]);

    cv::imwrite(filename, mat);
    return 0;
}

dehazing::dehazing(){}

/* 
	Constructor: dehazing constructor

	Parameters: 
		nW - width of input image
		nH - height of input image
		bPrevFlag - boolean for temporal cohenrence of video dehazing
		bPosFlag - boolean for postprocessing.
*/
dehazing::dehazing(int nW, int nH, bool bPrevFlag, bool bPosFlag)
{
	m_nWid			= nW;
	m_nHei			= nH;

	// Flags for temporal coherence & post processing
	m_bPreviousFlag	= bPrevFlag;
	m_bPostFlag		= bPosFlag;

	m_fLambda1		= 5.0f;
	m_fLambda2		= 1.0f;

	// Transmission estimation block size
	m_nTBlockSize		= 40;

	// Guided filter block size, step size(sampling step), & LookUpTable parameter
	m_nGBlockSize		= 40;
	m_nStepSize			= 2;
	m_fGSigma			= 10.0f;

	// to specify the region of atmospheric light estimation
	m_nTopLeftX		= 0;
	m_nTopLeftY		= 0;
	m_nBottomRightX	= m_nWid;
	m_nBottomRightY	= m_nHei;

	m_pfSmallTransP	= newdouble (320*240); // previous trans. (video only)
	m_pfSmallTrans	= newdouble (320*240); // init trans.
	m_pfSmallTransR	= newdouble (320*240); // refined trans.
	m_pnSmallYImg	= newint (320*240);   	
	m_pnSmallYImgP	= newint (320*240);
	m_pfSmallInteg	= newdouble(320*240);
	m_pfSmallDenom	= newdouble(320*240);
	m_pfSmallY		= newdouble(320*240);

	m_pfTransmission	= newdouble (m_nWid*m_nHei);
	m_pfTransmissionR	= newdouble(m_nWid*m_nHei);
	m_pfTransmissionP	= newdouble(m_nWid*m_nHei);
	m_pnYImg			= newint (m_nWid*m_nHei);
	m_pnYImgP			= newint (m_nWid*m_nHei);
	m_pfInteg			= newdouble(m_nWid*m_nHei);
	m_pfDenom			= newdouble(m_nWid*m_nHei);
	m_pfY				= newdouble(m_nWid*m_nHei);

	m_pfSmallPk_p		= newdouble(m_nGBlockSize*m_nGBlockSize);
	m_pfSmallNormPk		= newdouble(m_nGBlockSize*m_nGBlockSize);
	m_pfPk_p			= newdouble(m_nGBlockSize*m_nGBlockSize);
	m_pfNormPk			= newdouble(m_nGBlockSize*m_nGBlockSize);	
	m_pfGuidedLUT		= newdouble(m_nGBlockSize*m_nGBlockSize);	

    memset(m_pfTransmissionR, 0, sizeof(double) * m_nWid * m_nHei);
    memset(m_pfTransmissionP, 0, sizeof(double) * m_nWid * m_nHei);
    memset(m_pfTransmission, 0, sizeof(double) * m_nWid * m_nHei);
}


/* 
	Constructor: dehazing constructor using various options

	Parameters: 
		nW - width of input image
		nH - height of input image
		nTBLockSize - block size for transmission estimation 
		bPrevFlag - boolean for temporal cohenrence of video dehazing
		bPosFlag - boolean for postprocessing
		fL1 - information loss cost parameter (regulating)
		fL2 - temporal coherence paramter
		nGBlock - guided filter block size
*/
dehazing::dehazing(int nW, int nH, int nTBlockSize, bool bPrevFlag, bool bPosFlag, double fL1, double fL2, int nGBlock)
{
	m_nWid			= nW;
	m_nHei			= nH;

	// Flags for temporal coherence & post processing
	m_bPreviousFlag	= bPrevFlag;
	m_bPostFlag		= bPosFlag;

	// parameters for each cost (loss cost, temporal coherence cost)
	m_fLambda1		= fL1;
	m_fLambda2		= fL2;

	// block size for transmission estimation
	m_nTBlockSize		= nTBlockSize;

	// Guided filter block size, step size(sampling step), & LookUpTable parameter
	m_nGBlockSize		= nGBlock;
	m_nStepSize			= 2;	
	m_fGSigma			= 10.0f;

	// to specify the region of atmospheric light estimation
	m_nTopLeftX		= 0;
	m_nTopLeftY		= 0;
	m_nBottomRightX	= m_nWid;
	m_nBottomRightY	= m_nHei;

	m_pfSmallTransP		= newdouble (320*240);
	m_pfSmallTrans		= newdouble (320*240);
	m_pfSmallTransR		= newdouble (320*240);
	m_pnSmallYImg		= newint (320*240);
	m_pnSmallYImgP		= newint (320*240);

	m_pnSmallRImg		= newint (320*240);
	m_pnSmallRImgP		= newint (320*240);
	m_pnSmallGImg		= newint (320*240);
	m_pnSmallGImgP		= newint (320*240);
	m_pnSmallBImg		= newint (320*240);
	m_pnSmallBImgP		= newint (320*240);

	m_pfSmallInteg		= newdouble(320*240);
	m_pfSmallDenom		= newdouble(320*240);
	m_pfSmallY			= newdouble(320*240);

	m_pfTransmission	= newdouble (m_nWid*m_nHei);
	m_pfTransmissionR	= newdouble(m_nWid*m_nHei);
	m_pfTransmissionP	= newdouble(m_nWid*m_nHei);
	m_pnYImg			= newint (m_nWid*m_nHei);
	m_pnYImgP			= newint (m_nWid*m_nHei);

	m_pnRImg			= newint (m_nWid*m_nHei);
	m_pnGImg			= newint (m_nWid*m_nHei);
	m_pnBImg			= newint (m_nWid*m_nHei);

	m_pnRImgP			= newint (m_nWid*m_nHei);
	m_pnGImgP			= newint (m_nWid*m_nHei);
	m_pnBImgP			= newint (m_nWid*m_nHei);

	m_pfInteg			= newdouble(m_nWid*m_nHei);
	m_pfDenom			= newdouble(m_nWid*m_nHei);
	m_pfY				= newdouble(m_nWid*m_nHei);

	m_pfSmallPk_p		= newdouble(m_nGBlockSize*m_nGBlockSize);
	m_pfSmallNormPk		= newdouble(m_nGBlockSize*m_nGBlockSize);
	m_pfPk_p			= newdouble(m_nGBlockSize*m_nGBlockSize);
	m_pfNormPk			= newdouble(m_nGBlockSize*m_nGBlockSize);	
	m_pfGuidedLUT		= newdouble(m_nGBlockSize*m_nGBlockSize);	
}

dehazing::~dehazing(void)
{
	if(m_pfSmallTransP != NULL)
		delete [] m_pfSmallTransP;
	if(m_pfSmallTrans != NULL)
		delete [] m_pfSmallTrans;
	if(m_pfSmallTransR != NULL)
		delete [] m_pfSmallTransR;
	if(m_pnSmallYImg != NULL)
		delete [] m_pnSmallYImg;
	if(m_pfSmallY != NULL)
		delete [] m_pfSmallY;
	if(m_pnSmallYImgP != NULL)
		delete [] m_pnSmallYImgP;

	if(m_pnSmallRImg != NULL)
		delete [] m_pnSmallRImg;
	if(m_pnSmallRImgP != NULL)
		delete [] m_pnSmallRImgP;
	if(m_pnSmallGImg != NULL)
		delete [] m_pnSmallGImg;
	if(m_pnSmallGImgP != NULL)
		delete [] m_pnSmallGImgP;
	if(m_pnSmallBImg != NULL)
		delete [] m_pnSmallBImg;
	if(m_pnSmallBImgP != NULL)
		delete [] m_pnSmallBImgP;

	if(m_pfSmallInteg != NULL)
	delete [] m_pfSmallInteg;
	if(m_pfSmallDenom != NULL)
	delete [] m_pfSmallDenom;
	if(m_pfSmallNormPk != NULL)
	delete [] m_pfSmallNormPk;
	if(m_pfSmallPk_p != NULL)
	delete [] m_pfSmallPk_p;
	
	if(m_pnYImg != NULL)
		delete [] m_pnYImg; 
	if(m_pnYImgP != NULL)
		delete [] m_pnYImgP;
	if(m_pnRImg != NULL)
		delete [] m_pnRImg; 
	if(m_pnGImg != NULL)
		delete [] m_pnGImg; 
	if(m_pnBImg != NULL)
		delete [] m_pnBImg; 
	if(m_pnRImgP != NULL)
		delete [] m_pnRImgP; 
	if(m_pnGImgP != NULL)
		delete [] m_pnGImgP; 
	if(m_pnBImgP != NULL)
		delete [] m_pnBImgP; 

	if(m_pfTransmission != NULL)
		delete [] m_pfTransmission;
	if(m_pfTransmissionR != NULL)
		delete [] m_pfTransmissionR;
	if(m_pfTransmissionP != NULL)
		delete [] m_pfTransmissionP;

	if(m_pfInteg != NULL)
		delete [] m_pfInteg;
	if(m_pfDenom != NULL)
		delete [] m_pfDenom;
	if(m_pfY != NULL)
		delete [] m_pfY;
	if(m_pfPk_p != NULL)
		delete [] m_pfPk_p;
	if(m_pfNormPk != NULL)
		delete [] m_pfNormPk;
	if(m_pfGuidedLUT != NULL)
		delete [] m_pfGuidedLUT;

	m_pfSmallTransP		= NULL;
	m_pfSmallTrans		= NULL;
	m_pfSmallTransR		= NULL;
	m_pnSmallYImg		= NULL;
	m_pnSmallYImgP		= NULL;

	m_pnSmallRImg		= NULL;
	m_pnSmallRImgP		= NULL;
	m_pnSmallGImg		= NULL;
	m_pnSmallGImgP		= NULL;
	m_pnSmallBImg		= NULL;
	m_pnSmallBImgP		= NULL;

	m_pfSmallInteg		= NULL;
	m_pfSmallDenom		= NULL;
	m_pfSmallY			= NULL;

	m_pfTransmission	= NULL;
	m_pfTransmissionR	= NULL;
	m_pfTransmissionP	= NULL;
	m_pnYImg			= NULL;
	m_pnYImgP			= NULL;

	m_pnRImg			= NULL;
	m_pnGImg			= NULL;
	m_pnBImg			= NULL;

	m_pnRImgP			= NULL;
	m_pnGImgP			= NULL;
	m_pnBImgP			= NULL;

	m_pfInteg			= NULL;
	m_pfDenom			= NULL;
	m_pfY				= NULL;

	m_pfSmallPk_p		= NULL;
	m_pfSmallNormPk		= NULL;
	m_pfPk_p			= NULL;
	m_pfNormPk			= NULL;	
	m_pfGuidedLUT		= NULL;	
}

/*
	Function: AirlightEstimation
	Description: estimate the atmospheric light value in a hazy image.
			     it divides the hazy image into 4 sub-block and selects the optimal block, 
				 which has minimum std-dev and maximum average value.
				 *Repeat* the dividing process until the size of sub-block is smaller than 
				 pre-specified threshold value. Then, We select the most similar value to
				 the pure white.
				 IT IS A RECURSIVE FUNCTION.
	Parameter: 
		imInput - input image (cv iplimage)
	Return:
		m_anAirlight: estimated atmospheric light value
 */
void dehazing::AirlightEstimation(IplImage* imInput)
{
	int nMinDistance = 65536;
	int nDistance;

	int nX, nY;

	int nMaxIndex;
	double dpScore[3];
	double dpMean[3];
	double dpStds[3];

	double afMean[4] = {0};
	double afScore[4] = {0};
	double nMaxScore = 0;

	int nWid = imInput->width;
	int nHei = imInput->height;

	int nStep = imInput->widthStep;

	// 4 sub-block
	IplImage *iplUpperLeft = cvCreateImage(cvSize(nWid/2, nHei/2),IPL_DEPTH_8U, 3);
	IplImage *iplUpperRight = cvCreateImage(cvSize(nWid/2, nHei/2),IPL_DEPTH_8U, 3);
	IplImage *iplLowerLeft = cvCreateImage(cvSize(nWid/2, nHei/2),IPL_DEPTH_8U, 3);
	IplImage *iplLowerRight = cvCreateImage(cvSize(nWid/2, nHei/2),IPL_DEPTH_8U, 3);

	IplImage *iplR = cvCreateImage(cvSize(nWid/2, nHei/2),IPL_DEPTH_8U, 1);
	IplImage *iplG = cvCreateImage(cvSize(nWid/2, nHei/2),IPL_DEPTH_8U, 1);
	IplImage *iplB = cvCreateImage(cvSize(nWid/2, nHei/2),IPL_DEPTH_8U, 1);

	// divide 
	cvSetImageROI(imInput, cvRect(0, 0, nWid/2, nHei/2));
	cvCopyImage(imInput, iplUpperLeft);
	cvSetImageROI(imInput, cvRect(nWid/2+nWid%2, 0, nWid, nHei/2));
	cvCopyImage(imInput, iplUpperRight);
	cvSetImageROI(imInput, cvRect(0, nHei/2+nHei%2, nWid/2, nHei));
	cvCopyImage(imInput, iplLowerLeft);
	cvSetImageROI(imInput, cvRect(nWid/2+nWid%2, nHei/2+nHei%2, nWid, nHei));
	cvCopyImage(imInput, iplLowerRight);

	// compare to threshold(200) --> bigger than threshold, divide the block
	if(nHei*nWid > 200)
	{
		// compute the mean and std-dev in the sub-block
		// upper left sub-block
		cvCvtPixToPlane(iplUpperLeft, iplR, iplG, iplB, 0);

		cvMean_StdDev(iplR, dpMean, dpStds);
		cvMean_StdDev(iplG, dpMean+1, dpStds+1);
		cvMean_StdDev(iplB, dpMean+2, dpStds+2);
		// dpScore: mean - std-dev
		dpScore[0] = dpMean[0]-dpStds[0];
		dpScore[1] = dpMean[1]-dpStds[1];
		dpScore[2] = dpMean[2]-dpStds[2];

		afScore[0] = (double)(dpScore[0]+dpScore[1]+dpScore[2]);

		nMaxScore = afScore[0];
		nMaxIndex = 0;

		// upper right sub-block
		cvCvtPixToPlane(iplUpperRight, iplR, iplG, iplB, 0);

		cvMean_StdDev(iplR, dpMean, dpStds);
		cvMean_StdDev(iplG, dpMean+1, dpStds+1);
		cvMean_StdDev(iplB, dpMean+2, dpStds+2);
		
		dpScore[0] = dpMean[0]-dpStds[0];
		dpScore[1] = dpMean[1]-dpStds[1];
		dpScore[2] = dpMean[2]-dpStds[2];

		afScore[1] = (double)(dpScore[0]+dpScore[1]+dpScore[2]);
		
		if(afScore[1] > nMaxScore)
		{
			nMaxScore = afScore[1];
			nMaxIndex = 1;
		}

		// lower left sub-block
		cvCvtPixToPlane(iplLowerLeft, iplR, iplG, iplB, 0);

		cvMean_StdDev(iplR, dpMean, dpStds);
		cvMean_StdDev(iplG, dpMean+1, dpStds+1);
		cvMean_StdDev(iplB, dpMean+2, dpStds+2);
		
		dpScore[0] = dpMean[0]-dpStds[0];
		dpScore[1] = dpMean[1]-dpStds[1];
		dpScore[2] = dpMean[2]-dpStds[2];

		afScore[2] = (double)(dpScore[0]+dpScore[1]+dpScore[2]);
		
		if(afScore[2] > nMaxScore)
		{
			nMaxScore = afScore[2];
			nMaxIndex = 2;
		}

		// lower right sub-block
		cvCvtPixToPlane(iplLowerRight, iplR, iplG, iplB, 0);

		cvMean_StdDev(iplR, dpMean, dpStds);
		cvMean_StdDev(iplG, dpMean+1, dpStds+1);
		cvMean_StdDev(iplB, dpMean+2, dpStds+2);
		
		dpScore[0] = dpMean[0]-dpStds[0];
		dpScore[1] = dpMean[1]-dpStds[1];
		dpScore[2] = dpMean[2]-dpStds[2];

		afScore[3] = (double)(dpScore[0]+dpScore[1]+dpScore[2]);

		if(afScore[3] > nMaxScore)
		{
			nMaxScore = afScore[3];
			nMaxIndex = 3;
		}

		// select the sub-block, which has maximum score
		switch (nMaxIndex)
		{
		case 0:
			AirlightEstimation(iplUpperLeft); break;
		case 1:
			AirlightEstimation(iplUpperRight); break;
		case 2:
			AirlightEstimation(iplLowerLeft); break;
		case 3:
			AirlightEstimation(iplLowerRight); break;
		}
	}
	else
	{
		// select the atmospheric light value in the sub-block
		for(nY=0; nY<nHei; nY++)
		{
			for(nX=0; nX<nWid; nX++)
			{
				// 255-r, 255-g, 255-b
				nDistance = int(sqrt(double(255-(uchar)imInput->imageData[nY*nStep+nX*3])*double(255-(uchar)imInput->imageData[nY*nStep+nX*3])
					+double(255-(uchar)imInput->imageData[nY*nStep+nX*3+1])*double(255-(uchar)imInput->imageData[nY*nStep+nX*3+1])
					+double(255-(uchar)imInput->imageData[nY*nStep+nX*3+2])*double(255-(uchar)imInput->imageData[nY*nStep+nX*3+2])));
				if(nMinDistance > nDistance)
				{
					nMinDistance = nDistance;
					m_anAirlight[0] = (uchar)imInput->imageData[nY*nStep+nX*3];
					m_anAirlight[1] = (uchar)imInput->imageData[nY*nStep+nX*3+1];
					m_anAirlight[2] = (uchar)imInput->imageData[nY*nStep+nX*3+2];
				}
			}
		}
	}
	cvReleaseImage(&iplUpperLeft);
	cvReleaseImage(&iplUpperRight);
	cvReleaseImage(&iplLowerLeft);
	cvReleaseImage(&iplLowerRight);

	cvReleaseImage(&iplR);
	cvReleaseImage(&iplG);
	cvReleaseImage(&iplB);
}





/*
	Function: RestoreImage
	Description: Dehazed the image using estimated transmission and atmospheric light.
	Parameter: 
		imInput - Input hazy image.
	Return:
		imOutput - Dehazed image.
 */
void dehazing::RestoreImage(IplImage* imInput, IplImage* imOutput)
{
	int nStep = imInput->widthStep;

	int nX, nY;
	double fA_R, fA_G, fA_B;

	fA_B = (double)m_anAirlight[0];
	fA_G = (double)m_anAirlight[1];
	fA_R = (double)m_anAirlight[2];

	// post processing flag
	if(m_bPostFlag == true)
	{
		PostProcessing(imInput,imOutput);
	}
	else
	{
#pragma omp parallel for
		// (2) I' = (I - Airlight)/Transmission + Airlight
		for(nY=0; nY<m_nHei; nY++)
		{
			for(nX=0; nX<m_nWid; nX++)
			{	
                if (nY == m_nHei - 1 && nX == m_nWid -1) {
                    int a = 0;
                    a = 1;
                }

                int idxR = nY * nStep + nX * 3;
                int idxG = nY * nStep + nX * 3 + 1;
                int idxB = nY * nStep + nX * 3 + 2;
                double colorR = (uchar)imInput->imageData[idxR];
                double colorG = (uchar)imInput->imageData[idxG];
                double colorB = (uchar)imInput->imageData[idxB];

                double transmission = m_pfTransmissionR[nY * m_nWid + nX];
                transmission = CLIP_Z(transmission);

                colorR = (colorR - fA_B) / transmission + fA_B;
                colorG = (colorG - fA_G) / transmission + fA_G;
                colorB = (colorB - fA_R) / transmission + fA_R;
                
                int rstR = CLIP(colorR);
                int rstG = CLIP(colorG);
                int rstB = CLIP(colorB);

				// (3) Gamma correction using LUT
				imOutput->imageData[idxR] = (uchar)m_pucGammaLUT[(uchar)rstR];
				imOutput->imageData[idxG] = (uchar)m_pucGammaLUT[(uchar)rstG];
				imOutput->imageData[idxB] = (uchar)m_pucGammaLUT[(uchar)rstB];
			}
		}
	}
}

/*
	Function: PostProcessing
	Description: deblocking for blocking artifacts of mpeg video sequence.
	Parameter: 
		imInput - Input hazy frame.
	Return:
		imOutput - Dehazed frame.
 */
void dehazing::PostProcessing(IplImage* imInput, IplImage* imOutput)
{
	const int nStep = imInput->widthStep;
	const int nNumStep= 10;
	const int nDisPos= 20;

	double nAD0, nAD1, nAD2;
	int nS, nX, nY;
	double fA_R, fA_G, fA_B;

	fA_B = (double)m_anAirlight[0];
	fA_G = (double)m_anAirlight[1];
	fA_R = (double)m_anAirlight[2];

#pragma omp parallel for private(nAD0, nAD1, nAD2, nS)
	for(nY=0; nY<m_nHei; nY++)
	{
		for(nX=0; nX<m_nWid; nX++)
		{			
			// (1)  I' = (I - Airlight)/Transmission + Airlight
			imOutput->imageData[nY*nStep+nX*3]	 = (uchar)m_pucGammaLUT[(uchar)CLIP((((double)((uchar)imInput->imageData[nY*nStep+nX*3+0])-fA_B)/CLIP_Z(m_pfTransmissionR[nY*m_nWid+nX]) + fA_B))];
			imOutput->imageData[nY*nStep+nX*3+1] = (uchar)m_pucGammaLUT[(uchar)CLIP((((double)((uchar)imInput->imageData[nY*nStep+nX*3+1])-fA_G)/CLIP_Z(m_pfTransmissionR[nY*m_nWid+nX]) + fA_G))];
			imOutput->imageData[nY*nStep+nX*3+2] = (uchar)m_pucGammaLUT[(uchar)CLIP((((double)((uchar)imInput->imageData[nY*nStep+nX*3+2])-fA_R)/CLIP_Z(m_pfTransmissionR[nY*m_nWid+nX]) + fA_R))];

			// if transmission is less than 0.4, we apply post processing because more dehazed block yields more artifacts
			if( nX > nDisPos+nNumStep && m_pfTransmissionR[nY*m_nWid+nX-nDisPos] < 0.4 )
			{
				nAD0 = (double)((int)((uchar)imOutput->imageData[nY*nStep+(nX-nDisPos)*3])	- (int)((uchar)imOutput->imageData[nY*nStep+(nX-nDisPos-1)*3]));
				nAD1 = (double)((int)((uchar)imOutput->imageData[nY*nStep+(nX-nDisPos)*3+1])	- (int)((uchar)imOutput->imageData[nY*nStep+(nX-nDisPos-1)*3+1]));
				nAD2 = (double)((int)((uchar)imOutput->imageData[nY*nStep+(nX-nDisPos)*3+2])	- (int)((uchar)imOutput->imageData[nY*nStep+(nX-nDisPos-1)*3+2]));
				
				if(__max(__max(abs(nAD0),abs(nAD1)),abs(nAD2)) < 20 
					&&	 abs((uchar)imOutput->imageData[nY*nStep+(nX-nDisPos-1)*3+0]-(uchar)imOutput->imageData[nY*nStep+(nX-nDisPos-1-nNumStep)*3+0])
					+abs((uchar)imOutput->imageData[nY*nStep+(nX-nDisPos-1)*3+1]-(uchar)imOutput->imageData[nY*nStep+(nX-nDisPos-1-nNumStep)*3+1])
					+abs((uchar)imOutput->imageData[nY*nStep+(nX-nDisPos-1)*3+2]-(uchar)imOutput->imageData[nY*nStep+(nX-nDisPos-1-nNumStep)*3+2])
					+abs((uchar)imOutput->imageData[nY*nStep+(nX-nDisPos)*3+0]-(uchar)imOutput->imageData[nY*nStep+(nX-nDisPos-1+nNumStep)*3+0])
					+abs((uchar)imOutput->imageData[nY*nStep+(nX-nDisPos)*3+1]-(uchar)imOutput->imageData[nY*nStep+(nX-nDisPos-1+nNumStep)*3+1])
					+abs((uchar)imOutput->imageData[nY*nStep+(nX-nDisPos)*3+2]-(uchar)imOutput->imageData[nY*nStep+(nX-nDisPos-1+nNumStep)*3+2]) < 30 )
				{
					for( nS = 1 ; nS < nNumStep+1; nS++)
					{
						imOutput->imageData[nY*nStep+(nX-nDisPos-1+nS-nNumStep)*3+0]=(uchar)CLIP((double)((uchar)imOutput->imageData[nY*nStep+(nX-nDisPos-1+nS-nNumStep)*3+0])+(double)nS*nAD0/(double)nNumStep);
						imOutput->imageData[nY*nStep+(nX-nDisPos-1+nS-nNumStep)*3+1]=(uchar)CLIP((double)((uchar)imOutput->imageData[nY*nStep+(nX-nDisPos-1+nS-nNumStep)*3+1])+(double)nS*nAD1/(double)nNumStep);
						imOutput->imageData[nY*nStep+(nX-nDisPos-1+nS-nNumStep)*3+2]=(uchar)CLIP((double)((uchar)imOutput->imageData[nY*nStep+(nX-nDisPos-1+nS-nNumStep)*3+2])+(double)nS*nAD2/(double)nNumStep);
					}
				}
			}
		}
	}
}


/*
	Function: HazeRemoval
	Description: haze removal process

	Parameter:
		imInput - input image
		nFrame - frame number
	Return:
		imOutput - output image
 */
void dehazing::HazeRemoval(IplImage* imInput, IplImage* imOutput, int nFrame)
{
	if(nFrame == 0)
	{
		IplImage* imAir;
		// initializing
		MakeExpLUT(); 
		GuideLUTMaker(); 
		GammaLUTMaker(0.7f); 

		// specify the ROI region of atmospheric light estimation(optional)
		cvSetImageROI(imInput, cvRect(m_nTopLeftX, m_nTopLeftY, m_nBottomRightX-m_nTopLeftX, m_nBottomRightY-m_nTopLeftY));

		imAir = cvCreateImage(cvSize(m_nBottomRightX-m_nTopLeftX, m_nBottomRightY-m_nTopLeftY),IPL_DEPTH_8U, 3);
		cvCopyImage(imInput, imAir);
		
		AirlightEstimation(imAir);

		// Y value of atmosperic light
		m_nAirlight = (((int)(uchar)m_anAirlight[0]*25 + (int)(uchar)m_anAirlight[1]*129 + (int)(uchar)m_anAirlight[2]*66 +128) >>8) + 16;

		cvReleaseImage(&imAir);
		cvResetImageROI(imInput);
	}
	
	IplImageToInt(imInput);
	// down sampling to fast estimation
	DownsampleImage();
	
	// trnasmission estimation
	TransmissionEstimation(m_pnSmallYImg,m_pfSmallTrans, m_pnSmallYImgP, m_pfSmallTransP, nFrame, 320, 240);
		
	// store a data for temporal coherent processing
	memcpy(m_pfSmallTransP, m_pfSmallTrans, 320*240);
	memcpy(m_pnSmallYImgP, m_pnSmallYImg, 320*240);
	
	UpsampleTransmission();
	
	/*
	IplImage *test = cvCreateImage(cvSize(320, 240),IPL_DEPTH_8U, 1);
	for(int nK = 0; nK < 320*240; nK++)
		test->imageData[nK] = (uchar)(m_pnSmallYImg[nK]);
	cvNamedWindow("tests");
	cvShowImage("tests", test);
	cvWaitKey(-1);
	*/
	FastGuidedFilter();

	// (9) 영상 복원 수행
	RestoreImage(imInput, imOutput);
}

/*
	Function: ImageHazeRemoval
	Description: haze removal process for a single image

	Parameter:
		imInput - input image
	Return:
		imOutput - output image
 */
void dehazing::ImageHazeRemoval(IplImage* imInput, IplImage* imOutput)
{
	IplImage* imAir;
	IplImage* imSmallInput;
	// look up table creation
	MakeExpLUT(); 
	GuideLUTMaker(); 
	GammaLUTMaker(0.7f); 

    auto roirect = cvRect(m_nTopLeftX, m_nTopLeftY, m_nBottomRightX - m_nTopLeftX, m_nBottomRightY - m_nTopLeftY);
    auto roisize = cvSize(m_nBottomRightX - m_nTopLeftX, m_nBottomRightY - m_nTopLeftY);
	// specify the ROI region of atmospheric light estimation(optional)
	cvSetImageROI(imInput, roirect);

	imAir = cvCreateImage(roisize,IPL_DEPTH_8U, 3);
	imSmallInput = cvCreateImage(cvSize(320, 240), IPL_DEPTH_8U, 3);

    //imAir = cvCloneImage(imInput);
	cvCopyImage(imInput, imAir);
	
	AirlightEstimation(imAir);
	
	cvReleaseImage(&imAir);
	cvResetImageROI(imInput);
	
	// iplimage to int
	IplImageToIntColor(imInput);
		
	TransmissionEstimationColor(m_pnRImg, m_pnGImg, m_pnBImg, m_pfTransmission, m_pnRImg, m_pnGImg, m_pnBImg, m_pfTransmission, 0, m_nWid, m_nHei);
	
	GuidedFilter(m_nWid, m_nHei, 0.001f);
	//GuidedFilterShiftableWindow(0.001);
	/*
	IplImage *test = cvCreateImage(cvSize(m_nWid, m_nHei),IPL_DEPTH_8U, 1);
	for(int nK = 0; nK < m_nWid*m_nHei; nK++)
		test->imageData[nK] = (uchar)(m_pfTransmissionR[nK]*255);
	cvNamedWindow("tests");
	cvShowImage("tests", test);
	cvWaitKey(-1);
	*/
	// Restore image
	RestoreImage(imInput, imOutput);
	cvReleaseImage(&imSmallInput);
}

/*
	Function:GetAirlight
	Return: air light value 
 */
int* dehazing::GetAirlight()
{
	// (1) airlight값 주소 리턴
	return m_anAirlight;
}

/*
	Function:GetYImg
	Return: get y image array
 */
int* dehazing::GetYImg()
{
	// (1) Y영상 주소 리턴
	return m_pnYImg;
}

/*
	Function:GetTransmission
	Return: get refined transmission array
 */
double* dehazing::GetTransmission()
{
	// (1) 전달량 주소 리턴
	return m_pfTransmissionR;
}

/*
	Function:LambdaSetting
		chnage labmda values
	Parameter:
		fLambdaLoss - new weight coefficient for loss cost
		fLambdaTemp - new weight coefficient for temp cost
 */
void dehazing::LambdaSetting(double fLambdaLoss, double fLambdaTemp)
{
	// (1) 람다값을 재정의 할 때 사용하는 함수
	m_fLambda1 = fLambdaLoss;
	if(fLambdaTemp>0)
		m_fLambda2 = fLambdaTemp;
	else
		m_bPreviousFlag = false;
}

/*
	Function:PreviousFlag
		change boolean value
	Parameter
		bPrevFlag - flag
 */
void dehazing::PreviousFlag(bool bPrevFlag)
{
	// (1) 이전 데이터 사용 유무 결정
	m_bPreviousFlag = bPrevFlag;
}

/*
	Function:TransBlockSize
		change the block size of transmission estimation
	Parameter: 
		nBlockSize - new block size
 */
void dehazing::TransBlockSize(int nBlockSize)
{
	// (1) 전달량 계산의 블록 크기 결정
	m_nTBlockSize = nBlockSize;
}

/*
	Function:FilterBlockSize
		change the block size of guided filter
	Parameter: 
		nBlockSize - new block size
 */
void dehazing::FilterBlockSize(int nBlockSize)
{
	// (1) 전달량 계산의 블록 크기 결정
	m_nGBlockSize = nBlockSize;
}

/*
	Function:SetFilterStepSize
		change the step size of guided filter
	Parameter:
		nStepSize - new step size
 */
void dehazing::SetFilterStepSize(int nStepSize)
{
	// (1) guided filter의 step size 수정
	m_nStepSize = nStepSize;
}

/*
	Function: AirlightSearchRange
	Description: Specify searching range (block shape) by user
	Parameter:
		pointTopLeft - the top left point of searching block
		pointBottomRight - the bottom right point of searching block.
	Return:
		m_nTopLeftX - integer x point
		m_nTopLeftY - integer y point
		m_nBottomRightX - integer x point
		m_nBottomRightY - integer y point.
 */
void dehazing::AirlightSerachRange(POINT pointTopLeft, POINT pointBottomRight)
{
	m_nTopLeftX = pointTopLeft.x;
	m_nTopLeftY = pointTopLeft.y;
	m_nBottomRightX = pointBottomRight.x;
	m_nBottomRightY = pointBottomRight.y;
}

/*
	Function: FilterSigma
	Description: change the gaussian sigma value 
	Parameter:
		nSigma
 */
void dehazing::FilterSigma(double nSigma)
{
	m_fGSigma = nSigma;
}

/*
	Function: Decision
	Description: Decision function for re-estimation of atmospheric light
		in this file, we just implement the decision function and don't 
		apply the decision algorithm in the dehazing function.
	Parameter:
		imSrc1 - first frame 
		imSrc2 - second frame
		nThreshold - threshold value
	Return:
		boolean value 
 */
bool dehazing::Decision(IplImage* imSrc1, IplImage* imSrc2, int nThreshold)
{
	int nX, nY;
	int nStep;

	int nMAD;

	nMAD = 0;

	nStep = imSrc1->widthStep;
	if(imSrc2->widthStep != nStep)
	{
		std::cout << "Size of src1 and src2 is not the same" << std::endl;
		exit(1);
	}
		
	for(nY=0; nY<m_nHei; nY++)
	{
		for(nX=0; nX<m_nWid; nX++)
		{
			nMAD += abs((int)imSrc1->imageData[nY*nStep+nX*3+2]-(int)imSrc2->imageData[nY*nStep+nX*3+2]);
			nMAD += abs((int)imSrc1->imageData[nY*nStep+nX*3+1]-(int)imSrc2->imageData[nY*nStep+nX*3+1]);
			nMAD += abs((int)imSrc1->imageData[nY*nStep+nX*3+0]-(int)imSrc2->imageData[nY*nStep+nX*3+0]);
		}
	}

	nMAD /= (m_nWid*m_nHei);
	if(nMAD > nThreshold)
		return true;
	else
		return false; 
}
