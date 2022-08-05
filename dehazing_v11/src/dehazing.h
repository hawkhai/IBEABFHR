#include <Windows.h>
#include <opencv2\core.hpp>
#include <opencv2\highgui\highgui.hpp>
#include <opencv2\imgcodecs.hpp>
//#include <opencv2/opencv.hpp>
#include <opencv2/imgproc/imgproc_c.h>
#include <opencv2/imgproc/types_c.h>
#include <opencv2/highgui/highgui_c.h> 
#include <opencv2\imgproc.hpp>
//#include <omp.h>
//#include <emmintrin.h>

#define CLIP(x) ((x)<(0)?0:((x)>(255)?(255):(x)))
#define CLIP_Z(x) ((x)<(0)?0:((x)>(1.0f)?(1.0f):(x)))
#define CLIP_TRS(x) ((x)<(0.1f)?0.1f:((x)>(1.0f)?(1.0f):(x)))

#define cvCopyImage( src, dst )     cvCopy( src, dst, 0 )
#define cvCvtPixToPlane             cvSplit
#define cvCvtPlaneToPix             cvMerge
//#define cvMean_StdDev               cvAvgSdv
void cvMean_StdDev(const CvArr* image, double* mean, double* sdv, const CvArr* mask CV_DEFAULT(0));
enum {
    CV_LOAD_IMAGE_UNCHANGED = -1,
    CV_LOAD_IMAGE_GRAYSCALE = 0,
    CV_LOAD_IMAGE_COLOR = 1,
    CV_LOAD_IMAGE_ANYDEPTH = 2,
    CV_LOAD_IMAGE_ANYCOLOR = 4,
    CV_LOAD_IMAGE_IGNORE_ORIENTATION = 128
};
IplImage* cvLoadImage(const char* filename,
    int 	iscolor = CV_LOAD_IMAGE_COLOR
);
int cvSaveImage(const char* filename,
    const CvArr* image,
    const int* params = 0
);

double* newdouble(int size);
void checkdouble(double* arr, int size);
int* newint(int size);

using namespace std;

class dehazing
{
public:
	dehazing();
	dehazing(int nW, int nH, bool bPrevFlag, bool bPosFlag);
	dehazing(int nW, int nH, int nTBlockSize, bool bPrevFlag, bool bPosFlag, double fL1, double fL2, int nGBlock);
	~dehazing(void);
	
	void	HazeRemoval(IplImage* imInput, IplImage* imOutput, int nFrame);
	void	ImageHazeRemoval(IplImage* imInput, IplImage* imOutput);	
	void	LambdaSetting(double fLambdaLoss, double fLambdaTemp);
	void	DecisionUse(bool bChoice);
	void	TransBlockSize(int nBlockSize);
	void	FilterBlockSize(int nBlockSize);
	void	AirlightSerachRange(POINT pointTopLeft, POINT pointBottomRight);
	void	SetFilterStepSize(int nStepsize);
	void	PreviousFlag(bool bPrevFlag);
	void	FilterSigma(double nSigma);
	bool	Decision(IplImage* imInput, IplImage* imOutput, int nThreshold);

	int*	GetAirlight();
	int*	GetYImg();
	double*	GetTransmission();
	
private:

	//320*240 size
	double*	m_pfSmallY;			//Y image
	double*	m_pfSmallPk_p;		//(Y image) - (mean of Y image)
	double*	m_pfSmallNormPk;	//Normalize된 Y image
	double*	m_pfSmallInteg;		//Gaussian weight가 적용된 transmission 결과
	double*	m_pfSmallDenom;		//Gaussian weight가 저장된 행렬

	int*	m_pnSmallYImg;		//입력 영상의 Y채널

	int*	m_pnSmallRImg;		//입력 영상의 R채널
	int*	m_pnSmallGImg;		//입력 영상의 G채널
	int*	m_pnSmallBImg;		//입력 영상의 B채널

	double*	m_pfSmallTransP;	//이전 프레임의 transmission 영상
	double*	m_pfSmallTrans;		//초기 transmission 영상
	double*	m_pfSmallTransR;	//정련된 transmission 영상
	int*	m_pnSmallYImgP;		//이전 프레임의 Y채널

	int*	m_pnSmallRImgP;		//이전 프레임의 Y채널
	int*	m_pnSmallGImgP;		//이전 프레임의 Y채널
	int*	m_pnSmallBImgP;		//이전 프레임의 Y채널

	//Original size
	double*	m_pfY;				//Y image
	double*	m_pfPk_p;			//(Y image) - (mean of Y image)
	double*	m_pfNormPk;			//Normalize된 Y image
	double*	m_pfInteg;			//Gaussian weight가 적용된 transmission 결과
	double*	m_pfDenom;			//Gaussian weight가 저장된 행렬
	
	int*	m_pnYImg;			//입력 영상의 Y채널
	int*	m_pnYImgP;			//입력 영상의 Y채널

	int*	m_pnRImg;			//입력 영상의 Y채널
	int*	m_pnGImg;			//입력 영상의 Y채널
	int*	m_pnBImg;			//입력 영상의 Y채널

	int*	m_pnRImgP;			//입력 영상의 Y채널
	int*	m_pnGImgP;			//입력 영상의 Y채널
	int*	m_pnBImgP;			//입력 영상의 Y채널

	double*	m_pfTransmission;	//초기 transmission
	double*	m_pfTransmissionP;	//초기 transmission
	double*	m_pfTransmissionR;	//정련된 transmission 영상
	
	//////////////////////////////////////////////////////////////////////////
	int		m_nStepSize;		//Guided filter의 step size;
	double*	m_pfGuidedLUT;		//Guided filter 내의 gaussian weight를 위한 LUT
	double	m_fGSigma;			//Guided filter 내의 gaussian weight에 대한 sigma

	int		m_anAirlight[3];	// atmospheric light value
	uchar	m_pucGammaLUT[256];	//감마 보정을 위한 LUT
	double	m_pfExpLUT[256];	//Transmission 계산시, 픽셀 차이에 대한 weight용 LUT

	int		m_nAirlight;		//안개값(grey)
	
	bool	m_bPreviousFlag;	//이전 프레임 이용 여부
	double	m_fLambda1;			//Loss cost
	double	m_fLambda2;			//Temporal cost

	int		m_nWid;				//너비
	int		m_nHei;				//높이

	int		m_nTBlockSize;		// Block size for transmission estimation
	int		m_nGBlockSize;		// Block size for guided filter

	//Airlight search range
	int		m_nTopLeftX;				
	int		m_nTopLeftY;
	int		m_nBottomRightX;
	int		m_nBottomRightY;

	bool	m_bPostFlag;		// Flag for post processing(deblocking)
	// function.cpp
	
	void	DownsampleImage();
	void	DownsampleImageColor();
	void	UpsampleTransmission();
	void	MakeExpLUT();
	void	GuideLUTMaker();
	void	GammaLUTMaker(double fParameter);
	void	IplImageToInt(IplImage* imInput);
	void	IplImageToIntColor(IplImage* imInput);
	void	IplImageToIntYUV(IplImage* imInput);
	
	// dehazing.cpp
	void	AirlightEstimation(IplImage* imInput);
	void	RestoreImage(IplImage* imInput, IplImage* imOutput);
	void	PostProcessing(IplImage* imInput, IplImage* imOutput);

	// TransmissionRefinement.cpp
	void	TransmissionEstimation(int *pnImageY, double *pfTransmission, int *pnImageYP, double *pfTransmissionP, int nFrame, int nWid, int nHei);
	void	TransmissionEstimationColor(int *pnImageR, int *pnImageG, int *pnImageB, double *pfTransmission,int *pnImageRP, int *pnImageGP, int *pnImageBP, double *pfTransmissionP,int nFrame, int nWid, int nHei);

	double	NFTrsEstimation(int *pnImageY, int nStartX, int nStartY, int nWid, int nHei);
	double	NFTrsEstimationP(int *pnImageY, int *pnImageYP, double *pfTransmissionP, int nStartX, int nStartY, int nWid, int nHei);

	double	NFTrsEstimationColor(int *pnImageR, int *pnImageG, int *pnImageB, int nStartX, int nStartY, int nWid, int nHei);
	double	NFTrsEstimationPColor(int *pnImageR, int *pnImageG, int *pnImageB, int *pnImageRP, int *pnImageGP, int *pnImageBP, double *pfTransmissionP, int nStartX, int nStartY, int nWid, int nHei);	
	
	// guided filter.cpp
	void	CalcAcoeff(double* pfSigma, double* pfCov, double* pfA1, double* pfA2, double* pfA3, int nIdx);
	void	BoxFilter(double* pfInArray, int nR, int nWid, int nHei, double*& fOutArray);
	void	BoxFilter(double* pfInArray1, double* pfInArray2, double* pfInArray3, int nR, int m_nWid, int m_nHei, double*& pfOutArray1, double*& pfOutArray2, double*& pfOutArray3);

	void	GuidedFilterY(int nW, int nH, double fEps);
	void	GuidedFilter(int nW, int nH, double fEps);
	void	GuidedFilterShiftableWindow(double fEps);
	
	void	FastGuidedFilterS();
	void	FastGuidedFilter();

};