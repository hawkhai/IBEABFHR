/*
	The main function is an example of video dehazing 
	The core algorithm is in "dehazing.cpp," "guidedfilter.cpp," and "transmission.cpp". 
	You may modify the code to improve the results.

	The detailed description of the algorithm is presented
	in "http://mcl.korea.ac.kr/projects/dehazing". See also 
	J.-H. Kim, W.-D. Jang, Y. Park, D.-H. Lee, J.-Y. Sim, C.-S. Kim, "Temporally
	coherent real-time video dehazing," in Proc. IEEE ICIP, 2012.

	Last updated: 2013-02-14
	Author: Jin-Hwan, Kim.
 */
#include "dehazing.h"
#include <time.h>
#include <conio.h>

#include <shlwapi.h> 
#pragma comment(lib,"shlwapi.lib") 

int mainz(const char* input, const char* output)
{	
	IplImage *imInput = cvLoadImage(input, 1);
	
	int nWid = imInput->width;
	int nHei = imInput->height;
	
	IplImage *imOutput = cvCreateImage(cvSize(nWid, nHei),IPL_DEPTH_8U, 3);

	dehazing dehazingImg(nWid, nHei, 30, false, false, 5.0f, 1.0f, 40);

	//dehazingImg.ImageHazeRemoval(imInput, imOutput);
	//dehazingImg.ImageHazeRemovalYUV(imInput, imOutput);
	dehazingImg.ImageHazeRemoval(imInput, imOutput);
		
	cvSaveImage(output, imOutput);

	//cvReleaseImage(&imInput); 
 	cvReleaseImage(&imOutput);	
	return 0;
}

int main(int argc, char** argv) {

    //mainz("E:\\kSource\\pythonx\\cmyk\\IBEABFHR\\dehazing_v11\\src\\building3.png",
    //      "E:\\kSource\\pythonx\\cmyk\\IBEABFHR\\dehazing_v11\\src\\building3z.png");

    const char* rootdir = "E:\\kpdf\\pdfreader_image\\fastpdf-turbo\\image\\imagetest\\testdata\\contrast";
    const char* types[] = {
        ".png", ".jpg", ".jfif", ".webp", ".jpeg",
    };
    int size = sizeof(types) / sizeof(types[0]);

    for (int num = 0; num < 100; num++) {
        std::string fpath;
        for (int i = 0; i < size; i++) {
            char buffer[1024];
            sprintf_s(buffer, 1024, "%s\\%02d%s", rootdir, num, types[i]);
            if (PathFileExistsA(buffer)) {
                fpath = buffer;
                break;
            }
        }

        cv::Mat image = cv::imread(fpath.c_str());
        if (image.empty()) {
            continue;
        }
        std::string fpathO = fpath + ".dehazing.png";
        mainz(fpath.c_str(), fpathO.c_str());
    }
    getchar();
    return 0;
}