#include <iostream>
#include "opencv2/highgui/highgui.hpp"
#include <opencv2/imgcodecs/imgcodecs.hpp>
#include "image_brighten.h"
#include "Timer.h"

#include <shlwapi.h> 
#pragma comment(lib,"shlwapi.lib") 

int main()
{
    const char* rootdir = "D:\\kSource\\work\\pdfreader2\\fastpdf-turbo\\image\\imagetest\\testdata\\contrast";
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

        for (int i = 0; i < 10; i++) {

            cv::Mat image_brighten;
            std::string name;
            if (!brighten(image, image_brighten, i, name)) {
                break;
            }

            std::string fpathw = fpath + "." + name;
            fpathw.append(".brighten.png");
            cv::imwrite(fpathw, image_brighten);
        }
    }
    return 0;
}
