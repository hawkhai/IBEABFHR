#include <iostream>
#include "opencv2/highgui/highgui.hpp"
#include "image_brighten.h"
#include "Timer.h"

int main()
{
    const char* names[] = {
        "01.png",
            "02.png",
            "03.jpg",
            "04.jpg",
            "05.jpg",
            "06.jpg",
            "07.jpg",
            "08.png",
            "09.jpg",
            "10.jpg",
            "11.jpg",
            "12.jpg",
            "13.jpg",
            "14.jpg",
            "15.jpg",
            "16.jfif",
            "17.jfif",
            "18.png",
    "19.jpg",
    "20.png",
    "21.jpg",
    "22.jpg",

        "23.jpg",
    "24.webp",
    "25.jpg",
    "26.webp",
    "27.jpg",

    };
    int size = sizeof(names) / sizeof(names[0]);

    for (int i = 0; i < size; i++) {
        const char* name = names[i];
        std::string fpath("E:\\kpdf\\pdfreader_image\\fastpdf-turbo\\image\\imagetest\\testdata\\contrast\\");
        fpath.append(name);

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

            std::string fpathw = fpath + ".";
            fpathw.append(name);
            fpathw.append(".brighten.png");
            cv::imwrite(fpathw, image_brighten);
        }
    }
    return 0;
}
