#include "image_brighten.h"
#include <algorithm>
#include "opencv2/imgproc/imgproc.hpp"

void ImageBrighten::brighten(const cv::Mat& src, cv::Mat& dst)
{
	CV_Assert(src.type() == CV_8UC3 || src.type() == CV_8UC1);

	cv::Mat src_inverse = ~src;

	cv::Mat temp;
	fastHazeRemoval(src_inverse, temp);

	dst = ~temp;
}

void ImageBrighten::fastHazeRemoval(const cv::Mat& src, cv::Mat& dst)
{
	CV_Assert(src.type() == CV_8UC3 || src.type() == CV_8UC1);

	if (src.channels() == 1)
	{
		fastHazeRemoval_1Channel(src, dst);
	}
	else
	{
		fastHazeRemoval_3Channel(src, dst);
	}
}

void ImageBrighten::fastHazeRemoval_1Channel(const cv::Mat& src, cv::Mat& dst)
{
	CV_Assert(src.type() == CV_8UC1);

	// step 1. input H(x)
	const cv::Mat& H = src;

	// step 2. calc M(x)
	const cv::Mat& M = H;

	// step 3. calc M_ave(x)
	cv::Mat M_ave;
	const int radius = std::max(50, std::max(H.rows, H.cols) / 20); // should not be too small, or there will be a halo artifact 
	cv::boxFilter(M, M_ave, -1, cv::Size(2 * radius + 1, 2 * radius + 1));

	// step 4. calc m_av
	const float m_av = float(cv::mean(M)[0] / 255.0);

	// step 5. calc L(x)
	const float p = 1.0f - m_av + 0.9f; // a simple parameter selection strategy, for reference only
	const float coeff = std::min(p * m_av, 0.9f);
	cv::Mat L(H.size(), CV_32FC1);
	for (int y = 0; y < L.rows; ++y)
	{
		const uchar* M_line = M.ptr<uchar>(y);
		const uchar* M_ave_line = M_ave.ptr<uchar>(y);
		float* L_line = L.ptr<float>(y);
		for (int x = 0; x < L.cols; ++x)
		{
			L_line[x] = std::min(coeff * M_ave_line[x], float(M_line[x]));
		}
	}

	// step 6. calc A
	double max_H = 0.0;
	cv::minMaxLoc(H, nullptr, &max_H);
	double max_M_ave = 0.0;
	cv::minMaxLoc(M_ave, nullptr, &max_M_ave);
	const float A = 0.5f * float(max_H) + 0.5f * float(max_M_ave);

	// step 7. get F(x)
	cv::Mat F(H.size(), CV_8UC1);
	for (int y = 0; y < F.rows; ++y)
	{
		const uchar* H_line = H.ptr<uchar>(y);
		const float* L_line = L.ptr<float>(y);
		uchar* F_line = F.ptr<uchar>(y);
		for (int x = 0; x < F.cols; ++x)
		{
			const float l = L_line[x];
			const float factor = 1.0f / (1.0f - l / A);
			F_line[x] = cv::saturate_cast<uchar>((float(H_line[x]) - l) * factor);
		}
	}

	dst = F;
}

void ImageBrighten::fastHazeRemoval_3Channel(const cv::Mat& src, cv::Mat& dst)
{
	CV_Assert(src.type() == CV_8UC3);

	// step 1. input H(x)
	const cv::Mat& H = src;

	// step 2. calc M(x)
	cv::Mat M(H.size(), CV_8UC1);
	uchar max_H = 0; // used in step 6
	for (int y = 0; y < M.rows; ++y)
	{
		const cv::Vec3b* H_line = H.ptr<cv::Vec3b>(y);
		uchar* M_line = M.ptr<uchar>(y);
		for (int x = 0; x < M.cols; ++x)
		{
			const cv::Vec3b& h = H_line[x];
			M_line[x] = std::min(h[2], std::min(h[0], h[1]));
			max_H = std::max(std::max(h[0], h[1]), std::max(h[2], max_H));
		}
	}

	// step 3. calc M_ave(x)
	cv::Mat M_ave;
	const int radius = std::max(50, std::max(H.rows, H.cols) / 20); // should not be too small, or there will be a halo artifact 
	cv::boxFilter(M, M_ave, -1, cv::Size(2 * radius + 1, 2 * radius + 1));

	// step 4. calc m_av
	const float m_av = float(cv::mean(M)[0] / 255.0);

	// step 5. calc L(x)
	const float p = 1.0f - m_av + 0.9f; // a simple parameter selection strategy, for reference only
	const float coeff = std::min(p * m_av, 0.9f);
	cv::Mat L(H.size(), CV_32FC1);
	for (int y = 0; y < L.rows; ++y)
	{
		const uchar* M_line = M.ptr<uchar>(y);
		const uchar* M_ave_line = M_ave.ptr<uchar>(y);
		float* L_line = L.ptr<float>(y);
		for (int x = 0; x < L.cols; ++x)
		{
			L_line[x] = std::min(coeff * M_ave_line[x], float(M_line[x]));
		}
	}

	// step 6. calc A
	double max_M_ave = 0.0;
	cv::minMaxLoc(M_ave, nullptr, &max_M_ave);
	const float A = 0.5f * max_H + 0.5f * float(max_M_ave);

	// step 7. get F(x)
	cv::Mat F(H.size(), CV_8UC3);
	for (int y = 0; y < F.rows; ++y)
	{
		const cv::Vec3b* H_line = H.ptr<cv::Vec3b>(y);
		const float* L_line = L.ptr<float>(y);
		cv::Vec3b* F_line = F.ptr<cv::Vec3b>(y);
		for (int x = 0; x < F.cols; ++x)
		{
			const cv::Vec3b& h = H_line[x];
			const float l = L_line[x];
			cv::Vec3b& f = F_line[x];
			const float factor = 1.0f / (1.0f - l / A);
			f[0] = cv::saturate_cast<uchar>((float(h[0]) - l) * factor);
			f[1] = cv::saturate_cast<uchar>((float(h[1]) - l) * factor);
			f[2] = cv::saturate_cast<uchar>((float(h[2]) - l) * factor);
		}
	}

	dst = F;
}

void splitv(const cv::Mat& rgb, std::vector<cv::Mat>& hsvli, cv::Mat& hsvV) {
    cv::Mat hsv;
    cv::cvtColor(rgb, hsv, cv::COLOR_BGR2HSV_FULL);
    cv::split(hsv, hsvli);
    hsvV = hsvli[2];
}

void mergev(std::vector<cv::Mat>& hsvli, const cv::Mat& hsvV, cv::Mat& rgb) {
    hsvli[2] = hsvV;
    cv::Mat hsv;
    cv::merge(hsvli, hsv);
    cv::cvtColor(hsv, rgb, cv::COLOR_HSV2BGR_FULL);
}

void applyClahe(cv::Mat& hsvV) {
    cv::Ptr<cv::CLAHE> clahe = cv::createCLAHE(1.0);
    cv::Mat tempk;
    clahe->apply(hsvV, tempk);//equalizeHist(labL, labL);
    hsvV = tempk;
}

void applyBrighten(cv::Mat& hsvV) {
    cv::Mat tempk;
    ImageBrighten().brighten(hsvV, tempk);
    hsvV = tempk;
}

bool brighten(const cv::Mat& src, cv::Mat& dst, int arg, std::string& name)
{
    std::vector<cv::Mat> hsvli;
    cv::Mat hsvV;
    splitv(src, hsvli, hsvV);

    if (arg == 0) {
        applyClahe(hsvV);
        name = "CLAHE";
    }
    else if (arg == 1) {
        applyBrighten(hsvV);
        name = "brighten";
    }
    else if (arg == 2) {
        ImageBrighten().brighten(src, dst);
        name = "brighten3";
        return true;
    }
    else if (arg == 3) {
        applyClahe(hsvV);
        name = "CLAHE";
    }
    else {
        return false;
    }

    cv::Mat tempk;
    mergev(hsvli, hsvV, tempk);

    if (arg == 3) {
        name = "CLAHE+brighten";
        ImageBrighten().brighten(tempk, dst);
    }
    else {
        dst = tempk;
    }
    return true;
}
