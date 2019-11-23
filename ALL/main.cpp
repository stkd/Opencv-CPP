#include<iostream>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "Chapter3.h"
#include "Chapter4.h"
#include "Chapter5.h"
#include "Chapter6.h"
#include "Chapter7.h"
#include "Chapter8.h"
#include "Chapter9.h"
#include "Chapter10.h"
#include "Chapter11.h"
#include "Chapter12.h"
#include "Chapter13.h"
#include "Chapter14.h"
#include<opencv2\core\core.hpp>
#include<opencv2\opencv.hpp>
using namespace cv;
using namespace std;

int main()
{
	Mat image = imread("lenna.RGB.bmp", 0);
	Mat gradpic[] = { image.clone(),image.clone(),image.clone() };
	int i;

	gradient_difference(image, gradpic[0],3);
	gradient_roberts(image, gradpic[2], 3);

	Mat threshpic[] = { image.clone(),image.clone(),image.clone() };

	for (i = 0; i < 3; i++) {
		thresholds(gradpic[i], threshpic[i],100, 2);
	}

//	Mat histpic(256, 256, CV_8U, Scalar(255));
	Mat smpic[] = { image.clone(),image.clone(),image.clone() };
	for (i = 0; i < 3; i++) {
		thinning(threshpic[i], smpic[i]);
		imshow("grad" + i, gradpic[i]);
		imshow("thresh" + i, threshpic[i]);
		imshow("out"+i, smpic[i]);
		waitKey(0);
	}
	return 0;
}