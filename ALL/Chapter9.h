#include<iostream>
#include<math.h>
#include "Params.h"
#include<opencv2\opencv.hpp>
using namespace cv;
#define BIAS 128
void hist2_image(Mat image_in1, Mat image_in2, Mat image_hist)
{
	int i, j, kx, ky;
	int hx, hy, max, kk;

	for (i = 0; i < Y_SIZE; i++)
		for (j = 0; j < X_SIZE; j++)
			image_hist.at<uchar>(i, j) = 0;
	max = 0;
	ky = 256 / Y_SIZE;
	kx = 256 / X_SIZE;
	for (i = 0; i < Y_SIZE; i++) {
		for (j = 0; j < X_SIZE; j++) {
			hy = (HIGH - (int)image_in2.at<uchar>(i, j)) / ky;
			hx = ((int)image_in1.at<uchar>(i, j)) / kx;
			if (image_hist.at<uchar>(hy, hx) < HIGH) image_hist.at<uchar>(hy, hx)++;
			if (max < image_hist.at<uchar>(hy, hx)) max = image_hist.at<uchar>(hy, hx);
		}
	}
	for (i = 0; i < Y_SIZE; i++) {
		for (j = 0; j < X_SIZE; j++) {
			if (image_hist.at<uchar>(i, j) != 0) {
				kk = (long)image_hist.at<uchar>(i, j) * HIGH / max + BIAS;
				if (kk > HIGH) image_hist.at<uchar>(i, j) = HIGH;
				else		   image_hist.at<uchar>(i, j) = kk;
			}
		}
	}
	for (i = 0; i < Y_SIZE; i++)image_hist.at<uchar>(i,0) = HIGH;				/*X¶b*/
	for (j = 0; j = X_SIZE; j++)image_hist.at<uchar>(Y_SIZE - 1, j) = HIGH;		/*Y¶b*/
}
void thresh_rgb(Mat image_in_rgb[3], Mat image_out, int rmin, int rmax, int gmin, int gmax, int bmin, int bmax)
{
	int i, j;

	for (i = 0; i < Y_SIZE; i++) {
		for (j = 0; j < X_SIZE; j++) {
			if		(image_in_rgb[0].at<uchar>(i, j) < rmin)	image_out.at<uchar>(i, j) = LOW;
			else if (image_in_rgb[0].at<uchar>(i, j) > rmax)	image_out.at<uchar>(i, j) = LOW;
			else if (image_in_rgb[1].at<uchar>(i, j) < gmin)	image_out.at<uchar>(i, j) = LOW;
			else if (image_in_rgb[1].at<uchar>(i, j) > gmax)	image_out.at<uchar>(i, j) = LOW;
			else if (image_in_rgb[2].at<uchar>(i, j) < bmin)	image_out.at<uchar>(i, j) = LOW;
			else if (image_in_rgb[2].at<uchar>(i, j) > bmax)	image_out.at<uchar>(i, j) = LOW;
			else												image_out.at<uchar>(i, j) = HIGH;
		}
	}
}
void thresh_ysh(Mat image_in_ysh[3], Mat image_out, int ymin, int ymax, int smin, int smax, int hmin, int hmax)
{
	int i, j;

	if (hmax > 180) {
		hmax -= 360;
		for (i = 0; i < Y_SIZE; i++) {
			for (j = 0; j < X_SIZE; j++) {
				if		(image_in_ysh[0].at<uchar>(i, j) < ymin)	image_out.at<uchar>(i, j) = LOW;
				else if (image_in_ysh[0].at<uchar>(i, j) > ymax)	image_out.at<uchar>(i, j) = LOW;
				else if (image_in_ysh[1].at<uchar>(i, j) < smin)	image_out.at<uchar>(i, j) = LOW;
				else if (image_in_ysh[1].at<uchar>(i, j) > smax)	image_out.at<uchar>(i, j) = LOW;
				else if ((image_in_ysh[2].at<uchar>(i, j) < hmin) &&
						(image_in_ysh[2].at<uchar>(i, j) > hmax))	image_out.at<uchar>(i, j) = LOW;
				else												image_out.at<uchar>(i, j) = HIGH;
			}
		}
	}
	else {
		for (i = 0; i < Y_SIZE; i++) {
			for (j = 0; j < X_SIZE; j++) {
				if (image_in_ysh[0].at<uchar>(i, j) < ymin)	image_out.at<uchar>(i, j) = LOW;
				else if (image_in_ysh[0].at<uchar>(i, j) > ymax)	image_out.at<uchar>(i, j) = LOW;
				else if (image_in_ysh[1].at<uchar>(i, j) < smin)	image_out.at<uchar>(i, j) = LOW;
				else if (image_in_ysh[1].at<uchar>(i, j) > smax)	image_out.at<uchar>(i, j) = LOW;
				else if (image_in_ysh[2].at<uchar>(i, j) < hmin)	image_out.at<uchar>(i, j) = LOW;
				else if (image_in_ysh[2].at<uchar>(i, j) > hmax)	image_out.at<uchar>(i, j) = LOW;
				else												image_out.at<uchar>(i, j) = HIGH;
			}
		}
	}
}
void thresh_color_difference(Mat image_in_rgb[3], Mat image_out, int thresh, int type)
{
	int i, j, d;

	for (i = 0; i < Y_SIZE; i++) {
		for (j = 0; j < X_SIZE; j++) {
			switch (type)
			{
				case 1:
					d = ((int)image_in_rgb[1].at<uchar>(i, j)
						+ (int)image_in_rgb[2].at<uchar>(i, j)) / 2
						- (int)image_in_rgb[0].at<uchar>(i, j);
					if (d >= thresh)	image_out.at<uchar>(i, j) = 255;
					else				image_out.at<uchar>(i, j) = 0;
					break;
				case 2:
					d = ((int)image_in_rgb[2].at<uchar>(i, j)
						+ (int)image_in_rgb[0].at<uchar>(i, j)) / 2
						- (int)image_in_rgb[1].at<uchar>(i, j);
					if (d >= thresh)	image_out.at<uchar>(i, j) = 255;
					else				image_out.at<uchar>(i, j) = 0;
					break;
				case 3:
					d = ((int)image_in_rgb[0].at<uchar>(i, j)
						+ (int)image_in_rgb[1].at<uchar>(i, j)) / 2
						- (int)image_in_rgb[2].at<uchar>(i, j);
					if (d >= thresh)	image_out.at<uchar>(i, j) = 255;
					else				image_out.at<uchar>(i, j) = 0;
					break;
				default:
					break;
			}
		}
	}
}