#pragma once
#include<iostream>
#include<math.h>
#include "Params.h"
#include<opencv2\opencv.hpp>
using namespace cv;
int cx_difference[9] = { 0, 0, 0,
0, 1,-1,
0, 0, 0 };
int cy_difference[9] = { 0, 0, 0,
0, 1, 0,
0,-1, 0 };
int cx_roberts[9] = { 0, 0, 0,
0, 1, 0,
0, 0,-1 };
int cy_roberts[9] = { 0, 0, 0,
0, 0, 1,
0,-1, 0 };
int cx_sobel[9] = { -1, 0, 1,
-2, 0, 2,
-1, 0, 1 };
int cy_sobel[9] = { -1,-2,-1,
0, 0, 0,
1, 2, 1 };

void gradient(Mat image_in, Mat image_out, double amp, int cx[9], int cy[9]);
void gradient_difference(Mat image_in, Mat image_out, double amp)
{
	gradient(image_in, image_out, amp, cx_difference, cy_difference);
}
void gradient_roberts(Mat image_in, Mat image_out, double amp)
{
	gradient(image_in, image_out, amp, cx_roberts, cy_roberts);
}
void gradient_sobel(Mat image_in, Mat image_out, double amp)
{
	gradient(image_in, image_out, amp, cx_sobel, cy_sobel);
}
void gradient(Mat image_in, Mat image_out, double amp, int cx[9], int cy[9])
{
	int d[9];
	int i, j, dat;
	double xx, yy, zz;

	for (i = 1; i < image_in.rows - 1; i++) {
		for (j = 1; j < image_in.cols - 1; j++) {
			d[0] = image_in.at<uchar>(i - 1, j - 1);
			d[1] = image_in.at<uchar>(i - 1, j);
			d[2] = image_in.at<uchar>(i - 1, j + 1);
			d[3] = image_in.at<uchar>(i, j - 1);
			d[4] = image_in.at<uchar>(i, j);
			d[5] = image_in.at<uchar>(i, j + 1);
			d[6] = image_in.at<uchar>(i + 1, j - 1);
			d[7] = image_in.at<uchar>(i + 1, j);
			d[8] = image_in.at<uchar>(i + 1, j + 1);
			xx = (double)(cx[0] * d[0] + cx[1] * d[1] + cx[2] * d[2]
				+ cx[3] * d[3] + cx[4] * d[4] + cx[5] * d[5]
				+ cx[6] * d[6] + cx[7] * d[7] + cx[8] * d[8]);
			yy = (double)(cy[0] * d[0] + cy[1] * d[1] + cy[2] * d[2]
				+ cy[3] * d[3] + cy[4] * d[4] + cy[5] * d[5]
				+ cy[6] * d[6] + cy[7] * d[7] + cy[8] * d[8]);
			zz = (double)(amp*sqrt(xx*xx + yy*yy));
			dat = (int)zz;
			if (dat > 255)dat = 255;
			image_out.at<uchar>(i, j) = (char)dat;
		}
	}
}
void prewitt(Mat image_in, Mat image_out, double amp)
{
	int d0, d1, d2, d3, d4, d5, d6, d7, d8;
	int i, j, k, max, dat;
	int m[8];
	double zz;
	for (i = 1; i < image_in.rows - 1; i++) {
		for (j = 1; j = image_in.cols - 1; j++) {
			d0 = image_in.at<uchar>(i - 1, j - 1);
			d1 = image_in.at<uchar>(i - 1, j);
			d2 = image_in.at<uchar>(i - 1, j + 1);
			d3 = image_in.at<uchar>(i, j - 1);
			d4 = image_in.at<uchar>(i, j);
			d5 = image_in.at<uchar>(i, j + 1);
			d6 = image_in.at<uchar>(i + 1, j - 1);
			d7 = image_in.at<uchar>(i + 1, j);
			d8 = image_in.at<uchar>(i + 1, j + 1);
			m[0] = d0 + d1 + d2 + d3 - 2 * d4 + d5 - d6 - d7 - d8;
			m[1] = d0 + d1 + d2 + d3 - 2 * d4 - d5 + d6 - d7 - d8;
			m[2] = d0 + d1 - d2 + d3 - 2 * d4 - d5 + d6 + d7 - d8;
			m[3] = d0 - d1 - d2 + d3 - 2 * d4 - d5 + d6 + d7 + d8;
			m[4] = -d0 - d1 - d2 + d3 - 2 * d4 + d5 + d6 + d7 + d8;
			m[5] = -d0 - d1 + d2 - d3 - 2 * d4 + d5 + d6 + d7 + d8;
			m[6] = -d0 + d1 + d2 - d3 - 2 * d4 + d5 - d6 + d7 + d8;
			m[7] = d0 + d1 + d2 - d3 - 2 * d4 + d5 - d6 - d7 + d8;
			max = 0;
			for (k = 0; k < 8; k++)	if (max < m[k])max = m[k];
			zz = amp*(double)(max);
			dat = (int)(zz);
			if (dat > 255)dat = 255;
			image_out.at<uchar>(i, j) = (char)dat;
		}
	}
}
#define TMP 128
int ncon(int p[9]);
void thinning(Mat image_in, Mat image_out)
{
	int flg = 1;
	int i, j, k, n;
	int p[9];

	for (i = 0; i < image_in.rows; i++)
		for (j = 0; j < image_in.cols; j++)
			image_out.at<uchar>(i, j) = image_in.at<uchar>(i, j);
	while (flg != 0) {
		flg = 0;
		for (i = 1; i < image_in.rows - 1; i++) {
			for (j = 1; j < image_in.cols - 1; j++) {
				p[0] = image_out.at<uchar>(i, j);
				p[1] = image_out.at<uchar>(i, j + 1);
				p[2] = image_out.at<uchar>(i - 1, j + 1);
				p[3] = image_out.at<uchar>(i - 1, j);
				p[4] = image_out.at<uchar>(i - 1, j - 1);
				p[5] = image_out.at<uchar>(i, j - 1);
				p[6] = image_out.at<uchar>(i + 1, j - 1);
				p[7] = image_out.at<uchar>(i + 1, j);
				p[8] = image_out.at<uchar>(i + 1, j + 1);
				for (k = 0; k < 9; k++) {
					if (p[k] == 255)p[k] = 1;
					else if (p[k] == 0)p[k] = 0;
					else			   p[k] = -1;
				}
				/*條件1:圖形的一部份*/
				if (p[0] != 1)continue;
				/*條件2:境界像素(4個鄰近像素有1個以上是背景)*/
				if (p[1] * p[3] * p[5] * p[7] != 0)continue;
				/*條件3:保留端點(8個鄰近像素有2個以上是圖形)*/
				n = 0;
				for (k = 1; k < 9; k++)if (p[k] != 0)n++;
				if (n < 2)continue;
				/*條件4:保留獨立點(8個鄰近像素有1個以上是圖形)*/
				n = 0;
				for (k = 1; k < 9; k++)if (p[k] == 1)n++;
				if (n < 1)continue;
				/*條件5:保留連結性(8個連結數為1)*/
				if (ncon(p) != 1)continue;
				/*條件6:線條寬度為2時，只去除單邊(8個鄰近像素是否並非全部是-1，如果是-1時，其值為0，8個連結數為1)*/
				n = 0;
				for (k = 1; k < 9; k++) {
					if (p[k] != -1)n++;
					else if (p[k] == -1) {
						p[k] = 0;
						if (ncon(p) == 1)n++;
						p[k] = -1;
					}
				}
				if (n < 8)continue;
				/*條件1~6全部滿足時的削除對象*/
				image_out.at<uchar>(i, j) = TMP;
				flg++;
			}
		}
		for (i = 1; i < image_in.rows - 1; i++)
			for (j = 1; j < image_in.cols - 1; j++)
				if (image_out.at<uchar>(i, j) == TMP)image_out.at<uchar>(i, j) = 0;
	}
}
int ncon(int p[9])
{
	int i, i1, i2;
	int q[9];
	int n = 0;

	for (i = 0; i < 9; i++) {
		if ((p[i] == 1) || (p[i] == -1))q[i] = 0;
		else q[i] = 1;
	}
	for (i = 1; i < 9; i += 2) {
		i1 = i + 1;
		i2 = i + 2;
		if (i2 == 9)i2 = 1;
		n = n + q[i] - q[i] * q[i1] * q[i2];
	}
	return n;
}
/*OFFSET偏移量*/
void laplacian(Mat image_in, Mat image_out, double amp, int type)
{
	int i, j;
	int d;
	int c[3][9] = { 0,-1, 0,-1, 4,-1, 0,-1, 0,
		-1,-1,-1,-1, 8,-1,-1,-1,-1,
		1,-2, 1,-2, 4,-2, 1,-2, 1 };
	type = type - 1;
	if (type < 0)type = 0;
	if (type > 2)type = 2;
	for (i = 1; i < image_in.rows - 1; i++) {
		for (j = 1; j < image_in.cols - 1; j++) {
			d = c[type][0] * image_in.at<uchar>(i - 1, j - 1)
				+ c[type][1] * image_in.at<uchar>(i - 1, j)
				+ c[type][2] * image_in.at<uchar>(i - 1, j + 1)
				+ c[type][3] * image_in.at<uchar>(i, j - 1)
				+ c[type][4] * image_in.at<uchar>(i, j)
				+ c[type][5] * image_in.at<uchar>(i, j + 1)
				+ c[type][6] * image_in.at<uchar>(i + 1, j - 1)
				+ c[type][7] * image_in.at<uchar>(i + 1, j)
				+ c[type][8] * image_in.at<uchar>(i + 1, j + 1);
			d = (int)(d*amp) + OFFSET;
			if (d < 0)d = 0;
			if (d > 255)d = 255;
			image_out.at<uchar>(i, j) = (unsigned char)d;

		}
	}
}
void zero_cross(Mat image_in, Mat image_out)
{
	int i, j;

	for (i = 0; i < image_in.rows; i++)
		for (j = 0; j < image_in.cols; j++)
			image_out.at < uchar >(i, j) = 0;
	for (i = 1; i < image_in.rows - 1; i++) {
		for (j = 1; j < image_in.cols - 1; j++) {
			if ((int)image_in.at<uchar>(i, j) == OFFSET) {
				if (((int)image_in.at<uchar>(i, j + 1) - OFFSET)
					*((int)image_in.at<uchar>(i, j - 1) - OFFSET) < 0)image_out.at<uchar>(i, j) = 255;
				if (((int)image_in.at<uchar>(i + 1, j) - OFFSET)
					*((int)image_in.at<uchar>(i - 1, j) - OFFSET) < 0)image_out.at<uchar>(i, j) = 255;
				if (((int)image_in.at<uchar>(i + 1, j + 1) - OFFSET)
					*((int)image_in.at<uchar>(i - 1, j - 1) - OFFSET) < 0)image_out.at<uchar>(i, j) = 255;
				if (((int)image_in.at<uchar>(i + 1, j - 1) - OFFSET)
					*((int)image_in.at<uchar>(i - 1, j + 1) - OFFSET) < 0)image_out.at<uchar>(i, j) = 255;
			}
			else {
				if (((int)image_in.at<uchar>(i, j) - OFFSET)
					*((int)image_in.at<uchar>(i - 1, j - 1) - OFFSET) < 0)image_out.at<uchar>(i, j) = 255;
				if (((int)image_in.at<uchar>(i, j) - OFFSET)
					*((int)image_in.at<uchar>(i - 1, j) - OFFSET) < 0)image_out.at<uchar>(i, j) = 255;
				if (((int)image_in.at<uchar>(i, j) - OFFSET)
					*((int)image_in.at<uchar>(i - 1, j + 1) - OFFSET) < 0)image_out.at<uchar>(i, j) = 255;
				if (((int)image_in.at<uchar>(i, j) - OFFSET)
					*((int)image_in.at<uchar>(i, j - 1) - OFFSET) < 0)image_out.at<uchar>(i, j) = 255;
			}
		}
	}
}
#define PI 3.141592
#define DMAX 1000
void hough_line(Mat image_in, Mat image_out, double rho, double theta);
void hough(Mat image_in, Mat image_out, Mat image_hough, int thresh, char *buf)
{
	int i, j, u, v, n;
	double rho, theta, d, a;
	double p[DMAX][2];
	int posi, m;

	d = PI / image_in.cols;
	a = image_in.rows / 2 / sqrt(image_in.cols / 2 * image_in.cols / 2 + image_in.rows / 2 * image_in.rows / 2);
	for (i = 0; i < image_in.rows; i++)
		for (j = 0; j < image_in.cols; j++)
			image_hough.at<uchar>(i, j) = 0;
	for (i = 0; i < image_in.rows; i++) {
		for (j = 0; j < image_in.cols; j++) {
			if (image_in.at<uchar>(i, j) != 255)continue;
			for (u = 0; u < image_in.cols; u++) {
				theta = u * d;
				rho = (j - image_in.cols / 2)*cos(theta) + (image_in.rows / 2 - 1)*sin(theta);
				v = (int)(rho * a + image_in.rows / 2 + 0.5);
				if (v >= 0 && v < image_in.rows)
					if (image_hough.at<uchar>(v, u) < 255)image_hough.at<uchar>(v, u) += 1;
			}
		}
	}
	n = 0; posi = 0;
	for (u = 0; u < image_in.cols; u++)
		for (v = 0; v < image_in.rows; v++) {
			if (image_hough.at<uchar>(v, u) < thresh)continue;
			if (u != 0 && v != 0 &&
				image_hough.at<uchar>(v, u) < image_hough.at<uchar>(v - 1, u - 1))continue;
			if (v != 0 &&
				image_hough.at<uchar>(v, u) < image_hough.at<uchar>(v - 1, u))continue;
			if (u != image_in.cols - 1 && v != 0 &&
				image_hough.at<uchar>(v, u) < image_hough.at<uchar>(v - 1, u + 1))continue;
			if (u != 0 &&
				image_hough.at<uchar>(v, u) < image_hough.at<uchar>(v, u - 1))continue;
			if (u != image_in.cols - 1 &&
				image_hough.at<uchar>(v, u) < image_hough.at<uchar>(v, u + 1))continue;
			if (u != 0 && v != image_in.rows - 1 &&
				image_hough.at<uchar>(v, u) < image_hough.at<uchar>(v + 1, u))continue;
			if (u != image_in.cols - 1 && v != image_in.rows - 1 &&
				image_hough.at<uchar>(v, u) < image_hough.at<uchar>(v - 1, u - 1))continue;
			theta = u * d;
			rho = (v - image_in.rows / 2) / a;
			p[n][0] = rho; p[n][1] = theta; n++;
			m = printf(&buf[posi], "theta = %10.3lf, rho = %10.3lf, value = %5d\n", theta * 180 / PI, rho, image_hough.at<uchar>(v, u));
			posi += m;
			if (n == DMAX)return;
		};
	for (i = 0; i < image_in.rows; i++)
		for (j = 0; j < image_in.cols; j++)
			image_out.at<uchar>(i, j) = image_in.at<uchar>(i, j);
	for (i = 0; i < n; i++)
		hough_line(image_in, image_out, p[i][0], p[i][1]);
}
void hough_line(Mat image_in, Mat image_out, double rho, double theta)
{
	int i, j;

	if ((theta >= 0 && theta < PI / 4) || theta >= 3 * PI / 4) {
		for (i = 0; i < image_in.rows; i++) {
			j = (int)((rho - (image_in.rows / 2 - i)*sin(theta))
				/ cos(theta) + image_in.cols / 2 + 0.5);
			if (j > -0 && j < image_in.cols)image_out.at<uchar>(i, j) = 255;
		}
	}
	else {
		for (j = 0; j < image_in.cols; j++) {
			i = (int)((-rho + (j - image_in.cols / 2)*cos(theta))
				/ sin(theta) + image_in.rows / 2 + 0.5);
			if (j >= 0 && j < image_in.rows)image_out.at<uchar>(i, j) = 255;
		}
	}
}
void hough_cross(double theta1, double rho1, double theta2, double rho2, double *x, double *y)
{
	double d, t1, t2;

	t1 = theta1 * PI / 180.0;
	t2 = theta2 * PI / 180.0;
	d = sin(t1 - t2);
	if (d == 0)return;
	*x = (rho2 * sin(t1) - rho1 * sin(t2)) / d;
	*y = (rho1 * cos(t2) - rho2 * cos(t1)) / d;
}