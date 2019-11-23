#include<iostream>
#include<math.h>
#include "Params.h"
#include<opencv2\opencv.hpp>
using namespace cv;

void smooth(Mat image_in, Mat image_out, int type)
{
	int i, j, m, n, k, x, y;
	double  sum, num;

	k = type / 2;
	num = (double)type * type;
	for (i = 0; i < image_in.rows; i++) {
		for (j = 0; j < image_in.cols; j++) {
			sum = 0.0;
			for (m = -k; m <= k; m++) {
				for (n = -k; n <= k; n++) {
					y = i + m;
					x = j + n;
					if (y < 0) y = 0;
					if (x < 0) x = 0;
					if (y > image_in.rows)y = image_in.rows;
					if (x > image_in.cols)x = image_in.cols;
					sum += image_in.at<uchar>(i, j);
				}
			}
			sum = sum / num;
			if (sum < 0)sum = 0;
			if (sum > 255)sum = 255;
			image_out.at<uchar>(i, j) = (uchar)sum;
		}
	}
}
/*#define RAND_MAX 32767*/
/*有時定義為2147483647*/
void noise_rand(Mat image_in, Mat image_out, int level)
{
	int i, j;
	int data, noise;
	for (i = 0; i < image_in.rows; i++) {
		for (j = 0; j < image_in.cols; j++) {
			noise = (int)((rand() / (double)RAND_MAX - 0.5) * level * 2.0);
			data = image_out.at<uchar>(i, j) + noise;
			if (data > 255) image_out.at<uchar>(i, j) = 255;
			else if (data < 0) image_out.at<uchar>(i, j) = 0;
			else			   image_out.at<uchar>(i, j) = data;
		}
	}
}
void smooth_weighted(Mat image_in, Mat image_out, int type)
{
	int i, j, m, n;
	double f;
	double c[3][5][5] = { 0,		 0,		  0,	   0,		0,
		0,1.0 / 10,1.0 / 10,1.0 / 10,		0,
		0,1.0 / 10,2.0 / 10,1.0 / 10,		0,
		0,1.0 / 10,1.0 / 10,1.0 / 10,		0,
		0,		 0,		  0,	   0,		0,
		0,		 0,		  0,	   0,		0,
		0,1.0 / 16,2.0 / 16,1.0 / 16,		0,
		0,2.0 / 16,4.0 / 16,2.0 / 16,		0,
		0,1.0 / 16,2.0 / 16,1.0 / 16,		0,
		0,		 0,		  0,	   0,		0,
		0.0030,	0.0133,	 0.0219,  0.0133,  0.0030,
		0.0133,   0.0596,  0.0983,  0.0596,  0.0133,
		0.0219,	0.0983,  0.1621,  0.0983,  0.0219,
		0.0133,	0.0596,  0.0983,  0.0596,  0.0133,
		0.0030,	0.0133,  0.0219,  0.0133,  0.0030 };
	if (type < 1)type = 1;
	if (type > 3)type = 3;
	for (i = 1; i < image_in.rows - 1; i++) {
		for (j = 1; j < image_in.cols - 1; j++) {
			f = 0.0;
			for (m = -2; m <= 2; m++) {
				for (n = -2; n <= 2; n++) {
					if ((i + m >= 0) && (i + m <= image_in.rows - 1) && (j + n >= 0) && (j + n <= image_in.cols - 1))
						f += image_in.at<uchar>(i + m, j + n)*c[type - 1][2 + m][2 + n];
				}
			}
			if (f < 0)f = 0;
			if (f > 255)f = 255;
			image_out.at<uchar>(i, j) = (unsigned char)f;
		}
	}
}
int median_value(unsigned char c[9]);
void median(Mat image_in, Mat image_out)
{
	int i, j;
	unsigned char c[9];

	for (i = 1; i < image_in.rows - 1; i++) {
		for (j = 1; j < image_in.cols - 1; j++) {
			c[0] = image_in.at<uchar>(i - 1, j - 1);
			c[1] = image_in.at<uchar>(i - 1, j);
			c[2] = image_in.at<uchar>(i - 1, j + 1);
			c[3] = image_in.at<uchar>(i, j - 1);
			c[4] = image_in.at<uchar>(i, j);
			c[5] = image_in.at<uchar>(i, j + 1);
			c[6] = image_in.at<uchar>(i + 1, j - 1);
			c[7] = image_in.at<uchar>(i + 1, j);
			c[8] = image_in.at<uchar>(i + 1, j + 1);
			image_out.at<uchar>(i, j) = median_value(c);
		}
	}
}
int median_value(unsigned char c[9])
{
	int i, j, buf;

	for (j = 0; j < 8; j++) {
		for (i = 0; i < 8; i++) {
			if (c[i + 1] < c[i]) {
				buf = c[i + 1];
				c[i + 1] = c[i];
				c[i] = buf;
			}
		}
	}
	return c[4];
}
void noise_spike(Mat image_in, Mat image_out, int number, int level)
{
	int i, x, y;
	int data, noise;
	for (i = 0; i < number; i++) {
		x = (int)((rand() / (double)RAND_MAX) * image_in.cols);
		y = (int)((rand() / (double)RAND_MAX) * image_in.rows);
		noise = (int)((rand() / (double)RAND_MAX - 0.5) * level * 2.0);
		data = image_in.at<uchar>(y, x) + noise;
		if (data > 255)	image_out.at<uchar>(y, x) = 255;
		else if (data < 0)		image_out.at<uchar>(y, x) = 0;
		else					image_out.at<uchar>(y, x) = data;
	}
}
unsigned char average_minvar(unsigned char p[9][9]);
void smooth_edge_preserve(Mat image_in, Mat image_out)
{
	int i, j, k, m;
	unsigned char p[9][9];
	int patx[9][9] = { 0,-1, 0, 1,-1, 0, 1, 0, 0,
		0, 1, 2, 0, 1, 2, 1, 0, 0,
		0, 1, 2, 1, 2, 1, 2, 0, 0,
		0, 1, 0, 1, 2, 1, 2, 0, 0,
		0,-1, 0, 1,-1, 0, 1, 0, 0,
		0,-1,-2,-1, 0,-2,-1, 0, 0,
		0,-2,-1,-2,-1,-2,-1, 0, 0,
		0,-2,-1,-2,-1, 0,-1, 0, 0,
		-1, 0, 1,-1, 0, 1,-1, 0, 1 };
	int paty[9][9] = { 0,-2,-2,-2,-1,-1,-1, 0, 0,
		0,-2,-2,-1,-1,-1, 0, 0, 0,
		0,-1,-1, 0, 0, 1, 1, 0, 0,
		0, 0, 1, 1, 1, 2, 2, 0, 0,
		0, 1, 1, 1, 2, 2, 2, 0, 0,
		0, 0, 1, 1, 1, 2, 2, 0, 0,
		0,-1,-1, 0, 0, 1, 1, 0, 0,
		0,-2,-2,-1,-1,-1, 0, 0, 0,
		-1,-1,-1, 0, 0, 0, 1, 1, 1 };
	for (i = 0; i < image_in.rows; i++)
		for (j = 0; j < image_in.cols; j++)
			image_out.at<uchar>(i, j) = image_in.at<uchar>(i, j);
	for (i = 2; i < image_in.rows - 2; i++) {
		for (j = 2; j < image_in.cols - 2; j++) {
			for (k = 0; k < 9; k++) {
				for (m = 0; m < 9; m++) {
					p[k][m] = image_in.at<uchar>(i + paty[k][m], j + patx[k][m]);
				}
			}
			image_out.at<uchar>(i, j) = average_minvar(p);
		}
	}
}
unsigned char average_minvar(unsigned char p[9][9])
{
	int i, k, n;
	double ave[9], var[9], dmin;

	for (k = 0; k < 8; k++) {
		ave[k] = 0.0;
		for (i = 0; i < 7; i++)ave[k] += (double)p[k][i];
		ave[k] = ave[k] / 7.0;
		var[k] = 0.0;
		for (i = 0; i < 7; i++)
			var[k] += ((double)p[k][i] - ave[k]) * ((double)p[k][i] - ave[k]);
		var[k] = var[k] / 7.0;
	}
	ave[8] = 0.0;
	for (i = 0; i < 9; i++)ave[k] += (double)p[k][i];
	ave[8] = ave[k] / 9.0;
	var[8] = 0.0;
	for (i = 0; i < 9; i++)
		var[k] += ((double)p[k][i] - ave[k]) * ((double)p[k][i] - ave[k]);
	var[k] = var[k] / 9.0;
	dmin = var[0];
	n = 0;
	for (k = 1; k < 9; k++) {
		if (dmin > var[k]) {
			dmin = var[k];
			n = k;
		}
	}
	return (unsigned char)ave[n];
}
#define PI2 3.14
#define SQRT2 1.414
#define TAP 25
/*OFFSET偏移量*/
void laplacian_of_gaussian(Mat image_in, Mat image_out, double var, double amp);
void zero_cross(Mat image_in, Mat image_out);
void log_zero_cross(Mat image_in, Mat image_out, double var)
{
	Mat image_buf;
	laplacian_of_gaussian(image_in, image_buf, var, 1.0);
	zero_cross(image_buf, image_out);
}
void laplacian_of_gaussian(Mat image_in, Mat image_out, double var, double amp)
{
	double c[TAP * 2 + 1][TAP * 2 + 1];
	int i, j, k, m, n, x, y;
	double r2, v2, v4, d;

	k = (int)(3 * SQRT2 * var);
	if (k > TAP)k = TAP;
	v2 = var * var;
	v4 = 1 / (v2 * v2 * PI2);
	for (m = -k; m <= k; m++) {
		for (n = -k; n <= k; n++) {
			r2 = (m*m + n*n) / v2 / 2;
			c[TAP + m][TAP + n] = v4*(r2 - 1)*exp(-r2);
		}
	}
	for (i = 0; i < image_in.rows; i++) {
		for (j = 0; j < image_in.cols; j++) {
			d = OFFSET;
			for (m = -k; m <= k; m++) {
				for (n = -k; n <= k; n++) {
					y = i + m;
					x = j + n;
					if (y < 0)y = 0;
					if (x < 0)y = 0;
					if (y > image_in.rows - 1)y = image_in.rows - 1;
					if (x > image_in.cols - 1)x = image_in.cols - 1;
					d += c[TAP + m][TAP + n] * image_in.at<uchar>(y, x);
				}
			}
			if (d < 0)d = 0;
			if (d > 255)d = 255;
			image_in.at<uchar>(i, j) = (unsigned char)d;
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
void dilation(Mat image_in, Mat image_out)
{
	int i, j;

	for (i = 1; i < image_in.rows - 1; i++) {
		for (j = 1; j < image_in.cols - 1; j++) {
			image_out.at<uchar>(i, j) = image_in.at<uchar>(i, j);
			if (image_in.at<uchar>(i - 1, j - 1) == 255)image_out.at<uchar>(i, j) = 255;
			if (image_in.at<uchar>(i - 1, j) == 255)image_out.at<uchar>(i, j) = 255;
			if (image_in.at<uchar>(i - 1, j + 1) == 255)image_out.at<uchar>(i, j) = 255;
			if (image_in.at<uchar>(i, j - 1) == 255)image_out.at<uchar>(i, j) = 255;
			if (image_in.at<uchar>(i, j) == 255)image_out.at<uchar>(i, j) = 255;
			if (image_in.at<uchar>(i, j + 1) == 255)image_out.at<uchar>(i, j) = 255;
			if (image_in.at<uchar>(i + 1, j - 1) == 255)image_out.at<uchar>(i, j) = 255;
			if (image_in.at<uchar>(i + 1, j) == 255)image_out.at<uchar>(i, j) = 255;
			if (image_in.at<uchar>(i + 1, j + 1) == 255)image_out.at<uchar>(i, j) = 255;
		}
	}
}
void erosion(Mat image_in, Mat image_out)
{
	int i, j;

	for (i = 1; i < image_in.rows - 1; i++) {
		for (j = 1; j < image_in.cols - 1; j++) {
			image_out.at<uchar>(i, j) = image_in.at<uchar>(i, j);
			if (image_in.at<uchar>(i - 1, j - 1) == 0)image_out.at<uchar>(i, j) = 0;
			if (image_in.at<uchar>(i - 1, j) == 0)image_out.at<uchar>(i, j) = 0;
			if (image_in.at<uchar>(i - 1, j + 1) == 0)image_out.at<uchar>(i, j) = 0;
			if (image_in.at<uchar>(i, j - 1) == 0)image_out.at<uchar>(i, j) = 0;
			if (image_in.at<uchar>(i, j) == 0)image_out.at<uchar>(i, j) = 0;
			if (image_in.at<uchar>(i, j + 1) == 0)image_out.at<uchar>(i, j) = 0;
			if (image_in.at<uchar>(i + 1, j - 1) == 0)image_out.at<uchar>(i, j) = 0;
			if (image_in.at<uchar>(i + 1, j) == 0)image_out.at<uchar>(i, j) = 0;
			if (image_in.at<uchar>(i + 1, j + 1) == 0)image_out.at<uchar>(i, j) = 0;
		}
	}
}