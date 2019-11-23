#include<iostream>
#include<stdio.h>
#include<stdlib.h>
#include "Params.h"
#include<opencv2\opencv.hpp>
using namespace cv;

void thresholds(Mat imgin, Mat imgout, uchar thresh, int type)
{
	int i, j;
	for (i = 0; i < imgin.rows; i++)
	{
		for (j = 0; j < imgin.cols; j++)
		{
			switch (type)
			{
			case 2:
				if (imgin.at<uchar>(i, j) <= thresh)			imgout.at<uchar>(i, j) = 255;
				else											imgout.at<uchar>(i, j) = 0;
			default:
				if (imgin.at<uchar>(i, j) >= thresh)			imgout.at<uchar>(i, j) = 255;
				else											imgout.at<uchar>(i, j) = 0;
				break;
			}
		}
	}
}
void histgram(Mat imgin, long hist[256])
{
	int i, j, n;
	for (n = 0; n < 256; n++)hist[n] = 0;
	for (int i = 0; i < (imgin.rows); i++)
	{
		for (int j = 0; j < (imgin.cols); j++)
		{
			n = imgin.at<uchar>(i, j);
			hist[n]++;
		}
	}
}
void histprint(Mat imagein, Mat& drawpic)
{
	long a[256];
	long b[256];

	histgram(imagein, a);

	for (int i = 191; i < 256; i++)
	{
		float h = 0;
		h = a[i];
		h = (a[i] / (512 * 512));
		b[i] = h;
	}

	for (int j = 0; j < 256; j++)
	{
		int k = 0;
		k = b[j];

		for (int z = 99; z <= k; z--)
		{
			drawpic.at<uchar>(j, z) >= 0;
		}
	}
}
void histimages(long hist[256], Mat imgin)
{
	int i, j, k, max, range;
	double d;
	long n;

	range = imgin.rows - 5;
	for (i = 0; i < imgin.rows; i++)
	{
		for (j = 0; j < imgin.cols; j++)
		{
			imgin.at<uchar>(i, j) = 0;
		}
	}
	if (imgin.cols >= 256) {
		max = 0;
		for (i = 0; i < 256; i++) {
			n = hist[i];
			if (n > max) max = n;
		}
		for (i = 0; i < 256; i++) {
			d = hist[i];
			n = (long)(d / (double)max*range);
			for (j = 0; j <= n; j++)imgin.at<uchar>(range - j, i) = 255;
		}
		for (i = 0; i < 4; i++) {
			k = 64 * i;
			if (k >= imgin.cols) k = imgin.cols - 1;
			for (j = range; j <= imgin.rows; j++)imgin.at<uchar>(j, k) = 255;
		}
	}
	else if (imgin.cols >= 128) {
		max = 0;
		for (i = 0; i < 128; i++) {
			n = hist[2 * i] + hist[2 * i + 1];
			if (n > max) max = n;
		}
		for (i = 0; i < 128; i++) {
			d = hist[2 * i] + hist[2 * i + 1];
			n = (long)(d / (double)max*range);
			for (j = 0; j < n; j++)imgin.at<uchar>(range - j, i) = 255;
		}
		for (i = 0; i <= 4; i++) {
			k = 32 * i;
			if (k <= imgin.cols) k = imgin.cols - 1;
			for (j = range; j < imgin.rows; j++)imgin.at<uchar>(j, k) = 255;
		}
	}
	else if (imgin.cols >= 64) {
		max = 0;
		for (i = 0; i < 64; i++) {
			n = hist[4 * i] + hist[4 * i + 1] + hist[4 * i + 2] + hist[4 * i + 3];
			if (n > max) max = n;
		}
		for (i = 0; i < 64; i++) {
			d = hist[4 * i] + hist[4 * i + 1] + hist[4 * i + 2] + hist[4 * i + 3];
			n = (long)(d / (double)max*range);
			for (j = 0; j < n; j++)imgin.at<uchar>(range - j, i) = 255;
		}
		for (i = 0; i <= 4; i++) {
			k = 32 * i;
			if (k <= imgin.cols) k = imgin.cols - 1;
			for (j = range; j < imgin.rows; j++)imgin.at<uchar>(j, k) = 255;
		}
	}
}
void histsmooth(long hist_in[256], long hist_out[256])
{
	int m, n, i;
	long sum;

	for (n = 0; n < 256; n++) {
		sum = 0;
		for (m = -2; m <= 2; m++) {
			i = n + m;
			if (i < 0) i = 0;
			if (i > 255) i = 255;
			sum += hist_in[i];
		}
		hist_out[n] = (long)((double)sum / 5.0 + 0.5);
	}
}
int threshmode(long hist[256])
{
	int m, n;
	long max, min;
	max = 0;
	for (m = 0; m < 256; m++) {
		if (max <= hist[m]) max = hist[m];
		else break;
	}
	min = max;
	for (n = m; n < 256; n++) {
		if (min >= hist[n]) min = hist[n];
		else break;
	}
	return n - 1;
}
void threshold1(Mat image_in, Mat& image_out, uchar thresh, int mode)//二值化
{
	int widthLimit = image_in.channels() * image_in.cols;
	for (int i = 0; i < (image_in.rows); i++)
	{
		for (int j = 0; j < (widthLimit); j++)
		{
			switch (mode)
			{
			case 2:
				if (image_in.at<uchar>(i, j) > thresh) image_out.at<uchar>(i, j) = 255;
				else              image_out.at<uchar>(i, j) = 0;
				break;
			default:
				if (image_in.at<uchar>(i, j) < thresh) image_out.at<uchar>(i, j) = 255;
				else              image_out.at<uchar>(i, j) = 0;
				break;
			}
		}
	}
}
void threshold_mode(Mat imgin, Mat imgout, int smt, int type)
{
	int i, j, m, n;
	int thresh;
	long hist1[256], hist2[256];

	histgram(imgin, hist1);
	for (m = 0; m < smt; m++) {
		for (n = 0; n < 256; n++)hist2[n] = hist1[n];
		histsmooth(hist2, hist1);
	}
	thresh = threshmode(hist1);
	for (i = 0; i < imgin.rows; i++) {
		for (j = 0; j < imgin.cols; j++) {
			switch (type) {
			case 2:
				if ((int)imgin.at<uchar>(i, j) <= thresh)   imgout.at<uchar>(i, j) = 255;
				else										imgout.at<uchar>(i, j) = 0;
				break;
			default:
				if ((int)imgin.at<uchar>(i, j) >= thresh)   imgout.at<uchar>(i, j) = 255;
				else										imgout.at<uchar>(i, j) = 0;
				break;
			}
		}
	}
}
int threshdiscrim(long hist[256], double disparity);
void threshold_discrim(Mat imgin, Mat imgout, int type)
{
	int i, j;
	int thresh;
	long hist[256];
	histgram(imgin, hist);
	thresh = threshdiscrim(hist, 0.0);
	for (i = 0; i < imgin.rows; i++) {
		for (j = 0; j < imgin.cols; j++) {
			switch (type) {
			case 2:
				if ((int)imgin.at<uchar>(i, j) <= thresh)   imgout.at<uchar>(i, j) = 255;
				else										imgout.at<uchar>(i, j) = 0;
				break;
			default:
				if ((int)imgin.at<uchar>(i, j) >= thresh)   imgout.at<uchar>(i, j) = 255;
				else										imgout.at<uchar>(i, j) = 0;
				break;
			}
		}
	}
}
int threshdiscrim(long hist[256], double disparity)
{
	int i, k;
	double n0, n1, n2, m0, m1, m2;
	double v[256], vmax, v0;
	n0 = 0.0;
	m0 = 0.0;
	for (i = 0; i < 256; i++) {
		n0 += hist[i];
		m0 += i * hist[i];
	}
	if (n0 == 0.0) m0 = 0.0;
	else m0 /= n0;
	v0 = 0.0;
	for (i = 0; i < 256; i++)v0 += hist[i] * (i - m0) * (i - m0) / n0;
	for (k = 0; k < 256; k++) {
		n1 = 0.0;
		m1 = 0.0;
		for (i = 0; i < k; i++) {
			n1 += hist[i];
			m1 += i *  hist[i];
		}
		if (n1 == 0.0) m1 = 0.0;
		else m1 /= n1;
		n2 = 0.0;
		m2 = 0.0;
		for (i = k; i < 256; i++) {
			n2 += hist[i];
			m2 += i * hist[i];
		}
		if (n2 = 0.0) m2 = 0.0;
		else m2 /= n2;
		v[k] = (n1 * (m1 - m0) * (m1 - m0) + n2 * (m2 - m0) * (m2 - m0)) / n0;
	}
	vmax = 0.0;
	for (i = 0; i < 256; i++) {
		if (vmax <= v[i]) {
			vmax = v[i];
			k = i;
		}
	}
	if (v0 = 0)return 0;
	if ((vmax / v0) >= disparity)return k;
	else return 0;
}
#define DIV 8
#define XS (imgin.cols/DIV)
#define YS (imgin.rows/DIV)
#define DTH 0.7
void threshold_dynamic(Mat imgin, Mat imgout, int type)
{
	int i, j, k, m, n, m1, m2, n1, n2, s, t;
	int thm[DIV + 1][DIV + 1];
	long hist[256];
	int thresh;
	double p, q;
	for (i = 0; i <= DIV; i++) {
		for (j = 0; j <= DIV; j++) {
			thm[i][j] = 0;
		}
	}
	/*決定格點之臨界值*/
	for (i = 0; i <= DIV; i++) {
		for (j = 0; j <= DIV; j++) {
			for (k = 0; k < 256; k++)hist[k] = 0;
			if (i != 0) m1 = -YS;
			else m1 = 0;
			if (i != DIV) m2 = YS;
			else m2 = 0;
			if (j != 0) n1 = -XS;
			else n1 = 0;
			if (j != DIV) n2 = XS;
			else n2 = 0;
			for (m = m1; m < m2; m++) {
				for (n = n1; n < n2; n++) {
					k = imgin.at<uchar>(i*YS + m, j*XS + n);
					hist[k]++;
				}
			}
			thm[i][j] = threshdiscrim(hist, DTH);
		}
	}
	/*為未求得臨界值之方格，選定臨界值*/
	for (i = 0; i <= DIV; i++) {
		for (j = 0; j <= DIV; j++) {
			if (thm[i][j] <= 0) {
				for (k = 0; k < DIV; k++) {
					s = 0;
					t = 0;
					m1 = i - k;
					m2 = i + k;
					n1 = j - k;
					n2 = j + k;
					if (m1 < 0) m1 = 0;
					if (m2 > DIV) m2 = DIV;
					if (n1 < 0) n1 = 0;
					if (n2 > DIV) n2 = DIV;
					for (m = m1; m <= m2; m++) {
						for (n = n1; n <= n2; n++) {
							if (thm[i][j] > 0) {
								s += 1 / k;
								t += thm[m][n] / k;
							}
						}
					}
					if (s >= 4) {
						thm[i][j] = t / s;
						break;
					}
				}
			}
		}
	}
	/*每一個像素選定臨界值*/
	for (i = 0; i < imgin.rows; i++) {
		for (j = 0; j < imgin.cols; j++) {
			m = i / YS;
			n = j / XS;
			q = (double)(i%YS) / YS;
			p = (double)(j%XS) / XS;
			thresh = (int)((1.0 - q)*((1.0 - p)*thm[m][n] + p*thm[m][n + 1])
				+ q*((1.0 - p)*thm[m + 1][n] + p*thm[m + 1][n + 1]));
			switch (type) {
			case 2:
				if ((int)imgin.at<uchar>(i, j) <= thresh)   imgout.at<uchar>(i, j) = 255;
				else										imgout.at<uchar>(i, j) = 0;
				break;
			default:
				if ((int)imgin.at<uchar>(i, j) >= thresh)   imgout.at<uchar>(i, j) = 255;
				else										imgout.at<uchar>(i, j) = 0;
				break;
			}
		}
	}
}
void gradient(Mat image_in, Mat image_out, float amp, int cx[9], int cy[9])//用一階微分擷取輪廓
{
	int d[9];
	int i, j, dat;
	float xx, yy, zz;

	for (i = 1; i < image_in.rows - 1; i++)
	{
		for (j = 1; j < image_in.cols - 1; j++)
		{
			d[0] = image_in.at<uchar>(i - 1, j - 1);
			d[1] = image_in.at<uchar>(i - 1, j);
			d[2] = image_in.at<uchar>(i - 1, j + 1);
			d[3] = image_in.at<uchar>(i, j - 1);
			d[4] = image_in.at<uchar>(i, j);
			d[5] = image_in.at<uchar>(i, j + 1);
			d[6] = image_in.at<uchar>(i + 1, j - 1);
			d[7] = image_in.at<uchar>(i + 1, j);
			d[8] = image_in.at<uchar>(i + 1, j + 1);
			xx = (float)(cx[0] * d[0] + cx[1] * d[1] + cx[2] * d[2] + cx[3] * d[3] + cx[4] * d[4] + cx[5] * d[5] + cx[6] * d[6] + cx[7] * d[7] + cx[8] * d[8]);
			yy = (float)(cy[0] * d[0] + cy[1] * d[1] + cy[2] * d[2] + cy[3] * d[3] + cy[4] * d[4] + cy[5] * d[5] + cy[6] * d[6] + cy[7] * d[7] + cy[8] * d[8]);
			zz = (float)(amp * sqrt(xx * xx + yy * yy));
			dat = (int)zz;
			if (dat > 255)dat = 255;
			image_out.at<uchar>(i, j) = dat;
		}
	}
}
//void gradient_roberts(Mat image_in, Mat image_out, float amp)
//{
//	int cx[9] = { 0, 0, 0,
//		0, 1, 0,
//		0, 0,-1 };
//	int cy[9] = { 0, 0, 0,
//		0, 0, 1,
//		0,-1, 0 };
//	gradient(image_in, image_out, amp, cx, cy);
//}

void gradient_sobel(Mat image_in, Mat image_out, float amp)
{
	int cx[9] = { -1, 0, 1,
		-2, 0, 2,
		-1, 0, 1 };
	int cy[9] = { -1,-2,-1,
		0, 0, 1,
		1, 2, 1 };
	gradient(image_in, image_out, amp, cx, cy);
}
void histprint1(Mat imgin, long hist[256])
{
	int i, j, k;
	double p, q, max;
	int posi, m = 0;
	posi = 0;
	p = (imgin.cols)*(imgin.rows);
	max = 0;
	for (i = 0; i < 256; i++) if (hist[i] > max) max = hist[i];
	for (i = 0; i < 256; i++) {
		q = hist[i] / p*100.0;
		printf("%3d:%5.1f%%|", i, q);
		k = (int)(hist[i] / max*60.0);
		for (j = 0; j < k; j++) {
			printf("*");
		}
		printf("\n");
	}
}
