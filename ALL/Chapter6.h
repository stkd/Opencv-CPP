#include<iostream>
#include<math.h>
#include "Params.h"
#include<opencv2\opencv.hpp>
using namespace cv;

void amplify(Mat image_in, Mat image_out, double a, double b)
{
	int i, j, d;

	for (i = 0; i < image_in.rows; i++) {
		for (j = 0; j < image_in.cols; j++) {
			d = (int)(image_in.at<uchar>(i, j)*a + b);
			if (d < 0)d = 0;
			if (d > 255)d = 255;
			image_out.at<uchar>(i, j) = (unsigned char)d;
		}
	}
}
void range(Mat image_in, int *fmax, int *fmin)
{
	int i, j, n;
	*fmax = 0;
	*fmin = 255;
	for (i = 0; i < image_in.rows; i++) {
		for (j = 0; j < image_in.cols; j++) {
			n = (int)image_in.at<uchar>(i, j);
			if (n > *fmax)*fmax = n;
			if (n < *fmin)*fmin = n;
		}
	}
}
void expand(Mat image_in, Mat image_out, int fmax, int fmin)
{
	double a, b;
	a = 255.0 / (double)(fmax - fmin);
	b = -255.0 * fmin / (double)(fmax - fmin);
	amplify(image_in, image_out, a, b);
}
#define BUFF_MAX 10000
struct xyw {
	int x, y, w;
}buf[BUFF_MAX];

void sort(Mat image_in, struct xyw data[], int level);
void weight(Mat image_in, int i, int j, int *wt);
void plane(Mat image_in, Mat image_out, long hist[256])
{
	int i, j, iy, jx, sum;
	int delt;
	int low, high;
	int av;
	Mat image_buf;

	av = (int)((image_in.rows) * (image_in.cols) / 256);
	high = 255;
	low = 255;
	for (i = 0; i < image_in.rows; i++) {
		for (j = 0; j < image_in.cols; j++) {
			image_out.at<uchar>(i, j) = 0;
			image_buf.at<uchar>(i, j) = image_in.at<uchar>(i, j);
		}
	}
	for (i = 255; i > 0; i--) {
		for (sum = 0; sum < av; low--)sum = sum + hist[low];
		low++;
		delt = hist[low] - (sum - av);
		sort(image_buf, buf, low);
		if (low > high) {
			for (iy = 0; iy < image_in.rows; iy++) {
				for (jx = 0; jx < image_in.cols; jx++) {
					if (((int)image_buf.at<uchar>(iy, jx) >= low + 1) &&
						((int)image_buf.at<uchar>(iy, jx) <= high))
						image_out.at<uchar>(iy, jx) = (unsigned char)i;
				}
			}
		}
		for (j = 0; j < delt; j++) {
			image_out.at<uchar>(buf[j].y, buf[j].x) = (unsigned char)i;
			image_buf.at<uchar>(buf[j].y, buf[j].x) = (unsigned char)0;
		}
		hist[low] = hist[low] - delt;
		high = low;
	}
}
void sort(Mat image_in, struct xyw data[], int level)
{
	int i, j, inum, wt;
	struct xyw temp;

	inum = 0;
	for (i = 0; i < image_in.rows; i++) {
		for (j = 0; j < image_in.cols; j++) {
			if ((int)image_in.at<uchar>(i, j) == level) {
				weight(image_in, i, j, &wt);
				data[inum].y = i;
				data[inum].x = j;
				data[inum].w = wt;
				inum++;
			}
		}
	}

	for (i = 0; i < inum - 1; i++) {
		for (j = i + 1; j < inum; j++) {
			if (data[i].w <= data[j].w) {
				temp.y = data[i].y;
				temp.x = data[i].x;
				temp.w = data[i].w;
				data[i].y = data[j].y;
				data[i].x = data[j].x;
				data[i].w = data[j].w;
				data[j].y = temp.y;
				data[j].x = temp.x;
				data[j].w = temp.w;
			}
		}
	}
}
void weight(Mat image_in, int i, int j, int *wt)
{
	int dim, djm;
	int dip, djp;
	int k, d[8];

	dim = i - 1;
	djm = j - 1;
	dip = i + 1;
	djp = j + 1;
	if (dim < 0)dim = i;
	if (djm < 0)djm = j;
	if (dip > image_in.rows - 1)dip = i;
	if (djp > image_in.cols - 1)djp = j;
	d[0] = (int)image_in.at<uchar>(dim, djm);
	d[1] = (int)image_in.at<uchar>(dim, j);
	d[2] = (int)image_in.at<uchar>(dim, djp);
	d[3] = (int)image_in.at<uchar>(i, djm);
	d[4] = (int)image_in.at<uchar>(i, djp);
	d[5] = (int)image_in.at<uchar>(dip, djm);
	d[6] = (int)image_in.at<uchar>(dip, j);
	d[7] = (int)image_in.at<uchar>(dip, djp);
	for (k = 0; k < 8; k++)*wt = *wt + d[i];
}
void dither_ordered(Mat image_in, Mat image_out)
{
	static int thres[4][4] = { 0, 8, 2,10,
		12, 4,14, 6,
		3,11, 1, 9,
		15, 7,13, 5 };
	int i, j, m, n;

	for (i = 0; i < image_in.rows / 4; i++) {
		for (j = 0; j < image_in.cols / 4; j++) {
			for (m = 0; m < 4; m++) {
				for (n = 0; n < 4; n++) {
					if ((int)image_in.at<uchar>(i * 4 + m, j * 4 + n) > thres[m][n] * 16 + 8)
						image_out.at<uchar>(i * 4 + m, j * 4 + n) = 255;
					else image_out.at<uchar>(i * 4 + m, j * 4 + n) = 0;
				}
			}
		}
	}
}
void dither_minimized(Mat image_in, Mat image_out)
{
	int i, j;
	int t, d;
	int error[X_SIZE], error1[X_SIZE], error2[X_SIZE];
	t = (255 + 1) / 2;
	for (i = 0; i < image_in.cols; i++) {
		if (image_in.at<uchar>(0, i) > t)image_out.at<uchar>(0, i) = 255;
		else							 image_out.at<uchar>(0, i) = 0;
		error1[i] = image_in.at<uchar>(0, i) - image_out.at<uchar>(0, i);
	}
	for (i = 0; i < image_in.cols; i++) {
		if (image_in.at<uchar>(1, i) > t)image_out.at<uchar>(1, i) = 255;
		else							 image_out.at<uchar>(1, i) = 0;
		error2[i] = image_in.at<uchar>(1, i) - image_out.at<uchar>(1, i);
	}
	for (i = 2; i < image_in.rows; i++) {
		if (image_in.at<uchar>(i, 0) > t)image_out.at<uchar>(i, 0) = 255;
		else							 image_out.at<uchar>(i, 0) = 0;
		error[0] = image_in.at<uchar>(i, 0) - image_out.at<uchar>(i, 0);
		if (image_in.at<uchar>(i, 1) > t)image_out.at<uchar>(i, 1) = 255;
		else							 image_out.at<uchar>(i, 1) = 0;
		error[1] = image_in.at<uchar>(i, 1) - image_out.at<uchar>(i, 1);
		for (j = 2; j < image_in.cols - 2; j++) {
			d = (error1[j - 2] + error[j - 1] * 3 + error1[j] * 5
				+ error1[j + 1] * 3 + error1[j + 2]
				+ error2[j - 2] * 3 + error2[j - 1] * 5 + error2[j] * 7
				+ error2[j + 1] * 5 + error2[j + 2] * 3
				+ error[j - 2] * 5 + error[j - 1] * 7) / 48;
			if ((int)image_in.at<uchar>(i, j) + d > t)image_in.at<uchar>(i, j) = 255;
			else									  image_in.at<uchar>(i, j) = 0;
			error[j] = image_in.at<uchar>(i, j) + d - image_out.at<uchar>(i, j);
		}
		if (image_in.at<uchar>(i, image_in.cols - 2) > t)image_out.at<uchar>(i, image_in.cols - 2) = 255;
		else											 image_out.at<uchar>(i, image_in.cols - 2) = 0;
		error[image_in.cols - 2] = image_in.at<uchar>(i, image_in.cols - 2) - image_out.at<uchar>(i, image_in.cols - 2);
		if (image_in.at<uchar>(i, image_in.cols - 1) > t)image_out.at<uchar>(i, image_in.cols - 1) = 255;
		else											 image_out.at<uchar>(i, image_in.cols - 1) = 0;
		error[image_in.cols - 1] = image_in.at<uchar>(i, image_in.cols - 1) - image_out.at<uchar>(i, image_in.cols - 1);
		for (j = 0; j < image_in.cols; j++) {
			error1[j] = error2[j];
			error2[j] = error[j];
		}
	}
}
void d_quantize(int in, unsigned char *pout, int nq);
void dither_minimized_multi(Mat image_in, Mat image_out, int nq)
{
	int i, j;
	int d;
	int error[X_SIZE], error1[X_SIZE], error2[X_SIZE];

	for (i = 0; i < image_in.cols; i++) {
		d_quantize((int)image_in.at<uchar>(0, i), &image_out.at<uchar>(0, i), nq);
		error1[i] = (int)image_in.at<uchar>(0, i) - (int)image_out.at<uchar>(0, i);
	}
	for (i = 0; i < image_in.cols; i++) {
		d_quantize((int)image_in.at<uchar>(1, i), &image_out.at<uchar>(1, i), nq);
		error2[i] = (int)image_in.at<uchar>(1, i) - (int)image_out.at<uchar>(1, i);
	}
	for (i = 2; i < image_in.rows; i++) {
		d_quantize((int)image_in.at<uchar>(i, 0), &image_out.at<uchar>(i, 0), nq);
		error[i] = (int)image_in.at<uchar>(i, 0) - (int)image_out.at<uchar>(i, 0);
		d_quantize((int)image_in.at<uchar>(i, 1), &image_out.at<uchar>(i, 1), nq);
		error[i] = (int)image_in.at<uchar>(i, 1) - (int)image_out.at<uchar>(i, 1);
		for (j = 2; j < image_in.cols - 2; j++) {
			d = (error1[j - 2] + error[j - 1] * 3 + error1[j] * 5
				+ error1[j + 1] * 3 + error1[j + 2]
				+ error2[j - 2] * 3 + error2[j - 1] * 5 + error2[j] * 7
				+ error2[j + 1] * 5 + error2[j + 2] * 3
				+ error[j - 2] * 5 + error[j - 1] * 7) / 48;
			d_quantize((int)image_in.at<uchar>(i, j) + d, &image_out.at<uchar>(i, j), nq);
			error[j] = (int)image_in.at<uchar>(i, j) - (int)image_out.at<uchar>(i, j);
		}
		d_quantize((int)image_in.at<uchar>(i, image_in.cols - 2) + d, &image_out.at<uchar>(i, image_in.cols - 2), nq);
		error[image_in.cols - 2] = (int)image_in.at<uchar>(i, image_in.cols - 2)
			- (int)image_out.at<uchar>(i, image_in.cols - 2);
		d_quantize((int)image_in.at<uchar>(i, image_in.cols - 1) + d, &image_out.at<uchar>(i, image_in.cols - 1), nq);
		error[image_in.cols - 1] = (int)image_in.at<uchar>(i, image_in.cols - 1)
			- (int)image_out.at<uchar>(i, image_in.cols - 1);
		for (j = 0; j < image_in.cols; j++) {
			error1[j] = error2[j];
			error2[j] = error[j];
		}
	}
}
void quantize(Mat image_in, Mat image_out, int nq)
{
	int i, j;

	for (i = 0; i < image_in.rows; i++)
		for (j = 0; j < image_in.cols; j++)
			d_quantize(image_in.at<uchar>(i, j), &image_out.at<uchar>(i, j), nq);
}
void d_quantize(int in, unsigned char *pout, int nq)
{
	int t, i;

	t = (255 + 1) / (nq - 1);
	i = (int)((double)in / t + 0.5) * t;
	if (i > 255)*pout = 255;
	else *pout = i;
}
