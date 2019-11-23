#include<iostream>
#include<math.h>
#include "Params.h"
#include<opencv2\opencv.hpp>
using namespace cv;

#define DXMAX X_SIZE/8

void calc_ave_devi(int *pdata, int n, double *pave, double *pvei);
void calc_corre(int *pdata1, double ave1, double devi1, int *pdata2, double ave2, double devi2, int n, double *pcorre);

void stereo_diff(Mat image_l, Mat image_r, Mat image_d, int bsize)
{
	int i, j, dx, ddx;
	int x, y, xl, yl, xr, yr;
	double d, min;

	for (i = 0; i < Y_SIZE; i++) {
		for (j = 0; j < X_SIZE; j++) {
			for (dx = 0; dx <= DXMAX; dx++) {
				d = 0;
				for (y = -bsize / 2; y <= bsize / 2; y++) {
					yl = i + y; yr = yl;
					if (yl < 0)				yl = 0;
					if (yl > Y_SIZE - 1)	yl = Y_SIZE - 1;
					for (x = -bsize / 2; x <= bsize / 2; x++) {
						xl = j + x;
						if (xl < 0)				xl = 0;
						if (xl > X_SIZE - 1)	xl = X_SIZE - 1;
						xr = xl - dx;
						if (xr < 0)				xr = 0;
						if (xr > X_SIZE - 1)	xr = X_SIZE - 1;
						d += abs(image_r.at<uchar>(yr, xr) - image_l.at<uchar>(yl, xl));
					}
				}
				if (dx == 0) {
					min = d; ddx = dx;;
				}
				if (d < min) {
					min = d; ddx = dx; y;
				}
			}
			image_d.at<uchar>(i, j) = ddx;
		}
	}
}
void stereo_corre(Mat image_l, Mat image_r, Mat image_d, int bsize)
{
	int i, j, dx, ddx;
	int x, y, xl, yl, xr, yr;
	double max;
	double avel, aver, devil, devir, corre;
	int n;
	int *pdatal0, *pdatal, *pdatar0, *pdatar;

	pdatal0 = (int *)malloc((bsize + 1)*(bsize + 1) * sizeof(int));
	if (pdatal0 == 0)return;
	pdatar0 = (int *)malloc((bsize + 1)*(bsize + 1) * sizeof(int));
	if (pdatar0 == 0)return;
	for (i = 0; i < Y_SIZE; i++) {
		for (j = 0; j < X_SIZE; j++) {
			pdatal = pdatal0;
			n = 0;
			for (y = -bsize / 2; y <= bsize / 2; y++) {
				yl = i + y;
				if (yl < 0)				yl = 0;
				if (yl > Y_SIZE - 1)	yl = Y_SIZE - 1;
				for (x = -bsize / 2; x <= bsize / 2; x++) {
					xl = j + x;
					if (xl < 0)				xl = 0;
					if (xl > X_SIZE - 1)	xl = X_SIZE - 1;
					*pdatal++ = image_l.at<uchar>(yl, xl); n++;
				}
			}
			calc_ave_devi(pdatal0, n, &avel, &devil);
			for (dx = 0; dx <= DXMAX; dx++) {
				pdatar = pdatar0;
				for (y = -bsize / 2; y <= bsize / 2; y++) {
					yr = i + y;
					if (yr < 0)				yr = 0;
					if (yr > Y_SIZE - 1)	yr = Y_SIZE - 1;
					for (x = -bsize / 2; x <= bsize / 2; x++) {
						xr = j + x - dx;
						if (xr < 0)				xr = 0;
						if (xr > X_SIZE - 1)	xr = X_SIZE - 1;
						*pdatar++ = image_r.at<uchar>(yr, xr);
					}
				}
				calc_ave_devi(pdatar0, n, &aver, &devir);
				calc_corre(pdatal0, avel, devil, pdatar0, aver, devir, n, &corre);
				if (dx == 0) {
					max = corre; ddx = dx;
				}
				if (corre > max) {
					max = corre; ddx = dx;
				}
			}
			image_d.at<uchar>(i, j) = ddx;
		}
	}
}
void calc_ave_devi(int *pdata, int n, double *pave, double *pdevi)
{
	int i;
	int *p;
	double sum;
	
	sum = 0; p = pdata;
	for (i = 0; i < n; i++)
		sum += *p++;
	*pave = sum / n;
	sum = 0; p = pdata;
	for (i = 0; i < n; i++) {
		sum += (*p - *pave) * (*p - *pave);
		p++;
	}
	*pdevi = sqrt(sum);
}
void calc_corre(int *pdata1, double ave1, double devi1, int *pdata2, double ave2, double devi2, int n, double *pcorre)
{
	int i;
	int *p1, *p2;
	double sum;
	sum = 0;
	p1 = pdata1; p2 = pdata2;
	for (i = 0; i < n; i++)
		sum += (*p1++ - ave1) * (*p2++ - ave2);
	*pcorre = sum / (devi1 * devi2);
	if (*pcorre < 0)*pcorre = -*pcorre;
}