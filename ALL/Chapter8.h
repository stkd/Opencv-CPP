
#include<iostream>
#include<math.h>
#include "Params.h"
#include<opencv2\opencv.hpp>
using namespace cv;

void colorbar(Mat image_rgb[3], int level)
{
	int i, j, width;
	width = X_SIZE / 8;
	for (i = 0; i < Y_SIZE; i++) {
		for (j = 0; j < Y_SIZE; j++) {
			if (((j >= 0) && (j < 2 * width) ||
				(j >= 4 * width) && (j < 6 * width)))	/*R平面*/
				image_rgb[0].at<uchar>(i, j) = level;
			else image_rgb[0].at<uchar>(i, j) = 0;
			if ((j >= 0) && (j < 4 * width))			/*G平面*/
				image_rgb[1].at<uchar>(i, j) = level;
			else image_rgb[1].at<uchar>(i, j) = 0;
			if (((j >= 0) && (j < width)) ||
				((j >= 2 * width) && (j < 3 * width)) ||
				((j >= 4 * width) && (j < 5 * width)) ||
				((j >= 6 * width) && (j < 7 * width)))	/*B平面*/
				image_rgb[2].at<uchar>(i, j) = level;
			else image_rgb[2].at<uchar>(i, j) = 0;
		}
	}
}
#define PI6 3.141592
void rgb_to_ysh(Mat image_in_rgb[3], Mat image_out_ysh[3])
{
	int i, j;
	double r, g, b, y, cb, cr, s, h;
	for (i = 0; i < Y_SIZE; i++) {
		for (j = 0; j < X_SIZE; j++) {
			r = (double)image_in_rgb[0].at<uchar>(i, j);
			g = (double)image_in_rgb[1].at<uchar>(i, j);
			b = (double)image_in_rgb[2].at<uchar>(i, j);
			y = 0.2126 * r + 0.7152 * g + 0.0722 * b;
			cb = (b - y) / 1.8556;
			cr = (r - y) / 1.5748;
			s = sqrt(cb * cb + cr * cr);
			if (s != 0) h = atan2(cr, cb)*180.0 / PI6;
			else		h = 0;
			image_out_ysh[0].at<uchar>(i, j) = (int)y;
			image_out_ysh[1].at<uchar>(i, j) = (int)s;
			image_out_ysh[2].at<uchar>(i, j) = (int)h;
		}
	}
}
void ysh_ro_rgb(Mat image_in_ysh[3], Mat image_out_rgb[3])
{
	int i, j;
	double r, g, b, y, cb, cr, rad;
	for (i = 0; i < Y_SIZE; i++) {
		for (j = 0; j < X_SIZE; j++) {
			y = image_in_ysh[0].at<uchar>(i, j);
			rad = (double)(PI*image_in_ysh[2].at<uchar>(i, j) / 180.0);
			cb = image_in_ysh[1].at<uchar>(i, j)*cos(rad);
			cr = image_in_ysh[1].at<uchar>(i, j)*sin(rad);
			r = y + 1.5748 * cr;
			b = y + 1.8556 * cb;
			g = (y - 0.2126 * r - 0.0722 * b) / 0.7152;
			if (r < 0)		r = 0;
			if (r > 255)	r = 255;
			if (g < 0)		g = 0;
			if (g > 255)	g = 255;
			if (b < 0)		b = 0;
			if (b > 255)	b = 255;
			image_out_rgb[0].at<uchar>(i, j) = (unsigned char)r;
			image_out_rgb[1].at<uchar>(i, j) = (unsigned char)g;
			image_out_rgb[2].at<uchar>(i, j) = (unsigned char)b;

		}
	}
}
void y_image(Mat image_in_y, Mat image_out)
{
	int i, j, d;
	for (i = 0; i < Y_SIZE; i++) {
		for (j = 0; j < X_SIZE; j++) {
			d = (int)image_in_y.at<uchar>(i, j);
			if (d < 0)d = 0;
			if (d > 255)d = 255;
			image_out.at<uchar>(i, j) = (unsigned char)d;
		}
	}
}
void sat_image(Mat image_in_sat, Mat image_out)
{
	int i, j, d;
	for (i = 0; i < Y_SIZE; i++) {
		for (j = 0; j < X_SIZE; j++) {
			d = image_in_sat.at<uchar>(i, j);
			if (d < 0)d = 0;
			if (d > 255)d = 255;
			image_out.at<uchar>(i, j) = (unsigned char)d;
		}
	}
}
void hue_image(Mat image_in_sat, Mat image_in_hue, Mat image_out, int org)
{
	int i, j;
	double d;
	for (i = 0; i < Y_SIZE; i++) {
		for (j = 0; j < X_SIZE; j++) {
			if (image_in_sat.at<uchar>(i, j) > 0) {
				d = image_in_hue.at<uchar>(i, j);
				d = d - org;
				if (d < 0)d = 0;
				if (d > 180)d = 360.0 - d;
				d = 255.0 - d * 255.0 / 180.0;
				if (d < 0)d = 0;
				if (d > 255)d = 255;
				image_out.at<uchar>(i, j) = (unsigned char)d;
			}
			else image_out.at<uchar>(i, j) = 0;
		}
	}
}
void tran_ysh(Mat image_in_ysh[3], Mat image_out_ysh[3], double ya, double yb, double sa, double sb, double hb)
{
	int i, j;
	for (i = 0; i < Y_SIZE; i++) {
		for (j = 0; j < X_SIZE; j++) {
			image_out_ysh[0].at<uchar>(i, j) = (int)(image_in_ysh[0].at<uchar>(i, j) * ya + yb);
			image_out_ysh[1].at<uchar>(i, j) = (int)(image_in_ysh[1].at<uchar>(i, j) * sa + sb);
			image_out_ysh[2].at<uchar>(i, j) = (int)(image_in_ysh[2].at<uchar>(i, j) + hb);
			if (image_out_ysh[2].at<uchar>(i, j) > 180)image_out_ysh[2].at<uchar>(i, j) -= (uchar)360;
			if (image_out_ysh[2].at<uchar>(i, j) < -180)image_out_ysh[2].at<uchar>(i, j) += (uchar)360;
		}
	}
}
void range(Mat image_in, int *fmax, int *fmin);
void expand(Mat image_in, Mat iamge_out, int fmax, int fmin);
void expand_rgb(Mat image_in_rgb[3], Mat iamge_out_rgb[3], int type)
{
	int fmax[3], fmin[3], gmax, gmin;
	switch (type) {
		case 1:
			range(image_in_rgb[0], &fmax[0], &fmin[0]);
			range(image_in_rgb[1], &fmax[1], &fmin[1]);
			range(image_in_rgb[2], &fmax[2], &fmin[2]);
			gmax = fmax[0];
			if (fmax[1] > gmax)gmax = fmax[1];
			if (fmax[2] > gmax)gmax = fmax[2];
			gmin = fmin[0];
			if (fmin[1] < gmin)gmin = fmin[1];
			if (fmin[2] < gmin)gmin = fmin[2];
			expand(image_in_rgb[0], iamge_out_rgb[0], gmax, gmin);
			expand(image_in_rgb[1], iamge_out_rgb[1], gmax, gmin);
			expand(image_in_rgb[2], iamge_out_rgb[2], gmax, gmin);
		default:
			range(image_in_rgb[0], &fmax[0], &fmin[0]);
			range(image_in_rgb[1], &fmax[1], &fmin[1]);
			range(image_in_rgb[2], &fmax[2], &fmin[2]);
			expand(image_in_rgb[0], iamge_out_rgb[0], fmax[0], fmin[0]);
			expand(image_in_rgb[1], iamge_out_rgb[1], fmax[1], fmin[1]);
			expand(image_in_rgb[2], iamge_out_rgb[2], fmax[2], fmin[2]);
			break;
	}
}
void range_ys(Mat image_in_ysh[3], int *ymax, int *ymin, int *smax, int *smin);
void expand_ysh(Mat image_in_ysh[3], Mat image_out_ysh[3], int type)
{
	int ymax, ymin, smax, smin;
	double ya, yb, sa;
	switch (type)
	{
		case 1:
			range_ys(image_in_ysh, &ymax, &ymin, &smax, &smin);
			ya = 255.0 / (double)(ymax - ymin);
			yb = -255.0 * ymin / (double)(ymax - ymin);
			sa = 255.0 / (double)smax;
			tran_ysh(image_in_ysh, image_in_ysh, ya, yb, sa, 0.0, 0.0);
		default:
			range_ys(image_in_ysh, &ymax, &ymin, &smax, &smin);
			ya = 255.0 / (double)(ymax - ymin);
			yb = -255.0 * ymin / (double)(ymax - ymin);
			tran_ysh(image_in_ysh, image_in_ysh, ya, yb, 1.0, 0.0, 0.0);
		break;
	}
}
void range_ys(Mat image_in_ysh[3], int *ymax,int *ymin, int *smax, int *smin)
{
	int i, j, n;

	*ymax = 0;
	*ymin = 255;
	for (i = 0; i < Y_SIZE; i++) {
		for (j = 0; j < X_SIZE; j++) {
			n = image_in_ysh[0].at<uchar>(i, j);
			if (n > *ymax)*ymax = n;
			if (n < *ymin)*ymin = n;
		}
	}
	*smax = 0;
	*smin = 255;
	for (i = 0; i < Y_SIZE; i++) {
		for (j = 0; j < X_SIZE; j++) {
			n = image_in_ysh[1].at<uchar>(i, j);
			if (n > *smax)*smax = n;
			if (n < *smin)*smin = n;
		}
	}
}
void pseude_color(Mat image_in_m, Mat image_out_r, Mat image_out_g, Mat image_out_b, int type)
{
	int i, j;
	unsigned char lutr[256], lutg[256], lutb[256];
	switch (type)
	{
		case 1:
			for (i = 0; i < 256; i++)lutr[i] = i;
			for (i = 0; i < 256; i++)lutg[i] = i;
			for (i = 0; i < 256; i++)lutb[i] = 255;
		case 2:
			for (i = 0; i < 128; i++)lutr[i] = 0;
			for (i = 128; i < 256; i++)lutr[i] = 2 * i - 255;
			for (i = 0; i < 128; i++)lutg[i] = 2 * i;
			for (i = 128; i < 256; i++)lutg[i] = 510 - 2 * i;
			for (i = 0; i < 128; i++)lutb[i] = 255 - 2*i;
			for (i = 128; i < 256; i++)lutb[i] = 0;
		default:
			break;
	}
	for (i = 0; i < Y_SIZE; i++) {
		for (j = 0; j < X_SIZE; j++) {
			image_out_r.at<uchar>(i, j) = lutr[image_in_m.at<uchar>(i, j)];
			image_out_g.at<uchar>(i, j) = lutg[image_in_m.at<uchar>(i, j)];
			image_out_b.at<uchar>(i, j) = lutb[image_in_m.at<uchar>(i, j)];
		}
	}
}