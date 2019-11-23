#include<iostream>
#include<math.h>
#include "Params.h"
#include<opencv2\opencv.hpp>
using namespace cv;

#define L_BASE 100
void labelset(Mat image, int xs, int ys, int label);
int labeling(Mat image_in, Mat image_label, int *cnt, char *buf)
{
	int i, j, label;

	for (i = 0; i < image_in.rows; i++)
		for (j = 0; j < image_in.cols; j++)
			image_label.at<uchar>(i, j) = image_in.at<uchar>(i, j);
	label = L_BASE;
	for (i = 0; i < image_in.rows; i++)
		for (j = 0; j < image_in.cols; j++) {
			if (image_label.at<uchar>(i, j) == 255) {
				if (label >= 255) {
					printf(buf, "-Error! Too many labels. ");
					return -1;
				}
				labelset(image_label, j, i, label); label++;
			}
		}
	*cnt = label - L_BASE;
	printf(buf, "number of labels : %d ", *cnt);
	return 0;
}
void labelset(Mat image, int xs, int ys, int label)
{
	int i, j, cnt, im, ip, jm, jp;

	image.at<uchar>(ys, xs) = label;
	for (;;) {
		cnt = 0;
		for (i = 0; i < image.rows; i++)
			for (j = 0; j < image.cols; j++)
				if (image.at<uchar>(i, j) == label) {
					im = i - 1; ip = i + 1; jm = j - 1; jp = j + 1;
					if (im < 0)im = 0; if (jp >= image.rows)jp = image.rows - 1;
					if (im < 0)im = 0; if (jp >= image.cols)jp = image.cols - 1;
					if (image.at<uchar>(i, jp) == 255) {
						image.at<uchar>(i, jp) = label; cnt++;
					}
					if (image.at<uchar>(im, jp) == 255) {
						image.at<uchar>(im, jp) = label; cnt++;
					}
					if (image.at<uchar>(im, j) == 255) {
						image.at<uchar>(im, j) = label; cnt++;
					}
					if (image.at<uchar>(im, jm) == 255) {
						image.at<uchar>(im, jm) = label; cnt++;
					}
					if (image.at<uchar>(i, jm) == 255) {
						image.at<uchar>(i, jm) = label; cnt++;
					}
					if (image.at<uchar>(ip, jm) == 255) {
						image.at<uchar>(ip, jm) = label; cnt++;
					}
					if (image.at<uchar>(ip, j) == 255) {
						image.at<uchar>(ip, j) = label; cnt++;
					}
					if (image.at<uchar>(ip, jp) == 255) {
						image.at<uchar>(ip, jp) = label; cnt++;
					}
				}
		if (cnt == 0)break;
	}
}
#define PI8 (double)3.14159265
#define ROOT2 (double)1.41421356
double calc_size(Mat image_label, int label, int *cx, int *cy);
double calc_length(Mat image_label, int label);
double trace(Mat image_label, int xs, int ys);
void features(Mat image_label_in, Mat image_label_out, int cnt, double size[], double ratio[], char *buf)
{
	int i, j, center_x, center_y;
	double l;
	int posi, m;

	posi = 0;
	for (i = 0; i < image_label_in.rows; i++)
		for (j = 0; j < image_label_in.cols; j++)
			image_label_out.at<uchar>(i, j) = image_label_in.at<uchar>(i, j);
	m = printf(buf, " no area circum round grav(x,y)\n");
	posi += m;
	for (i = 0; i < cnt; i++) {
		size[i] = calc_size(image_label_out, i + L_BASE, &center_x, &center_y);
		l = calc_length(image_label_out, i + L_BASE);
		ratio[i] = 4 * PI8 * size[i] / (l*l);
		image_label_out.at<uchar>(center_y, center_x) = 255;	/*­«¤ß*/
		m = printf(&buf[posi], "%3d %6d %8.2f %8.4f (%3d,%3d)\n", i, (int)size[i], l, ratio[i], center_x, center_y);
		posi += m;
	}
}
double calc_size(Mat image_label, int label, int *cx, int *cy)
{
	int i, j;
	double tx, ty, total;

	tx = 0; ty = 0; total = 0;
	for (i = 0; i < image_label.rows; i++)
		for (j = 0; j < image_label.cols; j++)
			if (image_label.at<uchar>(i, j) == label) {
				tx += j; ty += i; total++;
			}
	if (total == 0.0)return 0.0;
	*cx = (int)(tx / total); *cy = (int)(ty / total);
	return total;
}
double calc_length(Mat image_label, int label)
{
	int i, j;
	double tarce();
	for (i = 0; i < image_label.rows; i++)
		for (j = 0; j < image_label.cols; j++)
			if (image_label.at<uchar>(i, j) == label)	return trace(image_label, j - 1, i);
	return 0;
}
double trace(Mat image_label, int xs, int ys)
{
	int x, y, no, vec;
	double l;
	l = 0; x = xs; y = ys; no = image_label.at<uchar>(y, x + 1); vec = 5;
	for (;;) {
		if (x == xs && y == ys && l != 0)return 1;
		image_label.at<uchar>(y, x) = 255;
		switch (vec) {
		case 3:
			if (image_label.at<uchar>(y, x + 1) != no && image_label.at<uchar>(y - 1, x + 1) == no)
			{
				x = x + 1; y = y; l++; vec = 0; continue;
			}
		case 4:
			if (image_label.at<uchar>(y - 1, x + 1) != no && image_label.at<uchar>(y - 1, x) == no)
			{
				x = x + 1; y = y - 1; l += ROOT2; vec = 1; continue;
			}
		case 5:
			if (image_label.at<uchar>(y - 1, x) != no && image_label.at<uchar>(y - 1, x - 1) == no)
			{
				x = x; y = y - 1; l++; vec = 2; continue;
			}
		case 6:
			if (image_label.at<uchar>(y - 1, x - 1) != no && image_label.at<uchar>(y, x - 1) == no)
			{
				x = x - 1; y = y - 1; l += ROOT2; vec = 3; continue;
			}
		case 7:
			if (image_label.at<uchar>(y, x - 1) != no && image_label.at<uchar>(y + 1, x - 1) == no)
			{
				x = x - 1; y = y; l++; vec = 4; continue;
			}
		case 0:
			if (image_label.at<uchar>(y + 1, x - 1) != no && image_label.at<uchar>(y + 1, x) == no)
			{
				x = x - 1; y = y + 1; l += ROOT2; vec = 5; continue;
			}
		case 1:
			if (image_label.at<uchar>(y + 1, x) != no && image_label.at<uchar>(y + 1, x + 1) == no)
			{
				x = x; y = y + 1; l++; vec = 6; continue;
			}
		case 2:
			if (image_label.at<uchar>(y + 1, x + 1) != no && image_label.at<uchar>(y, x + 1) == no)
			{
				x = x + 1; y = y + 1; l += ROOT2; vec = 7; continue;
			}
			vec = 3;
		}
	}
}
void extract_ratio(Mat image_label_in, Mat image_label_out, int cnt, double ratio[], double ratio_min, double ratio_max)
{
	int i, j, x, y;
	int lno[256];
	for (i = 0, j = 0; i < cnt; i++)
		if (ratio[i] >= ratio_min && ratio[i] <= ratio_max)
			lno[j++] = L_BASE + i;
	for (y = 0; y < image_label_in.rows; y++) {
		for (x = 0; x < image_label_in.cols; x++) {
			image_label_out.at<uchar>(y, x) = 0;
			for (i = 0; i < j; i++)
				if (image_label_in.at<uchar>(y, x) == lno[i])
					image_label_out.at<uchar>(y, x) = image_label_in.at<uchar>(y, x);
		}
	}
}
void masking(Mat image_in, Mat image_out, Mat image_mask)
{
	int i, j;

	for (i = 0; i < image_in.rows; i++) {
		for (j = 0; j < image_in.cols; j++) {
			if (image_mask.at<uchar>(i, j) == 255)		image_out.at<uchar>(i, j) = image_in.at<uchar>(i, j);
			else										image_out.at<uchar>(i, j) = 0;
		}
	}
}
void extract_size(Mat image_label_in, Mat image_label_out, int cnt, double size[], double size_min, double size_max)
{
	int i, j, x, y;
	int lno[256];

	for (i = 0, j = 0; i < cnt; i++)
		if (size[i] >= size_min && size[i] <= size_max)		lno[j++] = L_BASE + i;
	for (y = 0; y < image_label_in.rows; y++) {
		for (x = 0; x < image_label_in.cols; x++) {
			image_label_out.at<uchar>(y, x) = 0;
			for (i = 0; i < j; i++)
				if (image_label_in.at<uchar>(y, x) == lno[i])
					image_label_out.at<uchar>(y, x) = image_label_in.at<uchar>(y, x);
		}
	}
}
