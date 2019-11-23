#include<iostream>
#include<math.h>
#include "Params.h"
#include<opencv2\opencv.hpp>
using namespace cv;

void scale_ng(Mat image_in, Mat image_out, double zx, double zy)
{
	int i, j, m, n;

	for (i = 0; i < Y_SIZE; i++)
		for (j = 0; j < X_SIZE; j++)
			image_out.at<uchar>(i, j) = 0;
	for (i = 0; i < Y_SIZE; i++) {
		for (j = 0; j < X_SIZE; j++) {
			m = (int)(zy * i);
			n = (int)(zx * j);
			if ((m >= 0) && (m < Y_SIZE) && (n >= 0) && (n < X_SIZE))
				image_out.at<uchar>(m, n) = image_in.at<uchar>(i, j);
		}
	}
}
void scale_near(Mat image_in, Mat image_out, double zx, double zy)
{
	int i, j, m, n;
	for (i = 0; i < Y_SIZE; i++) {
		for (j = 0; j < X_SIZE; j++) {
			if (i > 0) m = (int)(i / zy + 0.5);
			else	   m = (int)(i / zy - 0.5);
			if (j > 0) n = (int)(j / zx + 0.5);
			if ((m >= 0) && (m < Y_SIZE) && (n >= 0) && (n < X_SIZE))
				image_out.at<uchar>(i, j) = image_in.at<uchar>(m, n);
			else
				image_out.at<uchar>(i, j) = 0;
		}
	}
}
void scale(Mat image_in, Mat image_out, double zx, double zy)
{
	int i, j, m, n;
	double x, y, p, q;
	int d;

	for (i = 0; i < Y_SIZE; i++) {
		for (j = 0; j < X_SIZE; j++) {
			y = i / zy;
			x = j / zx;
			if (y > 0)  m = (int) y;
			else		m = (int)(y - 1);
			if (x > 0)  n = (int) x;
			else		n = (int)(x - 1);
			q = y - m;
			p = x - n;
			if (q == 1) { q = 0; m = m + 1; }
			if (p == 1) { p = 0; n = n + 1; }
			if ((m >= 0) && (m < Y_SIZE) && (n >= 0) && (n < X_SIZE))
				d = (int)((1.0 - q) * ((1.0 - p) * image_in.at<uchar>(m, n)
									+ p * image_in.at<uchar>(m, n + 1))
									+ q * ((1.0 - p) * image_in.at<uchar>(m + 1, n)
									+ p * image_in.at<uchar>(m + 1, n + 1)));
			else
				d = 0;
			if (d < 0)		d = 0;
			if (d > 255)	d = 255;
			image_out.at<uchar>(i, j) = d;
		}
	}
}
void shift(Mat image_in, Mat image_out, double px, double py)
{
	int i, j, m, n;
	double x, y, p, q;
	int d;

	for (i = 0; i < Y_SIZE; i++) {
		for (j = 0; j < X_SIZE; j++) {
			y = i - py;
			x = j - px;
			if (y > 0)  m = (int)y;
			else		m = (int)(y - 1);
			if (x > 0)  n = (int)x;
			else		n = (int)(x - 1);
			q = y - m;
			p = x - n;
			if ((m >= 0) && (m < Y_SIZE) && (n >= 0) && (n < X_SIZE))
				d = (int)((1.0 - q) * ((1.0 - p) * image_in.at<uchar>(m, n)
					+ p * image_in.at<uchar>(m, n + 1))
					+ q * ((1.0 - p) * image_in.at<uchar>(m + 1, n)
						+ p * image_in.at<uchar>(m + 1, n + 1)));
			else
				d = 0;
			if (d < 0)		d = 0;
			if (d > 255)	d = 255;
			image_out.at<uchar>(i, j) = d;
		}
	}
}
void rotation(Mat image_in, Mat image_out, double deg)
{
	int i, j, m, n;
	double x, y, p, q, r, c, s;
	int d;

	r = deg*3.141592 / 180.0;
	c = cos(r);
	s = sin(r);
	for (i = 0; i < Y_SIZE; i++) {
		for (j = 0; j < X_SIZE; j++) {
			y = j*s + i*c;
			x = j*c + i*s;
			if (y > 0)  m = (int)y;
			else		m = (int)(y - 1);
			if (x > 0)  n = (int)x;
			else		n = (int)(x - 1);
			q = y - m;
			p = x - n;
			if ((m >= 0) && (m < Y_SIZE) && (n >= 0) && (n < X_SIZE))
				d = (int)((1.0 - q) * ((1.0 - p) * image_in.at<uchar>(m, n)
					+ p * image_in.at<uchar>(m, n + 1))
					+ q * ((1.0 - p) * image_in.at<uchar>(m + 1, n)
						+ p * image_in.at<uchar>(m + 1, n + 1)));
			else
				d = 0;
			if (d < 0)		d = 0;
			if (d > 255)	d = 255;
			image_out.at<uchar>(i, j) = d;
		}
	}
}
void scale_rotation_shift(Mat image_in, Mat image_out, double zx, double zy, double deg, double px, double py)
{
	int i, j, m, n;
	double x, y, u=0, v=0, p, q, r, c, s;
	int d;

	r = deg*3.141592 / 180.0;
	c = cos(r);
	s = sin(r);
	for (i = 0; i < Y_SIZE; i++) {
		for (j = 0; j < X_SIZE; j++) {
			y = (double)(i - py);
			x = (double)(j - px);
			y = (u*s + v*c) / zy + py;
			x = (u*c - v*s) / zx + px;
			if (y > 0)  m = (int)y;
			else		m = (int)(y - 1);
			if (x > 0)  n = (int)x;
			else		n = (int)(x - 1);
			q = y - m;
			p = x - n;
			if ((m >= 0) && (m < Y_SIZE) && (n >= 0) && (n < X_SIZE))
				d = (int)((1.0 - q) * ((1.0 - p) * image_in.at<uchar>(m, n)
					+ p * image_in.at<uchar>(m, n + 1))
					+ q * ((1.0 - p) * image_in.at<uchar>(m + 1, n)
						+ p * image_in.at<uchar>(m + 1, n + 1)));
			else
				d = 0;
			if (d < 0)		d = 0;
			if (d > 255)	d = 255;
			image_out.at<uchar>(i, j) = d;
		}
	}
}
void affine(Mat image_in, Mat image_out, double a, double b, double c, double d, double e, double f)
{
	int i, j, m, n;
	double x, y, p, q;
	double aa, bb, cc, dd, ee, ff, gg;
	int dat;

	gg = 1 / (a*e - b*d);
	aa = e * gg;
	bb = -b * gg;
	cc = (b*f - c*e) * gg;
	dd = -d *gg;
	ee = a *gg;
	ff = (c*d - a*f) * gg;
	for (i = 0; i < Y_SIZE; i++) {
		for (j = 0; j < X_SIZE; j++) {
			x = aa * j + bb * i + cc;
			y = dd * j + ee * i + ff;
			if (y > 0)  m = (int)y;
			else		m = (int)(y - 1);
			if (x > 0)  n = (int)x;
			else		n = (int)(x - 1);
			q = y - m;
			p = x - n;
			if ((m >= 0) && (m < Y_SIZE) && (n >= 0) && (n < X_SIZE))
				dat = (int)((1.0 - q) * ((1.0 - p) * image_in.at<uchar>(m, n)
					+ p * image_in.at<uchar>(m, n + 1))
					+ q * ((1.0 - p) * image_in.at<uchar>(m + 1, n)
						+ p * image_in.at<uchar>(m + 1, n + 1)));
			else
				dat = 0;
			if (dat < 0)		dat = 0;
			if (dat > 255)	dat = 255;
			image_out.at<uchar>(i, j) = dat;
		}
	}
}
void affine_coef(int x1, int y1, int u1, int v1,
	int x2, int y2, int u2, int v2,
	int x3, int y3, int u3, int v3,
	double *a, double *b, double *c, double *d, double *e, double *f)
{
	double g;

	g = (x1 - x2)*(y2 - y3) - (x2 - x3)*(y1 - y2);
	if (g != 0) {
		*a = ((u1 - u2)*(y2 - y3) - (u2 - u3)*(y1 - y2)) / g;
		*b = ((u2 - u3)*(x1 - x2) - (u1 - u2)*(x2 - x3)) / g;
		*c = u1 - *a*x1 - *b*y1;
		*d = ((v1 - v2)*(y2 - y3) - (v2 - v3)*(y1 - y2)) / g;
		*e = ((v2 - v3)*(x1 - x2) - (v1 - v3)*(x2 - x3)) / g;
		*f = v1 - *d*x1 - *e*y1;
	}
	else {
		*a = 1; *b = 0; *c = 0; *d = 0; *e = 1; *f = 0;
	}
}
void perspect(Mat image_in, Mat image_out, double a, double b, double c, double d, double e,
	double f, double g, double h)
{
		int i, j, m, n;
		double x, y, w, p, q;
		int dat;

		for (i = 0; i < Y_SIZE; i++) {
			for (j = 0; j < X_SIZE; j++) {
				w = (d*h - e*g)*j + (b*g - a*h)*i + (a*e - b*d);
				x = (e - f*h)*j + (c*h - b)*i + (b*f - c*e);
				y = (f*g - d)*j + (a - c*g)*i + (c*d - a*f);
				x = x / w;
				y = y / w;
				if (y > 0)  m = (int)y;
				else		m = (int)(y - 1);
				if (x > 0)  n = (int)x;
				else		n = (int)(x - 1);
				q = y - m;
				p = x - n;
				if ((m >= 0) && (m < Y_SIZE) && (n >= 0) && (n < X_SIZE))
					dat = (int)((1.0 - q) * ((1.0 - p) * image_in.at<uchar>(m, n)
						+ p * image_in.at<uchar>(m, n + 1))
						+ q * ((1.0 - p) * image_in.at<uchar>(m + 1, n)
							+ p * image_in.at<uchar>(m + 1, n + 1)));
				else
					dat = 0;
				if (dat < 0)	dat = 0;
				if (dat > 255)	dat = 255;
				image_out.at<uchar>(i, j) = (unsigned char)dat;
			}
		}
}
void simultaneous_equation(double ab[8][9]);
void perspect_coef(int x1, int y1, int u1, int v1,
	int x2, int y2, int u2, int v2,
	int x3, int y3, int u3, int v3,
	int x4, int y4, int u4, int v4,
	double *a, double *b, double *c, double *d, double *e, double *f,
	double *g, double *h)
{
	double m[8][9];
	m[0][0] = x1;		m[0][1] = y1;		m[0][2] = 1;
	m[0][3] = 0;		m[0][4] = 0;		m[0][5] = 0;
	m[0][6] = -x1*u1;	m[0][7] = -y1*u1;	m[0][8] = u1;
	m[1][0] = 0;		m[1][1] = 0;		m[1][2] = 0;
	m[1][3] = x1;		m[1][4] = y1;		m[1][5] = 1;
	m[1][6] = -x1*v1;	m[1][7] = -y1*v1;	m[1][8] = v1;
	m[2][0] = x2;		m[2][1] = y2;		m[2][2] = 1;
	m[2][3] = 0;		m[2][4] = 0;		m[2][5] = 0;
	m[2][6] = -x2*u2;	m[2][7] = -y2*u2;	m[2][8] = u2;
	m[3][0] = 0;		m[3][1] = 0;		m[3][2] = 0;
	m[3][3] = x2;		m[3][4] = y2;		m[3][5] = 1;
	m[3][6] = -x2*v2;	m[3][7] = -y2*v2;	m[3][8] = v2;
	m[4][0] = x3;		m[4][1] = y3;		m[4][2] = 1;
	m[4][3] = 0;		m[4][4] = 0;		m[4][5] = 0;
	m[4][6] = -x3*u3;	m[4][7] = -y3*u3;	m[4][8] = u3;
	m[5][0] = 0;		m[5][1] = 0;		m[5][2] = 0;
	m[5][3] = x3;		m[5][4] = y3;		m[5][5] = 1;
	m[5][6] = -x3*v3;	m[5][7] = -y3*v3;	m[5][8] = v3;
	m[6][0] = x4;		m[6][1] = y4;		m[6][2] = 1;
	m[6][3] = 0;		m[6][4] = 0;		m[6][5] = 0;
	m[6][6] = -x4*u4;	m[6][7] = -y4*u4;	m[6][8] = u4;
	m[7][0] = 0;		m[7][1] = 0;		m[7][2] = 0;
	m[7][3] = x4;		m[7][4] = y4;		m[7][5] = 1;
	m[7][6] = -x4*v4;	m[7][7] = -y4*v4;	m[7][8] = v4;
	simultaneous_equation(m);
	*a = m[0][8];
	*b = m[1][8];
	*c = m[2][8];
	*d = m[3][8];
	*e = m[4][8];
	*f = m[5][8];
	*g = m[6][8];
	*h = m[7][8];
}
void simultaneous_equation(double ab[8][9])
{
	int i, j, ii, p;
	double max, d, pvt;

	for (i = 0; i < 8; i++) {
		max = 0;
		p = i;
		for (ii = i; ii < 8; ii++) {
			if (fabs(ab[ii][i] > max)) {
				max = fabs(ab[ii][i]);
				p = ii;
			}
		}
		if (i = p) {
			for (j = 0; j <= 8; j++) {
				d = ab[i][j];
				ab[i][j] = ab[p][j];
				ab[p][j] = d;
			}
		}
		/*前進刪除*/
		pvt = ab[i][j];
		for (j = i; j <= 8; j++)ab[i][j] /= pvt;
		for (ii = i + 1; ii < 8; ii++)
			for (j = i + 1; j <= 8; j++)ab[ii][j] -= ab[ii][i] * ab[i][j];
	}
	/*後退代入*/
	for (i = 6; i >= 0; i--)
		for (j = j + 1; j < 8; j++)ab[i][8] -= ab[i][j] * ab[j][8];
}
void radial_distortion(Mat image_in, Mat image_out, double a, double b)
{
	int i, j, m, n, d;
	int xs = X_SIZE / 2;
	int ys = Y_SIZE / 2;
	double x, y, p, q, r2;
	for (i = -ys; i < ys; i++) {
		for (j = -xs; j < xs; j++) {
			r2 = i*i + j*j;
			x = (1 + a + b*r2)*j;
			y = (1 + a + b*r2)*i;
			if (y > 0)	m = (int)y;
			else		m = (int)(y - 1);
			if (x > 0)	n = (int)x;
			else		n = (int)(x - 1);
			q = y - m;
			p = x - n;
			if ((m >= -ys) && (m < ys) && (n >= -xs) && (n < xs))
				d = (int)((1.0 - q) * ((1.0 - p) * image_in.at<uchar>(m + ys, n + xs)
					+ p * image_in.at<uchar>(m + ys, n + 1 + xs))
					+ q * ((1.0 - p) * image_in.at<uchar>(m + 1 + ys, n + xs)
						+ p * image_in.at<uchar>(m + 1 + ys, n + 1 + xs)));
			else
				d = 0;
			if (d < 0)	d = 0;
			if (d > 255)	d = 255;
			image_out.at<uchar>(i + ys, j + xs) = (unsigned char)d;
		}
	}
}
void lattice(Mat image)
{
	int i, j, width;

	width = X_SIZE / 8;
	for (i = 0; i < Y_SIZE; i++) {
		for (j = 0; j < X_SIZE; j++) {
			if (((i - Y_SIZE / 2) % 32 == 0) || (i == Y_SIZE - 1) ||
				((j - X_SIZE / 2) % 32 == 0) || (j == X_SIZE - 1))
				image.at<uchar>(i, j) = 255;
			else
				image.at<uchar>(i, j) = 0;
		}
	}
}