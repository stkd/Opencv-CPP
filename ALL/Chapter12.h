#include<iostream>
#include<math.h>
#include "Params.h"
#include<opencv2\opencv.hpp>
using namespace cv;

#define B_VAL 128
void dpcm1(Mat image_in, int line, short data_out[X_SIZE])
{
	int pred, j;
	for (j = 0; j < X_SIZE; j++) {
		if (j == 0)	pred = B_VAL;
		else		pred = (int)image_in.at<uchar>(line, j - 1);
		data_out[j] = (int)image_in.at<uchar>(line, j) - pred;
	}
}
void dpcm2(Mat image_in, int line, short data_out[X_SIZE])
{
	int pred, j;
	if (line == 0) {
		for (j = 0; j < X_SIZE; j++) {
			if (j == 0)	pred = B_VAL;
			else		pred = (image_in.at<uchar>(line, j - 1) + B_VAL) / 2;
			data_out[j] = (int)image_in.at<uchar>(line, j) - pred;
		}
	}
	else {
		for (j = 0; j < X_SIZE; j++) {
			if (j == 0)	pred = (B_VAL + image_in.at<uchar>(line - 1, j)) / 2;
			else		pred = (image_in.at<uchar>(line, j - 1) + image_in.at<uchar>(line - 1, j)) / 2;
			data_out[j] = (int)image_in.at<uchar>(line, j) - pred;
		}
	}
}
void idpcm1(short data_in[X_SIZE], int line, Mat image_out)
{
	int pred;
	int j;
	
	for (j = 0; j < X_SIZE; j++) {
		if (j == 0) pred = B_VAL;
		else		pred = (int)image_out.at<uchar>(line, j - 1);
		image_out.at<uchar>(line, j) = (uchar)(pred+(int)data_in[j]);
	}
}
void idpcm2(short data_in[X_SIZE], int line, Mat image_out)
{
	int pred, j;

	if (line == 0) {
		for (j = 0; j < X_SIZE; j++) {
			if (j == 0) pred = B_VAL;
			else		pred = (image_out.at<uchar>(line, j - 1) + B_VAL) / 2;
			image_out.at<uchar>(line, j) = pred + (int)data_in[j];
		}
	}
	else {
		for (j = 0; j < X_SIZE; j++) {
			if (j == 0) pred = (B_VAL + image_out.at<uchar>(line - 1, j)) / 2;
			else		pred = (image_out.at<uchar>(line, j - 1) + image_out.at<uchar>(line-1, j)) / 2;
			image_out.at<uchar>(line, j) = pred + (int)data_in[j];
		}
	}
}
#define BYTESIZE 8
#define LEN 4
int vlcode(short int data_in[], int no, char vlc_out[])
{
	int i;
	int st = 0;
	int num = 0;
	int dl = BYTESIZE / LEN - 1;
	int mask = (1 << LEN) - 1;
	int dt, ms;

	vlc_out[num] = '\0';
	for (i = 0; i < no; i++) {
		dt = data_in[i];
		do {
			ms = dt >= mask ? mask : dt;
			vlc_out[num] |= (ms << LEN*(dl - st));
			dt -= mask;	st++;
			if (st > dl) {
				st = 0;	num++; vlc_out[num] = '\0';
			}
		} while (dt >= 0);
	}
		if (st != 0) {
			for (i = (dl - st); i >= 0; i++) {
				vlc_out[num] |= ms;
				ms <<= LEN;
			}
			num++;
		}
		return num;
}
void ivlcode(char vlc_in[], int no, short int data_out[]) {
	int i, j, k;
	int ino = 0;
	int num = 0;
	int dl = BYTESIZE / LEN - 1;
	int mask = (1 << LEN) - 1;
	
	for (i = 0; i < no; i++)data_out[i] = 0;
	do {
		for (j = dl; j >= 0; j--) {
			k = vlc_in[ino] & (mask << (LEN*j));
			k >>= (LEN*j);
			data_out[num] += k;
			if (k != mask) num++;
			if (num >= no) break;
		}
		ino++;
	} while (num < no);
}
int event(short dt) {
	int ev;

	if (dt <= 0) ev = -2 * dt;
	else		 ev = 2 * dt - 1;
	return ev;
}
int ievent(short ev)
{
	int dt;
	if (ev % 2 == 0) dt = -ev / 2;
	else			 dt = (ev + 1) / 2;
	return dt;
}
#define DPCM dpcm1
#define MAX_LENG X_SIZE*4
int dpcm_vlcode(Mat image_in, uchar image_buf[])
{
	int i, j, leng;
	long ptr, size;
	char vlc[MAX_LENG];
	short int data[X_SIZE];

	size = (long)X_SIZE*Y_SIZE;
	for (i = 0; i < size; i++)image_buf[i] = 0;
	ptr = 0;
	for (i = 0; i < Y_SIZE; i++) {
		DPCM(image_in, i, data);
		for (j = 0; j < X_SIZE; j++)
			data[j] = event(data[j]);
		leng = vlcode(data, X_SIZE, vlc);
		image_buf[ptr] = (uchar)((leng >> 8) & 0x00ff);
		ptr++;
		if (ptr > size) return -1;
		for (j = 0; j < leng; j++) {
			image_buf[ptr] = vlc[j];
			ptr++;
			if (ptr > size)return -1;
		}
	}
	return 0;
}
#define IDPCM idpcm1
int idpcm_vlcode(uchar image_buf[], Mat image_out)
{
	int i, j, leng;
	long ptr, size;
	char vlc[MAX_LENG];
	short int data[X_SIZE];

	size = (long)X_SIZE*Y_SIZE;
	ptr = 0;
	for (i = 0; i < Y_SIZE; i++) {
		leng = (int)image_buf[ptr];
		ptr++;
		if (ptr > size)return -1;
		leng = (leng >> 8) | image_buf[ptr];
		ptr++;
		for (j = 0; j < leng; j++) {
			vlc[j] = image_buf[ptr];
			ptr++;
			if (ptr > size)return -1;
		}
	ivlcode(vlc, X_SIZE, data);
	for (j = 0; j < X_SIZE; j++)
		data[j] = ievent(data[j]);
	IDPCM(data, i, image_out);
	}
	return 0;
}