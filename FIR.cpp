// FIR.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include <windows.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include<iomanip>


using namespace std;
int m;
//定义最大序列长度
 const int N = 2048;
int n;
const float PI =3.1415926 ;
/////快速傅里叶变换部分
inline void swap(float &a, float &b)
{
	float t;
	t = a;
	a = b;
	b = t;
}

void bitrp(float xreal[], float ximag[], int n)
{
	// 位反转置换
	int i, j, a, b, p;

	for (i = 1, p = 0; i < n; i *= 2)
	{
		p++;
	}
	for (i = 0; i < n; i++)
	{
		a = i;
		b = 0;
		for (j = 0; j < p; j++)
		{
			b = (b << 1) + (a & 1);    // b = b * 2 + a % 2;
			a >>= 1;        // a = a / 2;
		}
		if (b > i)
		{
			swap(xreal[i], xreal[b]);
			swap(ximag[i], ximag[b]);
		}
	}
}

void FFT(float xreal[], float ximag[], int n)
{
	// 快速傅立叶变换，将复数 x 变换后仍保存在 x 中，xreal, ximag 分别是 x 的实部和虚部
	float wreal[N / 2], wimag[N / 2], treal, timag, ureal, uimag, arg;
	int m, k, j, t, gIn, hIn;

	bitrp(xreal, ximag, n);

	// 计算 前 n / 2 个 n 次方根的共轭复数 W'j = wreal [j] + i * wimag [j] , j = 0, 1, ... , n / 2 - 1
	arg = -2 * PI / n;
	treal = cos(arg);
	timag = sin(arg);
	wreal[0] = 1.0;
	wimag[0] = 0.0;
	for (j = 1; j < n / 2; j++)
	{
		wreal[j] = wreal[j - 1] * treal - wimag[j - 1] * timag;
		wimag[j] = wreal[j - 1] * timag + wimag[j - 1] * treal;
	}

	for (m = 2; m <= n; m *= 2)
	{
		for (k = 0; k < n; k += m)
		{
			for (j = 0; j < m / 2; j++)
			{
				gIn = k + j;
				hIn = gIn + m / 2;
				t = n * j / m;    // 旋转因子 w 的实部在 wreal [] 中的下标为 t
				treal = wreal[t] * xreal[hIn] - wimag[t] * ximag[hIn];
				timag = wreal[t] * ximag[hIn] + wimag[t] * xreal[hIn];
				ureal = xreal[gIn];
				uimag = ximag[gIn];
				xreal[gIn] = ureal + treal;
				ximag[gIn] = uimag + timag;
				xreal[hIn] = ureal - treal;
				ximag[hIn] = uimag - timag;
			}
		}
	}
}

void  IFFT(float xreal[], float ximag[], int n)
{
	// 快速傅立叶反变换
	float wreal[N / 2], wimag[N / 2], treal, timag, ureal, uimag, arg;
	int m, k, j, t, gIn, hIn;

	bitrp(xreal, ximag, n);

	// 计算 1 的前 n / 2 个 n 次方根 Wj = wreal [j] + i * wimag [j] , j = 0, 1, ... , n / 2 - 1
	arg = 2 * PI / n;
	treal = cos(arg);
	timag = sin(arg);
	wreal[0] = 1.0;
	wimag[0] = 0.0;
	for (j = 1; j < n / 2; j++)
	{
		wreal[j] = wreal[j - 1] * treal - wimag[j - 1] * timag;
		wimag[j] = wreal[j - 1] * timag + wimag[j - 1] * treal;
	}

	for (m = 2; m <= n; m *= 2)
	{
		for (k = 0; k < n; k += m)
		{
			for (j = 0; j < m / 2; j++)
			{
				gIn = k + j;
				hIn = gIn + m / 2;
				t = n * j / m;    // 旋转因子 w 的实部在 wreal [] 中的下标为 t
				treal = wreal[t] * xreal[hIn] - wimag[t] * ximag[hIn];
				timag = wreal[t] * ximag[hIn] + wimag[t] * xreal[hIn];
				ureal = xreal[gIn];
				uimag = ximag[gIn];
				xreal[gIn] = ureal + treal;
				ximag[gIn] = uimag + timag;
				xreal[hIn] = ureal - treal;
				ximag[hIn] = uimag - timag;
			}
		}
	}

	for (j = 0; j < n; j++)
	{
		xreal[j] /= n;
		ximag[j] /= n;
	}
}

int main()
{
	//序列初始化
	float X1real[2048] = {0};
	float X1image[2048] = {0};
	float X2real[2048] = { 0 };
	float X2image[2048] = { 0 };
	float Hreal[2048] = {0};
	float Himage[2048] = { 0 };
	float HDreal[2048] = { 0 };
	float HDimage[2048] = { 0 };
	float Wreal[2048] = { 0 };
	float Wimage[2048] = { 0 };
	float Y1real[2048] = { 0 };
	float Y1image[2048] = { 0 };
	float Y2real[2048] = { 0 };
	float Y2image[2048] = { 0 };
	//序列（有值的）的长度
	int N1 =64;
	
	//截止频率
	float wp = 0.3*PI;
	float wst = 0.5*PI;
	float wc = (wp+ wst)*0.5;
	//查表得布莱克窗的过滤带宽与N之间的关系
	int N2 = 10 * PI / (wst - wp);
	float alpha = (N2 - 1) / 2;
	

	for (int i = 0;i < N1;i++)
	{
		X1real[i] = sin(0.275*PI*i )+ cos(0.655*PI*i);
		X2real[i]= cos(0.655*PI*i );
	}
	FFT(X1real, X1image, 128);
	FFT(X2real, X2image, 128);
	//BLACKMAN
	for (int i = 0;i < N2;i++)
	{
		/*Wreal[i] = 0.42 - 0.5*cos(2 * PI*i / (N2 - 1)) + 0.08*cos(4 * PI*i / (N2 - 1));*/
		Wreal[i] = 0.54 - 0.46*cos(2 * PI*i / (N2 - 1));
		if (i == alpha)
		{
			HDreal[i] = wc / PI;
		} else 
		HDreal[i] = sin((wc*(i - alpha)) )/ (PI*(i - alpha));
		Hreal[i] = Wreal[i] * HDreal[i];
	}
	
	FFT(Hreal, Himage, 128);
	for (int i = 0;i < N1 + N2 - 1;i++)
	{
		Y1real[i] = X1real[i] * Hreal[i] - X1image[i] * Himage[i];
		Y1image[i] = X1real[i] * Himage[i] + Hreal[i] * X1image[i];
		Y2real[i] = X2real[i] * Hreal[i] - X2image[i] * Himage[i];
		Y2image[i] = X2real[i] * Himage[i] + Hreal[i] * X2image[i];
	}
	IFFT(X1real, X1image, 256);
	IFFT(Y1real, Y1image, 256);
	IFFT(X2real, X2image, 256);
	IFFT(Y2real, Y2image, 256);

	for (int i = 0;i < 256;i++)
	{
		
		cout <<setw(15)<<X1real[i] << setw(15) << Y1real[i] << setw(15) << X2real[i] << setw(15) << Y2real[i] <<endl;
	}

	cin >> n;

    return 0;
}

