//============================================================================
// Name        : Sun.cpp
// Author      : Tomohiro Suzuki
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <cassert>

using namespace std;

// 閏年判定：閏年なら1、それ以外0
int leap( const int year )
{
	return (1 / (year % 4 + 1)) * (1 - 1 / (year % 100 + 1)) + (1 / (year % 400 + 1));
}

// 1月1日からの通日数
int GetDays( int y, int m, int d )
{
	// 各月の日数
	int days[12] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
	// 閏年なら2月を29日に
	days[1] += leap(y);

	int t = 0;
	for (int i=0; i<m-1; i++)
	{
		t += days[i];
	}
	return t + d;
}

int main( int argc, char* argv[] )
{
	if (argc < 4)
	{
		cerr << "Usage: a.out [year] [month] [day]\n";
		exit( EXIT_FAILURE );
	}
	const int year  = atoi(argv[1]);
	const int month = atoi(argv[2]);
	const int day   = atoi(argv[3]);

	// 三澤農場のセンサ位置
	const double longitude = 138.455; // 東経
	const double latitude = 35.786;   // 北緯

	cout << "東経: " << longitude << ", 北緯: " << latitude << endl;

//	θo=2π(dn-1)/365
	double th = 2.0 * M_PI * (double)(GetDays( year, month, day) -1) / (365.0 + leap(year));

// 赤緯
//	δ=0.006918-0.399912cos(θo)+0.070257sin(θo)-0.006758cos(2θo)+0.000907sin(2θo)-0.002697cos(3θo)+0.001480sin(3θo)
	double bt[7] = { 0.006918, -0.399912, 0.070257, -0.006758, 0.000907, -0.002697, 0.001480 };
	double delta = bt[0] + bt[1]*cos(th) + bt[2]*sin(th) + bt[3]*cos(2.0*th) + bt[4]*sin(2.0*th) + bt[5]*cos(3.0*th) + bt[6]*sin(3.0*th);

	cout << "太陽赤緯 = " << delta*180.0/M_PI << " [度]\n";

// 均時差
//	Eq=0.000075+0.001868cos(θo)-0.032077sin(θo)-0.014615cos(2θo)-0.040849sin(2θo)
	double ap[5] = { 0.000075, 0.001868, -0.032077, -0.014615, -0.040849 };
	double eq = ap[0] + ap[1]*cos(th) + ap[2]*sin(th) + ap[3]*cos(2.0*th) + ap[4]*sin(2.0*th);
	eq *= 12.0 / M_PI;

	cout << "均時差 = " << eq*60.0 << endl;

//// 太陽距離
//// r/r*=1/{1.000110+0.034221cos(θo)+0.001280sin(θo)+0.000719cos(2θo)+0.000077sin(2θo)}^0.5
//	double gm[5] = { 1.000110, 0.034221, 0.001280, 0.000719, 0.000077 };
//	double rr = 1.0 / sqrt( gm[0] + gm[1]*cos(th) + gm[2]*sin(th) + gm[3]*cos(2.0*th) + gm[4]*sin(2.0*th) );
//
//	cout << "太陽距離 = " << rr << " [天文単位]\n";

// 時角
//	h=(JST-12)π/12+標準子午線からの経度差+均時差(Eq)
}
