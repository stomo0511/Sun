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
//	double longitude = 138.455; // 東経
//	double latitude = 35.786;   // 北緯
	// 甲府市
	double longitude = 138.56833; // 東経
	double latitude = 35.66231;   // 北緯
	cout << "東経: " << longitude << ", 北緯: " << latitude << endl;

	const double rad = M_PI / 180.0;
	const double dlt0 = -23.4393*rad;  // 北半球の冬至の日赤緯（rad）（=-23.4393）

	int n = year - 1968;  // 年 Eq.(27)
	double d0 = 3.71 + 0.2596*(double)(n) - floor((double)(n) + 3.0) / 4.0;  // 平均軌道上の近日点通過日（日） Eq.(26)
	double M = 0.9856*( (int)(GetDays( year, month, day)) - d0 );  // 平均近点離角（度） Eq.(25)
	double eps = 12.3901 + 0.0172*( (double)(n) + M/360.0 );  // 近日点と冬至点の角度（度） Eq.(21)
	double v = M + 1.914*sin( M*rad ) + 0.02*sin( 2.0*M*rad );  // 真近点離角（度）Eq.(22)
	double ve = ( v + eps )*rad;  // v+eps （rad）
	double Et = (M-v) -  atan( 0.043*sin(2.0*ve) / (1.0 - 0.043*cos(2.0*ve)) ) / rad;  // 均時差（度） Eq.(20) (and (23), (24))
	double sindlt = cos( ve )*sin(dlt0); // Eq.(19)
	assert( fabs(sindlt) < 1.0 );
	double cosdlt = sqrt( 1.0 - sindlt*sindlt );

	double phirad = latitude*rad;
	for (double tm=7.0; tm<17.0; tm+=1.0)
	{
		double Tm = tm - 9.0; // 7:00 JST
		double t = 15.0*(Tm - 12.0) + (longitude - 135.0) + Et;  // 時角（度）Eq.(1)
		double trad = t*rad;

		double sinh = sin(phirad)*sindlt + cos(phirad)*cosdlt*cos(trad); // Eq.(2)
		assert( fabs(sinh) < 1.0 );
		double cosh = sqrt( 1.0 - sinh*sinh );
		double sinA = cosdlt*sin(trad) / sinh;  // Eq.(3)
		assert( fabs(sinA) < 1.0 );
		double cosA = (sinh*sin(phirad) - sindlt) / (cosh*cos(phirad)); // Eq.(4)
		assert( fabs(cosA) < 1.0 );
		double A = (sinA > 0) ? 90.0 + atan( cosA/sinA )/rad : -90 + atan( cosA/sinA )/rad;  // Eq.(5), (6)

		cout << "時間:" << Tm+9.0 << ", 高度：" << asin(sinh)/rad << ", 方位：" << A << " [度]\n";
	}
}
