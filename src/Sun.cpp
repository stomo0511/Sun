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

	latitude *= M_PI/180.0;  // [radian]

//	omega*J=2π(dn +0.5)/365
	double th = 2.0 * M_PI * (double)(GetDays( year, month, day) +0.5) / (365.0 + leap(year));

// 赤緯 delta [度]
	double bt[7] = { 0.33281, -22.984, -0.34990, -0.13980, 3.7872, 0.03250, 0.07187 };
	double delta = bt[0] + bt[1]*cos(th) + bt[2]*cos(2.0*th) + bt[3]*cos(3.0*th) + bt[4]*sin(th) + bt[5]*sin(2.0*th) + bt[6]*sin(3.0*th);
	cout << "太陽赤緯 = " << delta << " [度]\n";
	delta *= M_PI/180.0; // [radian]

// 均時差 eq [時]
	double ap[6] = { 0.0072, -0.0528, -0.0012, -0.1229, -0.1565, -0.0041 };
	double eq = ap[0]*cos(th) + ap[1]*cos(2.0*th) + ap[2]*cos(3.0*th) + ap[3]*sin(th) + ap[4]*sin(2.0*th) + ap[5]*sin(3.0*th);

	cout << "均時差 = " << eq*60.0 << " [分]\n";

//	日の出時刻：　t1　[単位： 時]
//	　　t0 = acos(-tan(δ)tan(φ))
//	　　T1 = (-t0 + 180)/15
//	　　t1 = T1 - (θ - 135)/15 - e
	double t0 = acos( -tan(delta)*tan(latitude) ); // [radian]
	t0 *= 180.0/M_PI;  // [度]
	cout << "（時角: " << t0 << "）\n";
	double T1 = (-t0 + 180.0)/15.0;
	double t1 = T1 - (longitude - 135.0)/15.0 - eq;

	cout << "日の出:" << floor(t1) << "時" << floor( (t1 - floor(t1))*60.0 ) << "分, ";

	double sA = cos(delta)*sin(t1-9.0);
	double cA = -sin(delta)/cos(latitude);
	double A = atan2(sA,cA) + M_PI;
	cout << "日の出方位角：" << A*180.0/M_PI << endl;

//	日の入時刻：　t2　[単位： 時]
//	　　t0 = acos(-tan(δ)tan(φ))
//	　　T2 = ( t0 + 180)/15
//	　　t2 = T2 - (θ - 135)/15 - e
	double T2 = (t0 + 180.0)/15.0;
	double t2 = T2 - (longitude - 135.0)/15.0 - eq;
	cout << "日の入:" << floor(t2) << "時" << floor( (t2 - floor(t2))*60.0 ) << "分, ";

	sA = cos(delta)*sin(t2-9.0);
	A = atan2(sA,cA) + M_PI;
	cout << "日の入方位角：" << A*180.0/M_PI << endl << endl;

	// 日の出から日の入りまでのループ
	int sr = (int)(floor(t1));
	int ss = (int)(floor(t2));

	for (int i=sr; i<=ss+1; i++)
	{
		cout << "時間: " << i << " [時], ";

// 時角
//		T = Ts + (θ - 135)/15 + e
//		　　t = 15T - 180
//		　　ここで、
//		　　　Ts：　時刻（中央標準時）
//		　　　θ：　東経
//		double T = ((double)(i) - 9.0) + (longitude - 135.0)/15.0 + eq;
		double T = ((double)(i) ) + (longitude - 135.0)/15.0 + eq;
		double t = 15.0*T - 180.0;
		cout << "時角 = " << t << ", ";
		t *= M_PI/180.0;           // [radian]

// 太陽高度 h [度]
//		h = asin(sin(φ)sin(δ) + cos(φ)cos(δ)cos(t))
		double h = asin( sin(latitude)*sin(delta) + cos(latitude)*cos(delta)*cos(t) ); // [radian]
		cout << "高度 = " << h*180.0/M_PI << ", ";
		h *= 180.0/M_PI;

//		方位角：　A（北 = 0, 東 = 90, 南 = 180, 西 = 270°）
//		　　sinA = cos(δ)sin(t)/cos(h)
//		　　cosA = (sin(h)sin(φ) - sin(δ))/cos(h)/cos(φ)
//		　　A = atan2(sinA, cosA) + π
//		double sA = cos(delta)*sin(t) / cos(h);
		double sA = cos(delta)*sin((double)(i) -9.0) / cos(h);
		double cA = (sin(h)*sin(latitude) - sin(delta)) / cos(h) / cos(latitude);
		double A = atan2(sA,cA) + M_PI;

//		double A;
//		if (sA > 0)
//			A = M_PI/2.0 - atan2(cA,sA);  // [radian]
//		else
//			A = -M_PI/2.0 - atan2(cA,sA);
		cout << "方位角 = " << A*180.0/M_PI << endl;
	}
}
