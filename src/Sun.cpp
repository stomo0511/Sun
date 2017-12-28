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

	cout << "均時差 = " << eq*60.0 << " [分]\n";

//// 太陽距離
//// r/r*=1/{1.000110+0.034221cos(θo)+0.001280sin(θo)+0.000719cos(2θo)+0.000077sin(2θo)}^0.5
//	double gm[5] = { 1.000110, 0.034221, 0.001280, 0.000719, 0.000077 };
//	double rr = 1.0 / sqrt( gm[0] + gm[1]*cos(th) + gm[2]*sin(th) + gm[3]*cos(2.0*th) + gm[4]*sin(2.0*th) );
//
//	cout << "太陽距離 = " << rr << " [天文単位]\n";

	//	日の出時刻：　t1　[単位： 時]
	//	　　t = acos(-tan(δ)tan(φ))
	//	　　T1 = (-t + 180)/15
	//	　　t1 = T1 - (θ - 135)/15 - e
	double t = acos( -tan(delta)*tan(latitude*M_PI/180.0) ); // [radian]
	t *= 180.0/M_PI;  // [度]
	double T1 = (-t + 180.0)/15.0;
	double t1 = T1 - (longitude - 135.0)/15.0 - eq;

	double m;
	cout << "日の出:" << floor(t1) << "時";
	m = t1 - floor(t1);
	cout << floor(m * 60.0) << "分, ";

	double sA = cos(delta)*sin(t1-9.0);
	double cA = -sin(delta)/cos(latitude*M_PI/180.0);
	double A = atan2(sA,cA) + M_PI;
	cout << "日の出方位角：" << A*180.0/M_PI << endl;

	//	日の入時刻：　t2　[単位： 時]
	//	　　T2 = ( t + 180)/15
	//	　　t2 = T2 - (θ - 135)/15 - e
	double T2 = (t + 180.0)/15.0;
	double t2 = T2 - (longitude - 135.0)/15.0 - eq;
	cout << "日の入:" << floor(t2) << "時";
	m = t2 - floor(t2);
	cout << floor(m * 60.0) << "分, ";

	sA = cos(delta)*sin(t2-9.0);
	A = atan2(sA,cA) + M_PI;
	cout << "日の入方位角：" << A*180.0/M_PI << endl << endl;

	// 日の出から日の入りまでのループ
	int sr = (int)(floor(t1));
	int ss = (int)(floor(t2));

	for (int i=sr-1; i<=ss+1; i++)
	{
		cout << "時間: " << i << " [時], ";

// 時角
//	h=(JST-12)π/12+標準子午線からの経度差+均時差(Eq)
		double h = ((double)(i) - 12.0)*M_PI/12.0 + (longitude - 135.0) + eq;
		cout << "時角 = " << h << ", ";
//		h *= M_PI/180.0;           // [radian]

// 太陽高度 α
// α=arcsin{sin(φ)sin(δ)+cos(φ)cos(δ)cos(h)}
		double a = acos( sin(latitude*M_PI/180.0)*sin(delta) + cos(latitude*M_PI/180.0)*cos(delta)*cos(h) ); // [radian]
		cout << "高度 = " << a*180.0/M_PI << ", \n";

//		方位角：　A（北 = 0, 東 = 90, 南 = 180, 西 = 270°）
//		　　sinA = cos(δ)sin(t)/cos(h)
//		　　cosA = (sin(h)sin(φ) - sin(δ))/cos(h)/cos(φ)
//		　　A = atan2(sinA, cosA) + π
//		double sA = cos(delta)*sin(t) / cos(h);
//		double cA = (sin(h)*sin(latitude*M_PI/180.0) - sin(delta)) / cos(h) / cos(latitude*M_PI/180.0);
//		double A;
//		if (sA > 0)
//			A = M_PI/2.0 - atan2(cA,sA);  // [radian]
//		else
//			A = -M_PI/2.0 - atan2(cA,sA);
//		cout << "方位角 = " << A*180.0/M_PI << endl;
	}
}
