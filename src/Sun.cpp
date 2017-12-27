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

//   太陽赤緯：　δ（太陽光線と地球の赤道面との角度、±23°27'の範囲で変化）　[単位： 度]
//　　δ = 0.33281 - 22.984 cos(ωJ ) - 0.34990 cos(2ωJ ) - 0.13980cos(3ωJ )
//　　　　　　　　　+ 3.7872 sin(ωJ ) + 0.03250 sin(2ωJ ) + 0.07187 sin(3ωJ )
//　　ここで、
//　　ω = 2π/365、閏年は ω = 2π/366、J： 元日からの通算日数 + 0.5
	double omega = 2.0*M_PI / (365.0 + leap(year));
	double J = GetDays( year, month, day ) + 0.5;
	double delta = 0.33281 - 22.984*cos(omega*J) - 0.34990*cos(2.0*omega*J) - 0.13980*cos(3.0*omega*J);
	              delta += + 3.7872*sin(omega*J) + 0.03250*sin(2.0*omega*J) + 0.07187*sin(3.0*omega*J);
	cout << "太陽赤緯 = " << delta << " [度]" << endl;
	delta *= M_PI/180.0; // radianに変換

//	均時差：　e（天球上を一定な速さで動くと考えた平均太陽と、実際の太陽との移動の差、17分未満）　[単位： 時間]
//	　　e = 0.0072 cos(ωJ ) - 0.0528 cos(2ωJ ) - 0.0012 cos(3ωJ )
//	　　　　- 0.1229 sin(ωJ ) - 0.1565 sin(2ωJ ) - 0.0041 sin(3ωJ )

// New
// Eq=0.000075+0.001868cos(θo)-0.032077sin(θo)-0.014615cos(2θo)-0.040849sin(2θo)

	double e = 0.0072*cos(omega*J) - 0.0528*cos(2.0*omega*J) - 0.0012*cos(3.0*omega*J);
	    e += - 0.1229*sin(omega*J) - 0.1565*sin(2.0*omega*J) - 0.0041*sin(3.0*omega*J);
	cout << "均時差 = " << e << " [時間] （" << e*60.0 << " [分]）" << endl;

//	日の出時刻：　t1　[単位： 時]
//	　　t = acos(-tan(δ)tan(φ))
//	　　T1 = (-t + 180)/15
//	　　t1 = T1 - (θ - 135)/15 - e
	double t = acos( -tan(delta)*tan(latitude*M_PI/180.0) ); // [radian]
	t *= 180.0/M_PI;  // [度]
	double T1 = (-t + 180.0)/15.0;
	double t1 = T1 - (longitude - 135.0)/15.0 - e;

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
	double t2 = T2 - (longitude - 135.0)/15.0 - e;
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

//		時角：　t　[単位： 度]
//		　　T = Ts + (θ - 135)/15 + e
//		　　t = 15T - 180
//		　　ここで、
//		　　　Ts：　時刻（中央標準時）
//		　　　θ：　東経
//		　　　φ：　北緯
		double ts = (double)(i) - 9.0;  // 中央標準時に変換
		double T = ts + (longitude - 135.0)/15.0 + e;
		double t = 15.0*T - 180.0; // [度]
		cout << "時角 = " << t << ", ";
		t *= M_PI/180.0;           // [radian]

//		高度（仰角）：　h
//		　　h = asin( sin(φ)sin(δ) + cos(φ)cos(δ)cos(t) )
		double h = asin( sin(latitude*M_PI/180.0)*sin(delta) + cos(latitude*M_PI/180.0)*cos(delta)*cos(t) ); // [radian]
		cout << "高度 = " << h*180.0/M_PI << ", ";

//		方位角：　A（北 = 0, 東 = 90, 南 = 180, 西 = 270°）
//		　　sinA = cos(δ)sin(t)/cos(h)
//		　　cosA = (sin(h)sin(φ) - sin(δ))/cos(h)/cos(φ)
//		　　A = atan2(sinA, cosA) + π
		double sA = cos(delta)*sin(t) / cos(h);
		double cA = (sin(h)*sin(latitude*M_PI/180.0) - sin(delta)) / cos(h) / cos(latitude*M_PI/180.0);
		double A;
		if (sA > 0)
			A = M_PI/2.0 - atan2(cA,sA);  // [radian]
		else
			A = -M_PI/2.0 - atan2(cA,sA);
		cout << "方位角 = " << A*180.0/M_PI << endl;
	}
}
