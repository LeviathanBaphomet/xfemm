/*
   This code is a modified version of an algorithm
   forming part of the software program Finite
   Element Method Magnetics (FEMM), authored by
   David Meeker. The original software code is
   subject to the Aladdin Free Public Licence
   version 8, November 18, 1999. For more information
   on FEMM see www.femm.info. This modified version
   is not endorsed in any way by the original
   authors of FEMM.

   This software has been modified to use the C++
   standard template libraries and remove all Microsoft (TM)
   MFC dependent code to allow easier reuse across
   multiple operating system platforms.

   Date Modified: 2011 - 11 - 10
   By: Richard Crozier
   Contact: richard.crozier@yahoo.co.uk
*/
#include <numeric>
#include <iostream>
#include <vector>
#include <algorithm>
#ifndef SPARS_H
#define SPARS_H

struct CEntry
{
public:
    int c;
    double x;
};

class CBigLinProb
{
public:

    std::vector<double> V;				// solution
    std::vector<double> P;				// search direction;
    std::vector<double> R;				// residual;
    std::vector<double> U;				// A * P;
    std::vector<double> Z;
    std::vector<double> b;				// RHS of linear equation
    std::vector<std::vector<CEntry>> M;				// pointer to list of matrix entries;
    int n = 0;					// dimensions of the matrix;
    int bdw;				// Optional matrix bandwidth parameter;
    double Precision;		// error tolerance for solution
    const double Lambda = 1.5;			// relaxation factor;

    int *Q; ///< Used by esolver and hsolver.

    int Create(int d, int bw);
    void Put(double v, int p, int q);
    double Get(int p, int q);
	void AddTo(double v, int p, int q);
	void SetValue(int i, double x);
    void Periodicity(int i, int j);
    void AntiPeriodicity(int i, int j);
    void Wipe();
	
    bool PCGSolve(int flag);
    void MultPC(const std::vector<double>& X, std::vector<double>& Y);
    void MultA(const std::vector<double>& X, std::vector<double>& Y);


private:

};

#endif
