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

#include "femmcomplex.h"
#include "spars.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <utility>

#define KLUDGE

int CBigLinProb::Create(int d, int bw) {

    bdw = bw;
    b.resize(d);
    V.resize(d);
    P.resize(d);
    R.resize(d);
    U.resize(d);
    Z.resize(d);
    M.resize(d);

    n = d;

    return 1;
}

void CBigLinProb::Put(double v, int p, int q) {
    if (q < p)
        std::swap(p, q);

    std::vector<CEntry>& entryList = M[p];

    auto it = std::lower_bound(entryList.begin(), entryList.end(), q,
        [](const CEntry& entry, int qVal) { return entry.c < qVal; });

    if (it != entryList.end() && it->c == q) {
        it->x = v;
    }
    else {
        CEntry newEntry{ q, v };
        entryList.insert(it, newEntry);
    }
}

double CBigLinProb::Get(int p, int q) {
    if (q < p)
        std::swap(p, q);

    const std::vector<CEntry>& entryList = M[p];

    auto it = std::lower_bound(entryList.begin(), entryList.end(), q,
        [](const CEntry& entry, int qVal) { return entry.c < qVal; });

    if (it != entryList.end() && it->c == q)
        return it->x;

    return 0;
}

void CBigLinProb::AddTo(double v, int p, int q)
{
	Put(Get(p,q)+v,p,q);
}

void CBigLinProb::MultA(const std::vector<double>& X, std::vector<double>& Y)
{
    for (int i = 0; i < n; i++) {
        Y[i] = 0.0;
    }

    std::fill(Y.begin(), Y.end(), 0.0);

    for (int i = 0; i < n; i++) {
        Y[i] += M[i][0].x * X[i];

        for (size_t j = 1; j < M[i].size(); j++) {
            const CEntry& entry = M[i][j];
            Y[i] += entry.x * X[entry.c];
            Y[entry.c] += entry.x * X[i];
        }
    }
}

void CBigLinProb::MultPC(const std::vector<double>& X, std::vector<double>& Y)
{
    double c = Lambda * (2.0 - Lambda);
    for (int i = 0; i < n; i++) {
        Y[i] = X[i] * c;
    }

    for (int i = 0; i < n; i++) {
        Y[i] /= M[i][0].x;

        for (size_t j = 1; j < M[i].size(); j++) {
            const CEntry& entry = M[i][j];
            Y[entry.c] -= entry.x * Y[i] * Lambda;
        }
    }

    for (int i = 0; i < n; i++) {
        Y[i] *= M[i][0].x;
    }

    for (int i = n - 1; i >= 0; i--) {
        for (size_t j = 1; j < M[i].size(); j++) {
            const CEntry& entry = M[i][j];
            Y[i] -= entry.x * Y[entry.c] * Lambda;
        }
        Y[i] /= M[i][0].x;
    }
}

bool CBigLinProb::PCGSolve(int flag)
{
    double res, res_o, res_new;
    double del, rho, pAp;
    double PrecisionSq = Precision * Precision;

    for (int i = 0; i < n; i++)
    {
        if (M[i].front().x == 0)
        {
            fprintf(stderr, "singular flag tripped at %i of %i\n", i, n);
            return false;
        }
    }

    MultPC(b, Z);

    res_o = std::transform_reduce(Z.begin(), Z.end(), b.begin(), 0.0);

    if (res_o == 0) return true;

    MultA(V, R);

    std::transform(b.begin(), b.end(), R.begin(), R.begin(), std::minus<double>());

    MultPC(R, Z);

    std::copy(Z.begin(), Z.end(), P.begin());

    res = std::transform_reduce(Z.begin(), Z.end(), R.begin(), 0.0);

    do
    {
        MultA(P, U);
        pAp = std::transform_reduce(P.begin(), P.end(), U.begin(), 0.0);
        del = res / pAp;

        std::transform(V.begin(), V.end(), P.begin(), V.begin(), [del](double v, double p) { return v + del * p; });
        std::transform(R.begin(), R.end(), U.begin(), R.begin(), [del](double r, double u) { return r - del * u; });

        MultPC(R, Z);

        res_new = std::transform_reduce(Z.begin(), Z.end(), R.begin(), 0.0);

        rho = res_new / res;
        res = res_new;

        std::transform(Z.begin(), Z.end(), P.begin(), P.begin(), [rho](double z, double p) { return z + rho * p; });

    } while (res / res_o > PrecisionSq);

    return true;
}

void CBigLinProb::SetValue(int i, double x)
{
    int k,fst,lst;
    double z;

    if(bdw==0)
    {
        fst=0;
        lst=n;
    }
    else
    {
        fst=i-bdw;
        if (fst<0) fst=0;
        lst=i+bdw;
        if (lst>n) lst=n;
    }

    for(k=fst; k<lst; k++)
    {
        z=Get(k,i);
        if(z!=0)
        {
            b[k]=b[k]-(z*x);
            if(i!=k) Put(0.,k,i);
        }
    }
    b[i]=Get(i,i)*x;
}

void CBigLinProb::Wipe() {
    for (int i = 0; i < n; i++) {
        b[i] = 0.0;
        M[i].clear();
    }
}

void CBigLinProb::AntiPeriodicity(int i, int j)
{
    int k,fst,lst;
    double v1,v2,c;

#ifdef KLUDGE
    int tmpbdw=bdw;
    bdw=0;
#endif

    if (j<i)
        std::swap(j,i);

    if(bdw==0)
    {
        fst=0;
        lst=n;
    }
    else
    {
        fst=i-bdw;
        if (fst<0) fst=0;
        lst=j+bdw;
        if (lst>n) lst=n;
    }

    for(k=fst; k<lst; k++)
    {
        if((k!=i) && (k!=j))
        {
            v1=Get(k,i);
            v2=Get(k,j);
            if ((v1!=0) || (v2!=0))
            {
                c=(v1-v2)/2.;
                Put(c,k,i);
                Put(-c,k,j);
            }
        }
        if((k==i+bdw) && (k<j-bdw) && (bdw!=0)) k=j-bdw;
    }

    c=0.5*(Get(i,i)+Get(j,j));
    Put(c,i,i);
    Put(c,j,j);

    c=0.5*(b[i]-b[j]);
    b[i]=c;
    b[j]=-c;

#ifdef KLUDGE
    bdw=tmpbdw;
#endif
}

void CBigLinProb::Periodicity(int i, int j)
{
    int k,fst,lst;
    double v1,v2,c;

#ifdef KLUDGE
    int tmpbdw=bdw;
    bdw=0;
#endif

    if (j<i)
        std::swap(j,i);

    if(bdw==0)
    {
        fst=0;
        lst=n;
    }
    else
    {
        fst=i-bdw;
        if (fst<0) fst=0;
        lst=j+bdw;
        if (lst>n) lst=n;
    }

    for(k=fst; k<lst; k++)
    {
        if((k!=i) && (k!=j))
        {
            v1=Get(k,i);
            v2=Get(k,j);
            if ((v1!=0) || (v2!=0))
            {
                c=(v1+v2)/2.;
                Put(c,k,i);
                Put(c,k,j);
            }
        }
        if((k==i+bdw) && (k<j-bdw) && (bdw!=0)) k=j-bdw;
    }

    c=(Get(i,i)+Get(j,j))/2.;
    Put(c,i,i);
    Put(c,j,j);

    c=0.5*(b[i]+b[j]);
    b[i]=c;
    b[j]=c;

#ifdef KLUDGE
    bdw=tmpbdw;
#endif
}