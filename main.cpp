
#include "./matrices.h"
#include <algorithm>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <stdio.h>
#include <sstream>
#include <ctime>
#include <boost/timer.hpp> 
#include <boost/progress.hpp>
#include <string>
#include <cstring>


using namespace std;
using namespace matrices;

void printFile(const array &results, int gS, double SMAX, char s[])
{

    double ds = SMAX/(double)gS;
    char rname [300];
    sprintf(rname, "%s",s);
    ofstream dataFile(rname, ios::out);
    dataFile.precision(10);
    cout.precision(10);

    if (!dataFile) {
        cout << " Error " << endl;
    }

        for (int i = 0; i <= results.size() - 1; i++) {

            dataFile << showpoint << setw(6) << ds*(i) << setw(25) << results[i] << endl;

        }

        dataFile << endl;

}

array hjbDiscretisation(int i, double dP, double V ,double dtow, double vol, double drift)
{
    array coeff(3);
    double alpha;
    double beta;
    double gamma;
    double alphacentral;
    double betacentral;
    double gammacentral;
    double alphaforward;
    double betaforward;
    double gammaforward;
    double gammaback;
    double alphaback;





    alphaforward = (vol*vol)*((i)*(i))/2.;
    gammaforward = (vol*vol) * i * i / (2.) + drift*i;




        gamma = gammaforward;
        alpha = alphaforward;
        beta = (alpha + gamma + V);

        coeff[0] = -alpha*dtow;

        coeff[1] = 1 + dtow*beta;
        coeff[2] = -gamma*dtow;

        return coeff;
    }

array triSolver(matrix &mat)
{

    double a1 = 0;
    double b1 = 0;
    double c1 = 0;
    double d1 = 0;

    double a2 = 0;
    double b2 = 0;
    double c2 = 0;
    double d2 = 0;

    array solution(mat.size());




    for (int i = 0; i != mat.size() - 1; i++) {

        a1 = mat[i][0];
        b1 = mat[i][1];
        c1 = mat[i][2];
        d1 = mat[i][3];
        a2 = mat[i+1][0];
        b2 = mat[i+1][1];
        c2 = mat[i+1][2];
        d2 = mat[i+1][3];

        if (i == 0)
        {
            mat[i+1][0] = 0;
            mat[i+1][1] = b2;
            mat[i+1][2] = c2;
            mat[i+1][3] = d2 - (a2 * d1 / a1);
        }
        else if (i > 0 && i < (mat.size() - 2))
        {
            mat[i+1][0] = 0;
            mat[i+1][1] = (b2) - (a2 * c1) / b1;
            mat[i+1][2] = (c2);
            mat[i+1][3] = d2 - (a2 * d1 / b1);
        }


    }
    solution[mat.size() - 1] = mat[mat.size() - 1][3] / mat[mat.size() - 1][2];

    for (int m = mat.size() - 2; m >= 1; m--)
    {
        solution[m] = (1. / mat[m][1])*(mat[m][3] - mat[m][2] * solution[m + 1]);

    }

    solution[0] = (mat[0][3]) / (mat[0][0]);

    return (solution);
}

double payoff(double strike, double S)
{
    double call;
    call = S - strike > 0 ? S - strike : 0 ;
    return call;
}

double barrier(double strike,double b,double s)
{

    if (s>=b)
        return 0;
    else if (s - strike > 0)
        return (s-strike);
    else
        return 0;

}
double abserror(double a, double b)
{
    return abs(a-b);
}


double rowValue(const array &row,const array &preSol, int i)
{
    double x = (row[0]*preSol[i-1] + row[1]*preSol[i] + row[2]*preSol[i+1]);
    return(x);
}

void policyIteration(double volmin,double volmax, double r, double K1, double K2, double K3, double SMAX, double T, int gT, int gS, int gC ,double tol, array &results,array &init, array & contret)
{

    double dc = (volmax-volmin)/double(gC);
    double dS = SMAX/double(gS);
    double dt = T/double(gT);
    array initsolution(gS+1);
    array solution(gS+1);
    array presolution(gS+1);
    array preapprox(gS+1);
    array approx(gS+1);
    array row(gC+1);
    matrix tri(gS+1, array(4));
    matrix preTri(gS+1, array(4));
    array control(gS+1);

    for(int s = 0; s !=gS+1; s++)
    {

        initsolution[s] = payoff(K1,s*dS) - 2.*payoff(K2,s*dS) + payoff(K3, s*dS);

    }

    presolution = initsolution;

    preapprox = initsolution;

    for(int t=0;t!=gT; t++)
    {
        cout << (t/double(gT))*100 << " %" <<endl;
        preapprox = presolution;

        for( int iter=0; iter!= 1000; iter++)
        {

                 array error(gS+1);

                 double controlvalue;


                 array controlloop(gC+1);
                 array hjb(3);

                for(int qs = 0 ; qs !=gS+1; qs++ )
                {
                    if(qs == 0)
                   {
                        control[qs] = 0;
                        continue;
                   }

                    if(qs!=0)
                    {
                    for (int controlMin = 0; controlMin != gC+1; controlMin++)
                    {
                        controlvalue = volmin + controlMin*dc;
                        controlloop[controlMin] = controlvalue;
                        hjb = hjbDiscretisation(qs,dS, r ,dt, controlvalue, r);
                        row[controlMin] = rowValue(hjb,preapprox,qs);
                    }

                    control[qs] = controlloop[(max_element(row.begin(),row.end())) - row.begin()];

                    continue;
                    }
                  if(qs == gS)
                    {
                        control[qs] = control[qs-1];
                        continue;
                    }

                }





                    for(int siter = 0; siter!=gS+1; siter++)
                    {
                         hjb = hjbDiscretisation(siter, dS , r ,dt, control[siter] , r) ;

                        if(siter == 0)
                        {

                            tri[0][0] = hjb[1];

                            tri[0][1]= 0;

                            tri[0][2]= 0;

                            tri[0][3]= presolution[0];

                            continue;

                        }

                        if(siter != 0 && siter != gS)
                        {
                            tri[siter][0] = hjb[0];
                            tri[siter][1] = hjb[1];
                            tri[siter][2] = hjb[2];
                            tri[siter][3] = presolution[siter];

                            continue;
                        }

                        if(siter == gS)
                        {
                           tri[siter][0] = 0;
                           tri[siter][1] = 0;
                           tri[siter][2] = 1.;
                           tri[siter][3] = 0;
                           continue;
                        }
                    }

                 solution = triSolver(tri);

                 transform(solution.begin(),solution.end(),preapprox.begin(),error.begin(),abserror);

               double maxerror = *(max_element(error.begin(), error.end()));


                 if(maxerror < tol)
                 {

                    break;

                 }else
                 {
                   preapprox = solution;
                 }


        }

        presolution = solution;
    }
    contret = control;
    results = presolution;
    init = initsolution;
    cout << endl;

}

double point(double istar, int maxsize, int minsize)
{

    int intistar = istar;
    double value = max(minsize,min(intistar,maxsize));
    return value;
}

double interpolate1d(double hS, double ds , int smax, array &vec)
{

    int s1;
    int s2;


    double ds1;


    ds1 = hS / ds;


    s1 = point(ds1,smax-1,0);
    s2 = s1+1;


    double step1 = ((s2 * ds - hS) / double(s2 * ds - s1 * ds)) * vec[s1] + ((hS - s1 * ds) / double(s2 * ds - s1 * ds)) * vec[s2];

    return step1;
}


int main()
{
double volmin = 0.3;
double volmax = 0.45;
double r = 0.04;
double K1 = 95;
double K2 =100;
double K3 = 105;
double SMAX = 500;
double T = 0.5;
int gT = 10000;
int gS = 5000;
int gC = 1;
double tol= 0.0000000001;
array results;
array init;
array control;
policyIteration(volmin, volmax, r, K1, K2, K3, SMAX, T, gT, gS, gC, tol,results,init,control);
cout << " value at 100 = [ " << interpolate1d(100.,SMAX/double(gS),gS,results) << " ]" << endl;

char name[] = "results.dat";
printFile(results, gS, SMAX,name);
char name1[] = "initresu.dat";
printFile(init,gS,SMAX,name1);
char name2[] = "control.dat";
printFile(control,gS,SMAX,name2);

return 1;

}
