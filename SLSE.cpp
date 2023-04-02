//
//  Created by Piet Groeneboom on May 1, 2022.
//  Copyright (c) 2022 Piet Groeneboom. All rights reserved.

// [[Rcpp::plugins("cpp11")]]

#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <chrono>
#include <random>
#include <Rcpp.h>

using namespace std;
using namespace Rcpp;

#define SQR(x) ((x)*(x))

typedef struct
{
    double x;
    double y;
}
data_object;

double  f0(double x);
int     compare(const void *a, const void *b);
void    order_data(int n, double **data);
void    convexminorant(int n, double cumw[], double cs[], double yy[]);
double  K(double x);
double  Kprime(double x);
double  KK(double x);
void    regression_estimate(int m, double tt[], int ngrid, double grid[], double h, double p[], double ff[], double fsmooth[]);
double  f0(double x);



// [[Rcpp::export]]

List SLSE(NumericMatrix X, int N)
{
    int     n,m,ngrid,i,j;
    double  **data;
    double  *cumw,*cs,*ff,*tt,*pp;
    double  h,*grid,*SLSE;
 
    n=(int)N;
    
    ngrid=100;
    
    h= 0.5*pow(n,-1.0/5);
    
    grid = new double[ngrid+1];
    
    for (i=0;i<=ngrid;i++)
        grid[i]=(double)i/ngrid;
    
    m=1;

    data = new double *[n+1];
    for (i=0;i<n+1;i++)
        data[i] = new double[2];
    
    for (i=0;i<n;i++)
    {
        for (j=0;j<2;j++)
              data[i+1][j]=(double)X(i,j);
    }
    
    order_data(n,data);
  
    cumw = new double[n+1];
    cs = new double[n+1];
    ff = new double[n+1];
    
    tt = new double[n+1];
    pp = new double[n+1];
    
    SLSE =new double[ngrid+1];
   
    cumw[0]=cs[0]=0;
    
    tt[0]=pp[0]=0;
    
    ff[0] = 0;
    
    for (i=1;i<=n;i++)
    {
        cumw[i]= (double)i;
        cs[i]=cs[i-1]+data[i][1];
    }
    
    convexminorant(n,cumw,cs,ff);

    j=0;
    ff[0]=ff[1];
            
    for (i=1;i<=n;i++)
    {
        if (ff[i] > ff[i-1])
        {
            j++;
            tt[j]=data[i][0];
            pp[j]=ff[i]-ff[i-1];
        }
    }
    
    m=j;
         
    for (i=1;i<=m;i++)
        ff[i]=ff[i-1]+pp[i];
    
    regression_estimate(m,tt,ngrid,grid,h,pp,ff,SLSE);
  
    NumericMatrix out1 = NumericMatrix(ngrid-1,2);
    
    for (i=0;i<ngrid-1;i++)
    {
        out1(i,0)=grid[i+1];
        out1(i,1)=SLSE[i+1];
    }
    
    NumericMatrix out2 = NumericMatrix(n,2);
    
    for (i=0;i<n;i++)
    {
        out2(i,0)=data[i+1][0];
        out2(i,1)=data[i+1][1];
    }
    
    // make the list for the output, containing the MLE, hazard, the bootstrap confidence intervals and -log likelihood
    
    List out = List::create(Rcpp::Named("SLSE")=out1,Rcpp::Named("data")=out2);
    // free memory
    
    for (i=0;i<n+1;i++)
        delete[] data[i];
    delete[] data;

    delete[] cumw; delete[] cs; delete[] ff;
    
    delete[] tt; delete[] pp;
    delete[] SLSE;;
    
    return out;
}

double f0(double x)
{
  return SQR(x)+x/5;
}

void order_data(int n, double **data)
{
  int    i;
  data_object *obs;
  
  obs = new data_object[n];
  
  for (i=0;i<n;i++)
  {
    obs[i].x = data[i+1][0];
    obs[i].y = data[i+1][1];
  }
  
  qsort(obs,n,sizeof(data_object),compare);
  
  data[0][0]=data[0][1]=0;
  for (i=1;i<=n;i++)
  {
    data[i][0] = obs[i-1].x;
    data[i][1] = obs[i-1].y;
  }
  
  delete[] obs;
}

int compare(const void *a, const void *b)
{
    if ((*(data_object *) a).x < (*(data_object *) b).x)
        return -1;
    if ((*(data_object *) a).x > (*(data_object *) b).x)
        return 1;
    return 0;
}



void convexminorant(int n, double cumw[], double cs[], double yy[])
{
    int    i,j,m;
    
    yy[1] = cs[1]/cumw[1];
    for (i=2;i<=n;i++)
    {
        yy[i] = (cs[i]-cs[i-1])/(cumw[i]-cumw[i-1]);
        if (yy[i-1]>yy[i])
        {
            j = i;
            while ((yy[j-1] > yy[i]) && (j>1))
            {
                j--;
                yy[i] = (cs[i]-cs[j-1])/(cumw[i]-cumw[j-1]);
                for (m=j;m<i;m++)
                    yy[m] = yy[i];
            }
        }
    }
}


double K(double x)
{
    double u,y;
    
    u=x*x;
    
    if (u<=1)
        y=(35.0/32)*pow(1-u,3);
    else
        y=0.0;
    
    return y;
}

double Kprime(double x)
{
    double u,y;
    
    u=x*x;
    
    if (u<=1) y = -(105.0/16)*x*SQR(1-u);
    else    y=0.0;
    
    return y;
}

void regression_estimate(int m, double tt[], int ngrid, double grid[], double h, double p[], double ff[], double fsmooth[])
{
    int i,k;
    double a,u,t;
    
    for (i=0;i<=ngrid;i++)
    {
        a=0;
        
        u=grid[i];
        
        if (u>=h && u <=1-h)
         {
             for (k=1;k<=m;k++)
             {
                 t=(u-tt[k])/h;
                 a += KK(t)*p[k];
             }
         }
        
        if (u<h)
        {
            for (k=1;k<=m;k++)
            {
                t=(h-tt[k])/h;
                a += (KK(t)+(u-h)*K(t)/h+0.5*SQR(u-h)*Kprime(t)/SQR(h))*p[k];
                //a += (KK(t)+(u-h)*K(t)/h)*p[k];
            }
        }
        
        if (u>1-h)
        {
            for (k=1;k<=m;k++)
            {
                t=(1-h-tt[k])/h;
                a += (KK(t)+(u-(1-h))*K(t)/h+0.5*SQR(u-(1-h))*Kprime(t)/SQR(h))*p[k];
                //a += (KK(t)+(u-(1-h))*K(t)/h)*p[k];
            }
        }
        fsmooth[i]=ff[0]+a;
    }
}

double KK(double x)
{
    double u,y;
    
    u=x*x;
    
    if (u<=1)
        y = (16.0 + 35*x - 35*pow(x,3) + 21*pow(x,5) - 5*pow(x,7))/32.0;
    else
    {
        if (x>1)
            y=1;
        else
            y=0;
        
    }
    
    return y;
}








