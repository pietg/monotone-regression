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

typedef struct
{
    double alpha, beta;
}
weight_t;

weight_t weight(double x);

double  f0(double x);
int     compare(const void *a, const void *b);
void    order_data(int n, double **data);
void    convexminorant(int n, double cumw[], double cs[], double yy[]);
double  K(double x);
double  Kprime(double x);
double  KK(double x);
double  regression_NW(int n, double **data, double h, double u);
double  f0(double x);

// [[Rcpp::export]]

List NW(NumericMatrix X, int N)
{
    int     n,ngrid,i,j;
    double  **data;
    double  h,*grid,*SLSE;
 
    n=(int)N;
    
    ngrid=100;
    
    h= 0.5*pow(n,-1.0/5);
    
    grid = new double[ngrid+1];
    
    for (i=0;i<=ngrid;i++)
        grid[i]=(double)i/ngrid;

    data = new double *[n+1];
    for (i=0;i<n+1;i++)
        data[i] = new double[2];
    
    for (i=0;i<n;i++)
    {
        for (j=0;j<2;j++)
              data[i+1][j]=(double)X(i,j);
    }
    
    order_data(n,data);
    
    SLSE =new double[ngrid+1];
    
    for (i=0;i<=ngrid;i++)
        SLSE[i] = regression_NW(n,data,h,grid[i]);
  
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
    
    delete[] SLSE;
    
    return out;
}

double f0(double x)
{
  return SQR(x)+x/5;
}

weight_t weight(double x)
{
    short i;
    double y1, y2, y3, y4;
    double *xx;
    weight_t temp;
    
    xx= new double[10];
    xx[1] = x;
    for (i=2;i<=9;i++)
        xx[i] = x * xx[i - 1];
    
    y1 = 0.5 + (35.0*xx[1])/32.0-(35.0*xx[3])/32.0+(21.0*xx[5])/32.0-(5.0*xx[7])/32.0;
    y2 = -35.0/256.0+(35.0 * xx[2])/64.0-(105.0*xx[4])/128.0+(35.0*xx[6])/64.0-(35.0*xx[8])/256.0;
    y3 = 1.0/18 + 25.0*xx[3]/96 - 21.0*xx[5]/32 + 15.0*xx[7]/32 - 35.0*xx[9]/288;
    
    
    y4 = y1 * y3 - (y2*y2);
    temp.alpha = y3 / y4;
    temp.beta = -y2 / y4;
    
    delete[] xx;
    
    return temp;
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

double regression_NW(int n, double **data, double h, double u)
{
    int             i;
    double          rho,sum1,sum2,x,value;
    weight_t        weights;
    
    sum1=sum2=value=0;
    if (u-h>=0 && u+h<=1)
    {
        for (i=1;i<=n;i++)
        {
            x=(u-data[i][0])/h;
            sum1+= data[i][1]*K(x)/h;
            sum2+= K(x)/h;
        }
        value = sum1/sum2;
    }
    else
    {
        if (u<h)
        {
            rho = u/h;
            weights=weight(rho);
            for (i=1;i<=n;i++)
            {
                x=(u-data[i][0])/h;
                sum1 += data[i][1]*(K(x)*weights.alpha + x*K(x)*weights.beta)/h;
                sum2 += (K(x)*weights.alpha + x*K(x)*weights.beta)/h;
            }
            value = sum1/sum2;
            
        }
        else
        {
            if (u>1-h)
            {
                rho = (1-u)/h;
                weights=weight(rho);
                for (i=1;i<=n;i++)
                {
                    x=(u-data[i][0])/h;
                    sum1 += data[i][1]*(K(x)*weights.alpha - x*K(x)*weights.beta)/h;
                    sum2 += (K(x)*weights.alpha - x*K(x)*weights.beta)/h;
                }
                value = sum1/sum2;
            }
        }
    }
     
    return value;
}









