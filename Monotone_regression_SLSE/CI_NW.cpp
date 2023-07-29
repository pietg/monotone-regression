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
void    data_bootstrap(int n, double **data, double SLSE_data[], double residu[],
                    double **bootstrap_data, int seed);
double  f0(double x);
double  regression_NW(int n, double **data, double h, double u);




// [[Rcpp::export]]

List CI_NW(NumericMatrix X, int N, int seed1)
{
    int     n,ngrid,i,j,seed;
    int     iter,NumIt;
    int     percentile1,percentile2,*percentage;
    double  *lowbound,*upbound,**f3,*f4;
    double  **data,**bootstrap_data,*residu;
    double  h,h2,*grid,*SLSE_data,*SLSE;
    double  *SLSE_bootstrap,*SLSE2;
    double  mean_residu;
    
    seed=(int)seed1;
    n=(int)N;
    
    ngrid=100;
    
    h= 0.5*pow(n,-1.0/5);
    h2= 0.7*pow(n,-1.0/9);
    
    NumIt=1000;
    
    grid = new double[ngrid+1];
    
    for (i=0;i<=ngrid;i++)
        grid[i]=(double)i/ngrid;
    
    percentile1=round(0.025*NumIt);
    percentile2=round(0.975*NumIt);
    
    data = new double *[n+1];
    for (i=0;i<n+1;i++)
        data[i] = new double[2];
    
    for (i=0;i<n;i++)
    {
        for (j=0;j<2;j++)
              data[i+1][j]=(double)X(i,j);
    }
    
    order_data(n,data);
    
    bootstrap_data = new double *[n+1];
    for (i=0;i<n+1;i++)
      bootstrap_data[i] = new double[2];
  
    SLSE =new double[ngrid+1];
    SLSE2 =new double[ngrid+1];
    SLSE_bootstrap = new double[ngrid+1];
    SLSE_data = new double[n+1];
    
    residu = new double[n+1];
    
    f3 = new double*[NumIt+1];
    
    for (iter=0;iter<NumIt+1;iter++)
        f3[iter] = new double[ngrid+1];
        
    f4 = new double[NumIt+1];
    
    percentage = new int[ngrid+1];
    
    for (i=0;i<=ngrid;i++)
        percentage[i]=0;
        
    lowbound = new double[ngrid+1];
    upbound  = new double[ngrid+1];
        
 
    data[0][0]=bootstrap_data[0][0]=0;
    data[0][1]=bootstrap_data[0][1]=0;
  
    for (i=0;i<=ngrid;i++)
        SLSE[i] = regression_NW(n,data,h,grid[i]);
    
    for (i=1;i<=n;i++)
        SLSE_data[i] = regression_NW(n,data,h2,data[i][0]);
    
    for (i=0;i<=ngrid;i++)
        SLSE2[i] = regression_NW(n,data,h2,grid[i]);

    for (i=1;i<=n;i++)
        residu[i]=data[i][1]-SLSE_data[i];
    
    mean_residu =0;
    for (i=1;i<=n;i++)
        mean_residu += residu[i];
    mean_residu /=n;
    
    for (i=1;i<=n;i++)
        residu[i] -= mean_residu;
    
    for (iter=0;iter<NumIt;iter++)
    {
        seed++;
        
        data_bootstrap(n,data,SLSE_data,residu,bootstrap_data,seed);
        
        // bootstrap SLSE at points of the grid
        
        for (i=0;i<=ngrid;i++)
            SLSE_bootstrap[i] = regression_NW(n,bootstrap_data,h,grid[i]);
        
        for (i=1;i<ngrid;i++)
            f3[iter][i] = SLSE_bootstrap[i]-SLSE2[i];
    }
    
    for (i=1;i<ngrid;i++)
    {
        for (iter=0;iter<NumIt;iter++)
            f4[iter]=f3[iter][i];
        
        qsort(f4,NumIt,sizeof(double),compare);
        
        lowbound[i] = SLSE[i]-f4[percentile2-1];
        upbound[i]  = SLSE[i]-f4[percentile1-1];
    }

    NumericMatrix out1 = NumericMatrix(ngrid-1,4);
    
    for (i=0;i<ngrid-1;i++)
    {
        out1(i,0)=grid[i+1];
        out1(i,1)=SLSE[i+1];
        out1(i,2)=lowbound[i+1];
        out1(i,3)=upbound[i+1];
    }
    
    NumericMatrix out2 = NumericMatrix(n,2);
    
    for (i=0;i<n;i++)
    {
        out2(i,0)=data[i+1][0];
        out2(i,1)=data[i+1][1];
    }
    
    // make the list for the output, containing the MLE, hazard, the bootstrap confidence intervals and -log likelihood
    
    List out = List::create(Rcpp::Named("CI_NW")=out1,Rcpp::Named("data")=out2);
    // free memory
    
    for (i=0;i<n+1;i++)
        delete[] data[i];
    delete[] data;
    
    for (i=0;i<n+1;i++)
        delete[] bootstrap_data[i];
    delete[] bootstrap_data;
    
    for (i=0;i<NumIt+1;i++)
        delete[] f3[i];
    delete[] f3;
    
    delete[] residu;
    
    delete[] SLSE; delete[] SLSE_data; delete[] SLSE2;
    
    delete[] lowbound; delete[] upbound; delete[] f4; delete[] percentage;
    
    return out;
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

void data_bootstrap(int n, double **data, double SLSE_data[], double residu[], double **bootstrap_data, int seed)
{
    int    i,j;
    
    std::mt19937_64 gen(seed);
    std::uniform_int_distribution<int> dis(1,n);
    
    for (i=1;i<=n;i++)
    {
        bootstrap_data[i][0] = data[i][0];
        j=dis(gen);
        bootstrap_data[i][1] = SLSE_data[i]+residu[j];
        //bootstrap_data[i][1] = f0(data[i][0])+residu[j];
    }
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

double smooth_LSE(double A, double B, int m, double t[], double p[], double u, double h)
{
    int       k;
    double    t1,t2,t3,sum;
    
    sum=0;
    
    for (k=1;k<=m;k++)
    {
        t1=(u-t[k])/h;
        t2=(u+t[k]-2*A)/h;
        t3=(2*B-u-t[k])/h;
        sum+= (KK(t1)+KK(t2)-KK(t3))*p[k];
        //sum+= KK(t1)*p[k];
    }

    return sum;
}










