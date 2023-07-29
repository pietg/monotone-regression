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

double  f(double x);
double  fprime(int m, double tt[], double pp[], double h, double u);
int     compare(const void *a, const void *b);
void    order_data(int n, double **data);
void    convexminorant(int n, double cumw[], double cs[], double yy[]);
double  K(double x);
double  Kprime(double x);
double  KK(double x);
void    data_bootstrap(int n, double **data, double SLSE_data[], double residual[],
                    double **bootstrap_data, int seed);
void    regression_estimate(int m, double tt[], int ngrid, double grid[], double h, double h0, double p[], double ff[], double fsmooth[]);
//double MSE_regression(int n,  double data0[], int B, double **data, double **bootstrap_data,  double SLSE_data[], double SLSE_bootstrap[], double tt[], double pp[], double tt_bootstrap[], double pp_bootstrap[], double ff_bootstrap[], double residual[], double h, double h0, int seed);

double MSE_regression(int n, int ngrid, double grid[], int B, double **data, double **bootstrap_data,  double SLSE_data[],
                      double SLSE2[], double SLSE_bootstrap[], double tt[], double pp[], double tt_bootstrap[], double pp_bootstrap[],
                       double ff_bootstrap[], double residual[], double h, double h0, int seed);



// [[Rcpp::export]]

List bandwidth(NumericMatrix X, int N, int seed1, int BB)
{
    int     n,m,B,ngrid,i,j,seed;
    int     iter,NumIt;
    double  **data,**bootstrap_data,*residual;
    double  *cumw,*cs,*ff,*tt,*pp,*tt_bootstrap,*pp_bootstrap,*ff_bootstrap,*data0;
    double  h,h0,hmin,*grid,*SLSE_data,*SLSE2,*MSE;
    double  *SLSE_bootstrap;
    double  mean_residual,min;
    
    seed=(int)seed1;
    n=(int)N;
    
    ngrid=100;
    
    h= 0.5*pow(n,-1.0/5);
    h0= 0.7*pow(n,-1.0/9);
    
    hmin=100;
    
    NumIt=100;
    B=(int)BB;
    
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
    
    data0 = new double[n+1];
    
    bootstrap_data = new double *[n+1];
    for (i=0;i<n+1;i++)
        bootstrap_data[i] = new double[2];
    
    
    cumw = new double[n+1];
    cs = new double[n+1];
    ff = new double[n+1];
    
    tt = new double[n+1];
    pp = new double[n+1];
    
    tt_bootstrap = new double[n+1];
    pp_bootstrap = new double[n+1];
    ff_bootstrap = new double[n+1];
    
    SLSE_bootstrap = new double[ngrid+1];
    SLSE_data = new double[n+1];
    SLSE2 = new double[ngrid+1];
    
    MSE = new double[NumIt+1];
    
    residual = new double[n+1];
        
    cumw[0]=cs[0]=0;
    
    tt[0]=pp[0]=0;
    tt_bootstrap[0]=pp_bootstrap[0]=0;
    
    data[0][0]=bootstrap_data[0][0]=0;
    data[0][1]=bootstrap_data[0][1]=0;
    
    ff[0] = ff_bootstrap[0]=0;
        
    for (i=1;i<=n;i++)
        data0[i]=data[i][0];
    
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
    
    regression_estimate(m,tt,n,data0,h0,h0,pp,ff,SLSE_data);
    regression_estimate(m,tt,ngrid,grid,h0,h0,pp,ff,SLSE2);
    
    for (i=1;i<=n;i++)
        residual[i]=data[i][1]-SLSE_data[i];
    
    mean_residual =0;
    for (i=1;i<=n;i++)
        mean_residual += residual[i];
    mean_residual /=n;
    
    for (i=1;i<=n;i++)
        residual[i] -= mean_residual;
    
    min=1.0e10;
    
    for (iter=1;iter<=NumIt;iter++)
    {
        seed +=B;
        h = iter*0.01*pow(n,-1.0/5);
                
        MSE[iter] =
        MSE_regression(n,ngrid,grid,B,data,bootstrap_data,SLSE_data,SLSE2, SLSE_bootstrap,tt,pp,tt_bootstrap,pp_bootstrap,ff_bootstrap,residual,h,h0,seed);
        //MSE_regression(n,data0,B,data,bootstrap_data,SLSE_data,SLSE_bootstrap,tt,pp, tt_bootstrap,pp_bootstrap,ff_bootstrap,residual,h,h0,seed);

        
        if (MSE[iter]<min)
        {
            min = MSE[iter];
            hmin = iter*0.01;
        }
        
        Rcout  << setw(10) << iter << setprecision(8) <<  setw(15) << h << setprecision(8) <<  setw(15) << MSE[iter] << std::endl;
    }
    

    NumericMatrix out1 = NumericMatrix(NumIt,2);
    for (i=0;i<NumIt;i++)
    {
        out1(i,0)=(i+1)*0.01;
        out1(i,1)=MSE[i+1];
    }
    

    NumericMatrix out2 = NumericMatrix(n,2);
    
    for (i=0;i<n;i++)
    {
        out2(i,0)=data[i+1][0];
        out2(i,1)=data[i+1][1];
    }
    
    double out3 = hmin;
    
    // make the list for the output, containing the MLE, hazard, the bootstrap confidence intervals and -log likelihood
    
    List out = List::create(Rcpp::Named("MSE")=out1,Rcpp::Named("data")=out2,Rcpp::Named("minimum_c")=out3);
    // free memory
    
    for (i=0;i<n+1;i++)
        delete[] data[i];
    delete[] data;
    
    for (i=0;i<n+1;i++)
        delete[] bootstrap_data[i];
    delete[] bootstrap_data;
    
    delete[] cumw; delete[] cs; delete[] ff;  delete[] residual;
    
    delete[] tt; delete[] pp; delete[] tt_bootstrap; delete[] pp_bootstrap;
    delete[] ff_bootstrap; delete[] SLSE_bootstrap; delete[] SLSE_data; delete[] SLSE2;
    
    return out;
}

double f(double x)
{
    return exp(4 * (x - 0.5))/(1+exp(4 * (x - 0.5)));
}

double MSE_regression(int n, int ngrid, double grid[], int B, double **data, double **bootstrap_data,  double SLSE_data[], double SLSE2[], double SLSE_bootstrap[], double tt[], double pp[], double tt_bootstrap[], double pp_bootstrap[], double ff_bootstrap[], double residual[], double h, double h0, int seed)
{
    int i,j,iter,m_bootstrap,seed1;
    double MSE,*MSE1,*cumw,*cs;
    
    cumw = new double[n+1];
    cs   = new double[n+1];
    MSE1   = new double[B+1];
    
    seed1=seed;
    
    cumw[0]=cs[0]=0;

    for (iter=1;iter<=B;iter++)
    {
        seed1++;
        data_bootstrap(n,data,SLSE_data,residual,bootstrap_data,seed1);
        
        for (i=1;i<=n;i++)
        {
            cumw[i]= (double)i;
            cs[i]=cs[i-1]+bootstrap_data[i][1];
        }
        
        convexminorant(n,cumw,cs,ff_bootstrap);
        
        ff_bootstrap[0]=ff_bootstrap[1];
                
        j=0;

        
        for (i=1;i<=n;i++)
        {
            if (ff_bootstrap[i] > ff_bootstrap[i-1])
            {
                j++;
                tt_bootstrap[j]=bootstrap_data[i][0];
                pp_bootstrap[j]=ff_bootstrap[i]-ff_bootstrap[i-1];
            }
        
        }
        
        m_bootstrap = j;
        
        for (i=1;i<=m_bootstrap;i++)
            ff_bootstrap[i]=ff_bootstrap[i-1]+pp_bootstrap[i];
        
        // bootstrap SLSE at observation points
                
        regression_estimate(m_bootstrap,tt_bootstrap,ngrid,grid,h,h0,pp_bootstrap,ff_bootstrap,SLSE_bootstrap);
        
        MSE1[iter] =0;
        //for (i=1;i<=n;i++)
            //MSE1[iter] += SQR(SLSE_bootstrap[i]-SLSE_data[i]);
        
        for (i=1;i<=ngrid;i++)
            MSE1[iter] += SQR(SLSE_bootstrap[i]-SLSE2[i]);
    }
    
    
    
    MSE=0;
    for (iter=1;iter<=B;iter++)
        MSE += MSE1[iter];
    
    delete[] cumw; delete[] cs; delete[] MSE1;
    
    return MSE/B;
}


double fprime(int m, double tt[], double pp[], double h, double u)
{
    int i;
    double t,sum=0;
    
    if (u>=h && u<=1-h)
    {
        for (i=1;i<=m;i++)
        {
            t=(u-tt[i])/h;
            sum += Kprime(t)*pp[i]/SQR(h);
        }
    }
    else
    {
        if (u<h)
        {
            for (i=1;i<=m;i++)
            {
                t=(h-tt[i])/h;
                sum += Kprime(t)*pp[i]/SQR(h);
            }
        }
        else
        {
            for (i=1;i<=m;i++)
            {
                t=(1-h-tt[i])/h;
                sum += Kprime(t)*pp[i]/SQR(h);
            }
        }
    }
   
    return sum;
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

void data_bootstrap(int n, double **data, double SLSE_data[], double residual[], double **bootstrap_data, int seed)
{
    int    i,j;
    
    std::mt19937_64 gen(seed);
    std::uniform_int_distribution<int> dis(1,n);
    
    for (i=1;i<=n;i++)
    {
        bootstrap_data[i][0] = data[i][0];
        j=dis(gen);
        bootstrap_data[i][1] = SLSE_data[i]+residual[j];
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

void regression_estimate(int m, double tt[], int ngrid, double grid[], double h, double h0, double p[], double ff[], double fsmooth[])
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
                a += (KK(t)+(u-h)*K(t)/h)*p[k];
            }
            a += 0.5*SQR(u-h)*fprime(m,tt,p,h0,h);
        }
        
        if (u>1-h)
        {
            for (k=1;k<=m;k++)
            {
                t=(1-h-tt[k])/h;
                a += (KK(t)+(u-(1-h))*K(t)/h)*p[k];
            }
            a += 0.5*SQR(u-(1-h))*fprime(m,tt,p,h0,1-h);
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










