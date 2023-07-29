//
//  Computation of bootstrap intervals for monotone regression
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
void    gen_data(int n, double **data, int seed);
void    convexminorant(int n, double cumw[], double cs[], double yy[]);
double  KK(double x);
double  K(double x);
double  Kprime(double x);
double  compute_sigma(int n, double **data, double ff[]);
void    data_bootstrap(int n, double **data, double SLSE_data[],
                       double residual[], double **bootstrap_data, int seed);
void    center_residuals(int n, double residual[]);
void    regression_estimate(int m, double tt[], int ngrid, double grid[], double h,
                            double h0, double p[], double ff[], double fsmooth[]);
double  fprime(int m, double tt[], double pp[], double h, double u);
double  f0(double x);


// [[Rcpp::export]]

List smooth_bootstrap(int N, double bandwidth)
{
    int     n,m,m_bootstrap,ngrid,i,j,seed;
    int     iter,iter2,NumIt,NumIt2,*index;
    int     percentile1,percentile2,*percentage;
    double  *lowbound,*upbound,**f3,*f4,*test;
    double  **data,**bootstrap_data,*LSE,*LSE_bootstrap,*residual,*residual_bootstrap;
    double  *cumw,*cs,*ff,*tt,*pp,*tt_bootstrap,*pp_bootstrap,*ff_bootstrap;
    double  h,h0,*grid,*SLSE_data,*SLSE,*LSE_data,*data0;
    
    
    Rcout << "For further information see:" << std::endl;
    Rcout << "Piet Groeneboom and Geurt Jongbloed: Credible intervals and bootstrap confidence intervals in monotone regression" << std::endl;
    Rcout << "The program produces 95% bootstrap confidence intervals for a monotone regression function, using a smoothed bootstrap methods" << std::endl;
    
    seed=1;
    
    n=(int)N;
    h= (double)bandwidth;
    h0= 0.7*pow(n,-1.0/7);
    
    ngrid=100;
    m=1;
    m_bootstrap =1;
    
    NumIt=1000;
    NumIt2=1000;
   
    grid = new double[ngrid+1];
    
    for (i=0;i<=ngrid;i++)
        grid[i]=(double)i/ngrid;
    
    percentile1=round(0.025*NumIt2);
    percentile2=round(0.975*NumIt2);
    
    printf("percentiles:\n\n");
    printf("%d  %d\n\n",percentile1,percentile2);
    
    data = new double *[n+1];
    for (i=0;i<n+1;i++)
        data[i] = new double[2];
    
    data0 = new double[n+1];
    
    bootstrap_data = new double *[n+1];
    for (i=0;i<n+1;i++)
        bootstrap_data[i] = new double[2];

    index = new int[n+2];
    cumw = new double[n+1];
    cs = new double[n+1];
    ff = new double[n+1];
    
    tt = new double[n+1];
    pp = new double[n+1];
    
    tt_bootstrap = new double[n+1];
    pp_bootstrap = new double[n+1];
    ff_bootstrap = new double[n+1];
    
    SLSE =new double[ngrid+1];
    SLSE_data = new double[n+1];
    
    LSE =new double[ngrid+1];
    LSE_bootstrap = new double[ngrid+1];
    LSE_data = new double[n+1];
    
    residual = new double[n+1];
    residual_bootstrap = new double[n+1];
    
    f3 = new double*[NumIt2+1];
    
    for (iter2=0;iter2<NumIt2+1;iter2++)
        f3[iter2] = new double[ngrid+1];
        
    f4 = new double[NumIt2+1];
    
    test = new double[NumIt+1];
    
    percentage = new int[ngrid+1];
    
    for (i=0;i<=ngrid;i++)
        percentage[i]=0;
        
    lowbound = new double[ngrid+1];
    upbound  = new double[ngrid+1];
        
    cumw[0]=cs[0]=0;
    
    tt[0]=pp[0]=0;
    tt_bootstrap[0]=pp_bootstrap[0]=0;
    
    data[0][0]=bootstrap_data[0][0]=0;
    data[0][1]=bootstrap_data[0][1]=0;
    
    ff[0] = ff_bootstrap[0]=0;
    
    for (iter=0;iter<NumIt;iter++)
    {
        seed++;
        
        Rcout  << setw(10) << iter+1 << std::endl;
        
        test[iter] = 0;
        
        //printf("%5d\n",iter+1);

        gen_data(n,data,seed);
        
        for (i=1;i<=n;i++)
            data0[i]=data[i][0];
        
        for (i=1;i<=n;i++)
        {
            cumw[i]= (double)i;
            cs[i]=cs[i-1]+data[i][1];
        }
        
        convexminorant(n,cumw,cs,ff);
        
        for (i=1;i<=n;i++)
            LSE_data[i]=ff[i];
        
        j=0;
        ff[0]=ff[1];
        
        for (i=1;i<=n;i++)
        {
            if (ff[i] > ff[i-1])
            {
                j++;
                tt[j]=data[i][0];
                pp[j]=ff[i]-ff[i-1];
                index[j]=i;
            }
        }
        
        m=j;
        
        index[0]=0;
             
        for (i=1;i<=m;i++)
            ff[i]=ff[i-1]+pp[i];
        
        for (i=0;i<=ngrid;i++)
        {
            if (grid[i]<tt[1])
                LSE[i]=ff[1];
            for (j=1;j<m;j++)
            {
                if (tt[j]<= grid[i] && grid[i]<tt[j+1])
                    LSE[i] = ff[j];
            }
            if (tt[m]<= grid[i])
                LSE[i]=ff[m];
        }
        
        regression_estimate(m,tt,ngrid,grid,h,h0,pp,ff,SLSE);
        regression_estimate(m,tt,n,data0,h,h0,pp,ff,SLSE_data);
        
        residual[0]=0;
        for (i=1;i<=n;i++)
            residual[i]=data[i][1]-SLSE_data[i];
        
        center_residuals(n,residual);
        
        for (iter2=0;iter2<NumIt2;iter2++)
        {
            seed++;
            
            data_bootstrap(n,data,SLSE_data,residual,bootstrap_data,seed);
            
            for (i=1;i<=n;i++)
            {
                cumw[i]= (double)i;
                cs[i]=cs[i-1]+bootstrap_data[i][1];
            }
            
            convexminorant(n,cumw,cs,ff_bootstrap);
            
            ff_bootstrap[0]=ff_bootstrap[1];
            
            //sigma_bootstrap = compute_sigma(n,bootstrap_data,ff_bootstrap);
            
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
            
            for (i=0;i<=ngrid;i++)
            {
                if (grid[i]<tt_bootstrap[1])
                    LSE_bootstrap[i]=ff_bootstrap[1];
                for (j=1;j<m_bootstrap;j++)
                {
                    if (tt_bootstrap[j]<= grid[i] && grid[i]<tt_bootstrap[j+1])
                        LSE_bootstrap[i] = ff_bootstrap[j];
                }
                if (tt_bootstrap[m]<= grid[i])
                    LSE_bootstrap[i]=ff_bootstrap[m_bootstrap];
            }
            
            if (LSE_bootstrap[ngrid/2]-SLSE[ngrid/2]<=LSE[ngrid/2]-f0(0.5))
                test[iter] += 1.0/NumIt2;
           
            for (i=1;i<ngrid;i++)
                f3[iter2][i] = LSE_bootstrap[i]-SLSE[i];
        }
        
        lowbound[0]=LSE[0]-f4[percentile2-1];
        upbound[0]=LSE[0]-f4[percentile1-1];
        
        for (i=1;i<ngrid;i++)
        {
            for (iter2=0;iter2<NumIt2;iter2++)
                f4[iter2]=f3[iter2][i];
            
            qsort(f4,NumIt2,sizeof(double),compare);
            
            lowbound[i] = LSE[i]-f4[percentile2-1];
            upbound[i]  = LSE[i]-f4[percentile1-1];
            
            lowbound[i]=fmax(lowbound[i-1],lowbound[i]);
            upbound[i]=fmax(upbound[i-1],upbound[i]);
        
        }
        
        // check whether real value is inside onterval
        
        for (i=1;i<ngrid;i++)
        {
            if (f0(grid[i])<lowbound[i] || f0(grid[i])>upbound[i])
                percentage[i]++;
        }
    }
    
    NumericMatrix out0 = NumericMatrix(m,2);
    
    for (i=0;i<m;i++)
    {
      out0(i,0)=tt[i+1];
      out0(i,1) = ff[i+1];
    }
    
    NumericMatrix out1 = NumericMatrix(ngrid,2);
    
    for (i=0;i<ngrid;i++)
    {
      out1(i,0)=grid[i+1];
      out1(i,1) = SLSE[i+1];
    }
    

    NumericMatrix out2 = NumericMatrix(ngrid-1,4);
    
    for (i=0;i<ngrid-1;i++)
    {
        out2(i,0)=grid[i+1];
        out2(i,1)=SLSE[i+1];
        out2(i,2)=lowbound[i+1];
        out2(i,3)=upbound[i+1];
    }
    
    NumericMatrix out3 = NumericMatrix(ngrid-1,2);
    
    for (i=0;i<ngrid-1;i++)
    {
        out3(i,0)=grid[i+1];
        out3(i,1)=(double)percentage[i+1]/NumIt;
    }
    
    NumericMatrix out4 = NumericMatrix(n,2);
    
    for (i=0;i<n;i++)
    {
        out4(i,0)=data[i+1][0];
        out4(i,1)=data[i+1][1];
    }
    
    NumericVector out5 = NumericVector(NumIt);
    
    for (i=0;i<NumIt;i++)
        out5(i)=test[i];
    
    Rcout << "Making output list" << std::endl;
    
    // make the list for the output, containing the MLE, hazard, the bootstrap confidence intervals and -log likelihood
    
    List out = List::create(Rcpp::Named("LSE")=out0,Rcpp::Named("SLSE")=out1,Rcpp::Named("CI_LSE")=out2,Rcpp::Named("percentages")=out3,Rcpp::Named("data")=out4,Rcpp::Named("test")=out5);

    
    ofstream file0_("LSE.txt");
    
    if (file0_.is_open())
    {
        for (i=1;i<=m;i++)
        {
            file0_ << setprecision(11) << setw(20) << tt[i];
            file0_ << setprecision(11) <<  setw(20) << ff[i];
            file0_ << "\n";
        }
        file0_.close();
    }
    
    ofstream file1_("SLSE.txt");
    
    if (file1_.is_open())
    {
        for (i=1;i<=ngrid;i++)
        {
            file1_ << setprecision(11) << setw(20) << grid[i];
            file1_ << setprecision(11) <<  setw(20) << SLSE[i];
            file1_ << "\n";
        }
        file1_.close();
    }
    
    
    ofstream file2_("CI_LSE.txt");
    
    if (file2_.is_open())
    {
        for (i=1;i<ngrid;i++)
        {
            file2_ << setprecision(10) << setw(20) << grid[i];
            file2_ << setprecision(10) << setw(20) << SLSE[i];
            file2_ << setprecision(11) <<  setw(20) << lowbound[i] << setprecision(11) <<  setw(20) << upbound[i];
            file2_ << "\n";
        }
        file2_.close();
    }
    
    ofstream file3_("percentages.txt");
    
    if (file3_.is_open())
    {
        for (i=1;i<ngrid;i++)
        {
            file3_ << setprecision(10) << setw(20) << grid[i];
            file3_ << setprecision(11) <<  setw(20) << (double)percentage[i]/NumIt;
            file3_ << "\n";
        }
        file3_.close();
    }
    
    ofstream file4_("data.txt");
    
    if (file4_.is_open())
    {
        for (i=1;i<=n;i++)
        {
            file4_ << setprecision(10) << setw(20) << data[i][0];
            file4_ << setprecision(11) <<  setw(20) << data[i][1];;
            file4_ << "\n";
        }
        file4_.close();
    }
    
    ofstream file5_("test.txt");
    
    if (file5_.is_open())
    {
        for (i=0;i<NumIt;i++)
        {
            file5_ << setprecision(11) <<  setw(20) << test[i];
            file5_ << "\n";
        }
        file5_.close();
    }

    
    Rcout << std::endl;

    // free memory
    
    for (i=0;i<n+1;i++)
        delete[] data[i];
    delete[] data;
    
    for (i=0;i<n+1;i++)
        delete[] bootstrap_data[i];
    delete[] bootstrap_data;
    
    delete[] cumw; delete[] cs; delete[] ff; delete[] data0;
    
    delete[] tt; delete[] pp; delete[] tt_bootstrap; delete[] pp_bootstrap;
    delete[] ff_bootstrap; delete[] SLSE; delete[] SLSE_data, delete[] test;
    delete[] LSE_data, delete[] residual; delete[] residual_bootstrap;
    
    return out;
}

double f0(double x)
{
  return SQR(x)+x/5;
}

void center_residuals(int n, double residual[])
{
  int i;
  double mean;
  
  mean=0;
  for (i=1;i<=n;i++)
    mean += residual[i];
  
  mean /=n;
  
  for (i=1;i<=n;i++) residual[i] -= mean;
}


void gen_data(int n, double **data, int seed)
{
    int    i;
    double x;
    data_object *obs;

    obs = new data_object[n];
    
    std::mt19937 gen(seed);
    std::uniform_real_distribution<> dis_unif(0,1);
    std::normal_distribution<double> dis_normal(0.0,0.1);
    
    for (i=0;i<n;i++)
    {
        x = obs[i].x = dis_unif(gen);
        obs[i].y = f0(x) + dis_normal(gen);
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

double compute_sigma(int n, double **data, double ff[])
{
    int i;
    double sum=0;
    
    for (i=1;i<=n;i++)
        sum += SQR(data[i][1]-ff[i]);
    
    return sqrt(sum/n);
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
                //a += (KK(t)+(u-h)*K(t)/h+0.5*SQR(u-h)*Kprime(t)/SQR(h))*p[k];
                a += (KK(t)+(u-h)*K(t)/h)*p[k];
            }
            a += 0.5*SQR(u-h)*fprime(m,tt,p,h0,h);
        }
        
        if (u>1-h)
        {
            for (k=1;k<=m;k++)
            {
                t=(1-h-tt[k])/h;
                //a += (KK(t)+(u-(1-h))*K(t)/h+0.5*SQR(u-(1-h))*Kprime(t)/SQR(h))*p[k];
                a += (KK(t)+(u-(1-h))*K(t)/h)*p[k];
            }
            a += 0.5*SQR(u-(1-h))*fprime(m,tt,p,h0,1-h);
        }
        fsmooth[i]=ff[0]+a;
    }
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
        //bootstrap_data[i][1] = f0(data[i][0])+residual[j];
    }
    
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










