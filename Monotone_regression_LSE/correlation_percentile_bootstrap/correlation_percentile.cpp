//
//  Computation of credible intervals for monotone regression
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
void    data_bootstrap(int n, double **data, double **bootstrap_data, int seed);
void    extract_nj(int n, int J, int nJ[], double **data, double ybar[]);
void    compute_LSE(int J, int nJ[], double cumw[], double cs[],
                    double theta[], double theta_reduced[], double theta_monotone[]);


// [[Rcpp::export]]

List correlation_percentile(int N)
{
    int     n,ngrid,i,k,J,*nJ,seed;
    int     iter,iter2,NumIt,NumIt2;
    int     percentile1,percentile2,*percentage;
    double  *lowbound,*upbound,**f3,*f4;
    double  **data,**bootstrap_data,*LSE,*ybar,*theta_reduced,*theta_monotone;
    double  *cumw,*cs,*test;
    double  *grid;

    Rcout << "For further information see:" << std::endl;
    Rcout << "Piet Groeneboom and Geurt Jongbloed:   Bootstrap confidence intervals and credible intervals in monotone regression" << std::endl;
    Rcout << "The program produces percentile bootstrap confidence intervals for a monotone regression function" << std::endl;
    
    seed=1;
    
    n=(int)N;
    
    ngrid=100;
    
    NumIt=1000;
    NumIt2=1000;
   
    J=(int)(pow(n,1.0/3)*log(n));
    //J=n;
    
    nJ =    new int[J+1];
    ybar =  new double[J+1];
    theta_reduced =  new double[J+1];
    theta_monotone =  new double[J+1];
        
    ngrid=100;
    
    LSE =  new double[ngrid+1];
    
    grid = new double[ngrid+1];
    
    for (i=0;i<=ngrid;i++)
        grid[i]=(double)i/ngrid;
    
    //percentile1=round(0.034*NumIt2);
    //percentile2=round(0.966*NumIt2);
    
    percentile1=round(0.025*NumIt2);
    percentile2=round(0.975*NumIt2);
    
    printf("percentiles:\n\n");
    printf("%d  %d\n\n",percentile1,percentile2);
    
    data = new double *[n+1];
    for (i=0;i<n+1;i++)
        data[i] = new double[2];
    
    bootstrap_data = new double *[n+1];
    for (i=0;i<n+1;i++)
        bootstrap_data[i] = new double[2];
    
    cumw = new double[J+1];
    cs = new double[J+1];
    
    test = new double[NumIt+1];
 
    f3 = new double*[NumIt2+1];
    
    for (iter2=0;iter2<NumIt2+1;iter2++)
        f3[iter2] = new double[ngrid+1];
        
    f4 = new double[NumIt2+1];
    
    percentage = new int[ngrid+1];
    
    for (i=0;i<=ngrid;i++)
        percentage[i]=0;
        
    lowbound = new double[ngrid+1];
    upbound  = new double[ngrid+1];
        
    cumw[0]=cs[0]=0;
    
    for (iter=0;iter<NumIt;iter++)
    {
        seed++;
        
        test[iter]=0;
        
        //printf("%5d\n",iter+1);

        gen_data(n,data,seed);
        
        for (iter2=0;iter2<NumIt2;iter2++)
        {
            seed++;
    
            data_bootstrap(n,data,bootstrap_data,seed);
            
            extract_nj(n,J,nJ,bootstrap_data,ybar);
                        
            compute_LSE(J,nJ,cumw,cs,ybar,theta_reduced,theta_monotone);
                        
            for (i=1;i<ngrid;i++)
            {
                for (k=1;k<=J;k++)
                {
                    if ((double)(k-1)/J<grid[i] && grid[i]<=(double)k/J)
                        LSE[i] = theta_monotone[k];
                }
            }
            if (LSE[50]-f0(0.5)<=0)
                test[iter] +=  1.0/NumIt2;
            
            for (i=1;i<=ngrid;i++)
                f3[iter2][i]= LSE[i];
            
        }
        
        for (i=1;i<ngrid;i++)
        {
            for (iter2=0;iter2<NumIt2;iter2++)
                f4[iter2]=f3[iter2][i];
            
            qsort(f4,NumIt2,sizeof(double),compare);
            
            lowbound[i] = f4[percentile1-1];
            upbound[i]  = f4[percentile2-1];
            
            if (f0(grid[i])<lowbound[i] || f0(grid[i])>upbound[i])
                percentage[i]++;
        }
        
        Rcout << setw(6) << iter+1 << std::endl;
    
    }

    NumericMatrix out0 = NumericMatrix(ngrid+1,2);
    
    out0(0,0)=0;
    out0(0,1) = 0;
    
    for (i=1;i<=ngrid;i++)
    {
      out0(i,0)=grid[i];
      out0(i,1) = LSE[i];
    }
    
    NumericMatrix out1 = NumericMatrix(ngrid+1,2);
    
    for (i=0;i<ngrid-1;i++)
    {
        out1(i,0)=grid[i+1];
        out1(i,1)=f0(grid[i+1]);
    }

    NumericMatrix out2 = NumericMatrix(ngrid-1,4);
    
    for (i=0;i<ngrid-1;i++)
    {
        out2(i,0)=grid[i+1];
        out2(i,1)=f0(grid[i+1]);
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
    
    Rcout << "Making output list" << std::endl;
    
    // make the list for the output, containing the MLE, hazard, the bootstrap confidence intervals and -log likelihood
    
    List out = List::create(Rcpp::Named("LSE")=out0,Rcpp::Named("f0")=out1,Rcpp::Named("CI_LSE")=out2,Rcpp::Named("percentages")=out3,Rcpp::Named("data")=out4);

    
    ofstream file0_("LSE.txt");
    
    if (file0_.is_open())
    {
        for (i=1;i<=ngrid;i++)
        {
            file0_ << setprecision(11) << setw(20) << grid[i];
            file0_ << setprecision(11) <<  setw(20) << LSE[i];
            file0_ << "\n";
        }
        file0_.close();
    }
    
    
    ofstream file2_("CI_LSE.txt");
    
    if (file2_.is_open())
    {
        for (i=1;i<ngrid;i++)
        {
            file2_ << setprecision(10) << setw(20) << grid[i];
            file2_ << setprecision(10) << setw(20) << LSE[i];
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
    
    delete[] cumw; delete[] cs; delete[] ybar; delete[] percentage;
    
    delete[] LSE; delete[] nJ;
    delete[] theta_monotone; delete[] theta_reduced; delete[] test;
    
    return out;
}

double f0(double x)
{
  return SQR(x)+x/5;
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
    obs[i].y = SQR(x) + x/5 + dis_normal(gen);
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

void data_bootstrap(int n, double **data, double **bootstrap_data, int seed)
{
  int    i,j;
  
  std::mt19937_64 gen(seed);
  std::uniform_int_distribution<int> dis(1,n);
  
  for (i=1;i<=n;i++)
  {
    j=dis(gen);
    bootstrap_data[i][0] = data[j][0];
    bootstrap_data[i][1] = data[j][1];
  }
  
}

void extract_nj(int n, int J, int nJ[], double **data, double ybar[])
{
  int i,j;
  
  for (j=1;j<=J;j++)
  {
    nJ[j]=0;
    ybar[j]=0;
  }
  
  for (i=1;i<=n;i++)
  {
    for (j=1;j<=J;j++)
    {
      if ((double)(j-1)/J< data[i][0] && data[i][0]<= (double)j/J)
      {
        ybar[j] += data[i][1];
        nJ[j]++;
      }
    }
  }
  
  for (j=1;j<=J;j++)
  {
    if (nJ[j]>0)
      ybar[j] /=nJ[j];
  }
}

void compute_LSE(int J, int nJ[], double cumw[], double cs[], double theta[],double theta_reduced[], double theta_monotone[])
{
  int i,J_reduced,*index;
  
  index = new int[J+1];
  
  J_reduced=0;
  index[0]=0;
  
  for (i=1;i<=J;i++)
  {
    if (nJ[i]>0)
    {
      J_reduced++;
      index[i] = J_reduced;
      cumw[index[i]]= cumw[index[i]-1] +nJ[i];
      cs[index[i]]= cs[index[i]-1]+nJ[i]*theta[i];
    }
    else
      index[i] = J_reduced;
  }
  
  convexminorant(J_reduced,cumw,cs,theta_reduced);
  
  for (i=1;i<=J;i++)
    theta_monotone[i] = theta_reduced[index[i]];
  
  delete[] index;
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

