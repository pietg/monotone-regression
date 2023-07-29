//
//  Computation of credible intervals for monotone regression
//
//  Created by Piet Groeneboom on May 1, 2022.
//  Copyright (c) 2022 Piet Groeneboom. All rights reserved.

// [[Rcpp::plugins("cpp11")]]


#include <iostream>
#include <sstream>
#include <fstream>
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

int     compare (const void * a, const void * b);
void    gen_data(int n, double **data, int seed);
void    extract_nj_sigmahat(int n, int J,  int nJ[], double **data, double zeta[], double lambda_sq[], double ybar[], double *sigma);
void    projection_credible_interval(int J, int nJ[], double zeta[], double lambda_sq[], double ybar[], double sigma, double post_mean[], double post_sigma[]);
void    data_post_sample(int J, double post_mean[], double post_sigma[], double theta[], int seed);
void    compute_LSE(int J, int nJ[], double cumw[], double cs[], double theta[],
                 double theta_reduced[], double theta_monotone[]);
void    convexminorant(int n1, double cumw[], double cs[], double yy[]);
double  f0(double x);
double  f1(double x);


// [[Rcpp::export]]

List credible(int N, int J1, NumericVector zeta0, NumericVector lambda_sq0)
{
    int     n,ngrid,i,k,seed,J,*nJ;
    int     iter,iter2,NumIt,NumIt2;
    int     percentile1,percentile2,*percentage;
    double   **data,*grid,sigma,*LSE;
    double  *zeta,*lambda_sq,*ybar,*post_mean,*post_sigma,*theta,*theta_monotone,*theta_reduced;
    double  *lowbound,*upbound,**f3,*f4,*test;
    double  *cumw,*cs;
    
    Rcout << "For further information see:" << std::endl;
    Rcout << "Moumita Chakraborty and Subhashis Ghosal   Coverage of Credible Intervals in Nonparammatric Monotone Regression" << std::endl;
    Rcout << "Annals of Statistics (2021), vol. 49, pp. 1011-1028" << std::endl << std::endl;
    Rcout << "The program produces the credible (pointwise) 95%  intervals for a monotone regression function" << std::endl;
    Rcout << std::endl << std::endl;
    
    seed=1;
    
    n=(int)N;
    
    ngrid=100;

    //J=(int)(pow(n,1.0/3)*log(n));
    
    J=(int)J1;
    
    nJ =    new int[J+1];
    zeta =  new double[J+1];
    lambda_sq = new double[J+1];
    ybar =  new double[J+1];
    post_mean =  new double[J+1];
    post_sigma =  new double[J+1];
    theta =  new double[J+1];
    theta_reduced =  new double[J+1];
    theta_monotone =  new double[J+1];
    cumw =  new double[J+1];
    cs =  new double[J+1];
    
    cumw[0]=cs[0]=0;
    
    LSE =  new double[ngrid+1];
    
    for (i=1;i<=J;i++)
    {
        zeta[i]=(double)zeta0[i-1];
        lambda_sq[i]=(double)lambda_sq0[i-1];
    }
    
    NumIt=1000;
    NumIt2=1000;
    
    test= new double[NumIt+1];
    
    Rcout << "Number of observations:" << std::setw(7) << n << std::endl << std::endl;
    Rcout << "Number of posterior samples:" << std::setw(10) << NumIt2 << std::endl;
    
    //percentile1=(int)(0.034*NumIt2);
    //percentile2=(int)(0.966*NumIt2);
    
    percentile1=(int)(0.025*NumIt2);
    percentile2=(int)(0.975*NumIt2);
    
    data = new double *[n+1];
    for (i=0;i<n+1;i++)
        data[i]= new double[2];

    grid = new double[ngrid+1];
                                  
    for (i=0;i<=ngrid;i++)
        grid[i] = (double)i/ngrid;
    
    percentage = new int[ngrid+1];
    
    for (i=0;i<=ngrid;i++)
        percentage[i]=0;
                                  
    f3 = new double*[NumIt2+1];
    for (iter2=0;iter2<NumIt2+1;iter2++)
         f3[iter2] = new double[ngrid+1];
                                      
    f4 = new double[NumIt2+1];
                                  
    lowbound = new double[ngrid+1];
    upbound  = new double[ngrid+1];
    
    Rcout << "f_0(x) = )x^2+x/5" << std::endl;
    Rcout << "Observation distribution is Uniform on [0,1]" << std::endl << std::endl;
    
    Rcout << "1000 samples with 1000 posterior samples from each sample:" << std::endl << std::endl;
    
    //Rcout << "     Iteration  " << "  f_0(1)  " << "     lower bound  "<< "  upper bound  " << "#{F_0(1) not in interval}  " << std::endl << std::endl;
    
    for (iter=0;iter<NumIt;iter++)
    {
        seed++;
        
        test[iter]=0;
        gen_data(n,data,seed);
        
        extract_nj_sigmahat(n,J,nJ,data,zeta,lambda_sq,ybar,&sigma);
        
        projection_credible_interval(J,nJ,zeta,lambda_sq,ybar,sigma,post_mean,post_sigma);
        
        for (iter2=0;iter2<NumIt2;iter2++)
        {
            //printf("%5d\n",iter2+1);
            seed++;
            
            data_post_sample(J,post_mean,post_sigma,theta,seed);

            compute_LSE(J,nJ, cumw,cs,theta,theta_reduced,theta_monotone);
                    
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
            
            for (i=1;i<ngrid;i++)
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
        
        Rcout  << setw(10) << iter+1 << std::endl;
    }
    
    NumericMatrix out1 = NumericMatrix(J,2);
    
    for (i=0;i<=J-1;i++)
    {
      out1(i,0)=(i+0.5)/J;
      out1(i,1) = post_mean[i+1];
    }
    

    NumericMatrix out2 = NumericMatrix(ngrid-1,4);
    
    for (i=0;i<ngrid-1;i++)
    {
        out2(i,0)=grid[i+1];
        out2(i,1)=LSE[i+1];
        out2(i,2)=lowbound[i+1];
        out2(i,3)=upbound[i+1];
    }
    
    NumericMatrix out3 = NumericMatrix(ngrid-1,2);
    
    for (i=0;i<ngrid-1;i++)
    {
        out3(i,0)=grid[i+1];
        out3(i,1)=(double)percentage[i+1]/NumIt;
    }
    
    double out4 = sigma;
    
    // make the list for the output, containing the MLE, hazard, the bootstrap confidence intervals and -log likelihood
    
    List out = List::create(Rcpp::Named("posterior_mean")=out1,Rcpp::Named("CI_LSE")=out2,Rcpp::Named("percentages")=out3,Rcpp::Named("sigma")=out4);

    Rcout << std::endl;
    
    ofstream file0_("posterior_mean.txt");
    
    if (file0_.is_open())
    {
        for (i=0;i<J;i++)
        {
            file0_ << setprecision(11) << setw(20) << (i+0.5)/J;
            file0_ << setprecision(11) <<  setw(20) << post_mean[i+1];
            file0_ << "\n";
        }
        file0_.close();
    }
    
    
    ofstream file1_("CI_LSE.txt");
    
    if (file1_.is_open())
    {
        for (i=1;i<ngrid;i++)
        {
            file1_ << setprecision(10) << setw(20) << grid[i];
            file1_ << setprecision(10) << setw(20) << LSE[i];
            file1_ << setprecision(11) <<  setw(20) << lowbound[i] << setprecision(11) <<  setw(20) << upbound[i];
            file1_ << "\n";
        }
        file1_.close();
    }
    
    ofstream file2_("percentages.txt");
    
    if (file2_.is_open())
    {
        for (i=1;i<ngrid;i++)
        {
            file2_ << setprecision(10) << setw(20) << grid[i];
            file2_ << setprecision(11) <<  setw(20) << (double)percentage[i]/NumIt;
            file2_ << "\n";
        }
        file2_.close();
    }
    
    ofstream file3_("test.txt");
    
    if (file3_.is_open())
    {
        for (i=0;i<NumIt;i++)
        {
            file3_ << setprecision(10) << setw(20) << test[i];
            file3_ << "\n";
        }
        file3_.close();
    }

    // free memory
    
    for (i=0;i<n+1;i++)
        delete data[i];
    delete[] data;
    
    delete[] cumw; delete[] cs;

    delete[] LSE; delete[] nJ; delete[] zeta; delete[]ybar;
    delete[] post_mean; delete[] post_sigma; delete[] theta;
    delete[] theta_monotone; delete[] theta_reduced; delete[] test;
    
    for (iter2 = 0;iter2 < NumIt2;iter2++)
        delete[] f3[iter2];
    delete[] f3;
    
    delete[] f4;
    
    return out;
}

double f0(double x)
{
  return SQR(x)+x/5;
}

double f1(double x)
{
    return exp(4*(x-0.5))/(1+exp(4* (x-0.5)));
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

void compute_LSE(int J, int nJ[], double cumw[], double cs[], double theta[],
                  double theta_reduced[], double theta_monotone[])
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
        {
            if (J_reduced==0)
                index[i]=1;
            else
                index[i] = J_reduced;
        }
    }
    
    convexminorant(J_reduced,cumw,cs,theta_reduced);
        
    for (i=1;i<=J;i++)
        theta_monotone[i] = theta_reduced[index[i]];
        
    delete[] index;
}


void extract_nj_sigmahat(int n, int J, int nJ[], double **data, double zeta[],
                         double lambda_sq[], double ybar[], double *sigma)
{
    int i,j;
    double sum;
    
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
        
    sum =0;
    for (i=1;i<=n;i++)
    {
        for (j=1;j<=J;j++)
        {
            if ((double)(j-1)/J< data[i][0] && data[i][0]<= (double)j/J)
                sum += (data[i][1]-zeta[j])*(data[i][1]-zeta[j]-nJ[j]*lambda_sq[j]*(ybar[j]-zeta[j])/(1+nJ[j]*lambda_sq[j]));
        }
    }
    *sigma = sqrt(sum/n);
}

void projection_credible_interval(int J, int nJ[], double zeta[], double lambda_sq[], double ybar[], double sigma, double post_mean[], double post_sigma[])
{
    int i;
    
    for (i=1;i<=J;i++)
    {
        post_mean[i] = (nJ[i]*ybar[i]+zeta[i]/lambda_sq[i])/(nJ[i]+1/lambda_sq[i]);
        post_sigma[i] = sigma/sqrt(nJ[i]+1/lambda_sq[i]);
    }
}

void data_post_sample(int J, double post_mean[], double post_sigma[], double theta[], int seed)
{
    int i;
    
    std::mt19937 gen(seed);

    for (i=1;i<=J;i++)
    {
        std::normal_distribution<double> dis_normal(post_mean[i],post_sigma[i]);
        theta[i] = dis_normal(gen);
    }
}

void convexminorant(int n, double cumw[], double cs[], double yy[])
{
    int    i,j,k;
    
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
                for (k=j;k<i;k++)
                    yy[k] = yy[i];
            }
        }
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


