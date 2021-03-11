#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <cmath>
#define numberofData 27
typedefstruct{
double r_t[numberofData];
double n_t[numberofData];
double h_t[numberofData];
double flat_portion[numberofData];
} f_parameters;
typedefstruct{
float mean;
float stdDev;
} array_properties;
typedefstruct{
double C_const =9.9;
float C_range =320;
float Cr_const =0.07;
float Cr_range =6;
float Cn_const =0.26;
float Cn_range =5;
float step_incr =0.01;
} trials;




// Function prototype declarations
f_parameters Defining_Parameters(float Web[],float Thickness[],float R,int
N[]);
void Coefficient_array(float Web[],float Thickness[], f_parameters
variables,int Fy,int sin_theta,float Ch,float FEA[],float uncertainty,
trials Coefficients);
double Equation_substitution(float Web[],float current_Thickness,float
C,float Cr,float Cn,int Fy,int sin_theta,float current_n_t,float
current_r_t,float current_h_t,float Ch,float FEA[],float uncertainty);
void Comparison(float FEA[],float results_array[],float C,float Cr,float
Cn,float uncertainty);
float Reliability_analysis(float FEA[],float results_array[]);
array_properties StdDev_equation(float results_array[]);
int main()
{
float R =1.2;
int N[numberofData]={};
float Web[numberofData]={};
float Thickness[numberofData]={};
float FEA[numberofData]={};
int sin_theta =1, Fy =284;
float Ch =0.031, uncertainty =4.4;
f_parameters variables;
trials Coefficients;
for(int i =0; i < numberofData; i++) 
{printf("%f ", FEA[i]);}
printf("\n\n");
variables = Defining_Parameters(Web, Thickness,R,N);
Coefficient_array(Web, Thickness, variables, Fy, sin_theta, Ch, FEA,
uncertainty, Coefficients);
return0;
}
f_parameters Defining_Parameters(float Web[],float Thickness[],float R,int
N[])
{
int i =0;
f_parameters Parameters;
for(i =0; i< numberofData; i++)
{
Parameters.flat_portion[i]= Web[i]-2* Thickness[i];
Parameters.r_t[i]= R / Thickness[i];
Parameters.n_t[i]= N[i]/ Thickness[i];
Parameters.h_t[i]= Parameters.flat_portion[i]/ Thickness[i];
};
return Parameters;
}
void Coefficient_array(float Web[],float Thickness[], f_parameters
variables,int Fy,int sin_theta,float Ch,float FEA[],float uncertainty,
trials Coefficients)
{
float results_array[numberofData]={0},Cr, Cn, C,current_n_t,
current_r_t,current_h_t, current_Thickness, y =0.00;
int i,j,k,l =0;
C = Coefficients.C_const;
for(i =0; i < Coefficients.C_range; i++)
{
Cr = Coefficients.Cr_const;
for(j =0; j < Coefficients.Cr_range; j++)
{ Cn = Coefficients.Cn_const;
for(k =0; k < Coefficients.Cn_range; k++)
{for(l =0; l < numberofData; l++)
{
current_r_t = variables.r_t[l];
current_n_t = variables.n_t[l];
current_h_t = variables.h_t[l];
current_Thickness = Thickness[l];
y = Equation_substitution(Web, current_Thickness, C,
Cr, Cn, Fy, sin_theta, current_n_t, current_r_t, current_h_t, Ch, FEA,
uncertainty);
results_array[l]= y;
}
Comparison(FEA, results_array, C, Cr, Cn, uncertainty);
Cn += Coefficients.step_incr;
}
Cr += Coefficients.step_incr;
}
C += Coefficients.step_incr;}
}
double Equation_substitution(float Web[],float current_Thickness,float
C,float Cr,float Cn,int Fy,int sin_theta,float current_n_t,float
current_r_t,float current_h_t,float Ch,float FEA[],float uncertainty)
{
double y =0;
y = C*current_Thickness*current_Thickness*Fy*sin_theta*(1-
Cr*sqrt(current_r_t))*(1+ Cn*sqrt(current_n_t))*(1- Ch*sqrt(current_h_t));
y /=1000;
return y;
}
void Comparison(float FEA[],float results_array[],float C,float Cr,float
Cn,float uncertainty)
{
int k, counter =0;
float R =0;
for(k =0; k < numberofData; k++)
{if(((results_array[k]- uncertainty)<= FEA[k])&&(FEA[k]<=(results_array[k]+
uncertainty)))
{counter++;}
}
if(counter >=19)
{R = Reliability_analysis(FEA, results_array);
if(R >=2.6)
{printf("%d Set C: %f Cr: %f Cn: %f R:%f \n", counter, C, Cr, Cn, R);}
}
}
float Reliability_analysis(float FEA[],float results_array[])
{
int n = numberofData, m = n –1, i =0;
float R_Rp[numberofData]={0}, Reliability =0;
array_properties numbers;
for(i =0; i < numberofData; i++)
{R_Rp[i]= FEA[i]/ results_array[i];}
numbers = StdDev_equation(R_Rp);
float Vp = numbers.stdDev / numbers.mean, Pm = numbers.mean, Cp =(((1+1/
n)*m)/(m -2)), Mm =1.1,Fm =1.0, phi =0.85, Vm =0.1, Vf =0.05;
printf("Mean: %f StdDev: %f\n", numbers.mean, numbers.stdDev);
Reliability = log(Mm*Fm*Pm /(0.657*phi))/ sqrt(Vm*Vm + Vf*Vf + Cp*Vp*Vp
+0.21*0.21);
return Reliability;
}
array_properties StdDev_equation(float results_array[])
{
float sum =0, mean, standardDeviation =0, stdDev =0;
int m,n =0;
array_properties numbers;
for(m =0; m < numberofData; m++)
{sum += results_array[m];}
mean = sum / numberofData;
for(n =0; n < numberofData; n++)
{standardDeviation += pow(results_array[n]- mean,2);}
stdDev = sqrt(standardDeviation / numberofData);
numbers.mean = mean;
numbers.stdDev = stdDev;
return numbers;
