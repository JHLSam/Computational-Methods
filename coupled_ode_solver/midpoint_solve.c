/*
Program that solves for coupled Ordinary Differential Equations(ODE), in
context of white dwarf star equations of state. 

Compile: gcc -Wall -lm <filename> -o <chosen_exe_name>

Input: No user input required. All feeding parameters are specified as advised.
User input scans removed. 

Output: Prints out logarithm product of linearly spaced central densities,
Radius and Mass of system. 
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//Initialization of constants
double m_p = 1.67e-27;
double m_e = 9.11e-31;
double c = 3.00e8;
double h = 6.626e-34;
double G = 6.67e-11;
double M_solar = 1.98892e30;

//Declaration of variables. Naming correspond to formulas. 
double rho_c, k_r, k_n, k_1;

//coupled ODE equation 1
double dtheta_dr (double m, double theta, double r) {
    double exponent = ((double)2)/3;
    double numerator, denominator;
    
    numerator = -3*pow(pow(theta,exponent)+1,1.5)*G*m*rho_c*theta;
    denominator = k_n*pow((k_r/k_n),5)*(4*pow(theta,exponent)+5)
    *pow(theta,exponent)*r*r;
    return numerator/denominator;
}

//coupled ODE equation 2
double dm_dr (double theta, double r) {
    double rho = rho_c*theta;
    return 4*M_PI*r*r*rho;
} 

int main () {
    //Declaration of variables
    double r,rho0,r_halfnew,r_new,delta_r,theta,theta_halfnew,theta_new,m,m_halfnew,m_new;
    //For ease of use, for the 'power' factors in the equations.
    double exponent = ((double)1)/3; 
    //Initialization of the previously declared corresponding global constants
    k_r = (h*c/8.0)*pow(3/M_PI,exponent)*1/pow(2*m_p,4*exponent);
    k_n = (pow(h,2)/(20*m_e))*pow(3/M_PI,2*exponent)*1/pow(2*m_p,5*exponent);
    rho_c = pow((k_r/k_n),3); //rho = rho_c*theta 

    k_1 = k_r * pow(rho_c,4*exponent);
    
    delta_r = 10; //step size

    /* 
    output from python print(numpy.linspace(math.log10(1e7),math.log10(1e16)).
    */
    double arr[] = {7.0000000,7.31034483,7.62068966,7.93103448,8.24137931,
        8.55172414,8.86206897,9.17241379,9.48275862,9.79310345,10.10344828,
        10.4137931,10.72413793,11.03448276,11.34482759,11.65517241,11.96551724,
        12.27586207,12.5862069,12.89655172,13.20689655,13.51724138,13.82758621,
        14.13793103,14.44827586,14.75862069,15.06896552,15.37931034,15.68965517,
        16.0000000};
        FILE *fp;
        fp = fopen("as05-lim-44589442-massvsradius.dat", "w");
    fprintf(fp, "%%log of Central Density(kg/m³) Radius(km) Mass(Solar Mass)\n");
    
    for(int i = 0; i<30; i++){ //iterate 30 times as per array index arr[0]:arr[29]
        rho0 = pow(10, arr[i]); //exponent power for each index element
        r = 100; //initial radius
        m = 0; //assumption/condition as advised
        theta = rho0/rho_c; //rho = rho_c*theta
    
        while(1){
        //As per the formula, basically midpoint method. Variables correspond.
            theta_halfnew = theta + (delta_r/2)*dtheta_dr(m,theta,r);
            m_halfnew = m + (delta_r/2)*dm_dr(theta,r);
            r_halfnew = r + delta_r/2;
    
            if(theta_halfnew < 0)//as per given exit condition
                break;
            
            theta_new = theta + delta_r*dtheta_dr(m_halfnew,theta_halfnew,r_halfnew);
            m_new = m + delta_r*dm_dr(theta_halfnew,r_halfnew);
            r_new = r + delta_r;
            
            //Update old variables for iteration
            theta = theta_new;
            m = m_new;
            r = r_new;
            
            if(theta < 0)//as per given exit condition
                break;
        }
        /*log_e(log_10(rho¹⁰))/log_e(10)) by log exponent and product rule simply cancels to log_10(rho)
        units defined in line 72 already*/
        fprintf(fp,"\t\t\t%lf \t\t %lf \t %lf\n", log(rho0)/log(10), r/1000, m/M_solar);
    }
    
    return 0;
}