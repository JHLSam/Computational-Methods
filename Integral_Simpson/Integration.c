/* 

Integration computation with Simpson's method 

Compile: gcc -Wall <filename> -o <chosen_exe_name>

Input: Queries user input for floating point numbers for the length of 
corrugated sheet(Lc),height of corrugated sheet(D) and wavelength(Lp) 

Output: Prints out the corresponding nth-step integral(N), the computed required
length of metal sheet(Lf), and the difference between successive computed 
integrals of order 2*N - N("Error").
*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define PI 3.1415927

double integrand(double Lc, double D, double Lp, double x){
    //Integrand to be integrated
    return ((double) sqrt(1+pow((2*PI*D/Lp),2)
        *pow(cos(2*PI*x/Lp),2))); //
}

double relativeError(double previous, double current) {
    /* Function for computing absolute relative difference between successive 
    terms. under iteration. */
    return fabs((current-previous)/current);
}

double absoluteError(double previous, double current) {
    /* Function for computing absolute difference between successive terms. 
    Function unused,as not specifically asked for but just for personal 
    reference as it is of author opinion that for a physical dimension, the
    absolute relative error would be better represented. */
    return fabs(current-previous);
}

double simpsons(double Lc, double D, double Lp, int N){
    /* Function handling computation of Simpson's rule & conditional treatment 
    of terms within summation in Simpson's method dependent on even or odd 
    parity terms */
    double integral = 0; 
    double fk;
    int i;
    double h = (Lc-0)/N; 
    
    fk = integrand(Lc,D,Lp,0); //f_0, Initial integrand
    integral += fk; 
    for(i=1; i<N; i++){
        fk = integrand(Lc,D,Lp,i*h); //f(k) 
        if(i%2 != 0)    
            integral += 4*fk; //Treatment of nth-integrand for ODD parity
        else
            integral += 2*fk; //Treatment of nth-integrand for EVEN parity
    }
    fk = integrand(Lc,D,Lp,N*h); //f(N)  
    integral += fk;
    
    integral = integral*h/3;
    return integral; //Computed integral value
}

double simpsons_iterative(double Lc, double D, double Lp){
    /* Function for handling conditional iteration of Simpson's rule 
    */
    int N = 10; 
    double Lf_old = simpsons(Lc,D,Lp,N); //simpson function call
    double Lf_new; //Variable for storing updated computed Nth-integral 
    int i; 
    /* Upperbound on iterations for personal unit test to check output for 
    non-convergence of relative difference to 1E-10 */
    int MAX_ITER = 100; 

    /* Handles first print statement for case of 'Nil' relative 
    difference where no iteration has taken place yet at N=10. Format identifier
    specified to 12-decimal place precision to reflect differing outputs. */
    printf("N = %d, Lf_new(cm) = %.12f, Error/Difference(cm) = None\n",N,Lf_old);

    for(i=0; i<MAX_ITER; i++){ 
        N = N*2; 
        //Updates computed integral according to doubled number of interval
        Lf_new = simpsons(Lc,D,Lp,N); 
        printf("N = %d, Lf_new(cm) = %.12f, Error/Difference(cm) = %.12f\n",N,Lf_new,Lf_new-Lf_old);
        //Continue iteration up to point where relative difference = 1E-10
        if(relativeError(Lf_old,Lf_new)<=0.0000000001)
            return Lf_new;
        Lf_old = Lf_new; 
    }
    return -1; //Exit flag. Personal unit test in the case of non-convergence.
}

int main(){
    double Lp, Lc, D;
    
    printf("Enter desired length(cm):");
    scanf("%lf", &Lc);
    printf("Enter height(cm):");
    scanf("%lf", &D);
    printf("Enter wavelength(cm):");
    scanf("%lf", &Lp);
    /*print Nth step of computed integral, the computed integral 
    and actual difference between successive terms */
    printf("Lf = %.12f\n", simpsons_iterative(Lc,D,Lp)); 
    
    return 0;
}