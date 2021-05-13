/*
About: Program that approximates solution of differential equation with Euler's method

Compile: gcc -Wall -lm <filename> -o <chosen_exe_name>

Input: Queries user input for floating point numbers for the initial values of  
water level(Litres),acid level(Litres) and the time step size for the Euler Integrator
and the time step size for each iteration of the output 

Output: Prints out approximated water level "w", acid level "a" and the time "t" for each
iteration step size
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//differential equation for rate of change of water
double dw_by_dt(double w, double a) { 
    /*edge cases not specifically asked for but in my cs courses
    it's usually done by default, if not needed just drop a feedback*/
    if(w+a==0)  
        return 2.5; 
    if((w+a)<3)
        return 2.5 - w;
    return (2.5 - 3*w/(w+a));
}
//differential equation for rate of change of acid
double da_by_dt(double w, double a) {
    if(w+a==0)
        return 2.5; 
    if((w+a)<3)
        return 2.5 - a;
    return (2.5 - 3*a/(w+a));
}

int main () {
 double w,a,t,delta_t,delta_t_out,wnew,anew,next_output_time; 
    FILE *fp; //file pointer
    fp = fopen("output2.dat", "w");
 
 printf("Initial water(Litres):"); //user input for water,acid, time step size
 scanf("%lf", &w);
 printf("Initial acid(Litres):");
 scanf("%lf", &a);
 if ((w+a < 0) || (w+a > 200)){ //edge case handler
    printf("Number cannot be negative and combined volume cannot be greater than 200");
    exit(0);
 }
 
 printf("Time step size(minutes):");
 scanf("%lf",&delta_t);
 printf("Time step size for output:");
 scanf("%lf",&delta_t_out);
 
 t = 0;
 next_output_time = t + delta_t_out; //first delta_t_out step size
 
 //Iteration handler, iterate all the way till tank is full as inflow > outflow
 while((a+w)<200){ 
    wnew = w + delta_t*dw_by_dt(w,a); //forward difference for each iteration
    anew = a + delta_t*da_by_dt(w,a); 
    t += delta_t; //forward difference for time for each iteration
    
    if(next_output_time <= t){
        fprintf(fp, "%lf\t%lf\t%lf\n",t,wnew,anew);
        next_output_time += delta_t_out;
    }
    w = wnew;
    a = anew;
 }
 fclose(fp);
}
