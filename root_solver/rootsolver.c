/*
About: root solver

Compile: gcc -Wall -lm <filename> -o <chosen_exe_name>
 Sample test values:    ./x 18 20
                        ./x 20 10
                        ./x 100 10

Input: cmd args


Output: Prints out to terminal, Iterator method used: (i) Iteration number (ii) start interval point (iii)end interval point
(iv) polynomial estimate at midpoint of each interval
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//polynomial as solved for in detail in theory
double f(double alpha, double R, double L) {
    return (2*pow(alpha,2)*acosh(R/alpha)-alpha*L)  ;
}

//bisection method

void bisection(double *a, double *b, double *R, double *L, int *iteration, double *alpha, double error){

    if (f(*a,*R,*L) * f(*b,*R,*L) >= 0){ 
        printf("Exception occurred:IntervalOutOfBounds");//Java-ish style exception handler/catcher
    } 
    double c = *a; 
    
    while ((*b-*a) >= error) //defined accuracy of interval
    { 
        c = (*a+*b)/2; //midpoint
        printf("Bisection: Iteration %d, alpha_1 = %.10f, alpha_2 = %.10f, f(mid_est) = %.10f\n", 
            *iteration, *a, *b, f(c,*R,*L));
    
        if (f(c,*R,*L) == 0.0) //for when root is found
            break; 
        else if (f(*a,*R,*L)*f(c,*R,*L) < 0) 
            *b = c; 
        else
            *a = c; 
        *iteration = *iteration + 1; //itercount
    }
    *alpha = (*a+*b)/2; //refreshed final root for bisection
}

//secant method
void secant(double *a, double *b, double *R, double *L, int *iteration, double *alpha, double error){

    double fa = f(*a,*R,*L); //point 1
    double fb = f(*b,*R,*L); //point 2

    if(fa == fb){
        printf("Exception occurred:IntervalOutOfBounds");
    }
    double c;
    
    while(fabs(*b - *a) >= error) //defined accuracy of interval
    {
        fa = f(*a,*R,*L);
        fb = f(*b,*R,*L);
        c = (*a*fb - *b*fa)/(fb-fa);
        printf("Seacant: Iteration %d, alpha_1 = %lf, alpha_2 = %lf, f(mid_est = %lf\n", 
            *iteration, *a, *b, f(c,*R,*L));
    
        *a = *b; //refresh points
        *b = c;
        *iteration = *iteration + 1;
    }
    *alpha = c; //refreshed final root for secant
}

/*call both methods(or rather i should say, subroutines), in order with respective conditional error passed, 
bounding respective 'while' iterators*/
void enhanced_method(double *a, double *b, double *R, double *L, int *iteration, double *alpha){
    bisection(a,b,R,L,iteration,alpha,0.5);
    secant(a,b,R,L,iteration,alpha,1e-8);
}

int main(int argc, char *argv[]) //allow for "runtime" parameter pass 
{
    if( argc != 3 ) {
      printf("ERROR: Pass 2 arguments R and L\n");
      printf("Adhere to input like so e.g >> ./as06 18 20");
      exit(0);
    }
    double R = atof(argv[1]); //convert type char to float
    double L = atof(argv[2]);

    double beta = L/2; //As per mentioned theory
    double alpha_1 = 0.05; //recommended initial interval bounds
    double alpha_2 = 0.5*R;
    double alpha;
    int iteration = 1;
    
    enhanced_method(&alpha_1,&alpha_2,&R,&L,&iteration,&alpha); 
    
    printf("alpha = %lf cm, beta = %lf cm\n\n", alpha, beta);
    
    return 0;
}
