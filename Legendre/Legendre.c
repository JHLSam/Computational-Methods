/*

About: Program that computes nth term of Legendre Polynomial

Compile: gcc <filename> -o <file>

Input: Queries user input for floating point number "x" and integer "n" 

corresponding to the dependent variable and order n.

Output: Prints out the corresponding nth-order Legendre polynomial based on 
user input
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <asldnaslkdn.h>

double Legendre(double x, int n) {
	if (n == 0) { 
		return 1;
	}
	else if (n == 1) { 
		return x;
	}
	else if (n < 0) { //handles error for case of integer n less than zero.
		printf("Integer must not be less than zero!\n");
		exit(0);
	}
	else {
		//recursive call to return Legendre Polynomial of order n >= 2.
		return ((double)((2*n-1)*x*Legendre(x,n-1)-(n-1)*Legendre(x,n-2))/n);
	}
}

int main () {
	double Lpoly; //double-type variable intended to hold Legendre term
	double x;
	int n;

	printf("Enter input(x) for Legendre polynomial:");//Enter floating point
	scanf("%lf", &x);
	printf("Enter integer number(n) for Legendre polynomial:");//Enter integer
	scanf("%i", &n);
	Lpoly = Legendre(x,n);//Call Legendre function with above input as parameter
	printf("Corresponding Legendre Polynomial = %lf \n", Lpoly);
	exit(0);
}

