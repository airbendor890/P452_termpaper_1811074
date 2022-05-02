#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include"termpaper.h"
void  main(){

	time_t t;
	srand((unsigned) time(&t));
	 int N_2=80;
	int MC=1000000;
	hubbard_1d_montecarlo(N_2,MC);
}