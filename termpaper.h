#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

#define PI 3.14159

 double mod_vandermonde_det_ratio(int N,int A_u[N],int B_u[N] ,int A_d[N],int B_d[N],int change_u,int change_d);
 void hubbard_1d_montecarlo(int N_2,int MC);