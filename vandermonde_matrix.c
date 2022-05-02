#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include"termpaper.h"


 double mod_vandermonde_det_ratio(int N,int A_u[N],int B_u[N] ,int A_d[N],int B_d[N],int change_u,int change_d){
	 //Number of sites=2N
	 //a=lattice spacing
	 //equal number of up and down spins
	 //Change of configuration from A_u(up spins),A_d(down spins) ---> B_u(up spins),B_d(down spins)
	 //By swapping of two opposite spins.
	 //cahnge_ud is the swapping index for swapping an up spin with down spin
	 long double temp1=1;
	 long double temp2=1;
	 
	 //long double temp3=1;
	 //long double temp4=1;
	 
	 //ratio of Mod of slatter determianant of up spins
	 for(int i=0;i<N;i++){
			if(i!=change_u){
			 temp1*=(pow(sin(PI*(1.0/(2*N-1))*(B_u[change_u]-B_u[i])),2))/(pow(sin(PI*(1.0/(2*N-1))*(A_u[change_u]-A_u[i])),2));
			 //printf("%d",A_u[change_u]-A_u[i]);
			}
		}
	 
	 //ratio of Mod of slatter determinant of down spins
	 for(int i=0;i<N;i++){
			if(i!=change_d)
			 temp2*=(pow(sin(PI*(1.0/(2*N-1))*(B_d[change_d]-B_d[i])),2))/(pow(sin(PI*(1.0/(2*N-1))*(A_d[change_d]-A_d[i])),2));
			 
		}
	//printf("\np=%Le\n",temp1*temp2);
	return temp1*temp2;  
	 
 }