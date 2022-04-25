#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include"termpaper.h"


//Hubbard model with 
//half filled ,No double occupancy
//calculation of average value of spin correlation function <s^z_i.s^z_i+1>
//Number of lattice sites 2N
//MC: number of montecarlo samples

double spin_correlation(int N,int A_u[N],int A_d[N]){
	
	int arr[2*N];
	for(int i=0;i<N;i++){
		//up spins
		int j=A_u[i];
		//printf(" j=%d",j);
		arr[j]=1;
	}
	
	for(int i=0;i<N;i++){
		//down spins
		int j=A_d[i];
		arr[j]=-1;
		
	}
	printf("\n");
	for(int i=0;i<2*N;i++){
		printf("%d",arr[i]);
	}
	printf("\n");
	//Sz spin correlation:<s^z_i.s^z_i+1>
		
		
	double sum=0;
	sum+=arr[0]*arr[2*N-1];	//due to periodic boundary condition
	for(int i=0;i<2*N-1;i++){
		sum+=arr[i]*arr[i+1];
	}
	//printf("\nsum=%lf\n",sum);
	
	return sum; 
	
}


void hubbard_1d_montecarlo(int N_2, int MC){
	int N=N_2/2;
	 int MC_t=1;
	double spin_corr;
	double sum_spin_corr=0;
	double spin_co[MC];
	FILE *file1;
	file1=fopen("spin_corr_vs_MC_step.txt","w");
	
	//configuraion A-->B,u:up spins d:down spins
	int A_u[N] , B_u[N] , A_d[N] , B_d[N];
	
	//initial state
	for(int i=0;i<N;i++){
			//up spins
			A_u[i]=i;
	}
	printf("\n");
	for(int i=0;i<N;i++){
			//down spins
			A_d[i]=(i+N);
			
	}
	spin_corr=spin_correlation(N,A_u,A_d)/4;
	spin_co[0]=spin_corr/N_2;
	sum_spin_corr+=spin_corr;
	fprintf(file1,"%d	%lf\n",MC_t,sum_spin_corr/MC_t);
	
	while(MC_t<=MC){
	//MC step
		//randomly excahnge oppposite spins
		int change_u=rand()%N;
		int change_d=rand()%N;
		double metropolis;
		for(int i=0;i<N;i++){
			if(i==change_u){
				B_u[i]=A_d[change_d];
			}		
			else
				B_u[i]=A_u[i];
		}
		for(int i=0;i<N;i++){
			if(i==change_d){
				B_d[i]=A_u[change_u];
			}		
			else
				B_d[i]=A_d[i];
		}
	/*	
			printf("\n");
			for(int i=0;i<N;i++)
				printf("%d",A_u[i]);
			printf("\n");
			for(int i=0;i<N;i++)
				printf("%d",A_d[i]);
			printf("\n\n");
	*/	

		//metropolis process

		metropolis=mod_vandermonde_det_ratio(N, A_u, B_u, A_d, B_d, change_u,change_d);
		//printf("\nmetro=%.3lf\n",metropolis);
			//accept with probability 
			double p=(1.0*rand())/RAND_MAX;
			//printf("%lf\t%lf",p,metropolis);
			if(p<=metropolis){
				//accept the new configuration;calculate with new one
				spin_corr=spin_correlation(N,B_u,B_d)/4.0;
				spin_co[MC_t]=spin_corr/N_2;
				//test print
				printf("\nsum=%lf",spin_corr*4.0);
				sum_spin_corr+=spin_corr;
				fprintf(file1,"%d	%lf\n",MC_t,sum_spin_corr/MC_t);
				
				for(int i=0;i<N;i++){
					//update
					A_u[i]=B_u[i];
					A_d[i]=B_d[i];
				}

			}
			else{
				//reject the new configuration ;calculate with old one
				printf("\nsum=%lf",spin_corr*4.0);
				sum_spin_corr+=spin_corr;
				spin_co[MC_t]=spin_corr/N_2;
				fprintf(file1,"%d	%lf\n",MC_t,sum_spin_corr/(N_2*MC_t));
				 
			}
	
		MC_t++;
	}
	
	printf("\nspin_corr/site=%lf\n",sum_spin_corr/(N_2*MC));
	
	//std deviation
	double std_dev=0;
	for( int i=0;i<MC;i++){
		std_dev+=pow((spin_co[i]-(sum_spin_corr/(N_2*MC))),2);
	}
	std_dev=pow(std_dev/MC,0.5);
	printf("std Dev:=%lf",std_dev);
	
}

