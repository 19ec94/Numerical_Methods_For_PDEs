#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main()
{

    int N,i,j,scheme;
    FILE *f = fopen("Edil.csv", "wb");
    printf("Enter the value of grid space\n");
    scanf("%d",&N);

    double h,U[N+1],X[N+1],mu,s,T,tau,E[N+1],E_temp,Exact[N+1],Z[N+1];
    h=1.0/(double)N;
    mu=0.9;
    tau =mu * h;

    int M=(int)(N/mu);

    T =M*tau;

    printf("Select the scheme you want to evaluate\n");
    printf(" For Upwind type 1, for Lax-wendroff type 2\n");
    scanf("%d",&scheme);
    if (scheme==1){
        s=mu;
    }
    else if (scheme==2){
        s= pow(mu,2);
    }
    for (i=0; i<N; i++){
        if (i<N/2){
            U[i]=1;
        }
        else
        {
            U[i]=0;
        }
        Z[i]=U[i];
        printf("%f\n",Z[i]);
    }
    printf("\n\n\n");

int p,q;
   for (j=0; j<=M; j++){                                                  //Time loop
     for (i=0; i<N; i++){
         p = i-1;
         q = i+1;

         if (p<0){
           p = N -abs(p);
         }
         if (i==N){
            q = q-N;
            }                                                            // approximate solution

           X[i]=U[i]- ((0.5*mu) *(U[abs(q)]-U[abs(p)])) + ((0.5*s) *(U[abs(q)]-2*U[i]+U[abs(p)]));
     }

        for (i=0; i<N; i++){
        U[i]=X[i];
        //printf("%lf",U[i]);
        }

    }

 for (i=0; i<N; i++){
        printf("%f\n",X[i]);
    fprintf(f,"%d\t, %f\t, %f\n ",i, Z[i],X[i]);
 }

  printf("\n\n");
 for (i=0; i<N; i++){
    //printf("%d\n",i);
 }


fclose(f);
}
