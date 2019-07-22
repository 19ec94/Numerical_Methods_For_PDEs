#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main()
{
    float c,a,b,s,tau,grid;
    float u[400],C[400],L[400][400],R[400][400], X[400], m[400];
    int N, t, n;
    printf("Enter the value of N:\n");
    scanf("%d", &N);
    n = 0.2*N;
    tau=1/(float)N;
    grid=1/(float)N;

    a=(tau/(grid * grid*2)) ;
    b=1+(tau/(grid*grid));
    c=1-(tau/(grid*grid));
    int j,i;


    for (i=1; i<=N-1; i++)
    {
        u[i]= sin(M_PI * i * grid);
       // printf("%f\n",u[i]);
    }

    //Left hand side matrix

    for (j=1; j<=N-1; j++){
        for(i=1; i<=N-1; i++){
            if(i==j)
            {
                L[i][j]=b;
            }
            else if (i==j+1)
            {
                L[i][j]=-a;
            }
            else if (i==j-1){
                L[i][j]=-a;
            }
            else{
                L[i][j]=0;
            }
            //printf("%f\t",L[i][j]);
        }
       // printf("\n");
    }
    //printf("\n");


    //Right hand side matrix
    for (j=1; j<=N-1; j++){
        for(i=1; i<=N-1; i++){
            if(i==j)
            {
                R[i][j]=c;
            }
            else if (i==j+1)
            {
                R[i][j]=a;
            }
            else if (i==j-1)
            {
                R[i][j]=a;
            }
            else
            {
                R[i][j]=0;
            }
           // printf("%f\t",R[i][j]);
        }
      //  printf("\n");
    }
   // printf("\n");

    for (i=2; i<N; i++) {
        m[i] = L[i][i-1]/L[i-1][i-1];
        L[i][i]= L[i][i]-(m[i]*L[i-1][i]);
        L[i][i-1]=0;
    }
    for (t=1; t<=n; t++) {          // Time Loop
    //matrix on rhs side
    float sum=0.0;
    int w;
    for (j=1; j<=N-1; j++)
    {
        for(i=1; i<=N-1; i++)
        {
            for(w=1; w<=N-1; w++)
            {
                sum=sum+R[j][w]*u[w];

            }

            C[j]=sum;
            sum=0.0;
        }

    }
   /* for (j=1; j<=N-1; j++){
        printf("%f\n",C[j]);
    }*/


   // printf("\n");

    for (i=2; i<N; i++) {
     //   m[i] = L[i][i-1]/L[i-1][i-1];
     //   L[i][i]= L[i][i]-(m[i]*L[i-1][i]);
        C[i] = C[i] - (m[i]*C[i-1]);
     //   L[i][i-1]=0;
    }


    X[N-1] = C[N-1]/L[N-1][N-1];
    u[N-1]= X[N-1];

    for (i=N-2; i>0; i--){
        X[i]= (C[i]-(L[i][i+1]*X[i+1]))/L[i][i];
        u[i]= X[i];
    }
}
    for (i=0; i<=N; i++){
        //printf("%f\n", X[i]);
    }

float o=0.2;
float U[400];
    for (i=0; i<=N; i++)
    {
        U[i]= sin(M_PI*i*grid) * exp(-pow(M_PI,2)*o);
      //  printf("%f\n",U[i]);
    }
    float E[400];
    float E_temp[400];
    printf("\n");
for (i=1; i<=N-1; i++)
{
    E_temp[i]=U[i]-u[i];
    E[i]=E_temp[i]*E_temp[i];
    E[i]=sqrt(E[i]);
    printf("%f\n",E[i]);
}
float location=0.0;
float maximum;
maximum=E[0];
for (i= 0; i <N; i++)
  {
    if (E[i] > maximum)
    {
       maximum  = E[i];
       location = i+1;
    }
  }
printf("Error value=%f",maximum);

printf("\n\n\n\n");

  return 0;
}
