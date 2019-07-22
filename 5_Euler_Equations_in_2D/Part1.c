#include <stdio.h>
#include <math.h>
#include <stdlib.h>
int main(){
int N;
printf("Enter the number of grid (N)=\n");
scanf("%d",&N);
FILE *fp= fopen("pressure","wb");
FILE *fd= fopen("density","wb");
FILE *fv= fopen("velocity","wb");
FILE *fe= fopen("energy","wb");
int i,j,p,q;
double U[3][N+1],velocity[N+1],P[N+1],f[3][N+1],C[N+1],Alpha_temp[N+1],Alpha_fwd[N+1],Alpha_bwd[N+1],maximum1,maximum2,f_fwd[3][N+1],f_bwd[3][N+1];
double lamda[N+1],maximum3,max_lamda,tau, T, del_x, U_num_update[3][N+1],U_ini[3][N+1];
del_x=1/(float)N;

for (i=0; i<N; i++){
    if (i < 0.5*N ){
        U[0][i]=1;
        U[2][i]=2.5;
    }
    else {
        U[0][i]=0.125;
        U[2][i]=0.25;
    }
}

for (i=0; i<3; i++)
for (j=0; j<N; j++)
        U_ini[i][j]=U[i][j];

for (i=0; i<N; i++)
    velocity[i]=0;

for (i=0; i<N; i++)
    U[1][i]= U[0][i]*velocity[i];

do {
for(i=0; i<N; i++){
        P[i] = 0.4 * ( U[2][i] - (0.5 * ( pow(U[1][i],2)/U[0][i]) ));
        f[0][i]= U[1][i];
        f[1][i]= ( pow(U[1][i],2)/U[0][i] )+ P[i];
        f[2][i]= ( U[1][i]/U[0][i]) * (U[2][i]+P[i]);
        C[i]=sqrt( 1.4*P[i]/U[0][i]);
        Alpha_temp[i] = fabs(U[1][i]/U[0][i])+C[i];
        }


maximum1=Alpha_temp[0];
for (i= 0; i <N; i++){
     q=i+1;
    if (i==N-1)
        q=N-1;

    if (Alpha_temp[q] > maximum1)
       maximum1  = Alpha_temp[i];

    Alpha_fwd[i]=maximum1;
}
maximum2=Alpha_temp[0];
for (i= 0; i <N; i++){
     p=i-1;
    if(p<0)
        p=0;

    if (Alpha_temp[p] > maximum2)
       maximum2  = Alpha_temp[i];

    Alpha_bwd[i]=maximum2;
  }





for (i=0; i<3; i++){
for (j=0; j<N; j++){
      p=j-1;
    if(p<0)
        p=0;
    q=j+1;
    if (j==N-1)
        q=N-1;

    f_fwd[i][j] = ( 0.5 * (f[i][q]+f[i][j])) - ( (Alpha_fwd[j]/2)*(U[i][q]-U[i][j]));
    f_bwd[i][j] = ( 0.5 * (f[i][j]+f[i][abs(p)])) - ( (Alpha_bwd[j]/2)*(U[i][j]-U[i][abs(p)]));
}
}


for(i=0; i<N; i++)
    lamda[i]=fabs(U[1][i]/U[0][i])+ C[i];

maximum3=lamda[0];
for (i= 0; i <N; i++){
       if (lamda[i] > maximum3)
       maximum3  = lamda[i];
}
max_lamda=maximum3;

tau =0.9* (del_x/max_lamda);
T=T+tau;



   for(i=0; i<3; i++){
        for (j=0; j<N; j++){
U_num_update[i][j] =U[i][j] - ( (tau/del_x) *( f_fwd[i][j]-f_bwd[i][j]));
        }
      }


for (i=0; i<3; i++){
        for (j=0; j<N; j++){
     U[i][j] = U_num_update[i][j];
}
}

}while(T < 0.2);





        printf("__________________________________________________________________________________________________\n");
        printf("Pressure values:\n");
        for (i=0; i<N; i++){
            printf("%f\n",P[i]);
             fprintf(fp, "%f\t, %f\n",(float)i/N,P[i]);
             fprintf(fd, "%f\t, %f\n",(float)i/N,U[0][i]);
             fprintf(fv, "%f\t, %f\n",(float)i/N,U[1][i]/U[0][i]);
             fprintf(fe, "%f\t, %f\n",(float)i/N,U[2][i]);
}

        printf("__________________________________________________________________________________________________\n");
        printf("Density values:\n");
        for (i=0; i<N; i++)
              printf("%f\t, %f\n",(float)i/N,U[0][i]);


        printf("__________________________________________________________________________________________________\n");
        printf("Velocity values:\n");
        for (i=0; i<N; i++)
              printf("%f\t, %f\n",(float)i/N,U[1][i]/U[0][i]);

        printf("__________________________________________________________________________________________________\n");
        printf("Energy values:\n");
        for (i=0; i<N; i++)
              printf("%f\t, %f\n",(float)i/N,U[2][i]);


         fclose(fp);
         fclose(fd);
         fclose(fv);
         fclose(fe);

}


