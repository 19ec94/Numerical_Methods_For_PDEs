#include <stdio.h>
#include <stdlib.h>
#include <math.h>


int main()
{   int i,j,k,p,q,N=10;
    double U[3][N+1],f[3][N+1],pressure[N+1],velocity[N+1],enthalphy[N+1],h1[N+1],v1[N+1],c1[N+1],R1[3][3],d1[3][3],RI1[3][3],a[N+1],gama1[N+1];
    double sum_mat1,Temp_mat1[3][3],sum_mat2,A_roe1[3][3],f_m_fwd[3][N+1],f_m_bwd[3][N+1];
    double v2[N+1],h2[N+1],c2[N+1],R2[3][3],d2[3][3],RI2[3][3],a2[N+1],gama2[N+1],sum_mat3,Temp_mat2[3][3],sum_mat4,A_roe2[3][3];
    double lamda[N+1],maximum3,tau,T=0,U_num_update[3][N+1],c[N+1],max_lamda,del_x=1/(float)N;

    for (i=0; i<N; i++){                                  //Initialization
        if(i < 0.5*N){
            U[0][i]=1;
            U[2][i]=2.5;
        }
        else{
            U[0][i]=0.125;
            U[2][i]=0.25;
        }
    }

    for(i=0; i<N; i++)
        velocity[i]=0;

    for(i=0; i<N; i++)
        U[1][i]= U[0][i] * velocity[i];

for(i=0; i<3; i++){
for(j=0; j<N; j++){
printf("U[%d][%d]=%f\n",i,j,U[i][j]);
}printf("\n");
}

printf("________________________________________________________________________________________\n");
while (T< 0.2){                                                                                        //Time loop
    for(i=0; i<N; i++){
        pressure[i]= 0.4 * (U[2][i] - ( 0.5 * (pow(U[1][i],2)/U[0][i]) ) );
        f[0][1]=U[1][i];
        f[1][i]=( pow(U[1][i],2) / U[0][i])+ pressure[i];
        f[2][i]=(U[1][i]/U[0][i])*(U[2][i]+pressure[i]);
        c[i]=sqrt(1.4 * pressure[i]/U[0][i]);
    }
for(i=0; i<3; i++){
for(j=0; j<N; j++){
printf("f[%d][%d]=%f\n",i,j,f[i][j]);
}printf("\n");
}

 for(i=0; i<N; i++)
        enthalphy[i] = ( U[2][i]+pressure[i])/U[0][i];


    for(i=0; i<N; i++){
            q=i+1;
            if(q==N)
            q=N-1;
            v1[i]= (sqrt(U[0][i])*(U[1][i]/U[0][i]))+(sqrt(U[0][q])*(U[1][q]/U[0][q]))/(sqrt(U[0][i])+ sqrt(U[0][q]));
            h1[i]= (sqrt(U[0][i])* enthalphy[i])+(sqrt(U[0][q])* enthalphy[q])/(sqrt(U[0][i])+ sqrt(U[0][q]));
            c1[i]=sqrt(( 0.4 * (h1[i] + (0.5* pow(v1[i],2)) )));

            R1[0][0]=1;
            R1[0][1]=1;
            R1[0][2]=1;
            R1[1][0]= v1[i]-c1[i];
            R1[1][1]=v1[i];
            R1[1][2]=v1[i]+c1[i];
            R1[2][0]=h1[i]-(v1[i]*c1[i]);
            R1[2][1]=0.5*(pow(v1[i],2));
            R1[2][2]=h1[i]+(v1[i]*c1[i]);

            d1[0][0]=fabs(v1[i]-c1[i]);
            d1[0][1]=0;
            d1[0][2]=0;
            d1[1][0]=0;
            d1[1][1]=fabs(v1[i]);
            d1[1][2]=0;
            d1[2][0]=0;
            d1[2][1]=0;
            d1[2][2]=fabs(v1[i]+c1[i]);

            a[i]=c1[i];
    for(i=0; i<3; i++){
        if(d1[i][i] < a[i]){
                d1[i][i]= 0.5 * (a[i] + (pow(d1[i][i],2)/a[i]));
        }
    }

            gama1[i]=0.4/pow(c1[i],2);
        RI1[0][0]= 0.5 * ( 1 + (gama1[i]*(pow(v1[i],2)-h1[i])) + (v1[i]/c1[i]) );
        RI1[0][1]= - 0.5 * ( (gama1[i]* v1[i])+(1/c1[i]) );
        RI1[0][2]= 0.5 * gama1[i];
        RI1[1][0]= -gama1[i]*(pow(v1[i],2)-pow(h1[i],2));
        RI1[1][1]= gama1[i]*v1[i];
        RI1[1][2]=-gama1[i];
        RI1[2][0]= 0.5 * ( 1 + (gama1[i]*(pow(v1[i],2)-h1[i])) - (v1[i]/c1[i]) );
        RI1[2][1]= - 0.5 * ( (gama1[i]* v1[i])-(1/c1[i]) );
        RI1[2][2]= 0.5 * gama1[i];

   for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) {
        for (k = 0; k < 3; k++) {
          sum_mat1 = sum_mat1 + d1[i][k]*RI1[k][j];
        }

        Temp_mat1[i][j] = sum_mat1;
        sum_mat1 = 0;
      }
    }

    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) {
        for (k = 0; k < 3; k++) {
          sum_mat2 = sum_mat2 + R1[i][k]*Temp_mat1[k][j];
        }

        A_roe1[i][j] = sum_mat2;                                                                // A_Roe (i+0.5)
        sum_mat2 = 0;
      }
    }
    }


for(i=0; i<3; i++){
for(j=0; j<N; j++){
printf("R[%d][%d]=%f\t",i,j,R1[i][j]);
}printf("\n");
}
printf("________________________________________________________________________________________\n");
for(i=0; i<3; i++){
for(j=0; j<N; j++){
printf("RI1[%d][%d]=%f\t",i,j,RI1[i][j]);
}printf("\n");
}
printf("________________________________________________________________________________________\n");
for(i=0; i<3; i++){
for(j=0; j<N; j++){
printf("d1[%d][%d]=%f\t",i,j,d1[i][j]);
}printf("\n");
}

 printf("________________________________________________________________________________________\n");

for(i=0; i<N; i++){
             p=i-1;
              if (p<0)
                p=0;
            v2[i]= (sqrt(U[0][abs(p)])*(U[1][abs(p)]/U[0][abs(p)]))+(sqrt(U[0][i])*(U[1][i]/U[0][i]))/(sqrt(U[0][abs(p)])+ sqrt(U[0][i]));
            h2[i]= (sqrt(U[0][abs(p)])* enthalphy[abs(p)])+(sqrt(U[0][i])* enthalphy[i])/(sqrt(U[0][abs(p)])+ sqrt(U[0][i]));
            c2[i]=sqrt(( 0.4 * (h2[i] + (0.5* pow(v2[i],2)) )));

            R2[0][0]=1;
            R2[0][1]=1;
            R2[0][2]=1;
            R2[1][0]= v2[i]-c2[i];
            R2[1][1]=v2[i];
            R2[1][2]=v2[i]+c2[i];
            R2[2][0]=h2[i]-(v2[i]*c2[i]);
            R2[2][1]=0.5*(pow(v2[i],2));
            R2[2][2]=h2[i]+(v2[i]*c2[i]);

            d2[0][0]=fabs(v2[i]-c2[i]);
            d2[0][1]=0;
            d2[0][2]=0;
            d2[1][0]=0;
            d2[1][1]=fabs(v2[i]);
            d2[1][2]=0;
            d2[2][0]=0;
            d2[2][1]=0;
            d2[2][2]=fabs(v2[i]+c2[i]);

            a2[i]=c2[i];
    for(i=0; i<3; i++){
        if(d1[i][i] < a2[i]){
                d2[i][i]= 0.5 * (a2[i] + (pow(d2[i][i],2)/a2[i]));
        }
    }


            gama2[i]=0.4/pow(c2[i],2);
        RI2[0][0]= 0.5 * ( 1 + (gama2[i]*(pow(v2[i],2)-h2[i])) + (v2[i]/c2[i]) );
        RI2[0][1]= - 0.5 * ( (gama2[i]* v2[i])+(1/c2[i]) );
        RI2[0][2]= 0.5 * gama2[i];
        RI2[1][0]= -gama2[i]*(pow(v2[i],2)-pow(h2[i],2));
        RI2[1][1]= gama2[i]*v2[i];
        RI2[1][2]=-gama2[i];
        RI2[2][0]= 0.5 * ( 1 + (gama2[i]*(pow(v2[i],2)-h2[i])) - (v2[i]/c2[i]) );
        RI2[2][1]= - 0.5 * ( (gama2[i]* v2[i])-(1/c2[i]) );
        RI2[2][2]= 0.5 * gama2[i];

        for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) {
        for (k = 0; k < 3; k++) {
          sum_mat3 = sum_mat3 + d2[i][k]*RI2[k][j];
        }

        Temp_mat2[i][j] = sum_mat3;
        sum_mat3 = 0;
      }
    }

    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) {
        for (k = 0; k < 3; k++) {
          sum_mat4 = sum_mat4 + R2[i][k]*Temp_mat2[k][j];
        }

        A_roe2[i][j] = sum_mat4;                                                                               //A_Roe (i-0.5)
        sum_mat4 = 0;
      }
    }
    }
for (i=0; i<3; i++){
for (j=0; j<N; j++){
      p=j-1;
    if(p<0)
        p=0;
    q=j+1;
    if (q==N)
        q=N-1;

    f_m_fwd[i][j] = ( 0.5 * (f[i][q]+f[i][j])) - ( (A_roe1[i][j]/2) * (U[i][q]-U[i][j]));                     // Flux averages
    f_m_bwd[i][j] = ( 0.5 * (f[i][j]+f[i][abs(p)])) - ( (A_roe2[i][j]/2)*(U[i][j]-U[i][abs(p)]));
}
}

for(i=0; i<N; i++)
    lamda[i]=fabs(U[1][i]/U[0][i])+c[i];

maximum3=lamda[0];
for (i= 0; i <N; i++){
       if (lamda[i] > maximum3)
       maximum3  = lamda[i];
}
max_lamda=maximum3;

tau = 0.9*(del_x/max_lamda);
T=T+tau;

for(i=0; i<3; i++){
        for (j=0; j<N; j++){
U_num_update[i][j] =U[i][j]- ( (tau/del_x) *( f_m_fwd[i][j]+f_m_bwd[i][j]));
        }
}
}

for (i=0; i<3; i++){
    for(j=0; j<N; j++){
      printf("%f\n",U_num_update[i][j]);
    }
    printf("\n");
}



}

