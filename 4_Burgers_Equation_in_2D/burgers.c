#include <stdio.h>
#include <stdlib.h>
#include<math.h>

int main()
{
int i,N,p,q,choice;
FILE *f = fopen("Edil.csv", "wb");
printf("Enter the number of mesh point\n");
scanf("%d",&N);
double U_ini[N+1],F_U[N+1],Alpha_bar_fwd[N+1],Alpha_temp_fwd[N+1],Alpha_bar_bwd[N+1],Alpha_temp_bwd[N+1],Alpha_fwd[N+1],Alpha_bwd[N+1],F_dash[N+1],F_m_fwd[N+1],F_m_bwd[N+1];
double C_fwd[N+1],C_bwd[N+1],U_approx_temp[N+1],tau[N+1];
double tau_min,T;
double h = 1/(float)N;
int condition;
printf("Enter type of initial condition for continuous type 1 for discontinuous 2\n");
scanf("%d",&condition);

printf("Enter the choice  number for alpha 1 or 2\n");
scanf("%d",&choice);


if (condition ==1){
for (i=0; i<N; i++){
    U_ini[i]= 1 + 0.5 * sin(2*M_PI*i*h);
    }
}


if (condition==2){
for (i=0; i<N; i++){
    if(i<=N*0.1)
    U_ini[i]=0;
    else if (i>=N*0.3)
    U_ini[i]=0;
    else
    U_ini[i]=1;
    }
}

while(T<=1){

for(i=0; i<N; i++){
    F_U[i]=0.5*(U_ini[i]*U_ini[i]);
    F_dash[i]=U_ini[i];
}

for (i=0; i<N; i++){
      q= i+1;
if (i==N-1){
    q = q-N;
    }

if(U_ini[q] != U_ini[i])
    Alpha_temp_fwd[i]= ( ( F_U[q]-F_U[i] ) / (U_ini[q]-U_ini[i]) );
else
    Alpha_temp_fwd[i]= F_dash[i];

Alpha_bar_fwd[i]= Alpha_temp_fwd[i];
}
for(i=0; i<N; i++){
if(choice==1)
    Alpha_fwd[i]=sqrt(Alpha_bar_fwd[i]*Alpha_bar_fwd[i]);
if (choice==2)
    Alpha_fwd[i]=2*fabs(Alpha_bar_fwd[i]);
}

for (i=0; i<N; i++){

       p=i-1;
        if (p<0){
           p = N -abs(p);
         }
if(U_ini[i] != U_ini[p])
    Alpha_temp_bwd[i]= ( ( F_U[i]-F_U[p] ) / (U_ini[i]-U_ini[p]) );
else
    Alpha_temp_bwd[i]= F_dash[i];

Alpha_bar_bwd[i]= Alpha_temp_bwd[i];
}
for(i=0; i<N; i++){
if(choice==1)
    Alpha_bwd[i]=sqrt(Alpha_bar_bwd[i]*Alpha_bar_bwd[i]);
if(choice==2)
    Alpha_bwd[i]=2*fabs(Alpha_bar_bwd[i]);

}

for (i=0; i<N; i++){
    C_fwd[i]= 0.5 * ( Alpha_fwd[i]-Alpha_bar_fwd[i] );
    C_bwd[i]=0.5 * (Alpha_bar_bwd[i]+ Alpha_bwd[i]);
}

for (i=0; i<N; i++){
        q= i+1;
if (i==N-1){
    q = q-N;
    }
        p=i-1;
if (p<0){
    p = N -abs(p);
    }
if(U_ini[q] != U_ini[i])
F_m_fwd[i]= F_U[i] - ( 0.5 * (Alpha_fwd[i] - Alpha_bar_fwd[i]) * (U_ini[q]-U_ini[i])  );
else
F_m_fwd[i]=F_U[i]+( 0.5 * (F_U[q]-F_U[i]) );

if(U_ini[i] != U_ini[p])
F_m_bwd[i]= F_U[i] - ( 0.5 * (Alpha_bwd[i] + Alpha_bar_bwd[i]) * (U_ini[i]-U_ini[p])  );
else
F_m_bwd[i]=F_U[i]-( 0.5 * (F_U[i]-F_U[p]) );

}

for (i=0; i<N; i++){
tau[i] = (h/ (C_fwd[i]+C_bwd[i]));
}

tau_min=tau[0];
for(i=0; i<N; i++){
    if(tau[i]<tau_min){
        tau_min= tau[i];
    }
}
  T=T+tau_min;

//Approximate sollution
for (i=0; i<N; i++){
    U_approx_temp[i]= U_ini[i] + (tau_min /h )*(-F_m_fwd[i]+F_m_bwd[i]);
    U_ini[i]=U_approx_temp[i];
    }


}

for(i=0; i<N; i++){
printf("%f\n",U_approx_temp[i]);
fprintf(f,"%f\t, %f\n ",(float)i/N, U_approx_temp[i]);
}

fclose(f);
}


