#include<stdio.h>
#include<stdlib.h>
#include<math.h>

int main(){
int i, j,N;
printf("Enter the Number of Grid:\n");
scanf("%d",&N);
double alpha=10,X[N+2],a[N+2][N+2],f[N+2],phi[N+2],Approx[N+2],Exact[N+2],
Error[N+2];

for(i=0; i<N+1; i++)
    X[i]=( exp((alpha*i)/(float)N)-1 )/( exp(alpha)-1 );

for(i=1; i<N; i++){
 for(j=1; j<N; j++){
  if(i==j){
     a[i][j]=(1/(X[i]-X[i-1]) )+(1/(X[i+1]-X[i]) );
  }
  else if (i-j==1 ) {
     a[i][j]= -1/(X[i]-X[i-1]);
  }
 else if (j-i==1 ) {
     a[i][j]= -1/(X[j]-X[j-1]);
  }
  else{
      a[i][j]=0;
 }
}
}

for (i=0; i<N+1; i++)
f[i]= pow(M_PI,2) * sin(M_PI*X[i]);

for(i=1; i<N+1; i++)
phi[i]= ( (X[i]-X[i-1])*0.5*f[i] ) + ( (X[i+1]-X[i])*0.5*f[i]);


double m;
for(i=2; i<N; i++){
       m =a[i][i-1]/a[i-1][i-1];
       a[i][i]=a[i][i]-(m * a[i-1][i]);
      // a[i][i-1]=0;
      phi[i]=phi[i]-(m*phi[i-1]);

}


Approx[N-1]=phi[N-1]/a[N-1][N-1];

for(i=N-2; i>=1; i--){
     Approx[i]=( phi[i]-a[i][i+1]*Approx[i+1] ) /a[i][i];
}


printf("______________________________________________________________\n");


for(i=0; i<N;i++)
Exact[i]=sin(M_PI*X[i]);

double grid,total_err;

 grid=X[1]-X[0] ;

 for (i=1;i<N;i++)
 {  if (X[i+1]-X[i]>grid)
    grid=X[i+1]-X[i] ;
 }


// error calculation
for (i=1;i<=N;i++)
total_err=total_err+0.5*(X[i]-X[i-1])*(fabs(Approx[i]-Exact[i])*fabs(Approx[i]-Exact[i])+fabs(Approx[i-1]-Exact[i-1])*fabs(Approx[i-1]-Exact[i-1]));
 total_err=sqrt(total_err);
//

printf("Error %f at location %f ",total_err,grid);

return 0;
}
