#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int factorial(int j)
{
    if (j == 0)
        return 1;

    else if(j==1)
        return 1;
    else
        return j * factorial(j - 1);

}
int main()


{
    int p,q,m;
    int i,j;

    printf("Enter order of derivative m\n");
    scanf("%d",&m);
    printf("Enter the value of lower stencil  p\n");
    scanf("%d",&p);
    printf("Enter the value of higher stencil q\n");
    scanf("%d",&q);

    double Matrix [100][100];

    for (i=0; i<=(p+q); i++){
        for (j=0; j<=(p+q+4);j++){
            Matrix [i][j]=(pow((-p+i),j))/factorial(j);
        }
    }


    double B[100][100];

    for (i =0; i<1; i++) {
        for (j=0; j<=(p+q); j++){
            if (j==m)
            {
                B[i][j]=1;
            }
            else
            {
                B[i][j]=0;
            }

        }

    }



    double C[100][100];


    printf("\n");

    for (i=0; i<=(p+q+1); i++){
        for (j=0; j<=(p+q);j++){

            if (i == p+q+1){
                C[i][j] = B[0][j];
            }
            else
            {
                C[i][j] = Matrix[i][j];
            }
        }
    }

printf("The augumented matrix\n");
    double A[100][100];

    for (i=0; i<=(p+q); i++){
        for (j=0; j<=(p+q+1);j++){
            A[i][j]=C[j][i];
            printf("%f\t",A[i][j]);
        }
        printf("\n");
    }
    printf("\n\n");

    float c;
    int k;
    for(j=0; j<=p+q+2; j++){
        for(i=0; i<=p+q+1; i++){
            if(i>j)
            {
                c=A[i][j]/A[j][j];
                for(k=0; k<=p+q+1; k++)
                {
                    A[i][k]=A[i][k]-c*A[j][k];
                }
            }
        }
    }


    double sum;
    double z[100];
    z [p+q]=A[p+q][p+q+1]/A[p+q][p+q];
    for(i=p+q-1; i>=0; i--){
        sum=0.0;
        for(j=i+1; j<=p+q+2; j++){
            sum=sum+A[i][j]*z[j];
        }
        z[i]=(A[i][p+q+1]-sum)/A[i][i];
    }
    printf("\nThe Coefficients are: \n");
    for(i=0; i<=p+q; i++){
        printf("\nz%d=%f\t",i,z[i]);
    }

    printf("\n\n");

printf("The Taylors table\n");

    double F[100][100];

    for (i=0; i<=(p+q); i++){
        for (j=0; j<=(p+q+4);j++){

            F[i][j]= z[i]*Matrix[i][j];
             printf("%f\t",F[i][j]);
        }
   printf("\n");
    }
    printf("\n\n\n");
   printf("The  summation of Error values for each terms in Taylors Table is\n");
    float slack=0;
    float E[100];
    for (j=0; j<=(p+q+4); j++){
        for (i=0; i<=(p+q);i++){

            slack=slack+ F[i][j];
        }
        E[j]=slack;
        printf("%d  Order term= %f\n",j,E[j]);
        slack=0;
    }


    printf("\n");

        if (p==2){
            if (q==2){
    for (j=p+q+2; j<p+q+3; j++){
        if (E[j] <= -1.1*pow(10,-17)){
        printf("Accuracy of approximation is %d\n",j-m);
         }
         }
          }
        }


printf("\n\n\n");
}



