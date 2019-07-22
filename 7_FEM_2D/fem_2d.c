//gcc Assignment8.c -o 	Assignment8 -I/usr/include  		-L/usr/lib -leng 	-lmx -lm
//		       Path for MATLAB-include  	Path for MATLAB/lib
#include<stdio.h>
#include<math.h>
#include<string.h>
#include<stdlib.h>
#include "engine.h"
#define PI 3.14159265358979323846

void main()
{
int i,j,k,N, Nint, Ncell,tvertices[201][3]={0.0},jg,kg,m=0;
double nodecoords[122][2]={0.0}, A[81][81]={0.0},C[81]={0.0}, actualsol[81] = {0.0},sum =0;

FILE *fp;
fp = fopen("MESH.dat", "r");
fscanf(fp,"%d %d", &N, &Nint);

for(i=0;i<N;i++)
{
fscanf(fp, "%lf %lf", &nodecoords[i][0],&nodecoords[i][1]);
}

fscanf(fp,"%d", &Ncell);
double* vol = (double*)malloc(Ncell*sizeof(double*));

for(i=0;i<Ncell;i++)
{
fscanf(fp,"%d %d %d", &tvertices[i][0], &tvertices[i][1], &tvertices[i][2]);
}
fclose(fp);
int p[2][3],irow,jcol;
p[0][0] = p[1][0] = -1;
p[0][1] = p[1][2] = 1;
p[0][2] = p[1][1] = 0;
double x1,x2,x3,y1,y2,y3,B1,B2,B3,B4,val1=0.0,val2 = 0.0,s,a,b,c;
for(i=0;i<Ncell;i++)
{

x1 = nodecoords[tvertices[i][0]-1][0];
x2 = nodecoords[tvertices[i][1]-1][0];
x3 = nodecoords[tvertices[i][2]-1][0];
y1 = nodecoords[tvertices[i][0]-1][1];
y2 = nodecoords[tvertices[i][1]-1][1];
y3 = nodecoords[tvertices[i][2]-1][1];
//printf("%lf\t%lf\t%lf\n",x1,x2,x3);
B1 = y3-y1;
B2 = y1-y2;
B3 = x1-x3;
B4 = x2-x1;


a = pow((pow((x1-x2),2)+pow((y1-y2),2)),0.5);
b = pow((pow((x2-x3),2)+pow((y2-y3),2)),0.5);
c = pow((pow((x1-x3),2)+pow((y1-y3),2)),0.5);

s = (a+b+c)/2.0;

vol[i] = pow(s*(s-a)*(s-b)*(s-c),0.5);

for(j=0;j<3;j++)
{

jg = tvertices[i][j]-1;
if(jg<Nint)
{
val1 =  2*pow(PI,2)*sin(PI*nodecoords[jg][0])*sin(PI*nodecoords[jg][1]) * vol[i]/3.0;
C[jg] += val1;

}
	for(k=j;k<3;k++)
	{
	kg = tvertices[i][k]-1;
	if((jg<Nint)&&(kg<Nint))
		{
		m++;
		irow = jg<kg?jg:kg;
		jcol = jg>kg?jg:kg;
		val2 = 0.25 *(( B1*p[0][j] + B2*p[1][j])*(B1*p[0][k] + B2*p[1][k]) + (B3*p[0][j] + B4*p[1][j]) * (B3*p[0][k] + B4*p[1][k]))/vol[i];

		A[irow][jcol] +=val2;
		A[jcol][irow] = A[irow][jcol];
		}
	}
}
actualsol[i] = sin(PI*nodecoords[i][0])*sin(PI*nodecoords[i][1]);

}

//MATLAB Solve
	int z = (Nint)*(Nint);
	Engine *ep;
	mxArray *A1 = NULL, *C1 =NULL, *Z = NULL,*res = NULL;
	double *ans = NULL;
	char buffer[BUFSIZ+1];


	if (!(ep = engOpen(""))) {
		fprintf(stderr, "\n Can't start MATLAB engine\n");
		exit;
	}
	A1 = mxCreateNumericMatrix(Nint, Nint, mxDOUBLE_CLASS, mxREAL);

	memcpy(mxGetData(A1),A,z*sizeof(double));
	engPutVariable(ep, "A1", A1);
	engEvalString(ep,"A1=A1';");

	C1 = mxCreateNumericMatrix(Nint, 1, mxDOUBLE_CLASS, mxREAL);
	memcpy(mxGetData(C1),C, (Nint)*sizeof(double));
	engPutVariable(ep, "C1", C1);
	Z = mxCreateNumericMatrix(Nint, 1, mxDOUBLE_CLASS, mxREAL);
	engEvalString(ep, "Z = mldivide(A1,C1);");
	res = engGetVariable(ep,"Z");
	mwSize nRow = mxGetM(res);
	ans = mxGetData(res);
	engClose;

for(i = 0; i < Ncell; i++)
{
	for(j = 0; j < 3; j++)
	{
	jg = tvertices[i][j]-1;
	if(jg<Nint)
		sum+= pow((actualsol[jg]-ans[jg]),2)*vol[i]/3.0;
	}
}
sum = pow(sum,0.5);

fp = fopen("Results.txt", "w");
fprintf(fp,"X\tY\tExact_solution\tNumerical_Solution\n");
for(i=0;i<Nint;i++)
{
fprintf(fp,"%f\t%f\t%f\t%f\n",nodecoords[i][0],nodecoords[i][1],actualsol[i],ans[i]);
}

printf("L2 Error = %f\n",sum);

}


