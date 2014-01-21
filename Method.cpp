	#include "Method.h"
#include "Thetraeder.h"
#include "Node.h"
#include "Mesh.h"
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <cstring>
#include <stdio.h>
#define R_I 1000000
#define R_D 1000000.0

double Method::axis[9]={0};
Method::Method()
{
}

Method::~Method()
{
}

int Method::init()
{
	srand(time(0));
	coeff[0] = coeff[1] = coeff[2] = 0.0;
	coeff[1] = 0.03;
	return 0;
}

double Method::scalar(double* _v1, double* _v2)
{
	return _v1[0]*_v2[0]+_v1[1]*_v2[1]+_v1[2]*_v2[2];
}

void Method::randomizeAxis(double* axes)
{
	for (int i=0; i<9; i++)
		axes[i]=0;

	axes[0] = (0.98 * (rand()%R_I)/R_D + 0.01) * (2.0*(rand()%2) - 1.0);
	axes[1] = (0.98 * (rand()%R_I)/R_D + 0.01) * (2.0*(rand()%2) - 1.0);
	axes[2] = (0.98 * (rand()%R_I)/R_D + 0.01) * (2.0*(rand()%2) - 1.0);
	double norma = sqrt(scalar(axes,axes));
	axes[0]/=norma; axes[1]/=norma; axes[2]/=norma;

	axes[3] = (0.98 * (rand()%R_I)/R_D + 0.01) * (2.0*(rand()%2) - 1.0);	
	axes[5] = (0.98 * (rand()%R_I)/R_D + 0.01) * (2.0*(rand()%2) - 1.0);
	axes[4] = -(axes[0]*axes[3]+axes[2]*axes[5])/axes[1];
	norma = sqrt(scalar(axes+3,axes+3));
	axes[3]/=norma; axes[4]/=norma; axes[5]/=norma;

	double znam = axes[0]*axes[4] - axes[3]*axes[1];
	if ( znam == 0.0 )
	{
		printf ("Warning! Recursion in axis random called!\n");
		randomizeAxis(axes);
		return;
	}
	
	double alfa = (axes[5]*axes[1] - axes[2]*axes[4])/znam,
		beta = (axes[3]*axes[2] - axes[0]*axes[5])/znam;
	axes[8] = 1.0/sqrt(alfa*alfa + beta*beta+ 1.0);
	axes[6] = alfa*axes[8]; axes[7] = beta*axes[8];
			
	//check!!!!!
	double check[6] = {scalar(axes,axes),scalar(axes+3,axes+3),scalar(axes+6,axes+6)
			,scalar(axes,axes+3),scalar(axes,axes+6),scalar(axes+3,axes+6)};
	if ((check[0] - 1.0)*(check[0] - 1.0) > 0.0001 || (check[1] - 1.0)*(check[1] - 1.0) > 0.0001 || (check[2] - 1.0)*(check[2] - 1.0) > 0.0001 || check[3]*check[3] > 0.0001 || check[4]*check[4] > 0.0001 || check[5]*check[5] > 0.0001)
	{
		printf ("Warning! Recursion in axis random called! %lf %lf %lf %lf %lf %lf   %lf %lf %lf\n",check[0],check[1],check[2],check[3],check[4],check[5], axes[0], axes[1], axes[2]);
		randomizeAxis(axes);
		return;
	}
}

void Method::randomizeAxis()
{
	for (int i=0; i<9; i++)
		axis[i]=0;
	//axis[0] = axis[4] = axis[8] = 1;
	//return;	//no randomization

	axis[0] = (0.98 * (rand()%R_I)/R_D + 0.01) * (2.0*(rand()%2) - 1.0);
	axis[1] = (0.98 * (rand()%R_I)/R_D + 0.01) * (2.0*(rand()%2) - 1.0);
	axis[2] = (0.98 * (rand()%R_I)/R_D + 0.01) * (2.0*(rand()%2) - 1.0);
	double norma = sqrt(scalar(axis,axis));
	axis[0]/=norma; axis[1]/=norma; axis[2]/=norma;

	axis[3] = (0.98 * (rand()%R_I)/R_D + 0.01) * (2.0*(rand()%2) - 1.0);	
	axis[5] = (0.98 * (rand()%R_I)/R_D + 0.01) * (2.0*(rand()%2) - 1.0);
	axis[4] = -(axis[0]*axis[3]+axis[2]*axis[5])/axis[1];
	norma = sqrt(scalar(axis+3,axis+3));
	axis[3]/=norma; axis[4]/=norma; axis[5]/=norma;

	double znam = axis[0]*axis[4] - axis[3]*axis[1];
	if ( znam == 0.0 )
	{
		printf ("Warning! Recursion in axis random called!\n");
		randomizeAxis();
		return;
	}
	
	double alfa = (axis[5]*axis[1] - axis[2]*axis[4])/znam,
		beta = (axis[3]*axis[2] - axis[0]*axis[5])/znam;
	axis[8] = 1.0/sqrt(alfa*alfa + beta*beta+ 1.0);
	axis[6] = alfa*axis[8]; axis[7] = beta*axis[8];
			
	//check!!!!!
	double check[6] = {scalar(axis,axis),scalar(axis+3,axis+3),scalar(axis+6,axis+6)
			,scalar(axis,axis+3),scalar(axis,axis+6),scalar(axis+3,axis+6)};
	if ((check[0] - 1.0)*(check[0] - 1.0) > 0.0001 || (check[1] - 1.0)*(check[1] - 1.0) > 0.0001 || (check[2] - 1.0)*(check[2] - 1.0) > 0.0001 || check[3]*check[3] > 0.0001 || check[4]*check[4] > 0.0001 || check[5]*check[5] > 0.0001)
	{
		printf ("Warning! Recursion in axis random called! %lf %lf %lf %lf %lf %lf   %lf %lf %lf\n",check[0],check[1],check[2],check[3],check[4],check[5], axis[0], axis[1], axis[2]);
		randomizeAxis();
		return;
	}

}

void Method::calculateCoeff(double* _c)
{
	_c[0] = coeff[0]*axis[0] + coeff[1]*axis[1] + coeff[2]*axis[2];
	_c[1] = coeff[0]*axis[3] + coeff[1]*axis[4] + coeff[2]*axis[5];
	_c[2] = coeff[0]*axis[6] + coeff[1]*axis[7] + coeff[2]*axis[8];
}

double Method::interpolate_1_order(Thetraeder* t, double* _crd, int val, Mesh* mesh)
{	
	Node* nodes[4] = {mesh->getNode(t->vert[0]), mesh->getNode(t->vert[1]), mesh->getNode(t->vert[2]), mesh->getNode(t->vert[3])};
	int min[3]={0},max[3]={0},diff=0;
	for (int m=0; m<3; m++)
	{
		for (int i=1; i<4; i++)
			if (nodes[i]->coords[m] < nodes[min[m]]->coords[m]) min[m] = i;
	}
	if (min[1]) {Node* tmp = nodes[0]; nodes[0] = nodes[min[1]]; nodes[min[1]] = tmp;};
	for (int m=0; m<3; m++)
	{
		for (int i=1; i<4; i++)
			if (nodes[i]->coords[m] > nodes[max[m]]->coords[m]) max[m] = i;
	}
	if (max[1] != 3) {Node* tmp = nodes[3]; nodes[3] = nodes[max[1]]; nodes[max[1]] = tmp;};

	double 	x0=nodes[0]->coords[0], x1=nodes[1]->coords[0], x2=nodes[2]->coords[0], x3=nodes[3]->coords[0], //
       		y0=nodes[0]->coords[1], y1=nodes[1]->coords[1], y2=nodes[2]->coords[1], y3=nodes[3]->coords[1],   //
       		z0=nodes[0]->coords[2], z1=nodes[1]->coords[2], z2=nodes[2]->coords[2], z3=nodes[3]->coords[2],   //
       		a,b,c,d;
    {                                                                
        double f0=nodes[0]->u[val],f1=nodes[1]->u[val],f2=nodes[2]->u[val],f3=nodes[3]->u[val]; //solve for 4x4
	double znam = x2*y1*z0 - x3*y1*z0 - x1*y2*z0 + x3*y2*z0 + x1*y3*z0 - x2*y3*z0 - 
                                    x2*y0*z1 + x3*y0*z1 + x0*y2*z1 - x3*y2*z1 - x0*y3*z1 + x2*y3*z1 + x1*y0*z2 - 
                                    x3*y0*z2 - x0*y1*z2 + x3*y1*z2 + x0*y3*z2 - x1*y3*z2 - 
                                    x1*y0*z3 + x2*y0*z3 + x0*y1*z3 - x2*y1*z3 - x0*y2*z3 + x1*y2*z3;
	
    	if (fabs(znam) < 0.000001) 
	{
		printf("%10lf<- skipping\n%10lf %10lf %10lf %10lf \n%10lf %10lf %10lf %10lf\n%10lf %10lf %10lf %10lf\n",znam, x0*2.0,x1*2.0,x2*2.0,x3*2.0,y0,y1,y2,y3,z0,z1,z2,z3); return 0.0;		
	}
        a = 	((-f1)*x3*y2*z0 + f1*x2*y3*z0 + f0*x3*y2*z1 - f0*x2*y3*z1 + 
    		f1*x3*y0*z2 - f0*x3*y1*z2 - f1*x0*y3*z2 + f0*x1*y3*z2 + 
    		f3*((-x2)*y1*z0 + x1*y2*z0 + x2*y0*z1 - x0*y2*z1 - x1*y0*z2 + 
      		x0*y1*z2) - f1*x2*y0*z3 + f0*x2*y1*z3 + f1*x0*y2*z3 - f0*x1*y2*z3 + 
    		f2*(x3*y1*z0 - x1*y3*z0 - x3*y0*z1 + x0*y3*z1 + x1*y0*z3 - x0*y1*z3))/
   		(x1*y2*z0 - x1*y3*z0 - x0*y2*z1 + x0*y3*z1 - x1*y0*z2 + x0*y1*z2 - 
    		x0*y3*z2 + x1*y3*z2 + x3*(y1*z0 - y2*z0 - y0*z1 + y2*z1 + y0*z2 - 
      		y1*z2) + x1*y0*z3 - x0*y1*z3 + x0*y2*z3 - x1*y2*z3 + 
    		x2*(y3*z0 + y0*z1 - y3*z1 - y0*z3 + y1*(z3 - z0)));
        b = 	(f1*y2*z0 - f1*y3*z0 - f0*y2*z1 + f0*y3*z1 - f1*y0*z2 + f0*y1*z2 - 
   	 	f0*y3*z2 + f1*y3*z2 + f3*(y1*z0 - y2*z0 - y0*z1 + y2*z1 + y0*z2 - 
      		y1*z2) + f1*y0*z3 - f0*y1*z3 + f0*y2*z3 - f1*y2*z3 + 
    		f2*(y3*z0 + y0*z1 - y3*z1 - y0*z3 + y1*(z3 - z0)))/
   		(x1*y2*z0 - x1*y3*z0 - x0*y2*z1 + x0*y3*z1 - x1*y0*z2 + x0*y1*z2 - 
    		x0*y3*z2 + x1*y3*z2 + x3*(y1*z0 - y2*z0 - y0*z1 + y2*z1 + y0*z2 - 
      		y1*z2) + x1*y0*z3 - x0*y1*z3 + x0*y2*z3 - x1*y2*z3 + 
    		x2*(y3*z0 + y0*z1 - y3*z1 - y0*z3 + y1*(z3 - z0)));
        c = 	((-f1)*x2*z0 + f1*x3*z0 + f0*x2*z1 - f0*x3*z1 + f1*x0*z2 - f0*x1*z2 + 
    		f0*x3*z2 - f1*x3*z2 + f3*((-x1)*z0 + x2*z0 + x0*z1 - x2*z1 - x0*z2 + 
      		x1*z2) - f1*x0*z3 + f0*x1*z3 - f0*x2*z3 + f1*x2*z3 + 
    		f2*((-x3)*z0 - x0*z1 + x3*z1 + x1*(z0 - z3) + x0*z3))/
   		(x1*y2*z0 - x1*y3*z0 - x0*y2*z1 + x0*y3*z1 - x1*y0*z2 + x0*y1*z2 - 
    		x0*y3*z2 + x1*y3*z2 + x3*(y1*z0 - y2*z0 - y0*z1 + y2*z1 + y0*z2 - 
      		y1*z2) + x1*y0*z3 - x0*y1*z3 + x0*y2*z3 - x1*y2*z3 + 
    		x2*(y3*z0 + y0*z1 - y3*z1 - y0*z3 + y1*(z3 - z0)));
        d =	(f1*x2*y0 - f1*x3*y0 - f0*x2*y1 + f0*x3*y1 - f1*x0*y2 + f0*x1*y2 - 
    		f0*x3*y2 + f1*x3*y2 + f3*(x1*y0 - x2*y0 - x0*y1 + x2*y1 + x0*y2 - 
      		x1*y2) + f1*x0*y3 - f0*x1*y3 + f0*x2*y3 - f1*x2*y3 + 
    		f2*(x3*y0 + x0*y1 - x3*y1 - x0*y3 + x1*(y3 - y0)))/
   		(x1*y2*z0 - x1*y3*z0 - x0*y2*z1 + x0*y3*z1 - x1*y0*z2 + x0*y1*z2 - 
    		x0*y3*z2 + x1*y3*z2 + x3*(y1*z0 - y2*z0 - y0*z1 + y2*z1 + y0*z2 - 
      		y1*z2) + x1*y0*z3 - x0*y1*z3 + x0*y2*z3 - x1*y2*z3 + 
    		x2*(y3*z0 + y0*z1 - y3*z1 - y0*z3 + y1*(z3 - z0)));
    }  
	return a+_crd[0]*b+_crd[1]*c+_crd[2]*d;
}
double fsq(double* coeff_num, double* _crd)
{
	return coeff_num[4]*_crd[0]*_crd[0] + coeff_num[5]*_crd[1]*_crd[1] + coeff_num[6]*_crd[2]*_crd[2] + 
			coeff_num[7]*_crd[0]*_crd[1] + coeff_num[8]*_crd[0]*_crd[2] + coeff_num[9]*_crd[1]*_crd[2] + 
         		coeff_num[1]*_crd[0] + coeff_num[2]*_crd[1] + coeff_num[3]*_crd[2] + coeff_num[0];
}
double fsq(double* coeff_num, double x, double y, double z)
{
	return coeff_num[4]*x*x + coeff_num[5]*y*y + coeff_num[6]*z*z + 
			coeff_num[7]*x*y + coeff_num[8]*x*z + coeff_num[9]*y*z + 
         		coeff_num[1]*x + coeff_num[2]*y + coeff_num[3]*z + coeff_num[0];
}

void Method::intoRandomAxes(double* x, double* y, double *z)
{
	double c[3]={	(*x)*axis[0]+(*y)*axis[1]+(*z)*axis[2],
			(*x)*axis[3]+(*y)*axis[4]+(*z)*axis[5],
			(*x)*axis[6]+(*y)*axis[7]+(*z)*axis[8]};
	*x=c[0];*y=c[1];*z=c[2];
}
void Method::intoRandomAxes(double* c)
{
	intoRandomAxes(c,c+1,c+2);
}

void Method::intoRandomAxes(double* x, double* y, double *z, double *axes)
{
	double c[3]={	(*x)*axes[0]+(*y)*axes[1]+(*z)*axes[2],
			(*x)*axes[3]+(*y)*axes[4]+(*z)*axes[5],
			(*x)*axes[6]+(*y)*axes[7]+(*z)*axes[8]};
	*x=c[0];*y=c[1];*z=c[2];
}
void Method::intoRandomAxes(double* c, double *axes)
{
	intoRandomAxes(c,c+1,c+2,axes);
}

void Method::fromRandomAxes(double* c, double *a)
{
	double det = -a[2]*a[4]*a[6] + a[1]*a[5]*a[6] + a[2]*a[3]*a[7] - a[0]*a[5]*a[7] - a[1]*a[3]*a[8] + a[0]*a[4]*a[8];
	double axesR[9]={(a[4]*a[8] - a[5]*a[7])/det,(a[2]*a[7] - a[1]*a[8])/det,(a[1]*a[5] - a[2]*a[4])/det,
			(a[5]*a[6] - a[3]*a[8])/det,(a[0]*a[8] - a[2]*a[6])/det,(a[2]*a[3] - a[0]*a[5])/det,
			(a[3]*a[7] - a[4]*a[6])/det,(a[1]*a[6] - a[0]*a[7])/det,(a[0]*a[4] - a[1]*a[3])/det};
	
	intoRandomAxes(c,axesR);
}

double Method::interpolate_2_order(Thetraeder* t, double* _crd, Mesh* mesh)
{	
	Node* nodes[4] = {mesh->getNode(t->vert[0]), mesh->getNode(t->vert[1]), mesh->getNode(t->vert[2]), mesh->getNode(t->vert[3])};
	//return (nodes[0]->v+nodes[1]->v+nodes[2]->v+nodes[3]->v)/4;

	//check if thetraeder is close to z-border
	/*int isZ=t->checkZborder(mesh);
	if (isZ)
	{
		for (int i=0; i<4; i++)
		{
			double tmp = nodes[i]->coords[2];
			nodes[i]->coords[2] = nodes[i]->coords[0];
			nodes[i]->coords[0] = tmp;
			tmp = nodes[i]->u[3];
			nodes[i]->u[3] = nodes[i]->u[1];
			nodes[i]->u[1] = tmp;
		}
		double tmp=_crd[2]; _crd[2]=_crd[0]; _crd[0]=tmp;
	}*/
	
	/*double 	x0=nodes[0]->coords[0], x1=nodes[1]->coords[0], x2=nodes[2]->coords[0], x3=nodes[3]->coords[0],   
       		y0=nodes[0]->coords[1], y1=nodes[1]->coords[1], y2=nodes[2]->coords[1], y3=nodes[3]->coords[1],   
		z0=nodes[0]->coords[2], z1=nodes[1]->coords[2], z2=nodes[2]->coords[2], z3=nodes[3]->coords[2];
	double	zn3y=y2*( z0 - z1) + y0*(z1 - z2) + y1*(-z0 + z2),
		zn2y=y3*(-z0 + z1) + y1*(z0 - z3) + y0*(-z1 + z3),
		zn1y=y3*(-z0 + z2) + y2*(z0 - z3) + y0*(-z2 + z3),
		zn  =y3*(z0 - z1) + y1*(z0 - z3) + y0*(-2.0*z0 + z1 + z3),
		znyz=(y0 - y3)*(y0 - y3)*(z0 - z1)*(z0 - z1) - (y0 - y1)*(y0 - y1)*(z0 - z3)*(z0 - z3);
	/*if (fabs(zn3y)<0.0000000001 || fabs(zn2y)<0.0000000001 || fabs(zn1y)<0.0000000001 || fabs(zn)<0.0000000001 || fabs(znyz)<0.0000000001)
	{
		isZ=1;
		//printf("fault znam: %lf %lf %lf %lf %lf\n",zn3y*10000,zn2y*10000,zn1y*10000,zn*10000,znyz*10000);
		for (int i=0; i<4; i++)
		{
			double tmp = nodes[i]->coords[0];
			nodes[i]->coords[0] = nodes[i]->coords[1];
			nodes[i]->coords[1] = tmp;
		}
		double tmp=_crd[0]; _crd[0]=_crd[1]; _crd[1]=tmp;
		x0=nodes[0]->coords[0], x1=nodes[1]->coords[0], x2=nodes[2]->coords[0], x3=nodes[3]->coords[0];   
       		y0=nodes[0]->coords[1], y1=nodes[1]->coords[1], y2=nodes[2]->coords[1], y3=nodes[3]->coords[1];   
		z0=nodes[0]->coords[2], z1=nodes[1]->coords[2], z2=nodes[2]->coords[2], z3=nodes[3]->coords[2];
		zn3y=y2*( z0 - z1) + y0*(z1 - z2) + y1*(-z0 + z2);
		zn2y=y3*(-z0 + z1) + y1*(z0 - z3) + y0*(-z1 + z3);
		zn1y=y3*(-z0 + z2) + y2*(z0 - z3) + y0*(-z2 + z3);
		zn  =y3*(z0 - z1) + y1*(z0 - z3) + y0*(-2.0*z0 + z1 + z3);
		znyz=(y0 - y3)*(y0 - y3)*(z0 - z1)*(z0 - z1) - (y0 - y1)*(y0 - y1)*(z0 - z3)*(z0 - z3);
		//printf("new znam: %lf %lf %lf %lf %lf\n\n",zn3y*10000,zn2y*10000,zn1y*10000,zn*10000,znyz*10000);
	}*/

	double 	r0,r1,r2,r3, coeff_num[10]={0.0},c[12]={
	 	nodes[0]->coords[0], nodes[1]->coords[0], nodes[2]->coords[0], nodes[3]->coords[0],   
       		nodes[0]->coords[1], nodes[1]->coords[1], nodes[2]->coords[1], nodes[3]->coords[1],   
		nodes[0]->coords[2], nodes[1]->coords[2], nodes[2]->coords[2], nodes[3]->coords[2]};

		//printf("%10lf %10lf %10lf %10lf <- before rand\n%10lf %10lf %10lf %10lf\n%10lf %10lf %10lf %10lf\n",y0,y1,y2,y3, x0,x1,x2,x3,z0,z1,z2,z3); 
	double axes[9]; randomizeAxis(axes);
	intoRandomAxes(c,c+4,c+8,axes);
	intoRandomAxes(c+1,c+5,c+9,axes);
	intoRandomAxes(c+2,c+6,c+10,axes);
	intoRandomAxes(c+3,c+7,c+11,axes);

		//printf("%10lf %10lf %10lf %10lf <- after rand\n%10lf %10lf %10lf %10lf\n%10lf %10lf %10lf %10lf\n\n",y0,y1,y2,y3, x0,x1,x2,x3,z0,z1,z2,z3); 
	int min=0;
	for (int i=1;i<4;i++) if (c[4+i]<c[4+min]) min=i;
	for (int i=0;i<3;i++) {double tmp=c[i*4+min]; c[i*4+min]=c[i*4]; c[i*4]=tmp;};
	int max=0;
	for (int i=1;i<4;i++) if (c[4+i]>c[4+max]) max=i;
	for (int i=0;i<3;i++) {double tmp=c[i*4+max]; c[i*4+max]=c[i*4+3]; c[i*4+3]=tmp;};
	double	x0=c[0],x1=c[1],x2=c[2],x3=c[3],
		y0=c[4],y1=c[5],y2=c[6],y3=c[7],
		z0=c[8],z1=c[9],z2=c[10],z3=c[11];
	/*int min[3]={0},max[3]={0},diff=0;
	//double diff;
	for (int m=0; m<3; m++)
	{
		for (int i=1; i<4; i++)
			if (nodes[i]->coords[m] < nodes[min[m]]->coords[m]) min[m] = i;
	}
	if (min[1]) {Node* tmp = nodes[0]; nodes[0] = nodes[min[1]]; nodes[min[1]] = tmp;};
	for (int m=0; m<3; m++)
	{
		for (int i=1; i<4; i++)
			if (nodes[i]->coords[m] > nodes[max[m]]->coords[m]) max[m] = i;
	}
	if (max[1] != 3) {Node* tmp = nodes[3]; nodes[3] = nodes[max[1]]; nodes[max[1]] = tmp;};*/

    	if (y3-y0 < 0.0001) 
	{
		printf("%10lf %10lf %10lf %10lf <- skipping by y diff\n%10lf %10lf %10lf %10lf\n%10lf %10lf %10lf %10lf\n",y0,y1,y2,y3, x0,x1,x2,x3,z0,z1,z2,z3); return 0.0;
	}
        double volume = (x2*y1*z0 - x3*y1*z0 - x1*y2*z0 + x3*y2*z0 + x1*y3*z0 - x2*y3*z0 - 
                                    x2*y0*z1 + x3*y0*z1 + x0*y2*z1 - x3*y2*z1 - x0*y3*z1 + x2*y3*z1 + x1*y0*z2 - 
                                    x3*y0*z2 - x0*y1*z2 + x3*y1*z2 + x0*y3*z2 - x1*y3*z2 - 
                                    x1*y0*z3 + x2*y0*z3 + x0*y1*z3 - x2*y1*z3 - x0*y2*z3 + x1*y2*z3);
	if (fabs(volume)<0.000001)
	{
		printf("%10lf %10lf %10lf %10lf <- skipping by zero volume\n%10lf %10lf %10lf %10lf\n%10lf %10lf %10lf %10lf\n",y0,y1,y2,y3, x0,x1,x2,x3,z0,z1,z2,z3); return 0.0;
	}
	double dy03=(y0 - y3)*(y0 - y3),dz03=(z0 - z3)*(z0 - z3),
		dy02=(y0 - y2)*(y0 - y2),dz02=(z0 - z2)*(z0 - z2),
		dy01=(y0 - y1)*(y0 - y1),dz01=(z0 - z1)*(z0 - z1);
		//f40=f0-f3-f4*y0+f4*y3-f5*z0+f5*z3,
		//f41=f1-f0+f4*y0-f4*y1+f5*z0-f5*z1,
		//f42=(dy03*dz01-dy01*dz03)*(dy03*(-f0+f2+f4*y0-f4*y2+f5*z0-f5*z2)+dy02*f40)
		//		-(dy03*dz02-dy02*dz03)*(dy03*f41+dy01*f40),
	double zn3y=y2*( z0 - z1) + y0*(z1 - z2) + y1*(-z0 + z2),
	zn2y=y3*(-z0 + z1) + y1*(z0 - z3) + y0*(-z1 + z3),
	zn1y=y3*(-z0 + z2) + y2*(z0 - z3) + y0*(-z2 + z3),
	zn  =y3*(z0 - z1) + y1*(z0 - z3) + y0*(-2.0*z0 + z1 + z3),
	znyz=dy03*dz01 - dy01*dz03;

	if (fabs(zn3y)<0.0000001 || fabs(zn2y)<0.0000001 || fabs(zn1y)<0.0000001 || fabs(zn)<0.0000001 || fabs(znyz)<0.0000001)
	{
		//printf("%10lf %10lf %10lf %10lf <- skipping by zero znam\n%10lf %10lf %10lf %10lf\n%10lf %10lf %10lf %10lf\n",y0,y1,y2,y3, x0,x1,x2,x3,z0,z1,z2,z3); 
		return interpolate_2_order(t, _crd, mesh);
	}

	double  grad0[3]={nodes[0]->u[1],nodes[0]->u[2],nodes[0]->u[3]},
		grad1[3]={nodes[1]->u[1],nodes[1]->u[2],nodes[1]->u[3]},
		grad2[3]={nodes[2]->u[1],nodes[2]->u[2],nodes[2]->u[3]},
		grad3[3]={nodes[3]->u[1],nodes[3]->u[2],nodes[3]->u[3]};
	fromRandomAxes(grad0,axes);	
	fromRandomAxes(grad1,axes);	
	fromRandomAxes(grad2,axes);	
	fromRandomAxes(grad3,axes);	
        double f0=grad0[0],f1=grad1[0],f2=grad2[0],f3=grad3[0];   //solve for grad_x: 4x4
        r0 = (f3*x2*y1*z0 - f2*x3*y1*z0 - f3*x1*y2*z0 + f1*x3*y2*z0 + f2*x1*y3*z0 - 
                                    f1*x2*y3*z0 - f3*x2*y0*z1 + f2*x3*y0*z1 + f3*x0*y2*z1 - f0*x3*y2*z1 - 
                                    f2*x0*y3*z1 + f0*x2*y3*z1 + f3*x1*y0*z2 - f1*x3*y0*z2 - f3*x0*y1*z2 + 
                                    f0*x3*y1*z2 + f1*x0*y3*z2 - f0*x1*y3*z2 - f2*x1*y0*z3 + 
                                    f1*x2*y0*z3 + f2*x0*y1*z3 - f0*x2*y1*z3 - f1*x0*y2*z3 + 
                                    f0*x1*y2*z3)/volume;
        r1 = -((f2*y1*z0 - f3*y1*z0 - f1*y2*z0 + f3*y2*z0 + f1*y3*z0 - f2*y3*z0 - f2*y0*z1 + 
                                    f3*y0*z1 + f0*y2*z1 - f3*y2*z1 - f0*y3*z1 + f2*y3*z1 + 
                                    f1*y0*z2 - f3*y0*z2 - f0*y1*z2 + f3*y1*z2 + f0*y3*z2 - 
                                    f1*y3*z2 - f1*y0*z3 + f2*y0*z3 + f0*y1*z3 - f2*y1*z3 - f0*y2*z3 + f1*y2*z3)/
                                    (2.0*volume));
        r2 = (f2*x1*z0 - f3*x1*z0 - f1*x2*z0 + f3*x2*z0 + f1*x3*z0 - f2*x3*z0 - f2*x0*z1 + 
                                    f3*x0*z1 + f0*x2*z1 - f3*x2*z1 - f0*x3*z1 + f2*x3*z1 + 
                                    f1*x0*z2 - f3*x0*z2 - f0*x1*z2 + f3*x1*z2 + f0*x3*z2 - f1*x3*z2 - 
                                    f1*x0*z3 + f2*x0*z3 + f0*x1*z3 - f2*x1*z3 - f0*x2*z3 + f1*x2*z3)/
                                    volume;
        r3 = ((-f2)*x1*y0 + f3*x1*y0 + f1*x2*y0 - f3*x2*y0 - f1*x3*y0 + f2*x3*y0 + f2*x0*y1 - f3*x0*y1 - 
                                    f0*x2*y1 + f3*x2*y1 + f0*x3*y1 - f2*x3*y1 - 
                                    f1*x0*y2 + f3*x0*y2 + f0*x1*y2 - f3*x1*y2 - f0*x3*y2 + f1*x3*y2 + f1*x0*y3 - 
                                    f2*x0*y3 - f0*x1*y3 + f2*x1*y3 + f0*x2*y3 - f1*x2*y3)/
                                    volume;
    
               f0=nodes[0]->u[0]-x0*r0 - x0*x0*r1 - x0*y0*r2 - x0*z0*r3;    //solve for the rest: 6x6
               f1=nodes[1]->u[0]-x1*r0 - x1*x1*r1 - x1*y1*r2 - x1*z1*r3;
               f2=nodes[2]->u[0]-x2*r0 - x2*x2*r1 - x2*y2*r2 - x2*z2*r3;
               f3=nodes[3]->u[0]-x3*r0 - x3*x3*r1 - x3*y3*r2 - x3*z3*r3;          
        double       f4=grad0[1]-x0*r2,f5=grad0[2]-x0*r3,
	
		fn=dy03*(-f0 + f1 + f4*y0 - f4*y1 + f5*z0 - f5*z1) - dy01*(-f0 + f3 + f4*y0 - f4*y3 + f5*z0 - f5*z3),
		fn0=dy03*(-f0 + f2 + f4*y0 - f4*y2 + f5*z0 - f5*z2) - dy02*(-f0 + f3 + f4*y0 - f4*y3 + f5*z0 - f5*z3);
        coeff_num[0]= /*f0 - f4*y0 - f5*z0 - y0*y0*f40/dy02 

		- ((2.0*y0*y3*z0*z1-y3*y3*z0*z0+y0*y0*z3*(z3-2.0*z0))*(dy03*f41+dy01*f40))
			/dy03/(dy03*dz01 - dy01*dz03)

		- ((y1*z0-y0*z1)*(y3*z0-y0*z3)*f42)
			/dy03
			/(y2*(z0-z1)+y0*(z1-z2)+y1*(z2-z0))
			/(y3*(z1-z0)+y1*(z0-z3)+y0*(z3-z1))
			/(y3*(z0-z1)+y1*(z0-z3)+y0*(z1+z3-2.0*z0))
			/(y3*(z2-z0)+y2*(z0-z3)+y0*(-z2+z3));*/
		f0 - f4*y0 - f5*z0 + y0*(-f0 + f3 + f4*y0 - f4*y3 + f5*z0 - f5*z3)/dy03*y0 - 

   	fn
	/dy03
	*(2.0*y0*y3*z0*z0 - y3*y3*z0*z0 + y0*y0*z3*(-2.0*z0 + z3))
	/znyz - 

   	1.0
	/zn
	*(-(dy03*dz02 - dy02*dz03)*fn + znyz*fn0)
	/dy03
	*(y3*z0 - y0*z3)
	/zn3y
	/zn2y
	*(y1*z0 - y0*z1)
	/zn1y;

        coeff_num[1]=r0;

        coeff_num[2]= /* f4 + 2.0*y0*f40/dy03 + 
   			((y3*z0*(z1-z0)+y1*z0*(z3-z0)+y0*(-2.0*z1*z3+z0*(z1+z3)))*f42)
			/dy03
			/(y2*(z0-z1)  f5 - (2.0*z0*(dy03*(-f0 + f1 + f4*y0 - f4*y1 + f5*z0 - f5*z1) - dy01*(-f0 + f3 + f4*y0 - f4*y3 + f5*z0 - f5*z3)))/
    (dy03*dz01 - dy01*dz03) - ((-2.0*y1*y3*z0 - y0*y0*(z1 + z3) + y0*(y3*(z0 + z1) + y1*(z0 + z3)))*
     ((-(dy03*dz02 - dy02*dz03))*(dy03*(-f0 + f1 + f4*y0 - f4*y1 + f5*z0 - f5*z1) - 
        dy01*(-f0 + f3 + f4*y0 - f4*y3 + f5*z0 - f5*z3)) + (dy03*dz01 - dy01*dz03)*
       (dy03*(-f0 + f2 + f4*y0 - f4*y2 + f5*z0 - f5*z2) - dy02*(-f0 + f3 + f4*y0 - f4*y3 + f5*z0 - f5*z3))))/
    (dy03*(y2*(z0 - z1) + y0*(z1 - z2) + y1*(-z0 + z2))*(y3*(-z0 + z1) + y1*(z0 - z3) + y0*(-z1 + z3))*
     (y3*(z0 - z1) + y1*(z0 - z3) + y0*(-2.0*z0 + z1 + z3))*(y3*(-z0 + z2) + y2*(z0 - z3) + y0*(-z2 + z3))), +y0*(z1-z2)+y1*(z2-z0))
			/(y3*(z1-z0)+y1*(z0-z3)+y0*(z3-z1))
			/(y3*(z0-z1)+y1*(z0-z3)+y0*(z1+z3-2.0*z0))
			/(y3*(z2-z0)+y2*(z0-z3)+y0*(-z2+z3));*/
		
	f4 - 2.0*(-f0 + f3 + f4*y0 - f4*y3 + f5*z0 - f5*z3)/dy03*y0 + 
   2.0*fn/dy03*y0/znyz*dz03 - 

	1.0
	/zn3y
   	*(y3*z0*(-z0 + z1) + y1*z0*(-z0 + z3) + y0*(-2.0*z1*z3 + z0*(z1 + z3)))
	/zn2y
	/dy03
	*(-(dy03*dz02 - dy02*dz03)*fn + znyz*fn0)
	/zn
	/zn1y;

        coeff_num[3]=/*f5 - (2.0*z0*(dy03*f41+dy01*f40))/dy03/(dy03*dz01 - dy01*dz03) - 

   ((-2.0*y1*y3*z0-y0*y0*(z1+z3)+y0*(y3*(z0+z1)+y1*(z0+z3)))*f42)
			/dy03
			/(y2*(z0-z1)+y0*(z1-z2)+y1*(z2-z0))
			/(y3*(z1-z0)+y1*(z0-z3)+y0*(z3-z1))
			/(y3*(z0-z1)+y1*(z0-z3)+y0*(z1+z3-2.0*z0))
			/(y3*(z2-z0)+y2*(z0-z3)+y0*(-z2+z3));*/
  f5 - 2.0*fn/znyz*z0
	- 
	1.0
	/zn3y
	*(-2.0*y1*y3*z0 - y0*y0*(z1 + z3) + y0*(y3*(z0 + z1) + y1*(z0 + z3)))
	/dy03
	/zn2y
	*(-(dy03*dz02 - dy02*dz03)*(dy03*(-f0 + f1 + f4*y0 - f4*y1 + f5*z0 - f5*z1) - dy01*(-f0 + f3 + f4*y0 - f4*y3 + f5*z0 - f5*z3)) + znyz*(dy03*(-f0 + f2 + f4*y0 - f4*y2 + f5*z0 - f5*z2) - dy02*(-f0 + f3 + f4*y0 - f4*y3 + f5*z0 - f5*z3)))
	/zn
	/zn1y;



        coeff_num[4]=r1;

        coeff_num[5]=/*-(f40
			+dz03*(dy03*f41+dy01*f40)/(dy03*dz01 - dy01*dz03) 
			+((z0-z1)*(z0-z3)*f42
				
			 )
			/(y2*(z0-z1)+y0*(z1-z2)+y1*(z2-z0))
			/(y3*(z1-z0)+y1*(z0-z3)+y0*(z3-z1))
			/(y3*(z0-z1)+y1*(z0-z3)+y0*(z1+z3-2.0*z0))
			/(y3*(z2-z0)+y2*(z0-z3)+y0*(-z2+z3))
			)/dy03;*/
  (-1.0/dy03)*(
	f0 - f3 - f4*y0 + f4*y3 - f5*z0 + f5*z3 + 
    	fn/znyz*dz03 + 
    	1.0
	/zn1y
	*(z0 - z1)
	/zn
	*(z0 - z3)
	/zn2y
	*(-(dy03*dz02 - dy02*dz03)*fn + znyz*fn0)
	/zn3y
	);

        coeff_num[6]=/*(dy03*f41+dy01*f40+
			((y0-y1)*f42
				)
				/((y0-y3)*(y2*(z0-z1)+y0*(z1-z2)+y1*(z2-z0))*
					(y3*(z2-z0)+y2*(z0-z3)+y0*(z3-z2)))			
			)/(dy03*dz01 - dy01*dz03) ;*/

  	(fn + 
    	1.0
	/zn3y
	*(y0 - y1)
	/(y0 - y3)
	*(-(dy03*dz02 - dy02*dz03)*fn + znyz*fn0)
	/zn1y)
	/znyz;

        coeff_num[7]=r2;
        coeff_num[8]=r3;
        coeff_num[9]= /*-f42
			/dy03
			/(y2*(z0-z1)+y0*(z1-z2)+y1*(z2-z0))
			/(y3*(z1-z0)+y1*(z0-z3)+y0*(z3-z1))
			/(y3*(z2-z0)+y2*(z0-z3)+y0*(z3-z2));*/
	1.0
	/dy03
	/zn3y
	*(-(dy03*dz02 - dy02*dz03)*fn + znyz*fn0)
	/(-zn2y)
	/zn1y;
        
    
	//printf("%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n", coeff_num[0], coeff_num[1], coeff_num[2], coeff_num[3], coeff_num[4], coeff_num[5], coeff_num[6], coeff_num[7], coeff_num[8], coeff_num[9]);
	//printf("%lf %lf %lf %lf %lf %lf  \n", coeff_num[0], coeff_num[1], coeff_num[2], coeff_num[4], coeff_num[7], coeff_num[8]);
	//printf("%lf %lf %lf %lf \n", coeff_num[3], coeff_num[5], coeff_num[6], coeff_num[9]);
	double _crdR[3]={_crd[0],_crd[1],_crd[2]};
	intoRandomAxes(_crdR);
	double res = fsq(coeff_num,_crdR), minU = nodes[0]->u[0], maxU = nodes[0]->u[0];
	for (int i=1; i<4; i++)
	{
		if (nodes[i]->u[0] < minU) minU = nodes[i]->u[0];
		if (nodes[i]->u[0] > maxU) maxU = nodes[i]->u[0];
	}
	if (fabs(res) > fabs(10000*(maxU+10))) 
	{
		printf("%lf %lf %lf %lf   %lf\n", fsq(coeff_num,c[0],c[4],c[8]), fsq(coeff_num,c[1],c[5],c[9]), fsq(coeff_num,c[2],c[6],c[10]), fsq(coeff_num,c[3],c[7],c[11]),/* nodes[0]->u[0], nodes[1]->u[0], nodes[2]->u[0], nodes[3]->u[0], nodes[0]->u[1], nodes[1]->u[1], nodes[2]->u[1], nodes[3]->u[1], nodes[0]->u[2], nodes[1]->u[2], nodes[2]->u[2], nodes[3]->u[2], nodes[0]->u[3], nodes[1]->u[3], nodes[2]->u[3], nodes[3]->u[3],/*coeff_num[0], coeff_num[1], coeff_num[2], coeff_num[3], coeff_num[4], coeff_num[5], coeff_num[6], coeff_num[7], coeff_num[8], coeff_num[9],*/ res);
		//printf("%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n\n", coeff_num[0], coeff_num[1], coeff_num[2], coeff_num[3], coeff_num[4], coeff_num[5], coeff_num[6], coeff_num[7], coeff_num[8], coeff_num[9]);
	printf("%30lf %30lf %30lf\n%30lf %30lf %30lf  \n", coeff_num[0], coeff_num[2], coeff_num[3], coeff_num[5], coeff_num[6], coeff_num[9]);
	printf("%lf %lf %lf %lf \n", coeff_num[1], coeff_num[4], coeff_num[7], coeff_num[8]);
	printf("%lf %lf %lf %lf %lf \n", zn3y,zn2y,zn1y,zn,znyz);
	printf("%lf %lf %lf %lf \n%lf %lf %lf %lf \n%lf %lf %lf %lf \n", c[0],c[1],c[2],c[3],c[4],c[5],c[6],c[7],c[8],c[9],c[10],c[11]);
	printf("%lf %lf %lf %lf \n%lf %lf %lf %lf \n%lf %lf %lf %lf \n\n--------------\n", nodes[0]->coords[0], nodes[1]->coords[0], nodes[2]->coords[0], nodes[3]->coords[0], nodes[0]->coords[1], nodes[1]->coords[1], nodes[2]->coords[1], nodes[3]->coords[1], nodes[0]->coords[2], nodes[1]->coords[2], nodes[2]->coords[2], nodes[3]->coords[2]);
		//printf("");
		//res = maxU;
	}
	if (res < minU) res = minU;
	if (res > maxU) res = maxU;
	//if (res!=res) res=0.0;
	
	//check if thetraeder is close to z-border
	/*if (isZ)
	{
		for (int i=0; i<4; i++)
		{
			double tmp = nodes[i]->coords[0];
			nodes[i]->coords[0] = nodes[i]->coords[2];
			nodes[i]->coords[2] = tmp;
			tmp = nodes[i]->u[3];
			nodes[i]->u[3] = nodes[i]->u[1];
			nodes[i]->u[1] = tmp;
		}
		double tmp=_crd[0]; _crd[0]=_crd[2]; _crd[2]=tmp;
	}*/

	return res;
}

double fcb(double* a, double* c)
{
       return a[0] + a[1]*c[0] + a[2]*c[1] + a[3]*c[2] + a[4]*c[0]*c[0] + a[5]*c[1]*c[1] + a[6]*c[2]*c[2] + a[7]*c[0]*c[1] + a[8]*c[0]*c[2] + a[9]*c[1]*c[2] + a[10]*c[0]*c[0]*c[0] + a[11]*c[1]*c[1]*c[1] + a[12]*c[2]*c[2]*c[2] + a[13]*c[0]*c[0]*c[1] + a[14]*c[0]*c[0]*c[2] + a[15]*c[1]*c[1]*c[0] + a[16]*c[1]*c[1]*c[2] + a[17]*c[2]*c[2]*c[0] + a[18]*c[2]*c[2]*c[1] + a[19]*c[0]*c[1]*c[2];
                        
}

double fcb(double* a, double x, double y, double z)
{
       return a[0] + a[1]*x + a[2]*y + a[3]*z + a[4]*x*x + a[5]*y*y + a[6]*z*z + a[7]*x*y + a[8]*x*z + a[9]*y*z + a[10]*x*x*x + a[11]*y*y*y + a[12]*z*z*z + a[13]*x*x*y + a[14]*x*x*z + a[15]*y*y*x + a[16]*y*y*z + a[17]*z*z*x + a[18]*z*z*y + a[19]*x*y*z;
                        
}
double Method::interpolate_3_order(Thetraeder* t, double* _crd, Mesh* mesh)
{
	Node* nodes[4] = {mesh->getNode(t->vert[0]), mesh->getNode(t->vert[1]), mesh->getNode(t->vert[2]), mesh->getNode(t->vert[3])};
	//int min[3]={0},max[3]={0},diff=0;
	double 	c[12]={
	 	nodes[0]->coords[0], nodes[1]->coords[0], nodes[2]->coords[0], nodes[3]->coords[0],   
       		nodes[0]->coords[1], nodes[1]->coords[1], nodes[2]->coords[1], nodes[3]->coords[1],   
		nodes[0]->coords[2], nodes[1]->coords[2], nodes[2]->coords[2], nodes[3]->coords[2]};

		//printf("%10lf %10lf %10lf %10lf <- before rand\n%10lf %10lf %10lf %10lf\n%10lf %10lf %10lf %10lf\n",y0,y1,y2,y3, x0,x1,x2,x3,z0,z1,z2,z3); 
	double axes[9]; randomizeAxis(axes);
	intoRandomAxes(c,c+4,c+8,axes);
	intoRandomAxes(c+1,c+5,c+9,axes);
	intoRandomAxes(c+2,c+6,c+10,axes);
	intoRandomAxes(c+3,c+7,c+11,axes);
	int min=0;
	for (int i=1;i<4;i++) if (c[4+i]<c[4+min]) min=i;
	for (int i=0;i<3;i++) {double tmp=c[i*4+min]; c[i*4+min]=c[i*4]; c[i*4]=tmp;};
	int max=0;
	for (int i=1;i<4;i++) if (c[4+i]>c[4+max]) max=i;
	for (int i=0;i<3;i++) {double tmp=c[i*4+max]; c[i*4+max]=c[i*4+3]; c[i*4+3]=tmp;};
	//for (int i=0;i<12;i++) c[i]+=10.0;
	double	x0=c[0],x1=c[1],x2=c[2],x3=c[3],
		y0=c[4],y1=c[5],y2=c[6],y3=c[7],
		z0=c[8],z1=c[9],z2=c[10],z3=c[11];

	/*for (int m=0; m<3; m++)
	{
		for (int i=1; i<4; i++)
			if (nodes[i]->coords[m] < nodes[min[m]]->coords[m]) min[m] = i;
	}
	if (min[1]) {Node* tmp = nodes[0]; nodes[0] = nodes[min[1]]; nodes[min[1]] = tmp;};
	for (int m=0; m<3; m++)
	{
		for (int i=1; i<4; i++)
			if (nodes[i]->coords[m] > nodes[max[m]]->coords[m]) max[m] = i;
	}
	if (max[1] != 3) {Node* tmp = nodes[3]; nodes[3] = nodes[max[1]]; nodes[max[1]] = tmp;};

	double 	x0=nodes[0]->coords[0], x1=nodes[1]->coords[0], x2=nodes[2]->coords[0], x3=nodes[3]->coords[0], //
       		y0=nodes[0]->coords[1], y1=nodes[1]->coords[1], y2=nodes[2]->coords[1], y3=nodes[3]->coords[1],   //
       		z0=nodes[0]->coords[2], z1=nodes[1]->coords[2], z2=nodes[2]->coords[2], z3=nodes[3]->coords[2],  */ //
       	double	coeff_num[20];
	double volume = x2*y1*z0 - x3*y1*z0 - x1*y2*z0 + x3*y2*z0 + x1*y3*z0 - x2*y3*z0 - x2*y0*z1 + x3*y0*z1 + x0*y2*z1 - x3*y2*z1 - x0*y3*z1 + x2*y3*z1 + x1*y0*z2 - x3*y0*z2 - x0*y1*z2 + x3*y1*z2 + x0*y3*z2 - x1*y3*z2 -x1*y0*z3 + x2*y0*z3 + x0*y1*z3 - x2*y1*z3 - x0*y2*z3 + x1*y2*z3,
		det=(y1*(y2*(z1 - z2)*(z0 - z3) - y3*(z0 - z2)*(z1 - z3)) + y0*(y3*(z1 - z2)*(z0 - z3) - y2*(z0 - z2)*(z1 - z3) + y1*(z0 - z1)*(z2 - z3)) + y2*y3*(z0 - z1)*(z2 - z3));
    //if (y3-y0 < 0.0001) {printf("%10f <- skipping\n",y3-y0); continue;};
    if (fabs(volume) < 0.00000001) {printf("%10f <- skipping by volume\n",volume*10000000); return 0.0;};
    if (fabs(det) < 0.0000000001) 
	{
		printf("%10f <- skipping by det\n",det*10000000); 
		//printf("%10lf %10lf %10lf %10lf \n%10lf %10lf %10lf %10lf\n%10lf %10lf %10lf %10lf\n\n", x0*2.0,x1*2.0,x2*2.0,x3*2.0,y0,y1,y2,y3,z0,z1,z2,z3);
		return interpolate_3_order(t,  _crd, mesh);
	};
	double  grad0[3]={nodes[0]->u[4],nodes[0]->u[5],nodes[0]->u[6]},
		grad1[3]={nodes[1]->u[4],nodes[1]->u[5],nodes[1]->u[6]},
		grad2[3]={nodes[2]->u[4],nodes[2]->u[5],nodes[2]->u[6]},
		grad3[3]={nodes[3]->u[4],nodes[3]->u[5],nodes[3]->u[6]};
	fromRandomAxes(grad0,axes);	
	fromRandomAxes(grad1,axes);	
	fromRandomAxes(grad2,axes);	
	fromRandomAxes(grad3,axes);	
	fromRandomAxes(grad0,axes);	
	fromRandomAxes(grad1,axes);	
	fromRandomAxes(grad2,axes);	
	fromRandomAxes(grad3,axes);	
    {                                                                
        double f0=grad0[0],f1=grad1[0],f2=grad2[0],f3=grad3[0];              
										//solve for grad_xx: 4x4
        coeff_num[4]  = ((-f1)*x3*y2*z0 + f1*x2*y3*z0 + f0*x3*y2*z1 - f0*x2*y3*z1 + f1*x3*y0*z2 - f0*x3*y1*z2 - f1*x0*y3*z2 + f0*x1*y3*z2 + f3*((-x2)*y1*z0 + x1*y2*z0 + x2*y0*z1 - x0*y2*z1 - x1*y0*z2 + x0*y1*z2) - f1*x2*y0*z3 + f0*x2*y1*z3 + f1*x0*y2*z3 - f0*x1*y2*z3 + f2*(x3*y1*z0 - x1*y3*z0 - x3*y0*z1 + x0*y3*z1 + x1*y0*z3 - x0*y1*z3))/(-2.0*volume);
        coeff_num[10] =   (f1*y2*z0 - f1*y3*z0 - f0*y2*z1 + f0*y3*z1 - f1*y0*z2 + f0*y1*z2 - f0*y3*z2 + f1*y3*z2 + f3*(y1*z0 - y2*z0 - y0*z1 + y2*z1 + y0*z2 - y1*z2) + f1*y0*z3 - f0*y1*z3 + f0*y2*z3 - f1*y2*z3 + f2*(y3*z0 + y0*z1 - y3*z1 - y0*z3 + y1*(-z0 + z3)))/(-6.0*volume);
        coeff_num[13] =   ((-f1)*x2*z0 + f1*x3*z0 + f0*x2*z1 - f0*x3*z1 + f1*x0*z2 - f0*x1*z2 + f0*x3*z2 - f1*x3*z2 + f3*((-x1)*z0 + x2*z0 + x0*z1 - x2*z1 - x0*z2 + x1*z2) - f1*x0*z3 + f0*x1*z3 - f0*x2*z3 + f1*x2*z3 + f2*((-x3)*z0 - x0*z1 + x3*z1 + x1*(z0 - z3) + x0*z3))/(-2.0*volume);
        coeff_num[14] = (f1*x2*y0 - f1*x3*y0 - f0*x2*y1 + f0*x3*y1 - f1*x0*y2 + f0*x1*y2 - f0*x3*y2 + f1*x3*y2 + f3*(x1*y0 - x2*y0 - x0*y1 + x2*y1 + x0*y2 - x1*y2) + f1*x0*y3 - f0*x1*y3 + f0*x2*y3 - f1*x2*y3 + f2*(x3*y0 + x0*y1 - x3*y1 - x0*y3 + x1*(-y0 + y3)))/(-2.0*volume);
    }  
    {                                                      
        double f4=grad0[1],f5=grad1[1],f6=grad2[1],f7=grad3[1];              
										//solve for grad_yy: 4x4
        coeff_num[5]  = ((-f5)*x3*y2*z0 + f5*x2*y3*z0 + f4*x3*y2*z1 - f4*x2*y3*z1 + f5*x3*y0*z2 - f4*x3*y1*z2 - f5*x0*y3*z2 + f4*x1*y3*z2 + f7*((-x2)*y1*z0 + x1*y2*z0 + x2*y0*z1 - x0*y2*z1 - x1*y0*z2 + x0*y1*z2) - f5*x2*y0*z3 + f4*x2*y1*z3 + f5*x0*y2*z3 - f4*x1*y2*z3 + f6*(x3*y1*z0 - x1*y3*z0 - x3*y0*z1 + x0*y3*z1 + x1*y0*z3 - x0*y1*z3))/(-2.0*volume);
        coeff_num[11] =   ((-f5)*x2*z0 + f5*x3*z0 + f4*x2*z1 - f4*x3*z1 + f5*x0*z2 - f4*x1*z2 + f4*x3*z2 - f5*x3*z2 + f7*((-x1)*z0 + x2*z0 + x0*z1 - x2*z1 - x0*z2 + x1*z2) - f5*x0*z3 + f4*x1*z3 - f4*x2*z3 + f5*x2*z3 + f6*((-x3)*z0 - x0*z1 + x3*z1 + x1*(z0 - z3) + x0*z3))/(-6.0*volume);
        coeff_num[15] =   (f5*y2*z0 - f5*y3*z0 - f4*y2*z1 + f4*y3*z1 - f5*y0*z2 + f4*y1*z2 - f4*y3*z2 + f5*y3*z2 + f7*(y1*z0 - y2*z0 - y0*z1 + y2*z1 + y0*z2 - y1*z2) + f5*y0*z3 - f4*y1*z3 + f4*y2*z3 - f5*y2*z3 + f6*(y3*z0 + y0*z1 - y3*z1 - y0*z3 + y1*(-z0 + z3)))/(-2.0*volume);
        coeff_num[16] =     (f5*x2*y0 - f5*x3*y0 - f4*x2*y1 + f4*x3*y1 - f5*x0*y2 + f4*x1*y2 - f4*x3*y2 + f5*x3*y2 + f7*(x1*y0 - x2*y0 - x0*y1 + x2*y1 + x0*y2 - x1*y2) + f5*x0*y3 - f4*x1*y3 + f4*x2*y3 - f5*x2*y3 + f6*(x3*y0 + x0*y1 - x3*y1 - x0*y3 + x1*(-y0 + y3)))/(-2.0*volume);
    } 
    {                                                      
        double f8=grad0[2],f9=grad1[2],f10=grad2[2],f11=grad3[2];              
										//solve for grad_zz: 4x4
        coeff_num[6]  = ((-f9)*x3*y2*z0 + f9*x2*y3*z0 + f8*x3*y2*z1 - f8*x2*y3*z1 + f9*x3*y0*z2 - f8*x3*y1*z2 - f9*x0*y3*z2 + f8*x1*y3*z2 + f11*((-x2)*y1*z0 + x1*y2*z0 + x2*y0*z1 - x0*y2*z1 - x1*y0*z2 + x0*y1*z2) - f9*x2*y0*z3 + f8*x2*y1*z3 + f9*x0*y2*z3 - f8*x1*y2*z3 + f10*(x3*y1*z0 - x1*y3*z0 - x3*y0*z1 + x0*y3*z1 + x1*y0*z3 - x0*y1*z3))/(-2.0*volume);
        coeff_num[12] =   (f9*x2*y0 - f9*x3*y0 - f8*x2*y1 + f8*x3*y1 - f9*x0*y2 + f8*x1*y2 - f8*x3*y2 + f9*x3*y2 + f11*(x1*y0 - x2*y0 - x0*y1 + x2*y1 + x0*y2 - x1*y2) + f9*x0*y3 - f8*x1*y3 + f8*x2*y3 - f9*x2*y3 + f10*(x3*y0 + x0*y1 - x3*y1 - x0*y3 + x1*(-y0 + y3)))/(-6.0*volume);
        coeff_num[17] =   (f9*y2*z0 - f9*y3*z0 - f8*y2*z1 + f8*y3*z1 - f9*y0*z2 + f8*y1*z2 - f8*y3*z2 + f9*y3*z2 + f11*(y1*z0 - y2*z0 - y0*z1 + y2*z1 + y0*z2 - y1*z2) + f9*y0*z3 - f8*y1*z3 + f8*y2*z3 - f9*y2*z3 + f10*(y3*z0 + y0*z1 - y3*z1 - y0*z3 + y1*(-z0 + z3)))/(-2.0*volume);
        coeff_num[18] =     ((-f9)*x2*z0 + f9*x3*z0 + f8*x2*z1 - f8*x3*z1 + f9*x0*z2 - f8*x1*z2 + f8*x3*z2 - f9*x3*z2 + f11*((-x1)*z0 + x2*z0 + x0*z1 - x2*z1 - x0*z2 + x1*z2) - f9*x0*z3 + f8*x1*z3 - f8*x2*z3 + f9*x2*z3 + f10*((-x3)*z0 - x0*z1 + x3*z1 + x1*(z0 - z3) + x0*z3))/(-2.0*volume);
    }     
    {              
	double  grad0[3]={nodes[0]->u[1],nodes[0]->u[2],nodes[0]->u[3]},
		grad1[3]={nodes[1]->u[1],nodes[1]->u[2],nodes[1]->u[3]},
		grad2[3]={nodes[2]->u[1],nodes[2]->u[2],nodes[2]->u[3]},
		grad3[3]={nodes[3]->u[1],nodes[3]->u[2],nodes[3]->u[3]};
	fromRandomAxes(grad0,axes);	
	fromRandomAxes(grad1,axes);	
	fromRandomAxes(grad2,axes);	
	fromRandomAxes(grad3,axes);	                                                 
        double 	f12=grad0[0]-2.0*x0*coeff_num[4]-3.0*x0*x0*coeff_num[10]-2.0*x0*y0*coeff_num[13]-2.0*x0*z0*coeff_num[14]-y0*y0*coeff_num[15]-z0*z0*coeff_num[17],
		f13=grad1[0]-2.0*x1*coeff_num[4]-3.0*x1*x1*coeff_num[10]-2.0*x1*y1*coeff_num[13]-2.0*x1*z1*coeff_num[14]-y1*y1*coeff_num[15]-z1*z1*coeff_num[17],
		f14=grad2[0]-2.0*x2*coeff_num[4]-3.0*x2*x2*coeff_num[10]-2.0*x2*y2*coeff_num[13]-2.0*x2*z2*coeff_num[14]-y2*y2*coeff_num[15]-z2*z2*coeff_num[17],
		f15=grad3[0]-2.0*x3*coeff_num[4]-3.0*x3*x3*coeff_num[10]-2.0*x3*y3*coeff_num[13]-2.0*x3*z3*coeff_num[14]-y3*y3*coeff_num[15]-z3*z3*coeff_num[17];
								            //solve for grad_x: 4x4
        coeff_num[1]  = ((-f13)*y0*y3*z0*z2 + f13*y2*y3*z0*z2 + f12*y1*y3*z1*z2 - f12*y2*y3*z1*z2 + f15*(y1*y2*z0*(z1 - z2) + y0*(y1*(z0 - z1)*z2 + y2*z1*(-z0 + z2))) + f13*y0*y2*z0*z3 - f13*y2*y3*z0*z3 - f12*y1*y2*z1*z3 + f12*y2*y3*z1*z3 - f13*y0*y2*z2*z3 + f12*y1*y2*z2*z3 + f13*y0*y3*z2*z3 - f12*y1*y3*z2*z3 + f14*(y1*y3*z0*(-z1 + z3) + y0*(y3*z1*(z0 - z3) + y1*(-z0 + z1)*z3)))/det;
        coeff_num[7] =     (f13*y0*z0*z2 - f13*y2*z0*z2 - f12*y1*z1*z2 + f12*y2*z1*z2 + f15*(y0*z0*(z1 - z2) + y2*(z0 - z1)*z2 + y1*z1*(-z0 + z2)) - f13*y0*z0*z3 + f13*y3*z0*z3 + f12*y1*z1*z3 - f12*y3*z1*z3 - f12*y2*z2*z3 + f13*y2*z2*z3 + f12*y3*z2*z3 - f13*y3*z2*z3 + f14*(y1*z1*(z0 - z3) + y3*(-z0 + z1)*z3 + y0*z0*(-z1 + z3)))/det;
        coeff_num[8] =   ((-f13)*y0*y2*z0 + f13*y0*y3*z0 + f12*y1*y2*z1 - f12*y1*y3*z1 + f13*y0*y2*z2 - f12*y1*y2*z2 + f12*y2*y3*z2 - f13*y2*y3*z2 + f15*(y1*y2*(-z1 + z2) + y0*((-y1)*z0 + y2*z0 + y1*z1 - y2*z2)) - f13*y0*y3*z3 + f12*y1*y3*z3 - f12*y2*y3*z3 + f13*y2*y3*z3 + f14*(y1*y3*(z1 - z3) + y0*(y1*z0 - y3*z0 - y1*z1 + y3*z3)))/det;
        coeff_num[19] =   (f13*y2*z0 - f13*y3*z0 - f12*y2*z1 + f12*y3*z1 - f13*y0*z2 + f12*y1*z2 - f12*y3*z2 + f13*y3*z2 + f15*(y1*z0 - y2*z0 - y0*z1 + y2*z1 + y0*z2 - y1*z2) + f13*y0*z3 - f12*y1*z3 + f12*y2*z3 - f13*y2*z3 + f14*(y3*z0 + y0*z1 - y3*z1 - y0*z3 + y1*(-z0 + z3)))/det;
    }  
    {                       
	double *a=coeff_num;                                         
        double 	f12=nodes[0]->u[0]-x0*a[1]-x0*x0*a[4]-y0*y0*a[5]-z0*z0*a[6]-x0*y0*a[7]-x0*z0*a[8]-x0*x0*x0*a[10]-y0*y0*y0*a[11]-z0*z0*z0*a[12]-x0*x0*y0*a[13]-x0*x0*z0*a[14]-y0*y0*x0*a[15]-y0*y0*z0*a[16]-z0*z0*x0*a[17]-z0*z0*y0*a[18]-x0*y0*z0*a[19],
		f13=nodes[1]->u[0]-x1*a[1]-x1*x1*a[4]-y1*y1*a[5]-z1*z1*a[6]-x1*y1*a[7]-x1*z1*a[8]-x1*x1*x1*a[10]-y1*y1*y1*a[11]-z1*z1*z1*a[12]-x1*x1*y1*a[13]-x1*x1*z1*a[14]-y1*y1*x1*a[15]-y1*y1*z1*a[16]-z1*z1*x1*a[17]-z1*z1*y1*a[18]-x1*y1*z1*a[19],
		f14=nodes[2]->u[0]-x2*a[1]-x2*x2*a[4]-y2*y2*a[5]-z2*z2*a[6]-x2*y2*a[7]-x2*z2*a[8]-x2*x2*x2*a[10]-y2*y2*y2*a[11]-z2*z2*z2*a[12]-x2*x2*y2*a[13]-x2*x2*z2*a[14]-y2*y2*x2*a[15]-y2*y2*z2*a[16]-z2*z2*x2*a[17]-z2*z2*y2*a[18]-x2*y2*z2*a[19],
		f15=nodes[3]->u[0]-x3*a[1]-x3*x3*a[4]-y3*y3*a[5]-z3*z3*a[6]-x3*y3*a[7]-x3*z3*a[8]-x3*x3*x3*a[10]-y3*y3*y3*a[11]-z3*z3*z3*a[12]-x3*x3*y3*a[13]-x3*x3*z3*a[14]-y3*y3*x3*a[15]-y3*y3*z3*a[16]-z3*z3*x3*a[17]-z3*z3*y3*a[18]-x3*y3*z3*a[19];
								              //solve for rest: 4x4
        coeff_num[0]  = ((-f13)*y0*y3*z0*z2 + f13*y2*y3*z0*z2 + f12*y1*y3*z1*z2 - f12*y2*y3*z1*z2 + f15*(y1*y2*z0*(z1 - z2) + y0*(y1*(z0 - z1)*z2 + y2*z1*(-z0 + z2))) + f13*y0*y2*z0*z3 - f13*y2*y3*z0*z3 - f12*y1*y2*z1*z3 + f12*y2*y3*z1*z3 - f13*y0*y2*z2*z3 + f12*y1*y2*z2*z3 + f13*y0*y3*z2*z3 - f12*y1*y3*z2*z3 + f14*(y1*y3*z0*(-z1 + z3) + y0*(y3*z1*(z0 - z3) + y1*(-z0 + z1)*z3)))/det;
        coeff_num[2] =   (f13*y0*z0*z2 - f13*y2*z0*z2 - f12*y1*z1*z2 + f12*y2*z1*z2 + f15*(y0*z0*(z1 - z2) + y2*(z0 - z1)*z2 + y1*z1*(-z0 + z2)) - f13*y0*z0*z3 + f13*y3*z0*z3 + f12*y1*z1*z3 - f12*y3*z1*z3 - f12*y2*z2*z3 + f13*y2*z2*z3 + f12*y3*z2*z3 - f13*y3*z2*z3 + f14*(y1*z1*(z0 - z3) + y3*(-z0 + z1)*z3 + y0*z0*(-z1 + z3)))/det;
        coeff_num[3] =   ((-f13)*y0*y2*z0 + f13*y0*y3*z0 + f12*y1*y2*z1 - f12*y1*y3*z1 + f13*y0*y2*z2 - f12*y1*y2*z2 + f12*y2*y3*z2 - f13*y2*y3*z2 + f15*(y1*y2*(-z1 + z2) + y0*((-y1)*z0 + y2*z0 + y1*z1 - y2*z2)) - f13*y0*y3*z3 + f12*y1*y3*z3 - f12*y2*y3*z3 + f13*y2*y3*z3 + f14*(y1*y3*(z1 - z3) + y0*(y1*z0 - y3*z0 - y1*z1 + y3*z3)))/det;
        coeff_num[9] =   (f13*y2*z0 - f13*y3*z0 - f12*y2*z1 + f12*y3*z1 - f13*y0*z2 + f12*y1*z2 - f12*y3*z2 + f13*y3*z2 + f15*(y1*z0 - y2*z0 - y0*z1 + y2*z1 + y0*z2 - y1*z2) + f13*y0*z3 - f12*y1*z3 + f12*y2*z3 - f13*y2*z3 + f14*(y3*z0 + y0*z1 - y3*z1 - y0*z3 + y1*(-z0 + z3)))/det;
    }  
	double _crdR[3]={_crd[0],_crd[1],_crd[2]};
	intoRandomAxes(_crdR);
	double res = fcb(coeff_num,_crdR), minU = nodes[0]->u[0], maxU = nodes[0]->u[0];
	for (int i=1; i<4; i++)
	{
		if (nodes[i]->u[0] < minU) minU = nodes[i]->u[0];
		if (nodes[i]->u[0] > maxU) maxU = nodes[i]->u[0];
	}
	if (res < minU) res = minU;
	if (res > maxU) res = maxU;
	
	if (fabs(res) > fabs(10000*(maxU+10))) 
	{
		printf("%lf %lf %lf %lf   %lf      det=%lf\n", fcb(coeff_num,c[0],c[4],c[8]), fcb(coeff_num,c[1],c[5],c[9]), fcb(coeff_num,c[2],c[6],c[10]), fcb(coeff_num,c[3],c[7],c[11]),/* nodes[0]->u[0], nodes[1]->u[0], nodes[2]->u[0], nodes[3]->u[0], nodes[0]->u[1], nodes[1]->u[1], nodes[2]->u[1], nodes[3]->u[1], nodes[0]->u[2], nodes[1]->u[2], nodes[2]->u[2], nodes[3]->u[2], nodes[0]->u[3], nodes[1]->u[3], nodes[2]->u[3], nodes[3]->u[3],/*coeff_num[0], coeff_num[1], coeff_num[2], coeff_num[3], coeff_num[4], coeff_num[5], coeff_num[6], coeff_num[7], coeff_num[8], coeff_num[9],*/ res,det);
		//printf("%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n\n", coeff_num[0], coeff_num[1], coeff_num[2], coeff_num[3], coeff_num[4], coeff_num[5], coeff_num[6], coeff_num[7], coeff_num[8], coeff_num[9]);
	//printf("%30lf %30lf %30lf\n%30lf %30lf %30lf  \n", coeff_num[0], coeff_num[2], coeff_num[3], coeff_num[5], coeff_num[6], coeff_num[9]);
	printf("%lf %lf %lf %lf \n", coeff_num[1], coeff_num[4], coeff_num[7], coeff_num[8]);
	//printf("%lf %lf %lf %lf %lf \n", zn3y,zn2y,zn1y,zn,znyz);
	printf("%lf %lf %lf %lf \n%lf %lf %lf %lf \n%lf %lf %lf %lf \n", c[0],c[1],c[2],c[3],c[4],c[5],c[6],c[7],c[8],c[9],c[10],c[11]);
	printf("%lf %lf %lf %lf \n%lf %lf %lf %lf \n%lf %lf %lf %lf \n\n--------------\n", nodes[0]->coords[0], nodes[1]->coords[0], nodes[2]->coords[0], nodes[3]->coords[0], nodes[0]->coords[1], nodes[1]->coords[1], nodes[2]->coords[1], nodes[3]->coords[1], nodes[0]->coords[2], nodes[1]->coords[2], nodes[2]->coords[2], nodes[3]->coords[2]);
		//printf("");
		//res = maxU;
	}
	return res;
}


void Method::count(Mesh* mesh, Node* node, double timeStep)
{
	randomizeAxis();
	int v_n = 4;//number of value components 		sizeof(node->values.name)/sizeof(double);
	Node* next = new Node(node); 	//we store new time step node in here
	double nextValues[7]={0.0}; 	//(!)4=v_n new time step values
	double a_r_coeff[3]={0.0}; 	//coeffs in xi-eta-teta coords (axis[])
	double coord_char[3]={0.0}; 	//coordinates of point in old time step, where characteristic falls
	Thetraeder* t = 0;		//thetr for interpolation
	//for (int v_c=0; v_c<v_n; v_c++) 	//for each value component
	//{
	calculateCoeff(a_r_coeff); 	//transform coefficients into random axes basis
	for (int ax=0; ax<3; ax++)	//axes splitting, temporarily additive
	{
		for (int i_crd=0; i_crd<3; i_crd++) //finding where characteristic for ax-th axis falls
			coord_char[i_crd] = node->coords[i_crd] - a_r_coeff[ax]*axis[3*ax+i_crd]*timeStep;
		t = mesh->findThetr(coord_char);
		if (!t) {printf("Fail! No thetr found for %lf %lf %lf\n",coord_char[0],coord_char[1],coord_char[2]); return;};

		nextValues[0] += interpolate_2_order(t, coord_char, mesh)/3.0; 
		nextValues[1] += interpolate_1_order(t, coord_char, 1, mesh)/3.0; 
		nextValues[2] += interpolate_1_order(t, coord_char, 2, mesh)/3.0; 
		nextValues[3] += interpolate_1_order(t, coord_char, 3, mesh)/3.0; 
		nextValues[4] += interpolate_1_order(t, coord_char, 4, mesh)/3.0; 
		nextValues[5] += interpolate_1_order(t, coord_char, 5, mesh)/3.0; 
		nextValues[6] += interpolate_1_order(t, coord_char, 6, mesh)/3.0; 
	}
	next->setValues(nextValues);	//copy new time step values into new time step node
	//memcpy(next->coords, node->coords, sizeof(double)*3); //new time step coords, no mesh movement while
	node->nextStep = next;		//add link from old node to the new one
}
