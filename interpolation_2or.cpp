#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<limits.h>
#include<math.h>
#define N 10
#define T 100000000	

struct vec3 
{
       double x,y,z;
};

void out(double* _a)
{
    putchar('\n'); 
    for (int i=0; i<N; i++)
        printf("%f ",_a[i]);
    putchar('\n');putchar('\n');
}
void out(vec3 v)
{
	printf("%10f %10f %10f\n",v.x,v.y,v.z);
}

/*double coeff_dif4(double* _a1, double * _a2)
{
    double res;
    for (int i=0; i<N; i++)
        printf("%f ",_a[i]);
    putchar('\n');putchar('\n');
}*/

double rand_1()
{
       return ((double)(rand()%5000000))/5000000.0;
}

double func(double* a, vec3 c)
{
       return a[4]*c.x*c.x + a[5]*c.y*c.y + a[6]*c.z*c.z + a[7]*c.x*c.y + a[8]*c.x*c.z + a[9]*c.y*c.z + 
                           a[1]*c.x + a[2]*c.y + a[3]*c.z + a[0];
}

double grad_x(double* a, vec3 c)
{
       return 2.0*a[4]*c.x + a[7]*c.y + a[8]*c.z + a[1];
}

double grad_y(double* a, vec3 c)
{
       return 2.0*a[5]*c.y + a[7]*c.x + a[9]*c.z + a[2];
}

double grad_z(double* a, vec3 c)
{
       return 2.0*a[6]*c.z + a[8]*c.x + a[9]*c.y + a[3];
}

int main()
{
    printf("Interpolation polynom:\na1*x^2+a2*y^2+a3*z^2+a4*x*y*+a5*x*z+a6*y*z+a7*x+a8*y+a9*z+a10\n");
    double dif4_max=0.0, dif4_average, d;
    double coeff_exact[N];                                  
    srand(time(NULL));                                       

    for (int all=0; all<T; all++)
    {                                          
				//exact polynom coefficients         //                                    
    for (int i=0; i<N; i++)                                          //
        coeff_exact[i] = rand_1()*10;                                //
        
    vec3 coord[4];                                                   //vertices
    double f[4], g_x[4], g_y[4], g_z[4];                             //exact function & gradient in vertices
    for (int i=0; i<4; i++)                                          //
    {                                                                //
        coord[i].x = rand_1();                                   //
        coord[i].y = rand_1();                                   //
        coord[i].z = rand_1();                                   //
        f[i]=func(coeff_exact,coord[i]);                             //
        g_x[i]=grad_x(coeff_exact,coord[i]);                         //
        g_y[i]=grad_y(coeff_exact,coord[i]);                         //
        g_z[i]=grad_z(coeff_exact,coord[i]);                         //
    }                                                                //
    
    
    double coeff_num[N];                                             //numerical coefficients
		//for (int i=0; i<4; i++) out(coord[i]);putchar('\n');
    int min=0,max=0;					//(y0-y3) appears in the following formulae,trying to 
    							//make the dif4erence as far from zero as possible
    vec3 temp; 
    double t;
    for (int i=1; i<4; i++)
	if (coord[i].y < coord[min].y) min=i;
    temp=coord[min]; coord[min]=coord[0]; coord[0]=temp;
    t=f[0]; f[0]=f[min]; f[min]=t;
    t=g_x[0]; g_x[0]=g_x[min]; g_x[min]=t;
    t=g_y[0]; g_y[0]=g_y[min]; g_y[min]=t;
    t=g_z[0]; g_z[0]=g_z[min]; g_z[min]=t;
    for (int i=1; i<4; i++)
	if (coord[i].y > coord[max].y) max=i;
    temp=coord[max]; coord[max]=coord[3]; coord[3]=temp;
    t=f[3]; f[3]=f[max]; f[max]=t;
    t=g_x[3]; g_x[3]=g_x[max]; g_x[max]=t;
    t=g_y[3]; g_y[3]=g_y[max]; g_y[max]=t;
    t=g_z[3]; g_z[3]=g_z[max]; g_z[max]=t;
		//for (int i=0; i<4; i++) out(coord[i]);putchar('\n');putchar('\n');
    double   	x0=coord[0].x, x1=coord[1].x, x2=coord[2].x, x3=coord[3].x,   //
       		y0=coord[0].y, y1=coord[1].y, y2=coord[2].y, y3=coord[3].y,   //
       		z0=coord[0].z, z1=coord[1].z, z2=coord[2].z, z3=coord[3].z,   //
       r0,r1,r2,r3;
    if (y3-y0 < 0.0001) {printf("%10f <- skipping\n",y3-y0); continue;};
    {                                                                
        double f0=g_x[0],f1=g_x[1],f2=g_x[2],f3=g_x[3];              //solve for grad_x: 4x4
        r0 = (f3*x2*y1*z0 - f2*x3*y1*z0 - f3*x1*y2*z0 + f1*x3*y2*z0 + f2*x1*y3*z0 - 
                                    f1*x2*y3*z0 - f3*x2*y0*z1 + f2*x3*y0*z1 + f3*x0*y2*z1 - f0*x3*y2*z1 - 
                                    f2*x0*y3*z1 + f0*x2*y3*z1 + f3*x1*y0*z2 - f1*x3*y0*z2 - f3*x0*y1*z2 + 
                                    f0*x3*y1*z2 + f1*x0*y3*z2 - f0*x1*y3*z2 - f2*x1*y0*z3 + 
                                    f1*x2*y0*z3 + f2*x0*y1*z3 - f0*x2*y1*z3 - f1*x0*y2*z3 + 
                                    f0*x1*y2*z3)/(x2*y1*z0 - x3*y1*z0 - x1*y2*z0 + x3*y2*z0 + x1*y3*z0 - x2*y3*z0 - 
                                    x2*y0*z1 + x3*y0*z1 + x0*y2*z1 - x3*y2*z1 - x0*y3*z1 + x2*y3*z1 + x1*y0*z2 - 
                                    x3*y0*z2 - x0*y1*z2 + x3*y1*z2 + x0*y3*z2 - x1*y3*z2 - 
                                    x1*y0*z3 + x2*y0*z3 + x0*y1*z3 - x2*y1*z3 - x0*y2*z3 + x1*y2*z3);
        r1 = -((f2*y1*z0 - f3*y1*z0 - f1*y2*z0 + f3*y2*z0 + f1*y3*z0 - f2*y3*z0 - f2*y0*z1 + 
                                    f3*y0*z1 + f0*y2*z1 - f3*y2*z1 - f0*y3*z1 + f2*y3*z1 + 
                                    f1*y0*z2 - f3*y0*z2 - f0*y1*z2 + f3*y1*z2 + f0*y3*z2 - 
                                    f1*y3*z2 - f1*y0*z3 + f2*y0*z3 + f0*y1*z3 - f2*y1*z3 - f0*y2*z3 + f1*y2*z3)/
                                    (2.0*((-x2)*y1*z0 + x3*y1*z0 + x1*y2*z0 - x3*y2*z0 - x1*y3*z0 + 
                                    x2*y3*z0 + x2*y0*z1 - x3*y0*z1 - x0*y2*z1 + x3*y2*z1 + x0*y3*z1 - x2*y3*z1 - 
                                    x1*y0*z2 + x3*y0*z2 + x0*y1*z2 - x3*y1*z2 - x0*y3*z2 + x1*y3*z2 + 
                                    x1*y0*z3 - x2*y0*z3 - x0*y1*z3 + x2*y1*z3 + x0*y2*z3 - x1*y2*z3)));
        r2 = (f2*x1*z0 - f3*x1*z0 - f1*x2*z0 + f3*x2*z0 + f1*x3*z0 - f2*x3*z0 - f2*x0*z1 + 
                                    f3*x0*z1 + f0*x2*z1 - f3*x2*z1 - f0*x3*z1 + f2*x3*z1 + 
                                    f1*x0*z2 - f3*x0*z2 - f0*x1*z2 + f3*x1*z2 + f0*x3*z2 - f1*x3*z2 - 
                                    f1*x0*z3 + f2*x0*z3 + f0*x1*z3 - f2*x1*z3 - f0*x2*z3 + f1*x2*z3)/
                                    ((-x2)*y1*z0 + x3*y1*z0 + x1*y2*z0 - x3*y2*z0 - x1*y3*z0 + x2*y3*z0 + 
                                    x2*y0*z1 - x3*y0*z1 - x0*y2*z1 + x3*y2*z1 + x0*y3*z1 - x2*y3*z1 - 
                                    x1*y0*z2 + x3*y0*z2 + x0*y1*z2 - x3*y1*z2 - x0*y3*z2 + x1*y3*z2 + x1*y0*z3 - 
                                    x2*y0*z3 - x0*y1*z3 + x2*y1*z3 + x0*y2*z3 - x1*y2*z3);
        r3 = ((-f2)*x1*y0 + f3*x1*y0 + f1*x2*y0 - f3*x2*y0 - f1*x3*y0 + f2*x3*y0 + f2*x0*y1 - f3*x0*y1 - 
                                    f0*x2*y1 + f3*x2*y1 + f0*x3*y1 - f2*x3*y1 - 
                                    f1*x0*y2 + f3*x0*y2 + f0*x1*y2 - f3*x1*y2 - f0*x3*y2 + f1*x3*y2 + f1*x0*y3 - 
                                    f2*x0*y3 - f0*x1*y3 + f2*x1*y3 + f0*x2*y3 - f1*x2*y3)/
                                    ((-x2)*y1*z0 + x3*y1*z0 + x1*y2*z0 - x3*y2*z0 - x1*y3*z0 + x2*y3*z0 + x2*y0*z1 - 
                                    x3*y0*z1 - x0*y2*z1 + x3*y2*z1 + x0*y3*z1 - x2*y3*z1 - 
                                    x1*y0*z2 + x3*y0*z2 + x0*y1*z2 - x3*y1*z2 - x0*y3*z2 + x1*y3*z2 + x1*y0*z3 - 
                                    x2*y0*z3 - x0*y1*z3 + x2*y1*z3 + x0*y2*z3 - x1*y2*z3);
    }  
    {                                                               
        double f0=f[0]-x0*r0 - x0*x0*r1 - x0*y0*r2 - x0*z0*r3,    //solve for the rest: 6x6
               f1=f[1]-x1*r0 - x1*x1*r1 - x1*y1*r2 - x1*z1*r3,
               f2=f[2]-x2*r0 - x2*x2*r1 - x2*y2*r2 - x2*z2*r3,
               f3=f[3]-x3*r0 - x3*x3*r1 - x3*y3*r2 - x3*z3*r3,          
               f4=g_y[0]-x0*r2,f5=g_z[0]-x0*r3;
	double dy03=(y0 - y3)*(y0 - y3),dz03=(z0 - z3)*(z0 - z3),
		dy02=(y0 - y2)*(y0 - y2),dz02=(z0 - z2)*(z0 - z2),
		dy01=(y0 - y1)*(y0 - y1),dz01=(z0 - z1)*(z0 - z1),
		//f40=f0-f3-f4*y0+f4*y3-f5*z0+f5*z3,
		//f41=f1-f0+f4*y0-f4*y1+f5*z0-f5*z1,
		//f42=(dy03*dz01-dy01*dz03)*(dy03*(-f0+f2+f4*y0-f4*y2+f5*z0-f5*z2)+dy02*f40)
		//		-(dy03*dz02-dy02*dz03)*(dy03*f41+dy01*f40),
		zn3y=y2*( z0 - z1) + y0*(z1 - z2) + y1*(-z0 + z2),
		zn2y=y3*(-z0 + z1) + y1*(z0 - z3) + y0*(-z1 + z3),
		zn1y=y3*(-z0 + z2) + y2*(z0 - z3) + y0*(-z2 + z3),
		zn  =y3*(z0 - z1) + y1*(z0 - z3) + y0*(-2.0*z0 + z1 + z3),
		znyz=dy03*dz01 - dy01*dz03,
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
        
    } 
    
    vec3 point,average={0,0,0};//,max=coord[0];                      // point inside the thetraedr
    for (int i=1; i<4; i++)                                          //
    {                                                                //
        //if (coord[i].x < min.x) min.x = coord[i].x;                  //
        //if (coord[i].x > max.x) max.x = coord[i].x;                  //
        //if (coord[i].y < min.y) min.y = coord[i].y;                  //
        //if (coord[i].y > max.y) max.y = coord[i].y;                  //
        //if (coord[i].x < min.z) min.x = coord[i].z;                  //
        //if (coord[i].x > max.z) max.x = coord[i].z;                  //
	average.x += coord[i].x;
	average.y += coord[i].y;
	average.z += coord[i].z;
    }                                                                //
    average.x /= 4.0;
    average.y /= 4.0;
    average.z /= 4.0;
	
    //point.x = (max.x+min.x)/2;//min.x + (max.x-min.x)*rand_1();                        //
    //point.y = (max.y+min.y)/2;//min.y + (max.y-min.y)*rand_1();                        //
    //point.z = (max.z+min.z)/2;//min.z + (max.z-min.z)*rand_1();                        //
    point = average;
    
    double f_p_exact,f_p_numer;
        f_p_exact = func(coeff_exact,point);                      //exact function in point
        f_p_numer = func(coeff_num  ,point);                      //numer function in point
        d=fabs((f_p_exact-f_p_numer)/f_p_exact);
        if (d > dif4_max) dif4_max = d;
	dif4_average += d;
	if (d > 1.000001)
	{
		printf("%6d::  E: %10f  N: %10f  D: %10f\n",all,f_p_exact,f_p_numer,d);
		for (int i=0; i<4; i++) out(coord[i]);
		out(coeff_exact); out(coeff_num);
	}
	if (!(all%100000)) printf("%10d|  Max: %10f  Average: %10f\n\n",all,dif4_max,dif4_average/T);
	//printf("%lf %lf %lf %lf\n",fabs(func(coeff_exact,coord[0])-func(coeff_num,coord[0])),fabs(func(coeff_exact,coord[1])-func(coeff_num,coord[1])),fabs(func(coeff_exact,coord[2])-func(coeff_num,coord[2])),fabs(func(coeff_exact,coord[3])-func(coeff_num,coord[3])));
    }
    ////double c_dif4=0.0;
    //for (int i=0; i<N; i++)
    //    c_dif4 += fabs(coeff_exact[i]-coeff_num[i])/fabs(coeff_exact[i]);
    dif4_average /= T;
    printf("\nMax: %10f  Average: %10f\n\n",dif4_max,dif4_average);
}
