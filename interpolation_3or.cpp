#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<limits.h>
#include<math.h>
#define N 20
#define T 100000000	

struct vec3 
{
       double x,y,z;
};

void out(double* _a)
{
    putchar('\n'); 
    for (int i=0; i<N; i++)
        printf("%20lf \n",_a[i]);
    putchar('\n');putchar('\n');putchar('\n');
}
void out(vec3 v)
{
	printf("%20f %20f %20f\n",v.x,v.y,v.z);
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
       return a[0] + a[1]*c.x + a[2]*c.y + a[3]*c.z + a[4]*c.x*c.x + a[5]*c.y*c.y + a[6]*c.z*c.z + a[7]*c.x*c.y + a[8]*c.x*c.z + a[9]*c.y*c.z + a[10]*c.x*c.x*c.x + a[11]*c.y*c.y*c.y + a[12]*c.z*c.z*c.z + a[13]*c.x*c.x*c.y + a[14]*c.x*c.x*c.z + a[15]*c.y*c.y*c.x + a[16]*c.y*c.y*c.z + a[17]*c.z*c.z*c.x + a[18]*c.z*c.z*c.y + a[19]*c.x*c.y*c.z;
                        
}

double grad_x(double* a, vec3 c)
{
       return a[1] + 2.0*a[4]*c.x + 3.0*a[10]*c.x*c.x + a[7]*c.y + 2.0*a[13]*c.x*c.y + a[15]*c.y*c.y + a[8]*c.z + 2.0*a[14]*c.x*c.z + a[19]*c.y*c.z + a[17]*c.z*c.z;
}

double grad_xx(double* a, vec3 c)
{
	return 2.0*a[4]+6.0*a[10]*c.x+2.0*a[13]*c.y+2.0*a[14]*c.z;
}

double grad_yy(double* a, vec3 c)
{
       return 2.0*a[5]+2.0*a[15]*c.x+6.0*a[11]*c.y+2.0*a[16]*c.z;//2.0*a[5]*c.y + a[7]*c.x + a[9]*c.z + a[2];
}

double grad_zz(double* a, vec3 c)
{
       return 2.0*a[6]+2.0*a[17]*c.x+2.0*a[18]*c.y+6.0*a[12]*c.z;//2.0*a[6]*c.z + a[8]*c.x + a[9]*c.y + a[3];
}

int main()
{
    printf("Interpolation polynom:\n3rd order\n");
    double dif4_max=0.0, dif4_average, d;
    double coeff_exact[N];                                  
    srand(time(NULL));                                       

    for (int all=0; all<T; all++)
    {                                          
				//exact polynom coefficients         //                                    
    for (int i=0; i<N; i++)                                          //
        coeff_exact[i] = rand_1()*10;                                //
        
    vec3 coord[4];                                                   //vertices
    double f[4], g_x[4], g_xx[4], g_yy[4], g_zz[4];                             //exact function & gradient & 2nd gradient in vertices
    for (int i=0; i<4; i++)                                          //
    {                                                                //
        coord[i].x = rand_1();                                   //
        coord[i].y = rand_1();                                   //
        coord[i].z = rand_1();                                   //
        f[i]=func(coeff_exact,coord[i]);                             //
        g_x[i]=grad_x(coeff_exact,coord[i]);                         //
        g_xx[i]=grad_xx(coeff_exact,coord[i]);                         //
        g_yy[i]=grad_yy(coeff_exact,coord[i]);                         //
        g_zz[i]=grad_zz(coeff_exact,coord[i]);                         //
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
    //t=f[0]; f[0]=f[min]; f[min]=t;
    //t=g_x[0]; g_x[0]=g_x[min]; g_x[min]=t;
    //t=g_y[0]; g_y[0]=g_y[min]; g_y[min]=t;
    //t=g_z[0]; g_z[0]=g_z[min]; g_z[min]=t;
    for (int i=1; i<4; i++)
	if (coord[i].y > coord[max].y) max=i;
    temp=coord[max]; coord[max]=coord[3]; coord[3]=temp;
    //t=f[3]; f[3]=f[max]; f[max]=t;
    //t=g_x[3]; g_x[3]=g_x[max]; g_x[max]=t;
    //t=g_y[3]; g_y[3]=g_y[max]; g_y[max]=t;
    //t=g_z[3]; g_z[3]=g_z[max]; g_z[max]=t;
    for (int i=0; i<4; i++)                                          //
    {                                                                //
        f[i]=func(coeff_exact,coord[i]);                             //
        g_x[i]=grad_x(coeff_exact,coord[i]);                         //
        g_xx[i]=grad_xx(coeff_exact,coord[i]);                         //
        g_yy[i]=grad_yy(coeff_exact,coord[i]);                         //
        g_zz[i]=grad_zz(coeff_exact,coord[i]);                         //
    }                                                                //
		//for (int i=0; i<4; i++) out(coord[i]);putchar('\n');putchar('\n');
    double   	x0=coord[0].x, x1=coord[1].x, x2=coord[2].x, x3=coord[3].x,   //
       		y0=coord[0].y, y1=coord[1].y, y2=coord[2].y, y3=coord[3].y,   //
       		z0=coord[0].z, z1=coord[1].z, z2=coord[2].z, z3=coord[3].z,   //
       r0,r1,r2,r3;
	double volume = x2*y1*z0 - x3*y1*z0 - x1*y2*z0 + x3*y2*z0 + x1*y3*z0 - x2*y3*z0 - x2*y0*z1 + x3*y0*z1 + x0*y2*z1 - x3*y2*z1 - x0*y3*z1 + x2*y3*z1 + x1*y0*z2 - x3*y0*z2 - x0*y1*z2 + x3*y1*z2 + x0*y3*z2 - x1*y3*z2 -x1*y0*z3 + x2*y0*z3 + x0*y1*z3 - x2*y1*z3 - x0*y2*z3 + x1*y2*z3,
		det=(y1*(y2*(z1 - z2)*(z0 - z3) - y3*(z0 - z2)*(z1 - z3)) + y0*(y3*(z1 - z2)*(z0 - z3) - y2*(z0 - z2)*(z1 - z3) + y1*(z0 - z1)*(z2 - z3)) + y2*y3*(z0 - z1)*(z2 - z3));
    if (fabs(volume) < 0.00000001) {printf("%10f <- skipping by volume\n",volume*10000000); continue;};
    if (fabs(det) < 0.00000001) {printf("%10f <- skipping by det\n",det*10000000); continue;};
    //if (y3-y0 < 0.0001) {printf("%10f <- skipping\n",y3-y0); continue;};
    {                                                                
        double f0=g_xx[0],f1=g_xx[1],f2=g_xx[2],f3=g_xx[3];              //solve for grad_xx: 4x4
        coeff_num[4]  = ((-f1)*x3*y2*z0 + f1*x2*y3*z0 + f0*x3*y2*z1 - f0*x2*y3*z1 + f1*x3*y0*z2 - f0*x3*y1*z2 - f1*x0*y3*z2 + f0*x1*y3*z2 + f3*((-x2)*y1*z0 + x1*y2*z0 + x2*y0*z1 - x0*y2*z1 - x1*y0*z2 + x0*y1*z2) - f1*x2*y0*z3 + f0*x2*y1*z3 + f1*x0*y2*z3 - f0*x1*y2*z3 + f2*(x3*y1*z0 - x1*y3*z0 - x3*y0*z1 + x0*y3*z1 + x1*y0*z3 - x0*y1*z3))/(-2*volume);
        coeff_num[10] =   (f1*y2*z0 - f1*y3*z0 - f0*y2*z1 + f0*y3*z1 - f1*y0*z2 + f0*y1*z2 - f0*y3*z2 + f1*y3*z2 + f3*(y1*z0 - y2*z0 - y0*z1 + y2*z1 + y0*z2 - y1*z2) + f1*y0*z3 - f0*y1*z3 + f0*y2*z3 - f1*y2*z3 + f2*(y3*z0 + y0*z1 - y3*z1 - y0*z3 + y1*(-z0 + z3)))/(-6*volume);
        coeff_num[13] =   ((-f1)*x2*z0 + f1*x3*z0 + f0*x2*z1 - f0*x3*z1 + f1*x0*z2 - f0*x1*z2 + f0*x3*z2 - f1*x3*z2 + f3*((-x1)*z0 + x2*z0 + x0*z1 - x2*z1 - x0*z2 + x1*z2) - f1*x0*z3 + f0*x1*z3 - f0*x2*z3 + f1*x2*z3 + f2*((-x3)*z0 - x0*z1 + x3*z1 + x1*(z0 - z3) + x0*z3))/(-2*volume);
        coeff_num[14] = (f1*x2*y0 - f1*x3*y0 - f0*x2*y1 + f0*x3*y1 - f1*x0*y2 + f0*x1*y2 - f0*x3*y2 + f1*x3*y2 + f3*(x1*y0 - x2*y0 - x0*y1 + x2*y1 + x0*y2 - x1*y2) + f1*x0*y3 - f0*x1*y3 + f0*x2*y3 - f1*x2*y3 + f2*(x3*y0 + x0*y1 - x3*y1 - x0*y3 + x1*(-y0 + y3)))/(-2*volume);
    }  
    {                                                      
        double f4=g_yy[0],f5=g_yy[1],f6=g_yy[2],f7=g_yy[3];              //solve for grad_yy: 4x4
        coeff_num[5]  = ((-f5)*x3*y2*z0 + f5*x2*y3*z0 + f4*x3*y2*z1 - f4*x2*y3*z1 + f5*x3*y0*z2 - f4*x3*y1*z2 - f5*x0*y3*z2 + f4*x1*y3*z2 + f7*((-x2)*y1*z0 + x1*y2*z0 + x2*y0*z1 - x0*y2*z1 - x1*y0*z2 + x0*y1*z2) - f5*x2*y0*z3 + f4*x2*y1*z3 + f5*x0*y2*z3 - f4*x1*y2*z3 + f6*(x3*y1*z0 - x1*y3*z0 - x3*y0*z1 + x0*y3*z1 + x1*y0*z3 - x0*y1*z3))/(-2*volume);
        coeff_num[11] =   ((-f5)*x2*z0 + f5*x3*z0 + f4*x2*z1 - f4*x3*z1 + f5*x0*z2 - f4*x1*z2 + f4*x3*z2 - f5*x3*z2 + f7*((-x1)*z0 + x2*z0 + x0*z1 - x2*z1 - x0*z2 + x1*z2) - f5*x0*z3 + f4*x1*z3 - f4*x2*z3 + f5*x2*z3 + f6*((-x3)*z0 - x0*z1 + x3*z1 + x1*(z0 - z3) + x0*z3))/(-6*volume);
        coeff_num[15] =   (f5*y2*z0 - f5*y3*z0 - f4*y2*z1 + f4*y3*z1 - f5*y0*z2 + f4*y1*z2 - f4*y3*z2 + f5*y3*z2 + f7*(y1*z0 - y2*z0 - y0*z1 + y2*z1 + y0*z2 - y1*z2) + f5*y0*z3 - f4*y1*z3 + f4*y2*z3 - f5*y2*z3 + f6*(y3*z0 + y0*z1 - y3*z1 - y0*z3 + y1*(-z0 + z3)))/(-2*volume);
        coeff_num[16] =     (f5*x2*y0 - f5*x3*y0 - f4*x2*y1 + f4*x3*y1 - f5*x0*y2 + f4*x1*y2 - f4*x3*y2 + f5*x3*y2 + f7*(x1*y0 - x2*y0 - x0*y1 + x2*y1 + x0*y2 - x1*y2) + f5*x0*y3 - f4*x1*y3 + f4*x2*y3 - f5*x2*y3 + f6*(x3*y0 + x0*y1 - x3*y1 - x0*y3 + x1*(-y0 + y3)))/(-2*volume);
    } 
    {                                                      
        double f8=g_zz[0],f9=g_zz[1],f10=g_zz[2],f11=g_zz[3];              //solve for grad_zz: 4x4
        coeff_num[6]  = ((-f9)*x3*y2*z0 + f9*x2*y3*z0 + f8*x3*y2*z1 - f8*x2*y3*z1 + f9*x3*y0*z2 - f8*x3*y1*z2 - f9*x0*y3*z2 + f8*x1*y3*z2 + f11*((-x2)*y1*z0 + x1*y2*z0 + x2*y0*z1 - x0*y2*z1 - x1*y0*z2 + x0*y1*z2) - f9*x2*y0*z3 + f8*x2*y1*z3 + f9*x0*y2*z3 - f8*x1*y2*z3 + f10*(x3*y1*z0 - x1*y3*z0 - x3*y0*z1 + x0*y3*z1 + x1*y0*z3 - x0*y1*z3))/(-2*volume);
        coeff_num[12] =   (f9*x2*y0 - f9*x3*y0 - f8*x2*y1 + f8*x3*y1 - f9*x0*y2 + f8*x1*y2 - f8*x3*y2 + f9*x3*y2 + f11*(x1*y0 - x2*y0 - x0*y1 + x2*y1 + x0*y2 - x1*y2) + f9*x0*y3 - f8*x1*y3 + f8*x2*y3 - f9*x2*y3 + f10*(x3*y0 + x0*y1 - x3*y1 - x0*y3 + x1*(-y0 + y3)))/(-6*volume);
        coeff_num[17] =   (f9*y2*z0 - f9*y3*z0 - f8*y2*z1 + f8*y3*z1 - f9*y0*z2 + f8*y1*z2 - f8*y3*z2 + f9*y3*z2 + f11*(y1*z0 - y2*z0 - y0*z1 + y2*z1 + y0*z2 - y1*z2) + f9*y0*z3 - f8*y1*z3 + f8*y2*z3 - f9*y2*z3 + f10*(y3*z0 + y0*z1 - y3*z1 - y0*z3 + y1*(-z0 + z3)))/(-2*volume);
        coeff_num[18] =     ((-f9)*x2*z0 + f9*x3*z0 + f8*x2*z1 - f8*x3*z1 + f9*x0*z2 - f8*x1*z2 + f8*x3*z2 - f9*x3*z2 + f11*((-x1)*z0 + x2*z0 + x0*z1 - x2*z1 - x0*z2 + x1*z2) - f9*x0*z3 + f8*x1*z3 - f8*x2*z3 + f9*x2*z3 + f10*((-x3)*z0 - x0*z1 + x3*z1 + x1*(z0 - z3) + x0*z3))/(-2*volume);
    }     
    {                                                                
        double 	f12=g_x[0]-2.0*x0*coeff_num[4]-3.0*x0*x0*coeff_num[10]-2.0*x0*y0*coeff_num[13]-2.0*x0*z0*coeff_num[14]-y0*y0*coeff_num[15]-z0*z0*coeff_num[17],
		f13=g_x[1]-2.0*x1*coeff_num[4]-3.0*x1*x1*coeff_num[10]-2.0*x1*y1*coeff_num[13]-2.0*x1*z1*coeff_num[14]-y1*y1*coeff_num[15]-z1*z1*coeff_num[17],
		f14=g_x[2]-2.0*x2*coeff_num[4]-3.0*x2*x2*coeff_num[10]-2.0*x2*y2*coeff_num[13]-2.0*x2*z2*coeff_num[14]-y2*y2*coeff_num[15]-z2*z2*coeff_num[17],
		f15=g_x[3]-2.0*x3*coeff_num[4]-3.0*x3*x3*coeff_num[10]-2.0*x3*y3*coeff_num[13]-2.0*x3*z3*coeff_num[14]-y3*y3*coeff_num[15]-z3*z3*coeff_num[17];
								            //solve for grad_x: 4x4
        coeff_num[1]  = ((-f13)*y0*y3*z0*z2 + f13*y2*y3*z0*z2 + f12*y1*y3*z1*z2 - f12*y2*y3*z1*z2 + f15*(y1*y2*z0*(z1 - z2) + y0*(y1*(z0 - z1)*z2 + y2*z1*(-z0 + z2))) + f13*y0*y2*z0*z3 - f13*y2*y3*z0*z3 - f12*y1*y2*z1*z3 + f12*y2*y3*z1*z3 - f13*y0*y2*z2*z3 + f12*y1*y2*z2*z3 + f13*y0*y3*z2*z3 - f12*y1*y3*z2*z3 + f14*(y1*y3*z0*(-z1 + z3) + y0*(y3*z1*(z0 - z3) + y1*(-z0 + z1)*z3)))/det;
        coeff_num[7] =     (f13*y0*z0*z2 - f13*y2*z0*z2 - f12*y1*z1*z2 + f12*y2*z1*z2 + f15*(y0*z0*(z1 - z2) + y2*(z0 - z1)*z2 + y1*z1*(-z0 + z2)) - f13*y0*z0*z3 + f13*y3*z0*z3 + f12*y1*z1*z3 - f12*y3*z1*z3 - f12*y2*z2*z3 + f13*y2*z2*z3 + f12*y3*z2*z3 - f13*y3*z2*z3 + f14*(y1*z1*(z0 - z3) + y3*(-z0 + z1)*z3 + y0*z0*(-z1 + z3)))/det;
        coeff_num[8] =   ((-f13)*y0*y2*z0 + f13*y0*y3*z0 + f12*y1*y2*z1 - f12*y1*y3*z1 + f13*y0*y2*z2 - f12*y1*y2*z2 + f12*y2*y3*z2 - f13*y2*y3*z2 + f15*(y1*y2*(-z1 + z2) + y0*((-y1)*z0 + y2*z0 + y1*z1 - y2*z2)) - f13*y0*y3*z3 + f12*y1*y3*z3 - f12*y2*y3*z3 + f13*y2*y3*z3 + f14*(y1*y3*(z1 - z3) + y0*(y1*z0 - y3*z0 - y1*z1 + y3*z3)))/det;
        coeff_num[19] =   (f13*y2*z0 - f13*y3*z0 - f12*y2*z1 + f12*y3*z1 - f13*y0*z2 + f12*y1*z2 - f12*y3*z2 + f13*y3*z2 + f15*(y1*z0 - y2*z0 - y0*z1 + y2*z1 + y0*z2 - y1*z2) + f13*y0*z3 - f12*y1*z3 + f12*y2*z3 - f13*y2*z3 + f14*(y3*z0 + y0*z1 - y3*z1 - y0*z3 + y1*(-z0 + z3)))/det;
    }  
    {                       
	double *a=coeff_num;                                         
        double 	f12=f[0]-x0*a[1]-x0*x0*a[4]-y0*y0*a[5]-z0*z0*a[6]-x0*y0*a[7]-x0*z0*a[8]-x0*x0*x0*a[10]-y0*y0*y0*a[11]-z0*z0*z0*a[12]-x0*x0*y0*a[13]-x0*x0*z0*a[14]-y0*y0*x0*a[15]-y0*y0*z0*a[16]-z0*z0*x0*a[17]-z0*z0*y0*a[18]-x0*y0*z0*a[19],
		f13=f[1]-x1*a[1]-x1*x1*a[4]-y1*y1*a[5]-z1*z1*a[6]-x1*y1*a[7]-x1*z1*a[8]-x1*x1*x1*a[10]-y1*y1*y1*a[11]-z1*z1*z1*a[12]-x1*x1*y1*a[13]-x1*x1*z1*a[14]-y1*y1*x1*a[15]-y1*y1*z1*a[16]-z1*z1*x1*a[17]-z1*z1*y1*a[18]-x1*y1*z1*a[19],
		f14=f[2]-x2*a[1]-x2*x2*a[4]-y2*y2*a[5]-z2*z2*a[6]-x2*y2*a[7]-x2*z2*a[8]-x2*x2*x2*a[10]-y2*y2*y2*a[11]-z2*z2*z2*a[12]-x2*x2*y2*a[13]-x2*x2*z2*a[14]-y2*y2*x2*a[15]-y2*y2*z2*a[16]-z2*z2*x2*a[17]-z2*z2*y2*a[18]-x2*y2*z2*a[19],
		f15=f[3]-x3*a[1]-x3*x3*a[4]-y3*y3*a[5]-z3*z3*a[6]-x3*y3*a[7]-x3*z3*a[8]-x3*x3*x3*a[10]-y3*y3*y3*a[11]-z3*z3*z3*a[12]-x3*x3*y3*a[13]-x3*x3*z3*a[14]-y3*y3*x3*a[15]-y3*y3*z3*a[16]-z3*z3*x3*a[17]-z3*z3*y3*a[18]-x3*y3*z3*a[19];
								              //solve for rest: 4x4
        coeff_num[0]  = ((-f13)*y0*y3*z0*z2 + f13*y2*y3*z0*z2 + f12*y1*y3*z1*z2 - f12*y2*y3*z1*z2 + f15*(y1*y2*z0*(z1 - z2) + y0*(y1*(z0 - z1)*z2 + y2*z1*(-z0 + z2))) + f13*y0*y2*z0*z3 - f13*y2*y3*z0*z3 - f12*y1*y2*z1*z3 + f12*y2*y3*z1*z3 - f13*y0*y2*z2*z3 + f12*y1*y2*z2*z3 + f13*y0*y3*z2*z3 - f12*y1*y3*z2*z3 + f14*(y1*y3*z0*(-z1 + z3) + y0*(y3*z1*(z0 - z3) + y1*(-z0 + z1)*z3)))/det;
        coeff_num[2] =   (f13*y0*z0*z2 - f13*y2*z0*z2 - f12*y1*z1*z2 + f12*y2*z1*z2 + f15*(y0*z0*(z1 - z2) + y2*(z0 - z1)*z2 + y1*z1*(-z0 + z2)) - f13*y0*z0*z3 + f13*y3*z0*z3 + f12*y1*z1*z3 - f12*y3*z1*z3 - f12*y2*z2*z3 + f13*y2*z2*z3 + f12*y3*z2*z3 - f13*y3*z2*z3 + f14*(y1*z1*(z0 - z3) + y3*(-z0 + z1)*z3 + y0*z0*(-z1 + z3)))/det;
        coeff_num[3] =   ((-f13)*y0*y2*z0 + f13*y0*y3*z0 + f12*y1*y2*z1 - f12*y1*y3*z1 + f13*y0*y2*z2 - f12*y1*y2*z2 + f12*y2*y3*z2 - f13*y2*y3*z2 + f15*(y1*y2*(-z1 + z2) + y0*((-y1)*z0 + y2*z0 + y1*z1 - y2*z2)) - f13*y0*y3*z3 + f12*y1*y3*z3 - f12*y2*y3*z3 + f13*y2*y3*z3 + f14*(y1*y3*(z1 - z3) + y0*(y1*z0 - y3*z0 - y1*z1 + y3*z3)))/det;
        coeff_num[9] =   (f13*y2*z0 - f13*y3*z0 - f12*y2*z1 + f12*y3*z1 - f13*y0*z2 + f12*y1*z2 - f12*y3*z2 + f13*y3*z2 + f15*(y1*z0 - y2*z0 - y0*z1 + y2*z1 + y0*z2 - y1*z2) + f13*y0*z3 - f12*y1*z3 + f12*y2*z3 - f13*y2*z3 + f14*(y3*z0 + y0*z1 - y3*z1 - y0*z3 + y1*(-z0 + z3)))/det;
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
		//out(coeff_exact); out(coeff_num);
	    	for (int i=0; i<N; i++)
	        	printf("%5d   %10lf %10lf\n",i,coeff_exact[i],coeff_num[i]);
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
