#ifndef METHOD
#define METHOD 1

class Node;
class Thetraeder;
class Mesh;
class Method
{
public:
	Method();
	~Method();
	int init();
	void count(Mesh* mesh, Node* n, double timeStep);//, int axis);

	static double scalar(double* _v1, double* _v2);
	static double interpolate_1_order(Thetraeder* t, double* _crd, int val, Mesh* mesh); 
	static double interpolate_2_order(Thetraeder* t, double* _crd, Mesh* mesh); 
	static double interpolate_3_order(Thetraeder* t, double* _crd, Mesh* mesh); 
private:
	static double axis[9]; 	// current local randomised axis
	double coeff[3]; 	// dv/dt + coeff[0]*dv/dx + coeff[1]*dv/vy + coeff[2]*dv/vz = 0
			 	// coeffs in global coordinates
	void randomizeAxis();
	static void randomizeAxis(double* axes);
	void calculateCoeff(double* _c); //transform coordinates in accordance with axis[]

	static void intoRandomAxes(double* x, double* y, double *z);
	static void intoRandomAxes(double* c);
	static void intoRandomAxes(double* x, double* y, double *z, double *axes);
	static void intoRandomAxes(double* c, double *axes);
	static void fromRandomAxes(double* c, double *axes);
};
#endif
