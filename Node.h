#ifndef NODE
#define NODE 1

class Thetraeder;
class Node
{
public:
	Node();
	Node(double* _c);
	Node(Node* _n);
	~Node();
	/*union crd
	{
		double c[3];
		struct nm
		{
			double x;
			double y;
			double z;
		} name; 
	}*/
	double coords[3];
	union 
	{
		double u[7];
		struct 
		{
			double v;
			double vx;
			double vy;
			double vz;
			double vxx;
			double vyy;
			double vzz;
		};
	};
	int local_num;
	void setValues(double* _v);
	Node* nextStep;
	Thetraeder* thetrs[30];//tethraeders
	int thetrsNum;//tethraeders
};
#endif
