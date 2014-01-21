#ifndef MESH
#define MESH 1

class Node;
class Thetraeder;
class Mesh
{
public:
	Mesh();
	~Mesh();
	int init(); //TODO: add filepath
	int getNodesNum() {return nn;};
	int getThetrsNum() {return nt;};
	Node* getNode(int num);
	Thetraeder* getThetr(int num);
	void transcend();
	Thetraeder* findThetr(double* _crd);
	int load_msh_file(char* file_name);
	void addNode(Node* newNode);
	void addThetr(Thetraeder* newThetr);
	void setInitialConditionsLinearX(double a);
	void setInitialConditionsLinearY(double a);
	void setInitialConditionsLinearZ(double a);
	void setInitialConditionsLinear(double a);
	void setInitialConditionsStep(double a, double w);
	void setInitialConditionsStepX(double a, double w);
	void setInitialConditionsStepY(double a, double w);
	void setInitialConditionsStepZ(double a, double w);
	void setInitialConditionsGradient();
	void setInitialConditionsGradientSecond();
	double getMinZ() {return borders[2];};
	double getMaxZ() {return borders[5];};
private:
	int nn,nt;
	Node** nodes;
	Thetraeder** thtrs;
	double borders[6];//0,1,2 - min values for x,y,z; 3,4,5 - max values
				//!!!Warning, borders only cyclic and for cubic mesh 
	double h; //space step
	//node list
	//tetr list
};
#endif
