#ifndef THETRAEDER
#define THETRAEDER 1
class Node;
class Mesh;
class Thetraeder
{
public:
	Thetraeder();
	Thetraeder(int* _v, Mesh* mesh);
	~Thetraeder();
	int check(double* _crd, Mesh* mesh);
	int vert[4];
	int local_num;
	int checkZborder(Mesh* mesh);
private:
	double aabb[6];
};

#endif
