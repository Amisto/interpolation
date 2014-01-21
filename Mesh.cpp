#include "Mesh.h"
#include "Thetraeder.h"
#include "Node.h"
#include "Method.h"
#include "stdlib.h"
#include "stdio.h"
#include <cstring>
#include <iostream>
#include <fstream>
#include <math.h>

using std::string;

Mesh::Mesh()
{
	nodes = 0;
	thtrs = 0;
	nn = nt = 0;
}

Mesh::~Mesh()
{
	if (nodes) 
	{
		for (int i=0; i<nn; i++)
			delete nodes[i];
		free(nodes);
	}
	if (thtrs) 
	{
		for (int i=0; i<nt; i++)
			delete thtrs[i];
		free(thtrs);
	}
}

int Mesh::init()
{
	/*FILE* f=fopen("mesh.node","r");
	if (!f)
	{
		printf("FNF\n");
		return 1;
	}
	int _t[4],_n;
	double _c[3];
	int _ta[4];
	fscanf(f,"%d",&nn);
	printf("%d nodes:\n",nn);
	nodes = (Node**)malloc(sizeof(Node*)*nn);
	for (int i=0; i<nn; i++)
	{
		fscanf(f,"%d%lf%lf%lf",&_n,_c,_c+1,_c+2);
		if (_n!=i) { printf("Node file damaged\n"); return 1;};
		nodes[i] = new Node(_c);
		printf(" node %d: %lf %lf %lf\n",i,nodes[i]->coords[0],nodes[i]->coords[1],nodes[i]->coords[2]);
	}
	fscanf(f,"%d",&nt);
	printf("%d thtrs:\n",nt);
	thtrs = (Thetraeder**)malloc(sizeof(Thetraeder*)*nt);
	for (int i=0; i<nt; i++)
	{
		fscanf(f,"%d%d%d%d%d",&_n,_t,_t+1,_t+2,_t+3);
		for (int j=0; j<4; j++)
		{		
			_ta[j]=_t[j];
		}
		if (_n!=i) { printf("Thtr file damaged\n"); return 1;}; 
		thtrs[i] = new Thetraeder(_ta);
		printf(" thtr %d: %d %d %d %d\n",i,_t[0],_t[1],_t[2],_t[3]);
		for (int j=0; j<4; j++)
		{	
			nodes[_t[j]]->thetrs[nodes[_t[j]]->thetrsNum++]=thtrs[i];
		}
	}*/
	load_msh_file("mesh0.045.msh");
	h=0.02;
	return 0;
} //TODO: add filepath

//int Mesh::getNodesNum()
//{
//	return 0;
//}

Node* Mesh::getNode(int num)
{
	return nodes[num];
}

Thetraeder* Mesh::getThetr(int num)
{
	return thtrs[num];
}

void Mesh::transcend()
{
	Node* tmp;
	for (int i=0; i<nn; i++)
	{
		tmp = nodes[i]->nextStep;
		if (tmp == 0) printf("NULL Node called\n");
		delete nodes[i];
		nodes[i] = tmp;
	}
}

Thetraeder* Mesh::findThetr(double* _crd)
{
	//border correction
	while (_crd[0] < borders[0])
		_crd[0] += borders[3] - borders[0];
	while (_crd[1] < borders[1])
		_crd[1] += borders[4] - borders[1];
	while (_crd[2] < borders[2])
		_crd[2] += borders[5] - borders[2];
	while (_crd[0] > borders[3])
		_crd[0] -= borders[3] - borders[0];
	while (_crd[1] > borders[4])
		_crd[1] -= borders[4] - borders[1];
	while (_crd[2] > borders[5])
		_crd[2] -= borders[5] - borders[2];

	for (int i=0; i<nt; i++)
		if (thtrs[i])
			if (thtrs[i]->check(_crd,this))
			{
				//printf("!!!Thetr %i found for %lf %lf %lf------------------------\n",i,_crd[0],_crd[1],_crd[2]);
				return thtrs[i];
			}
	return 0;
}


void Mesh::addNode(Node* newNode)
{
	//if (nn==2000) return 0;
	//printf("Adding node %5d\t",nn);
	nodes[nn] = newNode;
	nn++;
}

void Mesh::addThetr(Thetraeder* newThetr)
{
	//printf("Adding thetr %5d\t",nt);
	thtrs[nt] = newThetr;
	nt++;
}

int Mesh::load_msh_file(char* file_name)
{
	int fileVer;
	string str;
	int tmp_int;
	float tmp_float;
	int number_of_nodes;
	int number_of_elements;
	Node *new_node;

	std::ifstream infile;
	infile.open(file_name, std::ifstream::in);
	if(!infile.is_open())	printf("Can not open msh file");
	printf("Reading file...");

	infile >> str;
	if(strcmp(str.c_str(),"$MeshFormat") != 0)
		printf("Wrong file format");

	infile >> tmp_float >> tmp_int >> tmp_int;
	fileVer = (int)(tmp_float*10);

	infile >> str;
	if(strcmp(str.c_str(),"$EndMeshFormat") != 0)
		printf("Wrong file format");

	printf("INFO: Header Ok");
	infile >> str;
	if(strcmp(str.c_str(),"$Nodes") != 0)
		printf("Wrong file format");

	infile >> number_of_nodes;
	printf("the file contains %d nodes\n", number_of_nodes);
	nodes = (Node**)malloc(sizeof(Node*)*number_of_nodes);
	nn = 0;
	for(int i = 0; i < number_of_nodes; i++)
	{
		new_node = new Node();
		// Zero all values
		new_node->coords[0] = new_node->coords[1] = new_node->coords[2] = 0;
		new_node->u[0] = new_node->u[1] = new_node->u[2] = 0;
		new_node->u[3] = 0;

		infile >> new_node->local_num;
		if(new_node->local_num > 0)
		{
			new_node->local_num--;
			infile >> new_node->coords[0] >> new_node->coords[1] >> new_node->coords[2];
			if (!i) 
			{
				borders[0] = new_node->coords[0];
				borders[1] = new_node->coords[1];
				borders[2] = new_node->coords[2];
				borders[3] = new_node->coords[0];
				borders[4] = new_node->coords[1];
				borders[5] = new_node->coords[2];
			}
			else
			{
				if (borders[0] > new_node->coords[0])
					borders[0] = new_node->coords[0];
				if (borders[1] > new_node->coords[1])
					borders[1] = new_node->coords[1];
				if (borders[2] > new_node->coords[2])
					borders[2] = new_node->coords[2];
				if (borders[3] < new_node->coords[0])
					borders[3] = new_node->coords[0];
				if (borders[4] < new_node->coords[1])
					borders[4] = new_node->coords[1];
				if (borders[5] < new_node->coords[2])
					borders[5] = new_node->coords[2];
			}
		}
		else
			printf("Wrong file format");
		addNode(new_node);
	}
	printf("Finished reading nodes");
	
	infile >> str;
	if(strcmp(str.c_str(),"$EndNodes") != 0)
		printf("Wrong file format");

	printf("INFO: Nodes Ok");
	
	infile >> str;
	if(strcmp(str.c_str(),"$Elements") != 0)
		printf("Wrong file format");

	infile >> number_of_elements;
	thtrs = (Thetraeder**)malloc(sizeof(Thetraeder)*number_of_elements);
	nt = 0;
	printf("the file contains %d elements", number_of_elements);
	for(int i = 0; i < number_of_elements; i++)
	{
		int new_tetr_vert[4]={0}, local_num;
		infile >> tmp_int >> tmp_int;
		if(tmp_int != 4) {
			getline(infile, str);
			continue;
		} else if (tmp_int == 4) {
			local_num = nt;
			if( fileVer == 22 ) {
				infile >> tmp_int >> tmp_int >> tmp_int 
					>> new_tetr_vert[0] >> new_tetr_vert[1] >> new_tetr_vert[2] >> new_tetr_vert[3];
			} else {
				infile >> tmp_int >> tmp_int >> tmp_int >> tmp_int 
					>> new_tetr_vert[0] >> new_tetr_vert[1] >> new_tetr_vert[2] >> new_tetr_vert[3];
			}

			if( (new_tetr_vert[0] <= 0) || (new_tetr_vert[1] <= 0) || (new_tetr_vert[2] <= 0) || (new_tetr_vert[3] <= 0) )
				printf("Wrong file format");

			new_tetr_vert[0]--; new_tetr_vert[1]--; new_tetr_vert[2]--; new_tetr_vert[3]--;
			Thetraeder *new_tetr = new Thetraeder(new_tetr_vert,this);
			addThetr(new_tetr);
		}
	}
	printf("Finished reading elements");

	printf("INFO: Elements Ok");

	infile >> str;
	if(strcmp(str.c_str(),"$EndElements") != 0)
		printf("Wrong file format");

	printf("File successfully read.");

	infile.close();
	printf("\nMesh borders:\n%lf %lf %lf  %lf %lf %lf\n", borders[0], borders[1], borders[2], borders[3], borders[4], borders[5]);
	return 0;
};


void Mesh::setInitialConditionsStep(double a, double w)
{
	printf ("Step initial conditios:\n%lf %lf\n",a,w);
	double half[3] = {(borders[3]+borders[0])/2.0,(borders[4]+borders[1])/2.0,(borders[5]+borders[2])/2.0};
	for (int i=0; i<nn; i++)
	{
		if (	nodes[i]->coords[0] < half[0] + w && nodes[i]->coords[0] > half[0] - w &&
			nodes[i]->coords[1] < half[1] + w && nodes[i]->coords[1] > half[1] - w &&
			nodes[i]->coords[2] < half[2] + w && nodes[i]->coords[2] > half[2] - w)
		nodes[i]->v = a;
	}
	setInitialConditionsGradient();
}


void Mesh::setInitialConditionsStepX(double a, double w)
{
	printf ("Step initial conditios:\n%lf %lf\n",a,w);
	double half[3] = {(borders[3]+borders[0])/2.0,(borders[4]+borders[1])/2.0,(borders[5]+borders[2])/2.0};
	for (int i=0; i<nn; i++)
	{
		if (	nodes[i]->coords[0] < half[0] + w && nodes[i]->coords[0] > half[0] - w 
			//nodes[i]->coords[1] < half[1] + w && nodes[i]->coords[1] > half[1] - w &&
			//nodes[i]->coords[2] < half[2] + w && nodes[i]->coords[2] > half[2] - w)
		   )
		nodes[i]->v = a;
	}
	setInitialConditionsGradient();
}


void Mesh::setInitialConditionsStepY(double a, double w)
{
	printf ("Step initial conditios:\n%lf %lf\n",a,w);
	double half[3] = {(borders[3]+borders[0])/2.0,(borders[4]+borders[1])/2.0,(borders[5]+borders[2])/2.0};
	for (int i=0; i<nn; i++)
	{
		if (	//nodes[i]->coords[0] < half[0] + w && nodes[i]->coords[0] > half[0] - w &&
			nodes[i]->coords[1] < half[1] + w && nodes[i]->coords[1] > half[1] - w 
			//nodes[i]->coords[2] < half[2] + w && nodes[i]->coords[2] > half[2] - w)
		   )
		nodes[i]->v = a;
	}
	setInitialConditionsGradient();
}


void Mesh::setInitialConditionsStepZ(double a, double w)
{
	printf ("Step initial conditios:\n%lf %lf\n",a,w);
	double half[3] = {(borders[3]+borders[0])/2.0,(borders[4]+borders[1])/2.0,(borders[5]+borders[2])/2.0};
	for (int i=0; i<nn; i++)
	{
		if (	//nodes[i]->coords[0] < half[0] + w && nodes[i]->coords[0] > half[0] - w &&
			//nodes[i]->coords[1] < half[1] + w && nodes[i]->coords[1] > half[1] - w &&
			nodes[i]->coords[2] < half[2] + w && nodes[i]->coords[2] > half[2] - w)
		   
		nodes[i]->v = a;
	}
	setInitialConditionsGradient();
}


void Mesh::setInitialConditionsLinearZ(double a)
{
	printf ("Linear Z initial conditions: %lf\n",a);
	double half = (borders[5]-borders[2])/2.0;
	for (int i=0; i<nn; i++)
	{
		nodes[i]->v = a - a*(fabs(half-nodes[i]->coords[2]))/half;//fabs(half-nodes[i]->coords[0]) + fabs(half-nodes[i]->coords[1]) + 
	}
	setInitialConditionsGradient();
}


void Mesh::setInitialConditionsLinearY(double a)
{
	printf ("Linear Y initial conditios: %lf\n",a);
	double half = (borders[4]-borders[1])/2.0;
	for (int i=0; i<nn; i++)
	{
		nodes[i]->v = a - a*(fabs(half-nodes[i]->coords[1]))/half;//fabs(half-nodes[i]->coords[0]) + fabs(half-nodes[i]->coords[1]) + 
	}
	setInitialConditionsGradient();
}


void Mesh::setInitialConditionsLinearX(double a)
{
	printf ("Linear X initial conditios: %lf\n",a);
	double half = (borders[3]-borders[0])/2.0;
	for (int i=0; i<nn; i++)
	{
		nodes[i]->v = a - a*(fabs(half-nodes[i]->coords[0]))/half;//fabs(half-nodes[i]->coords[0]) + fabs(half-nodes[i]->coords[1]) + 
	}
	setInitialConditionsGradient();
}


void Mesh::setInitialConditionsLinear(double a)
{
	printf ("Linear initial conditios: %lf\n",a);
	double half = (borders[3]-borders[0])/2.0;
	for (int i=0; i<nn; i++)
	{
		nodes[i]->v = a - a*(fabs(half-nodes[i]->coords[2] + fabs(half-nodes[i]->coords[0]) + fabs(half-nodes[i]->coords[1])))/half;//
	}
	setInitialConditionsGradient();
}


void Mesh::setInitialConditionsGradient()
{
	double crd[3]={0.0},gP=0.0,gM=0.0;
	Thetraeder *tP=0, *tM=0;
	for (int i=0; i<nn; i++)
	{
		//h = ((nodes[i]->coords[0]-nodes[i]->coords[0])+(nodes[i]->coords[1]-nodes[i]->coords[1])+(nodes[i]->coords[2]-nodes[i]->coords[2]))/3.0;
		for (int ax=0; ax<3; ax++)
		{
			for (int c=0;c<3;c++) crd[c] = nodes[i]->coords[c];
			if (crd[ax]==borders[ax] || crd[ax]==borders[3+ax]) 
			{
				nodes[i]->u[1+ax] = 0.0;
				continue;
			}	
			crd[ax] += h/4.0;
			tP = findThetr(crd);
			gP = Method::interpolate_1_order(tP,crd,0,this);
			crd[ax] -= h/2.0;
			tM=findThetr(crd);
			gM = Method::interpolate_1_order(tM,crd,0,this);
			nodes[i]->u[1+ax] = 4.0*(gP-nodes[i]->u[0])/h;
		}
	}
	setInitialConditionsGradientSecond();
}

void Mesh::setInitialConditionsGradientSecond()
{
	double crd[3]={0.0},gP=0.0,gM=0.0;
	Thetraeder *tP=0, *tM=0;
	for (int i=0; i<nn; i++)
	{
		//h = ((nodes[i]->coords[0]-nodes[i]->coords[0])+(nodes[i]->coords[1]-nodes[i]->coords[1])+(nodes[i]->coords[2]-nodes[i]->coords[2]))/3.0;
		for (int ax=0; ax<3; ax++)
		{
			for (int c=0;c<3;c++) crd[c] = nodes[i]->coords[c];
			if (crd[ax]==borders[ax] || crd[ax]==borders[3+ax]) 
			{
				nodes[i]->u[4+ax] = 0.0;
				continue;
			}	
			crd[ax] += h/4.0;
			tP = findThetr(crd);
			gP = Method::interpolate_1_order(tP,crd,ax,this);
			crd[ax] -= h/2.0;
			tM=findThetr(crd);
			gM = Method::interpolate_1_order(tM,crd,ax,this);
			nodes[i]->u[4+ax] = 4.0*(gP-nodes[i]->u[ax])/h;
		}
	}
}
