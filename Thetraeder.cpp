#include "Thetraeder.h"
#include "Mesh.h"
#include "Node.h"
#include <math.h>
#include <stdio.h>


Thetraeder::Thetraeder()
{
}

Thetraeder::Thetraeder(int* _v, Mesh* mesh)
{
	for (int i=0; i<4; i++)
		vert[i]=_v[i];
	Node* nodes[4] = {mesh->getNode(vert[0]), mesh->getNode(vert[1]), mesh->getNode(vert[2]), mesh->getNode(vert[3])};
	aabb[0] = aabb[3] = nodes[0]->coords[0];
	aabb[1] = aabb[4] = nodes[0]->coords[1];
	aabb[2] = aabb[5] = nodes[0]->coords[2];
	for (int i=1; i<4; i++)
	{
		if (aabb[0] > nodes[i]->coords[0])
			aabb[0] = nodes[i]->coords[0];
		if (aabb[1] > nodes[i]->coords[1])
			aabb[1] = nodes[i]->coords[1];
		if (aabb[2] > nodes[i]->coords[2])
			aabb[2] = nodes[i]->coords[2];
		if (aabb[3] < nodes[i]->coords[0])
			aabb[3] = nodes[i]->coords[0];
		if (aabb[4] < nodes[i]->coords[1])
			aabb[4] = nodes[i]->coords[1];
		if (aabb[5] < nodes[i]->coords[2])
			aabb[5] = nodes[i]->coords[2];
	}
}

Thetraeder::~Thetraeder()
{
}

int Thetraeder::check(double* _crd, Mesh* mesh)
{
	if (_crd[0] < aabb[0] || _crd[1] < aabb[1] || _crd[2] < aabb[2] || _crd[0] > aabb[3] || _crd[1] > aabb[4] || _crd[2] > aabb[5]) return 0; //fast cut of far thetrs

	Node* nodes[4] = {mesh->getNode(vert[0]), mesh->getNode(vert[1]), mesh->getNode(vert[2]), mesh->getNode(vert[3])};
	double axis[9];
	//check for a plane against node 0
	axis[0] = nodes[1]->coords[0] - nodes[3]->coords[0];
	axis[1] = nodes[1]->coords[1] - nodes[3]->coords[1];
	axis[2] = nodes[1]->coords[2] - nodes[3]->coords[2];
	axis[3] = nodes[2]->coords[0] - nodes[3]->coords[0];
	axis[4] = nodes[2]->coords[1] - nodes[3]->coords[1];
	axis[5] = nodes[2]->coords[2] - nodes[3]->coords[2];
	double znam=0.0, alfa=0.0, beta=0.0;
	if (znam = axis[0]*axis[4] - axis[3]*axis[1])
	{
		alfa = (axis[5]*axis[1] - axis[2]*axis[4])/znam;
		beta = (axis[3]*axis[2] - axis[0]*axis[5])/znam;
		axis[8] = 1.0/sqrt(alfa*alfa + beta*beta+ 1.0);
		axis[6] = alfa*axis[8]; axis[7] = beta*axis[8];
	}
	else if (znam = axis[0]*axis[5] - axis[3]*axis[2])
	{
		alfa = (axis[4]*axis[2] - axis[1]*axis[5])/znam;
		beta = (axis[3]*axis[1] - axis[0]*axis[4])/znam;
		axis[7] = 1.0/sqrt(alfa*alfa + beta*beta+ 1.0);
		axis[6] = alfa*axis[7]; axis[8] = beta*axis[7];
	}
	else if (znam = axis[2]*axis[4] - axis[5]*axis[1])
	{
		alfa = (axis[3]*axis[1] - axis[0]*axis[4])/znam;
		beta = (axis[5]*axis[0] - axis[2]*axis[3])/znam;
		axis[6] = 1.0/sqrt(alfa*alfa + beta*beta+ 1.0);
		axis[8] = alfa*axis[6]; axis[7] = beta*axis[6];
	}
	else 
	{
		printf("Fail in determining axis direction");
		return 0;
	}
	//double 	c1 = axis[0]*axis[6]+axis[1]*axis[7]+axis[2]*axis[8] ,
	//	c2 = axis[3]*axis[6]+axis[4]*axis[7]+axis[5]*axis[8];
	//if (c1*c1 > 0.0000000001 || c2*c2 > 0.000000001)
	//printf("CHECK %lf %lf\n", c1,c2);
	if (   ((_crd[0]-nodes[3]->coords[0])*axis[6]+
		(_crd[1]-nodes[3]->coords[1])*axis[7]+
		(_crd[2]-nodes[3]->coords[2])*axis[8])*
	       ((nodes[0]->coords[0]-nodes[3]->coords[0])*axis[6]+
		(nodes[0]->coords[1]-nodes[3]->coords[1])*axis[7]+
		(nodes[0]->coords[2]-nodes[3]->coords[2])*axis[8])
									<-0.0000001)
	{
		//printf ("%lf %lf %lf node 0\n%lf %lf %lf node 1\n%lf %lf %lf node 2\n%lf %lf %lf node 3\n%lf %lf %lf crd\n",nodes[0]->coords[0],nodes[0]->coords[1],nodes[0]->coords[2],nodes[1]->coords[0],nodes[1]->coords[1],nodes[1]->coords[2],nodes[2]->coords[0],nodes[2]->coords[1],nodes[2]->coords[2],nodes[3]->coords[0],nodes[3]->coords[1],nodes[3]->coords[2],_crd[0],_crd[1],_crd[2]);

		return 0;
	}

	//check for a plane against node 1
	axis[0] = nodes[0]->coords[0] - nodes[3]->coords[0];
	axis[1] = nodes[0]->coords[1] - nodes[3]->coords[1];
	axis[2] = nodes[0]->coords[2] - nodes[3]->coords[2];
	axis[3] = nodes[2]->coords[0] - nodes[3]->coords[0];
	axis[4] = nodes[2]->coords[1] - nodes[3]->coords[1];
	axis[5] = nodes[2]->coords[2] - nodes[3]->coords[2];
	znam=0.0; alfa=0.0; beta=0.0;
	if (znam = axis[0]*axis[4] - axis[3]*axis[1])
	{
		alfa = (axis[5]*axis[1] - axis[2]*axis[4])/znam;
		beta = (axis[3]*axis[2] - axis[0]*axis[5])/znam;
		axis[8] = 1.0/sqrt(alfa*alfa + beta*beta+ 1.0);
		axis[6] = alfa*axis[8]; axis[7] = beta*axis[8];
	}
	else if (znam = axis[0]*axis[5] - axis[3]*axis[2])
	{
		alfa = (axis[4]*axis[2] - axis[1]*axis[5])/znam;
		beta = (axis[3]*axis[1] - axis[0]*axis[4])/znam;
		axis[7] = 1.0/sqrt(alfa*alfa + beta*beta+ 1.0);
		axis[6] = alfa*axis[7]; axis[8] = beta*axis[7];
	}
	else if (znam = axis[2]*axis[4] - axis[5]*axis[1])
	{
		alfa = (axis[3]*axis[1] - axis[0]*axis[4])/znam;
		beta = (axis[5]*axis[0] - axis[2]*axis[3])/znam;
		axis[6] = 1.0/sqrt(alfa*alfa + beta*beta+ 1.0);
		axis[8] = alfa*axis[6]; axis[7] = beta*axis[6];
	}
	else 
	{
		printf("Fail in determining axis direction");
		return 0;
	}
	if (   ((_crd[0]-nodes[3]->coords[0])*axis[6]+
		(_crd[1]-nodes[3]->coords[1])*axis[7]+
		(_crd[2]-nodes[3]->coords[2])*axis[8])*
	       ((nodes[1]->coords[0]-nodes[3]->coords[0])*axis[6]+
		(nodes[1]->coords[1]-nodes[3]->coords[1])*axis[7]+
		(nodes[1]->coords[2]-nodes[3]->coords[2])*axis[8])
									<-0.0000001)
		return 0;
	
	//check for a plane against node 2
	axis[0] = nodes[0]->coords[0] - nodes[3]->coords[0];
	axis[1] = nodes[0]->coords[1] - nodes[3]->coords[1];
	axis[2] = nodes[0]->coords[2] - nodes[3]->coords[2];
	axis[3] = nodes[1]->coords[0] - nodes[3]->coords[0];
	axis[4] = nodes[1]->coords[1] - nodes[3]->coords[1];
	axis[5] = nodes[1]->coords[2] - nodes[3]->coords[2];
	znam=0.0; alfa=0.0; beta=0.0;
	if (znam = axis[0]*axis[4] - axis[3]*axis[1])
	{
		alfa = (axis[5]*axis[1] - axis[2]*axis[4])/znam;
		beta = (axis[3]*axis[2] - axis[0]*axis[5])/znam;
		axis[8] = 1.0/sqrt(alfa*alfa + beta*beta+ 1.0);
		axis[6] = alfa*axis[8]; axis[7] = beta*axis[8];
	}
	else if (znam = axis[0]*axis[5] - axis[3]*axis[2])
	{
		alfa = (axis[4]*axis[2] - axis[1]*axis[5])/znam;
		beta = (axis[3]*axis[1] - axis[0]*axis[4])/znam;
		axis[7] = 1.0/sqrt(alfa*alfa + beta*beta+ 1.0);
		axis[6] = alfa*axis[7]; axis[8] = beta*axis[7];
	}
	else if (znam = axis[2]*axis[4] - axis[5]*axis[1])
	{
		alfa = (axis[3]*axis[1] - axis[0]*axis[4])/znam;
		beta = (axis[5]*axis[0] - axis[2]*axis[3])/znam;
		axis[6] = 1.0/sqrt(alfa*alfa + beta*beta+ 1.0);
		axis[8] = alfa*axis[6]; axis[7] = beta*axis[6];
	}
	else 
	{
		printf("Fail in determining axis direction");
		return 0;
	}
	if (   ((_crd[0]-nodes[3]->coords[0])*axis[6]+
		(_crd[1]-nodes[3]->coords[1])*axis[7]+
		(_crd[2]-nodes[3]->coords[2])*axis[8])*
	       ((nodes[2]->coords[0]-nodes[3]->coords[0])*axis[6]+
		(nodes[2]->coords[1]-nodes[3]->coords[1])*axis[7]+
		(nodes[2]->coords[2]-nodes[3]->coords[2])*axis[8])
									<-0.0000001)
		return 0;

	//check for a plane against node 3
	axis[0] = nodes[0]->coords[0] - nodes[1]->coords[0];
	axis[1] = nodes[0]->coords[1] - nodes[1]->coords[1];
	axis[2] = nodes[0]->coords[2] - nodes[1]->coords[2];
	axis[3] = nodes[2]->coords[0] - nodes[1]->coords[0];
	axis[4] = nodes[2]->coords[1] - nodes[1]->coords[1];
	axis[5] = nodes[2]->coords[2] - nodes[1]->coords[2];
	znam=0.0; alfa=0.0; beta=0.0;
	if (znam = axis[0]*axis[4] - axis[3]*axis[1])
	{
		alfa = (axis[5]*axis[1] - axis[2]*axis[4])/znam;
		beta = (axis[3]*axis[2] - axis[0]*axis[5])/znam;
		axis[8] = 1.0/sqrt(alfa*alfa + beta*beta+ 1.0);
		axis[6] = alfa*axis[8]; axis[7] = beta*axis[8];
	}
	else if (znam = axis[0]*axis[5] - axis[3]*axis[2])
	{
		alfa = (axis[4]*axis[2] - axis[1]*axis[5])/znam;
		beta = (axis[3]*axis[1] - axis[0]*axis[4])/znam;
		axis[7] = 1.0/sqrt(alfa*alfa + beta*beta+ 1.0);
		axis[6] = alfa*axis[7]; axis[8] = beta*axis[7];
	}
	else if (znam = axis[2]*axis[4] - axis[5]*axis[1])
	{
		alfa = (axis[3]*axis[1] - axis[0]*axis[4])/znam;
		beta = (axis[5]*axis[0] - axis[2]*axis[3])/znam;
		axis[6] = 1.0/sqrt(alfa*alfa + beta*beta+ 1.0);
		axis[8] = alfa*axis[6]; axis[7] = beta*axis[6];
	}
	else 
	{
		printf("Fail in determining axis direction");
		return 0;
	}
	if (   ((_crd[0]-nodes[1]->coords[0])*axis[6]+
		(_crd[1]-nodes[1]->coords[1])*axis[7]+
		(_crd[2]-nodes[1]->coords[2])*axis[8])*
	       ((nodes[3]->coords[0]-nodes[1]->coords[0])*axis[6]+
		(nodes[3]->coords[1]-nodes[1]->coords[1])*axis[7]+
		(nodes[3]->coords[2]-nodes[1]->coords[2])*axis[8])
									<-0.0000001)
		return 0;
//printf ("%lf %lf %lf node 0\n%lf %lf %lf node 1\n%lf %lf %lf node 2\n%lf %lf %lf node 3\n%lf %lf %lf crd\n",nodes[0]->coords[0],nodes[0]->coords[1],nodes[0]->coords[2],nodes[1]->coords[0],nodes[1]->coords[1],nodes[1]->coords[2],nodes[2]->coords[0],nodes[2]->coords[1],nodes[2]->coords[2],nodes[3]->coords[0],nodes[3]->coords[1],nodes[3]->coords[2],_crd[0],_crd[1],_crd[2]);
	return 1;
}


int Thetraeder::checkZborder( Mesh* mesh)
{	
	Node* nodes[4] = {mesh->getNode(vert[0]), mesh->getNode(vert[1]), mesh->getNode(vert[2]), mesh->getNode(vert[3])};
	double min=mesh->getMinZ(),
		max=mesh->getMaxZ();	
	if (nodes[0]->coords[2] == min && nodes[1]->coords[2] == min && nodes[2]->coords[2] == min) return 1;
	if (nodes[0]->coords[2] == min && nodes[1]->coords[2] == min && nodes[3]->coords[2] == min) return 1;
	if (nodes[0]->coords[2] == min && nodes[2]->coords[2] == min && nodes[3]->coords[2] == min) return 1;
	if (nodes[1]->coords[2] == min && nodes[2]->coords[2] == min && nodes[3]->coords[2] == min) return 1;
	if (nodes[0]->coords[2] == max && nodes[1]->coords[2] == max && nodes[2]->coords[2] == max) return 1;
	if (nodes[0]->coords[2] == max && nodes[1]->coords[2] == max && nodes[3]->coords[2] == max) return 1;
	if (nodes[0]->coords[2] == max && nodes[2]->coords[2] == max && nodes[3]->coords[2] == max) return 1;
	if (nodes[1]->coords[2] == max && nodes[2]->coords[2] == max && nodes[3]->coords[2] == max) return 1;
	return 0;
}
