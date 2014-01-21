#include "Node.h"

Node::Node()
{
	nextStep = 0;
	thetrsNum = 0;
}

Node::Node(double* _c)
{
	for (int i=0; i<3; i++)
	{
		coords[i]=_c[i];
	}
	nextStep = 0;
	thetrsNum = 0;
}

Node::Node(Node* _n)
{
	for (int i=0; i<3; i++)
		coords[i]=_n->coords[i];
	setValues(_n->u);
	local_num = _n->local_num;
	nextStep = 0;
	thetrsNum = _n->thetrsNum;	
	for (int i=0; i<thetrsNum; i++)
		thetrs[i]=_n->thetrs[i];
}

Node::~Node()
{
	//do NOT delete nextStep!
}

void Node::setValues(double* _v)
{
	if (!_v) return; 
	for (int i=0;i<4;i++) 
		u[i]=_v[i];
}
