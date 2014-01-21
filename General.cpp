#include "General.h"
#include "Mesh.h"
#include "Method.h"
#include "Node.h"
#include "VTKSnapshotWriter.h"

General::General()
{
	mesh = 0;
	method = 0;
}

General::~General()
{
	if (mesh) delete mesh;
	if (method) delete method;
}

int General::init()
{
	finalStep = 100;
	time = 0;
	timeStep = 0.001;
	snapStep = 1;			//TODO:load parameters from file
	if (mesh = new Mesh())
		{if (mesh->init()) return 1;} 		//TODO: add filepath
	else return 1;
	if (method = new Method())
		{if (method->init()) return 1;}		//TODO: when added different methods, must be chosen from xml
	else return 1;
	if (sw = new VTKSnapshotWriter())
		{if (sw->init()) return 1;}		//TODO: when added different methods, must be chosen from xml
	else return 1;

	mesh->setInitialConditionsStepY(10.0,0.15);
	return 0;
}

void General::setTimeStep()
{	
	timeStep = 1;
}

int General::step(int currentStep)
{
	if (!currentStep)
		sw->dump_vtk(mesh,0);
	setTimeStep();
	time += timeStep;
	//for (int axis=0; axis<3; axis++)
		for (int i=0; i<mesh->getNodesNum(); i++)
		{
			Node* n = mesh->getNode(i);
			method->count(mesh, n, timeStep);//, axis);
		}
	mesh->transcend();
	if (!(currentStep%snapStep))
		sw->dump_vtk(mesh,currentStep+1);
		
}
