#include <iostream>
#include "Point.H"
#include "FlowVariables.H"

#define segnum 3

struct wallface
{
	int ci;
	int dir; 
	int side;
	Pointxyz flowdir;
	Pointxyz segcenter[segnum][segnum];
	double segdistowall[segnum][segnum];
	FlowVec newflux;
};