#ifndef _GhostPair_H_
#define _GhostPair_H_

#include "Boxloc.H"

struct CtoFPair
{
	int ci;
	Boxson<int> fi;

	CtoFPair()
	{}

	CtoFPair(const int & i0)
	{
		ci = i0;
	}
};

struct Domainghost
{
	int cell;
	int refcell;
	int nbface;
	Pointxyz nmv;
	bool nmvflag;

	Domainghost()
	{}

	Domainghost(const int & abox, const int & aface, const int & rc)
	{
		cell = abox;
		nbface = aface;
		refcell = rc;
		nmvflag = false;
	}
};

struct NodePair
{
	int innode;
	Boxloc outnode;
#ifdef PASSAGE_ANGLE	
	double theta;
#endif
	NodePair()
	{}

	NodePair(const int & i0, const int & j0, const int & k0)
	{		
		innode = i0;
		outnode.node = j0;
		outnode.index = k0;
#ifdef PASSAGE_ANGLE		
		theta = 0.0;
#endif		
	}

	bool insamenode()
	{
		if (node == outnode.node) return true;
		else return false;
	}
};

struct Facepair
{
	int cface;
	int fface[2][2];

	FlowVec cvflux;
	FlowVec cinvflux;
	FlowVec fvflux[2][2];
	FlowVec finvflux[2][2];

	Facepair()
	{}

	Facepair(const int & ci0)
	{
		cface = ci0;
	}
};

#endif