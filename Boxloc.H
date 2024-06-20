#ifndef _Boxloc_H_
#define _Boxloc_H_

#include <cmath>
#include "Point_iterator.H"
#include "AMRmpi.H"

struct Boxloc
{	
	int node;
	int index;

	Boxloc()
	{
		node = -1;
		index = -1;
	}

	void init()
	{
		node = -1;
		index = -1;
	}

	Boxloc(const int & ilevel, const int & myindex0)
	{
		node = ilevel;
		index = myindex0;
	}
	inline int & operator[](const int & i0)
	{
		return *(&node+i0);
	}

};

template <class T>
struct Boxson
{
	T   son[2][2][2];
	Boxson()
	{}
	Boxson(const T & a0)
	{
		for (int i = 0; i < 2; ++i)
		{
			for (int j = 0; j < 2; ++j)
			{
				for (int k = 0; k < 2; ++k)
				{
					son[i][j][k] = a0;
				}
			}
		}
	}
};
#endif
