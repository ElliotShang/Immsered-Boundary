#ifndef _Face_H_
#define _Face_H_

#include "Point.H"
/*keisa is an unit normal vector of the face*/
/*which is different from the keisa of the box geom*/
/*the keisa in box geom is partial derivative of the coordinates*/

struct Face
{
static Point unitnv[3];
	   int bi[2];
	   int fnv;
	   double area;
	   Pointxyz keisa;

	   Pointxyz keisa1;
	   Pointxyz keisa2;
	   Pointxyz fc1;
	   Pointxyz fc2;

	   Face()
	   {
	   	bi[0] = -1;
	   	bi[1] = -1;
	   	fnv = -1;
	   }

	   inline int & operator[](const int & i0)
	   {
	   	return bi[i0];
	   }

	   inline bool hastwobox()
	   {
	   	if (bi[0] > -1 && bi[1] > -1)
	   	{
	   		return true;
	   	}
	   	else return false;
	   }
	   inline bool haszerobox()
	   {
	   	if (bi[0] == -1 && bi[1] == -1)
	   	{
	   		return true;
	   	}
	   	else return false;
	   }
};

#endif