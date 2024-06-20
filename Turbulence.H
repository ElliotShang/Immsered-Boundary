#ifndef _Turbulence_H_
#define _Turbulence_H_

#include "Box.H"
#include "Point.H"
#include "Constants.H"
#include "FlowVariables.H"

extern const double cb1;
extern const double segma;
extern const double cb2;
extern const double kappa;
extern const double cw2;
extern const double cw3;
extern const double cv1;
extern const double ct3;
extern const double ct4;
extern const double cw1;
extern const double cv1three;
extern const double far_vis_ratio;
extern const double max_viseddy;
extern const double min_viseddy;
extern const double min_viseddy_wall;

double Viseddy2mut(double & cellvis, double & cellroe, double & cellmu);

void ImageViseddy(IBCell & a_ibcell, Pointxyz & patv, Pointxyz & nmv, BoxtoWall & ibboxtopatch);

void ImageViseddy_Domain(FlowVariables & refv, FlowVariables & bcvar, Pointxyz & patv, Pointxyz & nmv, 
	double & refdis, double & targetdis);

void uplus_as_function_of_yplus(double & uplus, double & yplus);

void yplus_as_function_of_re(double & hg_yplus, double & hg_re);

void Patchut(Surfpatch & mypatch, Pointxyz & patv, Pointxyz & nmv, Pointxyz & vtdir);

void VelocityDecom(Pointxyz & inputvel, Pointxyz & nmv,
				   Pointxyz & tangvel, double & tangspeed, Pointxyz & tangdir,
				   Pointxyz & nmvel, double & nmspeed);

void Cellut(FlowVariables & cellfv, double & celldis_to_wall, Pointxyz & patv, Pointxyz & nmv, double & patut, Pointxyz & vtdir);				   

#endif