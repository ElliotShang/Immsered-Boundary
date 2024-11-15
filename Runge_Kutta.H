#ifndef _Runge_Kutta_H_
#define _Runge_Kutta_H_

#include <iostream>
#include "math.h"
#include <vector>
#include <cmath>

using namespace std;

void forthRK(double & y, double & ydot,
             double & mass, double & damp, double & stif,
             double & force, double & timestep);

void computeDerivative(double h,double & dh,double & d2h,double theta,double & dtheta,double & d2theta,
                       const double & mass, const double & intertia, const double & damp,
                       const double & rotdamp, const double & stif, const double & rotstif,
                       const double & force, const double & moment, const double & b);

void computeDerivative_Anal(double h,double & dh,double & d2h, double theta,double & dtheta,double & d2theta,
                            const double & mass, const double & intertia, const double & damp,
                            const double & rotdamp, const double & stif, const double & rotstif,
                            const double & force, const double & moment, const double & b);

void forthRK2DOF(double & y, double & ydot, double & theta, double & thetadot,
                 const double & mass, const double & intertia, const double & damp,
                 const double & rotdamp, const double & stif, const double & rotstif,
                 const double & force, const double & moment, const double & b, double & timestep);

void forthRK2DOF_Anal(double & y, double & ydot, double & theta, double & thetadot,
                      const double & mass, const double & intertia, const double & damp,
                      const double & rotdamp, const double & stif, const double & rotstif,
                      const double & force, const double & moment, const double & b, double & timestep, double & theta0);

#endif