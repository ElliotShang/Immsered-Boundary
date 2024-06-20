#include "Runge_Kutta.H"
#include <vector>
#include "math.h"
#include <cmath>
#include "Constants.H"

void forthRK(double & y, double & ydot,
             double & mass, double & damp, double & stif,
             double & force, double & timestep)
{
    double yddot_1 = (force - damp*ydot - stif*y)/mass; // k21
    double ydot_1 = ydot;  // k11
    double y_1 = y;

    double ydot_2 = ydot_1 + yddot_1 * timestep / 2.0;  // k12
    double yddot_2 = (force - damp*(ydot_1+0.5*timestep*yddot_1) - stif*(y_1+0.5*timestep*ydot_1))/mass;   // k22

    double ydot_3 = ydot_1 + yddot_2 * timestep / 2.0;       // k23
    double yddot_3 = (force - damp*(ydot_1+0.5*timestep*yddot_2) - stif*(y_1+0.5*timestep*ydot_2))/mass;  // k13

    double ydot_4 = ydot_1 + yddot_3 * timestep;           //k24
    double yddot_4 = (force - damp*(ydot_1+timestep*yddot_3) - stif*(y_1+timestep*ydot_3))/mass;   // k14

    y = y_1 + timestep/6.0*(ydot_1 + 2.0*ydot_2 + 2.0*ydot_3 + ydot_4);
    ydot = ydot_1 + timestep/6.0*(yddot_1 + 2.0*yddot_2 + 2.0*yddot_3 + yddot_4);
}

void computeDerivative(double h,double & dh,double & d2h, double theta,double & dtheta,double & d2theta,
                       const double & mass, const double & intertia, const double & damp,
                       const double & rotdamp, const double & stif, const double & rotstif,
                       const double & force, const double & moment, const double & b)
{
    // Initial guess for d2h and d2theta
    d2h = (force-damp*dh-stif*h)/mass;
    d2theta = (moment-rotdamp*dtheta-rotstif*theta)/intertia;
    // fixed-point iteration
    for (int i = 0; i < 40; ++i) {
        double new_d2h = (force - damp* dh - stif * h + b * mass * cos(theta) * d2theta - b * mass * sin(theta) * pow(dtheta, 2)) / mass;
        double new_d2theta = (moment - rotdamp * dtheta - rotstif * theta + b * mass * cos(theta) * new_d2h) / intertia;
        // Convergence check (can be adjusted)
        if (fabs(new_d2h - d2h) < 1e-6 && fabs(new_d2theta - d2theta) < 1e-6) {
            break;
        }
        d2h = new_d2h;
        d2theta = new_d2theta;
    }
}

// 不动点迭代求解
void forthRK2DOF(double & y, double & ydot, double & theta, double & thetadot,
                 const double & mass, const double & intertia, const double & damp,
                 const double & rotdamp, const double & stif, const double & rotstif,
                 const double & force, const double & moment, const double & b, double & timestep)
{
    double y_1 = y;
    double ydot_1 = ydot;
    double theta_1 = theta;
    double thetadot_1 = thetadot;
    double ydotk1,yddotk1,thetadotk1,thetaddotk1;
    double ydotk2,yddotk2,thetadotk2,thetaddotk2;
    double ydotk3,yddotk3,thetadotk3,thetaddotk3;
    double ydotk4,yddotk4,thetadotk4,thetaddotk4;

    ydotk1 = ydot_1;
    thetadotk1 = thetadot_1;
    computeDerivative(y_1,ydotk1,yddotk1, theta_1,thetadotk1,thetaddotk1,
                      mass,intertia,damp,rotdamp,stif,rotstif,force,moment,b);

    ydotk2 = ydot_1+0.5*timestep*yddotk1;
    thetadotk2 = thetadot_1+0.5*timestep*thetaddotk1;
    computeDerivative(y_1+0.5*timestep*ydotk1, ydotk2, yddotk2,
                      theta_1+0.5*timestep*thetadotk1,thetadotk2,thetaddotk2,
                      mass,intertia,damp,rotdamp,stif,rotstif,force,moment,b);

    ydotk3 = ydot_1+0.5*timestep*yddotk2;
    thetadotk3 = thetadot_1+0.5*timestep*thetaddotk2;
    computeDerivative(y_1+0.5*timestep*ydotk2, ydotk3, yddotk3,
                      theta_1+0.5*timestep*thetadotk2,thetadotk3,thetaddotk3,
                      mass,intertia,damp,rotdamp,stif,rotstif,force,moment,b);

    ydotk4 = ydot_1+timestep*yddotk3;
    thetadotk4 = thetadot_1+timestep*thetaddotk3;
    computeDerivative(y_1+timestep*ydotk3, ydotk4, yddotk4,
                      theta_1+timestep*thetadotk3,thetadotk4,thetaddotk4,
                      mass,intertia,damp,rotdamp,stif,rotstif,force,moment,b);

    y = y_1 + timestep/6.0*(ydotk1 + 2.0*ydotk2 + 2.0*ydotk3 + ydotk4);
    ydot = ydot_1 + timestep/6.0*(yddotk1 + 2.0*yddotk2 + 2.0*yddotk3 + yddotk4);
    theta = theta_1 + timestep/6.0*(thetadotk1 + 2.0*thetadotk2 + 2.0*thetadotk3 + thetadotk4);
    thetadot = thetadot_1 + timestep/6.0*(thetaddotk1 + 2.0*thetaddotk2 + 2.0*thetaddotk3 + thetaddotk4);
}

void computeDerivative_Anal(double h,double & dh,double & d2h, double theta,double & dtheta,double & d2theta,
                            const double & mass, const double & intertia, const double & damp,
                            const double & rotdamp, const double & stif, const double & rotstif,
                            const double & force, const double & moment, const double & b){
    double numerator = intertia*mass-b*b*cos(theta)*cos(theta)*mass*mass;
    d2h = -1.0*(intertia*stif/numerator)*h+(-1.0)*(b*cos(theta)*mass*rotstif/numerator)*(theta)
          +(-1.0)*(damp*intertia/numerator)*dh+(-1.0)*((b*cos(theta)*rotdamp*mass/numerator)+(b*intertia*mass*sin(theta)*dtheta/numerator))*dtheta+
          (force*intertia/numerator)+(b*cos(theta)*mass*moment/numerator);

    d2theta = -1.0*(b*cos(theta)*stif*mass/numerator)*h+(-1.0)*(rotstif*mass/numerator)*(theta)
              +(-1.0)*(b*cos(theta)*damp*mass/intertia)*dh+(-1.0)*((rotdamp*mass/numerator)+(b*b*cos(theta)*mass*mass*sin(theta)*dtheta/numerator))*dtheta+
              (b*cos(theta)*mass*force/numerator)+(mass*moment/numerator);
}

// 解析求解耦合方程组
void forthRK2DOF_Anal(double & y, double & ydot, double & theta, double & thetadot,
                      const double & mass, const double & intertia, const double & damp,
                      const double & rotdamp, const double & stif, const double & rotstif,
                      const double & force, const double & moment, const double & b, double & timestep)
{
    double y_1 = y;
    double ydot_1 = ydot;
    double theta_1 = theta;
    double thetadot_1 = thetadot;
    double ydotk1,yddotk1,thetadotk1,thetaddotk1;
    double ydotk2,yddotk2,thetadotk2,thetaddotk2;
    double ydotk3,yddotk3,thetadotk3,thetaddotk3;
    double ydotk4,yddotk4,thetadotk4,thetaddotk4;

    ydotk1 = ydot_1;
    thetadotk1 = thetadot_1;
    computeDerivative_Anal(y_1,ydotk1,yddotk1, theta_1,thetadotk1,thetaddotk1,
                           mass,intertia,damp,rotdamp,stif,rotstif,force,moment,b);

    ydotk2 = ydot_1+0.5*timestep*yddotk1;
    thetadotk2 = thetadot_1+0.5*timestep*thetaddotk1;
    computeDerivative_Anal(y_1+0.5*timestep*ydotk1, ydotk2, yddotk2,
                           theta_1+0.5*timestep*thetadotk1,thetadotk2,thetaddotk2,
                           mass,intertia,damp,rotdamp,stif,rotstif,force,moment,b);

    ydotk3 = ydot_1+0.5*timestep*yddotk2;
    thetadotk3 = thetadot_1+0.5*timestep*thetaddotk2;
    computeDerivative_Anal(y_1+0.5*timestep*ydotk2, ydotk3, yddotk3,
                           theta_1+0.5*timestep*thetadotk2,thetadotk3,thetaddotk3,
                           mass,intertia,damp,rotdamp,stif,rotstif,force,moment,b);

    ydotk4 = ydot_1+timestep*yddotk3;
    thetadotk4 = thetadot_1+timestep*thetaddotk3;
    computeDerivative_Anal(y_1+timestep*ydotk3, ydotk4, yddotk4,
                           theta_1+timestep*thetadotk3,thetadotk4,thetaddotk4,
                           mass,intertia,damp,rotdamp,stif,rotstif,force,moment,b);

    y = y_1 + timestep/6.0*(ydotk1 + 2.0*ydotk2 + 2.0*ydotk3 + ydotk4);
    ydot = ydot_1 + timestep/6.0*(yddotk1 + 2.0*yddotk2 + 2.0*yddotk3 + yddotk4);
    theta = theta_1 + timestep/6.0*(thetadotk1 + 2.0*thetadotk2 + 2.0*thetadotk3 + thetadotk4);
    thetadot = thetadot_1 + timestep/6.0*(thetaddotk1 + 2.0*thetaddotk2 + 2.0*thetaddotk3 + thetaddotk4);

}