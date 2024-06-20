#include "Constants.H"
#include <cmath>
#include "FlowVariables.H"

double gama = 1.4;
double Rg_ref = 287.4;
double Re = roe_ref*u_ref*L_ref/mu_ref;
const double Ma = u_ref/sqrt(gama*Rg_ref*T_ref);

double Pr = 0.72;
double Pr_tur = 0.9;
double Rg = 1.0/(gama*pow(Ma, 2));

double Cv = Rg/(gama-1.0);
double Cp = Cv + Rg;

double mu = 1.0;
double qk = mu/(Re*Pr*(gama-1)*pow(Ma,2));
double viseddy = 1.0;

//double qk = gama*mu/(Re*Pr*(gama-1.0));

FlowVariables initvar;
double inlet_flow_angle;

double Cv_ref = Rg_ref/(gama-1.0);
double Cp_ref = Cv_ref + Rg_ref;
double qk_ref = mu_ref*Cp_ref/Pr;


double p_ref = roe_ref*pow(u_ref, 2);
double vis_ref = mu_ref/roe_ref;
double e_ref = pow(u_ref, 2);
//double e_ref = Rg_ref*T_ref;

double t_ref = L_ref/u_ref;
double omiga_ref = u_ref/L_ref;

double shear_ref = mu_ref*u_ref/L_ref;
double area_ref = L_ref*L_ref;

double S_over_T_ref = 110.33333333/T_ref;

double Pr_cubicroot = pow(Pr, 1.0/3.0);
