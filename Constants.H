#ifndef _Constants_H_
#define _Constants_H_

// #define bcny 1
// #define bcnz 1

#define ighost 2
//Data transfer layer number
#define dtln 3

#define pai 3.141592654
//#define OTHER_VARS 1

extern double gama;
extern double Rg;
extern double Cv;
extern double Cp;
extern const double Ma;
extern double Re;
extern double mu;
extern double Pr;
extern double Pr_tur;
extern double qk;
extern double viseddy;


/*------------Reference value-------*/
extern double mu_ref;
extern double Rg_ref;
extern double Cv_ref;
extern double Cp_ref;
extern double qk_ref;

extern double u_ref;
extern double roe_ref;
extern double T_ref;
extern double p_ref;
extern double p_inlet;
extern double e_ref;

extern double L_ref;
extern double t_ref;
extern double omiga_ref;

extern double shear_ref;
extern double area_ref;
extern double vis_ref;
extern double S_over_T_ref;
extern double Pr_cubicroot;
#endif