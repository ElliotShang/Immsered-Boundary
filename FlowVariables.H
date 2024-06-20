#ifndef _FlowVariables_H_
#define _FlowVariables_H_

#include <cmath>
#include "Constants.H"
#include "AMRmpi.H"
#include <string>

extern double roe_ref;
extern double L_ref;
extern double T_ref;
extern double mu_ref;
extern double u_ref;

//#define OTHER_VARS 2

//#define TURBULENCE

#define kineticenergy(fv) (0.5*(fv.u*fv.u+fv.v*fv.v+fv.w*fv.w))

#define Get_p_ideal_gas(fv) fv.p = fv.roe*Rg*fv.T

#define Get_T_ideal_gas(fv) fv.T = fv.p/(fv.roe*Rg)

#define Get_roe_ideal_gas(fv) fv.roe = fv.p/(Rg*fv.T)

#define speed(fv) sqrt(fv.u*fv.u+fv.v*fv.v+fv.w*fv.w)

#define soundspeed(fv) sqrt(gama*Rg*fv.T)

#define Manumber(fv) (speed(fv)/soundspeed(fv))

#define Get_E(fv) fv.e = Cv*fv.T+kineticenergy(fv)

#define Get_p_energy(fv) fv.p = (fv.e-kineticenergy(fv))*(gama-1.0)*fv.roe

#define normalvel(fv, nv) (fv.u*nv[0]+fv.v*nv[1]+fv.w*nv[2]);

#define Total_Pressure(fv) fv.p*pow(1.0+0.2*(fv.u*fv.u+fv.v*fv.v+fv.w*fv.w)/(gama*Rg*fv.T), 1.4/0.4)

#define Total_Temperature(fv) fv.T + kineticenergy(fv)/Cp

#define Get_Cons_Vars(fv,cv)\
		cv[0] = fv.roe;\
		cv[1] = fv.roe*fv.u;\
		cv[2] = fv.roe*fv.v;\
		cv[3] = fv.roe*fv.w;\
		cv[4] = fv.roe*fv.e

#define Get_Prim_Vars(fv,cv)\
		fv.roe = cv[0];\
		fv.u = cv[1]/cv[0];\
		fv.v = cv[2]/cv[0];\
		fv.w = cv[3]/cv[0];\
		fv.e = cv[4]/cv[0];\
		Get_p_energy(fv);\
		Get_T_ideal_gas(fv)
	
#define Get_inv_flux(fv, flowvelocity, a0, kx, ky, kz)\
		a0[0] = fv.roe*flowvelocity;\
		a0[1] = fv.roe*fv.u*flowvelocity+fv.p*kx;\
		a0[2] = fv.roe*fv.v*flowvelocity+fv.p*ky;\
		a0[3] = fv.roe*fv.w*flowvelocity+fv.p*kz;\
		a0[4] = (fv.p+fv.e*fv.roe)*flowvelocity

struct FlowVec
{
	double fv[5];

	FlowVec()
	{
		for (int i = 0; i < 5; ++i)
		{
			fv[i] = 0.0;
		}
	}

	FlowVec(const double & a0)
	{
		fv[0] = a0;
		fv[1] = a0;
		fv[2] = a0;
		fv[3] = a0;
		fv[4] = a0;
	}

	inline void setvalue(const int & a0)
	{
		fv[0] = a0;
		fv[1] = a0;
		fv[2] = a0;
		fv[3] = a0;
		fv[4] = a0;
	}

	inline FlowVec operator+(FlowVec & a0)
	{
		FlowVec b0;
		for (int i = 0; i < 5; ++i)
		{
			b0[i] = this->fv[i] + a0.fv[i];
		}
		return b0;
	}

	inline FlowVec operator-(const FlowVec & a0)
	{
		FlowVec b0;
		for (int i = 0; i < 5; ++i)
		{
			b0[i] = this->fv[i] - a0.fv[i];
		}
		return b0;
	}

	inline FlowVec operator*(const double & a0)
	{
		FlowVec b0;
		for (int i = 0; i < 5; ++i)
		{
			b0[i] = this->fv[i] * a0;
		}
		return b0;
	}
	inline FlowVec operator/(const double & a0)
	{
		FlowVec b0;
		for (int i = 0; i < 5; ++i)
		{
			b0[i] = this->fv[i]/a0;
		}
		return b0;
	}
	inline double & operator[](const int & n0)
	{
		return fv[n0];
	}

	inline void operator+=(const FlowVec & a0)
	{
		fv[0] += a0.fv[0];
		fv[1] += a0.fv[1];
		fv[2] += a0.fv[2];
		fv[3] += a0.fv[3];
		fv[4] += a0.fv[4];
	}

	inline void operator-=(const FlowVec & a0)
	{
		fv[0] -= a0.fv[0];
		fv[1] -= a0.fv[1];
		fv[2] -= a0.fv[2];
		fv[3] -= a0.fv[3];
		fv[4] -= a0.fv[4];
	}

	inline void operator/=(const double & a0)
	{
		fv[0] /= a0;
		fv[1] /= a0;
		fv[2] /= a0;
		fv[3] /= a0;
		fv[4] /= a0;
	}

	inline void operator*=(const double & a0)
	{
		fv[0] *= a0;
		fv[1] *= a0;
		fv[2] *= a0;
		fv[3] *= a0;
		fv[4] *= a0;
	}

	bool hasnan(int & l0, int & b0, const int & dir)
	{
		bool nanflag = false;
		for (int i = 0; i < 5; ++i)
		{
			if (isnan(fv[i]))
			{
				printf("P%dL%dF%d FlowVec [%d] is nan!!!Face direction is %d\n", 
					srank, l0, b0, i, dir);
				nanflag = true; 
			}
		}
		return nanflag;
	}

	bool hasnan(const string & ast)
	{
		bool nanflag = false;
		for (int i = 0; i < 5; ++i)
		{
			if (isnan(fv[i]))
			{
				printf("FlowVec [%d] is nan flag is %s!!!\n", 
					i, ast.c_str());
				nanflag = true; 
			}
		}
		return nanflag;
	}

	// FlowVec & operator-()
	// {
	// 	for (int i = 0; i < 5; ++i)
	// 	{
	// 		fv[i] = 0.0-fv[i];
	// 	}
	// 	return *this;
	// }

	void showdata(const int & l0, const int & b0, const int & dir)
	{
		printf("P%dL%dF%d is (%f,%f,%f,%f,%f) direction is %d\n", 
			srank,l0,b0,fv[0],fv[1],fv[2],fv[3],fv[4],dir);
	}

	void showdata(const string & tag0)
	{
		printf("%s (%f,%f,%f,%f,%f)\n", 
			tag0.c_str(),fv[0],fv[1],fv[2],fv[3],fv[4]);
	}

	void initzero()
	{
		for (int i = 0; i < 5; ++i)
		{
			fv[i] = 0.0;
		}
	}
};

struct FlowVariables
{
	FlowVariables()
	{}

	inline void zero()
	{
		roe = 0.0;
		u = 0.0;
		v = 0.0;
		w = 0.0;
		p = 0.0;
		e = 0.0;
		T = 0.0;
#ifdef TURBULENCE
		viseddy = 0.0;
#endif
#ifdef OTHER_VARS
		for (int i = 0; i < OTHER_VARS; ++i)
		{
			var[i] = 0.0;
		}
#endif			
	}

	inline double & operator[](const int & i)
	{
		double * aimvar = &roe;
		return aimvar[i];
	}
	inline double & at(const int & i)
	{
		double * aimvar = &roe;
		return aimvar[i];
	}
	bool hasnan(const int & l0, 
							const int & b0, 
							const string & arrname)
	{
		bool nanflag = false;
		double * initaddr = &roe;
		for (int i = 0; i < 7; ++i)
		{
			if (isnan(initaddr[i]))
			{
				printf("P%dL%dB%d FlowVariables %s[%d] is nan!!!\n", nrank, l0, b0,
					arrname.c_str(), i);
				nanflag = true; 
			}
		}
		if (isnan(Manumber((*this))))
		{
			printf("P%dL%dB%d Ma is nan T is %f!!!\n", nrank, l0, b0, T);
			nanflag = true;
		}
		return nanflag;
	}

	bool roeisright(const string & roec) const
	{
		bool right = true;
		if (roe < 0.0000001)
		{
			right = false;
			printf("Flow density in array %s is %f, please check it!!!\n", roec.c_str(), roe);
		}
		return right;
	}
#ifdef TURBULENCE
	void showdata(const string & tag0) const
	{
		printf("roe %f u %f v %f w %f p %f T %f e %f viseddy %f tag is %s\n", roe,u,v,w,p,T,e,viseddy,tag0.c_str());
	}
#else
	void showdata(const string & tag0) const
	{
		printf("roe %f u %f v %f w %f p %f T %f e %f tag is %s\n", roe,u,v,w,p,T,e,tag0.c_str());
	}
#endif	
	void showdata(const int & i0, const int & j0, const int & k0, const double & dis0,const double & dis1) const
	{
		printf("roe %f u %f v %f w %f p %f T %f e %f i0 %d j0 %d k0 %d f1 %f f2 %f\n", 
			roe,u,v,w,p,T,e,i0,j0,k0,dis0,dis1);
	}
	/*the order of variables can not be modified*/
	double roe;
	double u;
	double v;
	double w;
	double p;
	double T;	
	double e;
#ifdef TURBULENCE
	double viseddy;
#endif		
	//FlowVec his[3];
#ifdef OTHER_VARS
	double var[OTHER_VARS];
#endif
};

#define DEFINE_MPI_FV() \
        MPI_Type_contiguous(sizeof(FlowVariables)/sizeof(double), MPI_DOUBLE, &MPI_FV);\
        MPI_Type_commit(&MPI_FV)

extern FlowVariables initvar;
extern double inlet_flow_angle;

#endif