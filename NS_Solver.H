#ifndef _NS_Solver_H
#define _NS_Solver_H

#include "AMR.H"
#include "Mesh.H"
#include "FlowVariables.H"
#include "AMRmpi.H"
#include "AMRSpaceTime.H"
#include "BCValues.H"
#include "Body.H"
#include "Slopelimiter.H"
#include "GhostPair.H"
#include "Turbulence.H"

#define VISCOUSITY
#define SN 2
// #define MUSCL 0	
#define ROE_ALLSPEED
//#define STEADY

#define comptqright(rightvalue, q_i, qi, qi1)\
		SlopeCorrection_Right(q_i, qi, qi1);\
		rightvalue = qi - 0.5*dqii*limitvalue;

#define comptqleft(leftvalue, q_i, qi, qi1)\
		SlopeCorrection_Left(q_i, qi, qi1);\
		leftvalue = qi + 0.5*dqi*limitvalue;

// #endif				

#define FaceLeftRight() \
		comptqleft(left.roe, t_flowptr[0]->roe, t_flowptr[1]->roe, t_flowptr[2]->roe);\
		comptqright(right.roe, t_flowptr[1]->roe, t_flowptr[2]->roe, t_flowptr[3]->roe);\
		comptqleft(left.u, t_flowptr[0]->u, t_flowptr[1]->u, t_flowptr[2]->u);\
		comptqright(right.u, t_flowptr[1]->u, t_flowptr[2]->u, t_flowptr[3]->u);\
		comptqleft(left.v, t_flowptr[0]->v, t_flowptr[1]->v, t_flowptr[2]->v);\
		comptqright(right.v, t_flowptr[1]->v, t_flowptr[2]->v, t_flowptr[3]->v);\
		comptqleft(left.w, t_flowptr[0]->w, t_flowptr[1]->w, t_flowptr[2]->w);\
		comptqright(right.w, t_flowptr[1]->w, t_flowptr[2]->w, t_flowptr[3]->w);\
		comptqleft(left.p, t_flowptr[0]->p, t_flowptr[1]->p, t_flowptr[2]->p);\
		comptqright(right.p, t_flowptr[1]->p, t_flowptr[2]->p, t_flowptr[3]->p);\
		Get_T_ideal_gas(left);\
		Get_E(left);\
		Get_T_ideal_gas(right);\
		Get_E(right)	

#define FaceLeftRight_2ndcd()\
		left = *t_flowptr[1];\
		right = *t_flowptr[2]	

#define Roeaverage()\
		double r2lroe = sqrt(right.roe/left.roe);\
		double r2l0 = 1.0/(1.0+r2lroe);\
		double h_left = left.e+left.p/left.roe;\
		double h_right = right.e+right.p/right.roe;\
		roe0.roe = sqrt(left.roe*right.roe);\
		roe0.u = (left.u+right.u*r2lroe)*r2l0;\
		roe0.v = (left.v+right.v*r2lroe)*r2l0;\
		roe0.w = (left.w+right.w*r2lroe)*r2l0;\
		roe0.h = (h_left+h_right*r2lroe)*r2l0;\
		roe0.c = sqrt((gama-1.0)*(roe0.h-(roe0.u*roe0.u+roe0.v*roe0.v+roe0.w*roe0.w)/2.0));\
		roe0.ubar = roe0.u*kx+roe0.v*ky+roe0.w*kz;



#ifdef ROE_ALLSPEED
#define roe0c_coef()\
		double T_a = (left.T+right.T)*0.5;\
		double V_a = (lv+rv)*0.5;\
		double Ma_a = abs(V_a/sqrt(gama*Rg*T_a));\
		double f_Ma = min(max(Ma_a, Ma), 1.0);\
		double s0 = (1.0-f_Ma)*(1.0-f_Ma);\
		double namda0 = sqrt(4.0*roe0.c*roe0.c*f_Ma+s0*(roe0.ubar-frame_vel)*(roe0.ubar-frame_vel));\
		double namda1 = abs(0.5*((1.0+f_Ma)*(roe0.ubar-frame_vel)+namda0));\
		double namda2 = abs(0.5*((1.0+f_Ma)*(roe0.ubar-frame_vel)-namda0))
		//f_Ma = min(Ma_a*sqrt(4.0+(1.0-Ma_a*Ma_a)*(1.0-Ma_a*Ma_a))/(1.0+Ma_a*Ma_a), 1.0);\
		//f_Ma = max(f_Ma, 0.1)
#else
#define roe0c_coef()\
		double namda1 = abs(roe0.ubar-frame_vel+roe0.c);\
		double namda2 = abs(roe0.ubar-frame_vel-roe0.c);\
		double critcinamda = 0.05*roe0.c;\
		if (namda1 < critcinamda) namda1 = critcinamda;\
		if (namda2 < critcinamda) namda2 = critcinamda; 
#endif				

#define ComptFaceFlux(a_flux, facearea)\
		const double csquare = roe0.c*roe0.c;\
		const double namda3 = abs(roe0.ubar-frame_vel);\
		const double coef12 = 1.0/(2.0*csquare)*facearea;\
		const double roe0cxdubar = roe0.roe*roe0.c*dubar;\
		const double faceareaxabsubar = facearea*namda3;\
		const double halffacearea = 0.5*facearea;\
		alpha[0] = faceareaxabsubar*(droe-dp/csquare);\
		roe0c_coef();\
		alpha[1] = coef12*namda1*(dp+roe0cxdubar);\
		alpha[2] = coef12*namda2*(dp-roe0cxdubar);\
		alpha[3] = alpha[0]+alpha[1]+alpha[2];\
		alpha[4] = roe0.c*(alpha[1]-alpha[2]);\
		alpha[5] = faceareaxabsubar*(roe0.roe*(du[0]-kx*dubar));\
		alpha[6] = faceareaxabsubar*(roe0.roe*(du[1]-ky*dubar));\
		alpha[7] = faceareaxabsubar*(roe0.roe*(du[2]-kz*dubar));\
		a_flux[0] = -alpha[3];\
		a_flux[1] = -(roe0.u*alpha[3]+(kx)*alpha[4]+alpha[5]);\
		a_flux[2] = -(roe0.v*alpha[3]+(ky)*alpha[4]+alpha[6]);\
		a_flux[3] = -(roe0.w*alpha[3]+(kz)*alpha[4]+alpha[7]);\
		a_flux[4] = -(roe0.h*alpha[3]+roe0.ubar*alpha[4]+roe0.u*alpha[5]+\
			roe0.v*alpha[6]+roe0.w*alpha[7]-(csquare/(gama-1.0))*alpha[0]);\
		a_flux[0] += (primflux[0][0]+primflux[1][0])*halffacearea;\
		a_flux[1] += (primflux[0][1]+primflux[1][1])*halffacearea;\
		a_flux[2] += (primflux[0][2]+primflux[1][2])*halffacearea;\
		a_flux[3] += (primflux[0][3]+primflux[1][3])*halffacearea;\
		a_flux[4] += (primflux[0][4]+primflux[1][4])*halffacearea;

#define VarReflect(pr, fmc, dfk0, dfk1, dfk2, dir)\
		pr = dfk0*fmc[0][dir]+dfk1*fmc[1][dir]+dfk2*fmc[2][dir]
		
#ifdef TURBULENCE
#define VarFromKesaitoX()\
				VarReflect(dux, fk0, duk0, duk1, duk2, 0);\
				VarReflect(duy, fk0, duk0, duk1, duk2, 1);\
				VarReflect(duz, fk0, duk0, duk1, duk2, 2);\
				VarReflect(dvx, fk0, dvk0, dvk1, dvk2, 0);\
				VarReflect(dvy, fk0, dvk0, dvk1, dvk2, 1);\
				VarReflect(dvz, fk0, dvk0, dvk1, dvk2, 2);\
				VarReflect(dwx, fk0, dwk0, dwk1, dwk2, 0);\
				VarReflect(dwy, fk0, dwk0, dwk1, dwk2, 1);\
				VarReflect(dwz, fk0, dwk0, dwk1, dwk2, 2);\
				VarReflect(dTx, fk0, dTk0, dTk1, dTk2, 0);\
				VarReflect(dTy, fk0, dTk0, dTk1, dTk2, 1);\
				VarReflect(dTz, fk0, dTk0, dTk1, dTk2, 2);\
				VarReflect(dvisx, fk0, dvisk0, dvisk1, dvisk2, 0);\
				VarReflect(dvisy, fk0, dvisk0, dvisk1, dvisk2, 1);\
				VarReflect(dvisz, fk0, dvisk0, dvisk1, dvisk2, 2)
#else						
#define VarFromKesaitoX()\
				VarReflect(dux, fk0, duk0, duk1, duk2, 0);\
				VarReflect(duy, fk0, duk0, duk1, duk2, 1);\
				VarReflect(duz, fk0, duk0, duk1, duk2, 2);\
				VarReflect(dvx, fk0, dvk0, dvk1, dvk2, 0);\
				VarReflect(dvy, fk0, dvk0, dvk1, dvk2, 1);\
				VarReflect(dvz, fk0, dvk0, dvk1, dvk2, 2);\
				VarReflect(dwx, fk0, dwk0, dwk1, dwk2, 0);\
				VarReflect(dwy, fk0, dwk0, dwk1, dwk2, 1);\
				VarReflect(dwz, fk0, dwk0, dwk1, dwk2, 2);\
				VarReflect(dTx, fk0, dTk0, dTk1, dTk2, 0);\
				VarReflect(dTy, fk0, dTk0, dTk1, dTk2, 1);\
				VarReflect(dTz, fk0, dTk0, dTk1, dTk2, 2)
#endif
#define pfpkesai0(pr, fleft, fright)\
		pr = fright - fleft;\

#define pfpkesai1(pr, up1, up2, down1, down2)\
		pr = 0.25*(up1+up2-down1-down2);\

#define ComptFaceVel(ff, fv, k1, k2, k3)\
		fv = ff.u*k1+ff.v*k2+ff.w*k3;

struct Roeavedata
{
	double roe;
	double u;
	double v;
	double w;
	double h;
	double c;
	double ubar;
	bool hasnan()
	{
		bool nanflag = false;
		double * initaddr = &roe;
		for (int i = 0; i < 7; ++i)
		{
			if (isnan(initaddr[i]))
			{
				printf("Roe0[%d] is nan!!!\n", i);
				nanflag = true; 
			}
		}
		return nanflag;
	}

	void showdata()
	{
		printf("roe0 is roe %f, u %f, v %f, w %f, h %f, c %f, ubar %f\n", 
			roe, u, v, w, h, c, ubar);
	}
};

class NS_Solver
{
	friend class Pltoutput;
public:

	NS_Solver(Mesh & the_mesh, Domain & dm, const int & abodynum):
	a_mesh(the_mesh), bodynum(abodynum)
	{
		int fcnum = dm.TellFaceNum();
		meshbc = BCValues(fcnum);
		sub_time_step = 3;
		tp.resize(sub_time_step);
		if (sub_time_step == 4)
		{
#ifndef STEADY			
			tp[0] = 0.1084;
			tp[1] = 0.2602;
			tp[2] = 0.5052;
			tp[3] = 1.0;
#else 
			tp[0] = 0.0833;
			tp[1] = 0.2069;
			tp[2] = 0.4265;
			tp[3] = 1.0;
#endif						
		}
		else if (sub_time_step == 3)
		{
#ifndef STEADY			
			tp[0] = 0.1918;
			tp[1] = 0.4929;
			tp[2] = 1.0;
#else
			tp[0] = 0.1481;
			tp[1] = 0.4000;
			tp[2] = 1.0;
#endif						
		}
		else if (sub_time_step == 1)
		{
			tp[0] = 1.0;
		}
		bj[0] = 2.30255809;
		bj[1] = -2.49100760;
		bj[2] = 1.57434093;
		bj[3] = -0.38589142;
		
	}

	~NS_Solver()
	{}

	void testperiodic()
	{
		int level0 = a_mesh.MyCurNum();
		for (int level_n = 0; level_n < level0; ++level_n)
		{				
			int start0 = a_mesh.m_level[level_n].m_box.ps();
			int end0 = a_mesh.m_level[level_n].m_box.pe();
			for (int i = start0; i < end0; ++i)
			{
				if (a_mesh.m_level[level_n].m_box[i].iy() < 2 && a_mesh.m_level[level_n].m_box[i].iy() > -1)
				{
					a_mesh.m_level[level_n].m_data[i].roe = 3.0;
				}
				else if (a_mesh.m_level[level_n].m_box[i].iy() < highpt.xy[1]*pow(2, level_n) && a_mesh.m_level[level_n].m_box[i].iy() > highpt.xy[1]*pow(2, level_n)-3)
				{
					a_mesh.m_level[level_n].m_data[i].roe = 6.0;
				}
			}
		}
		//a_mesh.DataExchange_alllevels(0);
	}

	void ResinDomain()
	{
		FlowVec res_sum(0);
#ifdef TURBULENCE
		double res_vis_sum = 0.0;
#endif				
		int boxnum = 0;
		int level0 = 0;
		if (bodynum == 0)
		{
			level0 = a_mesh.MyCurNum();
		}
		else
		{
			level0 = a_mesh.MyCurNum()-1;
		}
		for (int level_n = 0; level_n < level0; ++level_n)
		{
#ifdef TEMPORAL_REFINE
			if (ts%marching_step[level_n] == marching_left_step[level_n])
			{
#endif					
			boxnum += a_mesh.ShowLevelBox(level_n).size();
			int start0 = a_mesh.m_level[level_n].m_box.ps();
			int end0 = a_mesh.m_level[level_n].m_box.pe();
			for (int box_b = start0; box_b < end0; ++box_b)
			{
				int l_b0 = box_b - start0;
				for (int i = 0; i < 5; ++i)
				{
					res_sum[i] = max(abs(source[0][level_n][l_b0][i]), res_sum[i]);										
				}
#ifdef TURBULENCE
				res_vis_sum = max(abs(mutsource[level_n][l_b0]), res_vis_sum);
#endif				
			}
#ifdef TEMPORAL_REFINE
			}
#endif												
		}
		if (bodynum > 0)
		{
			for (int level_n = a_mesh.MyCurNum()-1; level_n < a_mesh.MyCurNum(); ++level_n)
			{							
				boxnum += a_mesh.ShowLevelBox(level_n).size();
				int start0 = a_mesh.m_level[level_n].m_box.ps();
				int end0 = a_mesh.m_level[level_n].m_box.pe();
				for (int box_b = start0; box_b < end0; ++box_b)
				{
					if (a_mesh.infectbox[box_b] == -1 && a_mesh.m_level[level_n].m_box[box_b].pair.signdis > 0.0)
					{
						int l_b0 = box_b - start0;
						for (int i = 0; i < 5; ++i)
						{
							res_sum[i] = max(abs(source[0][level_n][l_b0][i]), res_sum[i]);
						}
#ifdef TURBULENCE
						res_vis_sum = max(abs(mutsource[level_n][l_b0]), res_vis_sum);
#endif						
					}
				}
			}

		}
		MPI_Allreduce(MPI_IN_PLACE, &res_sum[0], 5, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#ifdef TURBULENCE
		MPI_Allreduce(MPI_IN_PLACE, &res_vis_sum, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif	
		if (ts == 0)
		{
			init_res = res_sum;
#ifdef TURBULENCE
			init_res_mut = res_vis_sum;
#endif						
		}	
		for (int i = 0; i < 5; ++i)
		{			
			res_sum[i] /=init_res[i];			
		}
#if DIM == 2
		res_sum[3] = 0.0;
#endif				
#ifdef TURBULENCE
		res_vis_sum /= init_res_mut;
#endif				
		//MPI_Barrier(MPI_COMM_WORLD);
		if (nrank == 0)
		{
#ifndef TURBULENCE			
			printf("Time step %d roe %e u %e v %e w %e e %e\n",
				ts, res_sum[0], res_sum[1], res_sum[2], res_sum[3], res_sum[4]);
			for (int i = 0; i < 5; ++i)
			{
				if (isnan(res_sum[i]))
				{
					MPI_Abort(MPI_COMM_WORLD, 268);
				}
			}
#else
			printf("Time step %d roe %e u %e v %e w %e e %e mut %e\n", ts,
				res_sum[0], res_sum[1], res_sum[2], res_sum[3], res_sum[4], res_vis_sum);
			for (int i = 0; i < 5; ++i)
			{
				if (isnan(res_sum[i]))
				{
					MPI_Abort(MPI_COMM_WORLD, 268);
				}
			}
			if (isnan(res_vis_sum))
			{
				MPI_Abort(MPI_COMM_WORLD, 268);
			}
#endif
		}
	}

	inline BCValues & BCinfo()
	{
		return meshbc;
	}

	void Comptdt();

	inline void PrintDt()
	{
		if (nrank == 0) printf("dt based on CFL is %16.10f, next step is %16.10f present time is %f\n\n", CFL_dt, dt, t);
	}

	void ComptOldConsVar()
	{
		for (int i = 0; i < a_mesh.MyCurNum(); ++i)
		{
#ifdef TEMPORAL_REFINE
			if (ts%marching_step[i] == marching_left_step[i])
			{
			//if (nrank == 0) printf("ts %d marching level %d!!!\n", i);
#endif			
			int start0 = a_mesh.LevelBoxStart(i);
			int end0 = a_mesh.LevelBoxEnd(i);
			for (int bn = start0; bn < end0; ++bn)
			{
				int b0 = bn - start0;
				Get_Cons_Vars(a_mesh.m_level[i].m_data[bn], cv_old[i][b0]);
				//cv_old[i][b0].showdata(i, bn, -1);
#ifdef TURBULENCE
				c_vis_old[i][b0] = a_mesh.m_level[i].m_data[bn].viseddy;				
#endif
#ifdef DEBUG							
				if (cv_old[i][b0].hasnan(i, bn, 0))
				{
					cv_old[i][b0].showdata(i, bn, -1);
					a_mesh.m_level[i].m_data[bn].showdata("old cv computation");
				}
				Assert(!cv_old[i][b0].hasnan(i, bn, 0), "old cv check before update", 121);
#endif				
			}
#ifdef TEMPORAL_REFINE
			}
#endif						
		}
	}

	void UpdateConsVars(const double & dt0);

	void RestoreConsVar();

	void ComptNewPrimVar();

	void InitFlow()
	{
		for (int i = 0; i < a_mesh.MyCurNum(); ++i)
		{
			for (int b0 = a_mesh.m_level[i].m_box.ps(); b0 < a_mesh.m_level[i].m_box.pe(); ++b0)
			{
				meshbc.FlowFieldInit(a_mesh.m_level[i].m_data[b0]);
			}
			for (int b0 = a_mesh.m_level[i].m_box.gps(); b0 < a_mesh.m_level[i].m_box.gpe(); ++b0)
			{
				meshbc.FlowFieldInit(a_mesh.m_level[i].m_data[b0]);
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}

	void InitFlow_User()
	{
		for (int i = 0; i < a_mesh.MyCurNum(); ++i)
		{
			for (int b0 = a_mesh.m_level[i].m_box.ps(); b0 < a_mesh.m_level[i].m_box.pe(); ++b0)
			{
				meshbc.UserInit(a_mesh, i, b0);
			}
			MPI_Barrier(MPI_COMM_WORLD);
			a_mesh.DataExchange(a_mesh.m_level[i].nd, i);
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}

	template <class T>
	void CheckDoubleVectSizeBox(vector<vector<T> > & acontainer)
	{
		if (acontainer.size() < a_mesh.MyCurNum())
		{
			acontainer.resize(a_mesh.MyCurNum());
		}
		// else if (acontainer.capacity() > a_mesh.MyCurNum())
		// {
		// 	//printf("[%d]acontainer capacity is %d size is %d\n", nrank, (int)acontainer.capacity(), (int)acontainer.size());
		// 	vector<vector<T> > acttemp(a_mesh.MyCurNum());
		// 	acontainer.swap(acttemp);
		// }
		for (int i = 0; i < a_mesh.MyCurNum(); ++i)
		{
			int etnum = a_mesh.LevelBoxEnd(i) - a_mesh.LevelBoxStart(i);
			acontainer[i].resize(etnum);
			// else if (acontainer[i].capacity() > etnum + 20)
			// {
			// 	vector<T> ctemp(etnum);
			// 	acontainer[i].swap(ctemp);
			// }
		}
	}

	template <class T>
	void CheckDoubleVectSizeFace(vector<vector<T> > & acontainer)
	{
		if (acontainer.size() < a_mesh.MyCurNum())
		{
			acontainer.resize(a_mesh.MyCurNum());
		}
		// else if (acontainer.capacity() > a_mesh.MyCurNum())
		// {
		// 	vector<vector<T> > acttemp(a_mesh.MyCurNum());
		// 	acontainer.swap(acttemp);
		// }
		for (int i = 0; i < a_mesh.MyCurNum(); ++i)
		{
			int etnum = a_mesh.Facestart(i) - a_mesh.Faceend(i);
			acontainer[i].resize(etnum);
			// if (acontainer[i].size() < etnum)
			// {
				
			// }
			// else if (acontainer[i].capacity() > etnum + 20)
			// {
			// 	vector<T> ctemp(etnum);
			// 	acontainer[i].swap(ctemp);
			// }
		}
	}

	template <class T>
	void CheckVectDAFace(vector<DataArray<T> > & acontainer)
	{
		if (acontainer.size() < a_mesh.MyCurNum())
		{
			acontainer.resize(a_mesh.MyCurNum());
			for (int i = 0; i < a_mesh.MyCurNum(); ++i)
			{
				if (acontainer[i].size() == 0) acontainer[i] = DataArray<T>();
			}
		}
		for (int i = 0; i < a_mesh.MyCurNum(); ++i)
		{
			acontainer[i].setnum_nocopy(a_mesh.Meshface(i).size(), 0);
		}
	}

	template <class T>
	void CheckVectDABox(vector<DataArray<T> > & acontainer)
	{
		if (acontainer.size() < a_mesh.MyCurNum())
		{
			acontainer.resize(a_mesh.MyCurNum());
			for (int i = 0; i < a_mesh.MyCurNum(); ++i)
			{
				if (acontainer[i].size() == 0) acontainer[i] = DataArray<T>();
			}
		}
		for (int i = 0; i < a_mesh.MyCurNum(); ++i)
		{
			acontainer[i].setnum_nocopy(a_mesh.m_level[i].m_box.size(), 0);
		}
	}

	void CheckDataContainer()
	{
		CheckDoubleVectSizeBox(cv_old);
		CheckDoubleVectSizeBox(cv_new);
		CheckVectDAFace(faceflux);
		CheckVectDAFace(vsflux);
		for (int i = 0; i < SN; ++i)
		{
			CheckDoubleVectSizeBox(source[i]);
		}
#ifdef LOCAL_TIME_STEPPING		
		CheckVectDABox(local_dt);
		CheckVectDABox(local_cfl);
		CheckVectDABox(local_cfl_ismax);
#endif		
#ifdef TURBULENCE
		CheckDoubleVectSizeBox(c_vis_old);
		CheckVectDAFace(mutflux);
		CheckDoubleVectSizeBox(mutsource);
#endif				
	}

	void ComptSource();

	void ComptViscousStress();

/*------------------------------------------------------------------------------------*/
#ifdef TIMECOUNTING
	void TimeMarching(AMR & amr)
	{

		if (ts%dt_var_step == 0) Comptdt();
		//PrintDt();
		double tm[20], tsm;
		tm[0] = MPI_Wtime();
		CheckDataContainer();
		GiveAFlag("finish CheckDataContainer", 5);
		tm[1] = MPI_Wtime();
		ComptOldConsVar();
		GiveAFlag("finish ComptOldConsVar", 5);
		tm[2] = MPI_Wtime();
		meshbc.BCBoxTreat(a_mesh);
		tm[3] = MPI_Wtime();
#ifdef VISCOUSITY
		ComptViscousStress();			
		GiveAFlag("finish compute VISCOUSITY", 5);
		tm[4] = MPI_Wtime();
#endif
		step_inv_time = 0.0;
		step_exchange_time = 0.0;
#ifdef TURBULENCE
		step_tur_time = 0.0;
#endif
		double ivst[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
		for (int subt = 0; subt < sub_time_step; ++subt)
		{
			double dt0 = tp[subt]*dt;
			tm[5] = MPI_Wtime();
			ComputeFaceFlux();
			GiveAFlag("start to compute source!!!", 5);
			tm[6] = MPI_Wtime();
			ivst[0] += tm[6] - tm[5];
			ComptSource();
			GiveAFlag("finish to compute source!!!", 5);
			tm[7] = MPI_Wtime();
			ivst[1] += tm[7] - tm[6];
#ifdef TURBULENCE
			ComptMutSource();
			tm[8] = MPI_Wtime();
			ivst[5] += tm[8] - tm[7];
#endif			
			UpdateConsVars(dt0);
			GiveAFlag("finish UpdateConsVars!!!", 5);
			tm[9] = MPI_Wtime();
			ivst[2] += tm[9] - tm[8];			
			MPI_Barrier(MPI_COMM_WORLD);
			ComptNewPrimVar();
			tm[10] = MPI_Wtime();
			ivst[3] += tm[10] - tm[9];
			if (subt < 3) a_mesh.DataExchange_alllevels(ts);
			tm[11] = MPI_Wtime();
			ivst[4] += tm[11]-tm[10];
		}
		tm[12] = MPI_Wtime();
		amr.InterfaceExchange();
		tm[13] = MPI_Wtime();
		ResinDomain();
		tsm = MPI_Wtime();
		double atsm = tsm - tm[0];
		if (nrank == 0)
		{
			printf("Time for CheckDataContainer is  %6.2f%\n", (tm[1]-tm[0])/atsm*100.0);
			printf("Time for ComptOldConsVar is     %6.2f%\n", (tm[2]-tm[1])/atsm*100.0);
			printf("Time for BCBoxTreat is          %6.2f%\n", (tm[3]-tm[2])/atsm*100.0);
			printf("Time for ComptViscousStress is  %6.2f%\n", (tm[4]-tm[3])/atsm*100.0);
			printf("Time for ComptFaceFlux is       %6.2f%\n", ivst[0]/atsm*100.0);
			printf("Time for ComptSource is         %6.2f%\n", ivst[1]/atsm*100.0);
#ifdef TURBULENCE			
			printf("Time for ComptmutSource is      %6.2f%\n", ivst[5]/atsm*100.0);
#endif			
			printf("Time for UpdateConsVars is      %6.2f%\n", ivst[2]/atsm*100.0);
			printf("Time for ComptNewPrimVar is     %6.2f%\n", ivst[3]/atsm*100.0);
			printf("Time for DataExchange_alllev is %6.2f%\n", ivst[4]/atsm*100.0);
			printf("Time for InterfaceExchange is   %6.2f%\n", (tm[13]-tm[12])/atsm*100.0);
			printf("Time for ResinDomain is         %6.2f%\n", (tsm-tm[13])/atsm*100.0);
		}		
	}
#else
	void TimeMarching(AMR & amr)
	{
		CheckDataContainer();
		GiveAFlag("finish CheckDataContainer", 5);
		if (ts%dt_var_step == 0) Comptdt();
		GiveAFlag("Finish compute dt...", 5);
		//PrintDt();
		ComptOldConsVar();
		GiveAFlag("finish ComptOldConsVar", 5);
		meshbc.BCBoxTreat(a_mesh);
		GiveAFlag("finish BCBoxTreat", 5);
		//amr.InterfaceExchange();
#ifdef VISCOUSITY
		ComptViscousStress();			
		GiveAFlag("finish compute VISCOUSITY", 5);
		//BalanceFaceFlux(a_mesh.mdyface, vsflux);
#endif
#ifdef SHOWTIME
		step_inv_time = 0.0;
		step_exchange_time = 0.0;
#ifdef TURBULENCE
		step_tur_time = 0.0;
#endif				
#endif
		for (subt = 0; subt < sub_time_step; ++subt)
		{
			double dt0 = tp[subt]*dt;
			ComputeFaceFlux();
			GiveAFlag("start to compute source!!!", 5);
			//if (subt == sub_time_step-1) BalanceFaceFlux(a_mesh.mdyface, faceflux);
			ComptSource();
			GiveAFlag("finish to compute source!!!", 5);
#ifdef TURBULENCE
			ComptMutSource();
			GiveAFlag("Finish compute mutsource!!!", 5);
#endif						
			UpdateConsVars(dt0);
			GiveAFlag("finish UpdateConsVars!!!", 5);
			ComptNewPrimVar();
			GiveAFlag("finish compute new prim var!!!", 5);						
			//if (subt < 3) meshbc.BCBoxTreat(a_mesh);
			if (subt < sub_time_step-1) a_mesh.DataExchange_alllevels(ts);
		}
		GiveAFlag("finish InterfaceExchange!!!", 5);
		ResinDomain();
	}
#endif	
	void ComputeFaceFlux()
	{
		int dr[2] = {-1,1};
		int dn[3] = {0,0,0};
		FlowVec primflux[2];
		FlowVariables left, right;
		Roeavedata  roe0;
		double alpha[8];
		int * dnptr = &dn[0];
#ifdef SHOWTIME
		double starttime = MPI_Wtime();
#endif
		double limitvalue;
		double dqi, dqii, qr;	
		for (int level_n = 0; level_n < a_mesh.MyCurNum(); ++level_n)
		{
#ifdef TEMPORAL_REFINE
			if (ts%marching_step[level_n] == marching_left_step[level_n])
			{
#endif			
			int end0 = a_mesh.m_level[level_n].m_face.pe();
			for (int bn = a_mesh.m_level[level_n].m_face.ps(); bn < end0; ++bn)
			{
				Face & theface = a_mesh.m_level[level_n].m_face[bn];
#if DIM == 2
				if (theface.fnv != 2)
				{
#endif								
				dnptr[theface.fnv] = 1;
				int ndx, ndy, ndz;
				int fn1, fn2;
				for (int fn = 0; fn < 2; ++fn)
				{
					ndx = 1+dn[0]*dr[fn];
					ndy = 1+dn[1]*dr[fn];
					ndz = 1+dn[2]*dr[fn];
					int aneib = a_mesh.m_level[level_n].m_box[theface[fn]].neib[ndx][ndy][ndz];
					Assert(theface[fn] > -1 && theface[fn] < a_mesh.m_level[level_n].m_box.realsize(), "Negative index for a face side box!!!", 208);
					// Assert(aneib > -1 && aneib < a_mesh.m_level[level_n].m_box.realsize(), "Negative index for a second face side box!!!", 209);
					if (aneib < 0) 
					{
						int rfn = 1-fn;
						int ndx0 = 1 + dn[0]*dr[rfn];
						int ndy0 = 1 + dn[1]*dr[rfn];
						int ndz0 = 1 + dn[2]*dr[rfn];
						aneib = a_mesh.m_level[level_n].m_box[theface[rfn]].neib[ndx0][ndy0][ndz0];
					}	
					fn1 = fn+1;
					fn2 = fn+1+dr[fn];
					t_flowptr[fn1] = &a_mesh.m_level[level_n].m_data[theface[fn]];
					t_flowptr[fn2] = &a_mesh.m_level[level_n].m_data[aneib];
				}				
				// t_flowptr[0]->showdata("0");
				// t_flowptr[1]->showdata("1");
				// t_flowptr[2]->showdata("2");
				// t_flowptr[3]->showdata("3");
				// if (ts == 0)
				// {
				// 	printf("L%dF%d dir is (%f,%f,%f)\n", level_n, bn, theface.keisa[0], theface.keisa[1], theface.keisa[2]);
				// }
#ifdef SECOND_CENTRALDIFF
				FaceLeftRight_2ndcd();
#else			
				FaceLeftRight();
#endif			
				const double & kx = a_mesh.m_level[level_n].m_face[bn].keisa[0];
				const double & ky = a_mesh.m_level[level_n].m_face[bn].keisa[1];
				const double & kz = a_mesh.m_level[level_n].m_face[bn].keisa[2];

				Pointxyz rot_vel1 = frame_rotate_vel.cross(theface.fc1 - frame_rotate_origin) + frame_trans_vel;
				Pointxyz rot_vel2 = frame_rotate_vel.cross(theface.fc2 - frame_rotate_origin) + frame_trans_vel;
				double frame_vel = rot_vel1.dot(theface.keisa1) + rot_vel2.dot(theface.keisa2);

				const double lv = left.u*kx + left.v*ky + left.w*kz - frame_vel;
				const double rv = right.u*kx + right.v*ky + right.w*kz - frame_vel;

				Get_inv_flux(left, lv, primflux[0], kx, ky, kz);
				Get_inv_flux(right, rv, primflux[1], kx, ky, kz);

				primflux[0][4] += frame_vel*left.p;
				primflux[1][4] += frame_vel*right.p;

#ifdef SECOND_CENTRALDIFF
				faceflux[level_n][bn] = (primflux[0]+primflux[1])*0.5*a_mesh.m_level[level_n].m_face[bn].area;
#else
				const double droe = right.roe-left.roe;
				const double dp = right.p-left.p;
				const double du[3] = {right.u-left.u, right.v-left.v, right.w-left.w};
				const double dubar = rv - lv;
				Roeaverage();								
				ComptFaceFlux(faceflux[level_n][bn], a_mesh.m_level[level_n].m_face[bn].area);						
#ifdef DEBUG
				if (roe0.hasnan())
				{
					roe0.showdata();
					left.showdata("roe0 left");
					right.showdata("roe0 right");
					MPI_Abort(MPI_COMM_WORLD, 33);
				}
				if (left.hasnan(level_n, bn, "left") || right.hasnan(level_n, bn, "right"))
				{
					MPI_Abort(MPI_COMM_WORLD, 279);
				}
				if (faceflux[level_n][bn].hasnan(level_n, bn, a_mesh.Levelface(level_n, bn).fnv))
				{
					MPI_Abort(MPI_COMM_WORLD, 283);
				}
#endif

#endif
				dnptr[theface.fnv] = 0;
#if DIM == 2
				}
#endif								
			}
#ifdef TEMPORAL_REFINE
			}
#endif						
			//MPI_Win_unlock_all(faceflux[level_n].arraywin());
			//MPI_Win_fence(0, faceflux[level_n].arraywin());
		}
		//GiveAFlag("finish flux computation!!!", 5);
		MPI_Barrier(MPI_COMM_WORLD);
#ifdef SHOWTIME		
		double endtime = MPI_Wtime();
		step_inv_time += endtime - starttime;
		total_inv_time += step_inv_time;
#endif		
	}

	void BalanceFaceFlux(vector<DataArray<Facepair> > & mdyface, vector<DataArray<FlowVec> > & fluxarray);

	void MiddleStepFlux(vector<DataArray<Facepair> > & mdyface);

	void ConfineWallface();

#ifdef TURBULENCE
	void ComptMutSource();
#endif
/*----------------------------------------------------------------------------------------*/
protected:
	Mesh & a_mesh;
	BCValues meshbc;
	vector<vector<FlowVec> > cv_old;
	vector<vector<FlowVec> > cv_new;
	vector<vector<FlowVec> > source[SN];
	vector<DataArray<FlowVec> > faceflux;
	vector<DataArray<FlowVec> > vsflux;
	vector<DataArray<double> > local_dt;
	vector<DataArray<double> > local_cfl;
	vector<DataArray<bool> > local_cfl_ismax;
	const FlowVariables * vsflowptr[3][3][3];
	const FlowVariables * t_flowptr[4];
	/*Temporal scheme*/
	int sub_time_step;
	int subt;
	vector<double> tp;
	double bj[4];
	const int bodynum;
	/*For limiter*/
	//Slopelimiter roelimit;
	FlowVec init_res;
#ifdef TURBULENCE
	vector<DataArray<double> > mutflux;
	vector<vector<double> > mutsource;
	vector<vector<double> > c_vis_old;
	double init_res_mut;
#endif		
};

#endif