#ifndef _BCValues_H_
#define _BCValues_H_

#include <vector>
#include "Box.H"
#include "FlowVariables.H"
#include "Boxloc.H"
#include "Domain.H"
#include "Mesh.H"
#include "Turbulence.H"

using namespace std;

class BCValues
{
private:
	vector<FlowVariables> facevar;
	vector<FlowVariables> inletparams;
	double backpressure;
	vector<Pointxyz> face_force_rotate;
	vector<Pointxyz> face_force_move;
	vector<Pointxyz> face_force_rot_center;
public:

	friend class NS_Solver;
	friend class Pltoutput;
	BCValues()
	{}

	BCValues(const int & fnum)
	{
		facevar.resize(fnum);
		face_force_move.resize(fnum);
		face_force_rotate.resize(fnum);
		face_force_rot_center.resize(fnum);
		SetInitValues();
		SetBDFaceVars();
	}

	void SetInitValues();

	void FlowFieldInit(FlowVariables & abox)
	{
		abox = initvar;
	}

	void UserInit(Mesh & amesh, const int & ilevel, const int & bn);

	void SetBackPressure(Mesh & amesh);

	void BCBoxTreat(Mesh & amesh)
	{
		SetBackPressure(amesh);
		for (int nlevel = 0; nlevel < amesh.cur_level_num; ++nlevel)
		{
#ifdef TEMPORAL_REFINE
			if (ts%marching_step[nlevel] == marching_left_step[nlevel])
			{
#endif			
			DataArray<Domainghost> & dmg = amesh.MeshGhostArray(nlevel);
			int ds0 = dmg.ps(); int de0 = dmg.pe();
			for (int i = ds0; i < de0; ++i)
			{
				int f0 = dmg[i].nbface;
				int bdtype = amesh.m_dm.facetype[f0];
#if DIM == 2
				if (amesh.m_dm.facedir[f0][2] == 0)
				{
#endif

				if (dmg[i].refcell > -1)
				{															
					if (bdtype == -99) {} //Periodic Boundary Condition
					else if (bdtype == -50)	SetBCType_50(amesh.m_level[nlevel].m_data[dmg[i].cell], 
																										 amesh.m_level[nlevel].m_data[dmg[i].refcell], dmg[i].nbface);
					else if (bdtype == -51)	
					{
						Boxloc & refface = amesh.m_dm.refboxface[dmg[i].nbface];
						Assert((refface.node > -1 && refface.node < 3 &&
										refface.index > -1 && refface.index < 2), 
										"Miss a right ref box face for bc type 51!!!", 242);
						Assert(dmg[i].nmvflag, "The type 51 bc cell does not have a normal vector!!!",95);
						int f0 = amesh.m_level[nlevel].m_box[dmg[i].refcell].faces[refface.node][refface.index];
						SetBCType_51(amesh.m_level[nlevel].m_data[dmg[i].cell], 
												 amesh.m_level[nlevel].m_data[dmg[i].refcell], 
												 dmg[i].nbface,
												 dmg[i].nmv);
					}
#ifdef TURBULENCE					
					else if (bdtype == -81)	
					{
						Boxloc & refface = amesh.m_dm.refboxface[dmg[i].nbface];						
						Assert((refface.node > -1 && refface.node < 3 &&
										refface.index > -1 && refface.index < 2), 
										"Miss a right ref box face for bc type 81!!!", 242);
						Assert(dmg[i].nmvflag, "The type 81 bc cell does not have a normal vector!!!",95);
						int f0 = amesh.m_level[nlevel].m_box[dmg[i].refcell].faces[refface.node][refface.index];
						Pointxyz myfc = (amesh.m_level[nlevel].m_face[f0].fc1+amesh.m_level[nlevel].m_face[f0].fc2)*0.5;
						Pointxyz fcvel = face_force_rotate[dmg[i].nbface].cross(myfc - face_force_rot_center[dmg[i].nbface]) + face_force_move[dmg[i].nbface];

						Point neib0(1,1,1);
						*(&neib0[0]+refface.node) = 2-2*refface.index;
						int refneib = amesh.m_level[nlevel].m_box[dmg[i].refcell].neib[neib0[0]][neib0[1]][neib0[2]];						
						Assert(amesh.m_level[nlevel].m_box[refneib].type == Normalcell, "the bd cell ref neib must a normal cell!!!", 101);
						if (amesh.m_level[nlevel].m_box[refneib].pair.distance_to_dm > 1.0)
						{
							amesh.m_level[nlevel].m_data[dmg[i].cell] = amesh.m_level[nlevel].m_data[dmg[i].refcell];
							// printf("L%dB(%d,%d,%d) distance to domain is %f which can not a ref cell for a domain ghost!!!!\n",
							// 	nlevel,
							// 	amesh.m_level[nlevel].m_box[refneib].ix(),
							// 	amesh.m_level[nlevel].m_box[refneib].iy(),
							// 	amesh.m_level[nlevel].m_box[refneib].iz(),
							// 	amesh.m_level[nlevel].m_box[refneib].pair.distance_to_dm);
							// MPI_Abort(MPI_COMM_WORLD, 112);
						}				
						else
						{
							ImageViseddy_Domain(amesh.m_level[nlevel].m_data[refneib], 
											amesh.m_level[nlevel].m_data[dmg[i].cell], 
											fcvel, 
											dmg[i].nmv, 
											amesh.m_level[nlevel].m_box[refneib].pair.distance_to_dm, 
											amesh.m_level[nlevel].m_box[dmg[i].refcell].pair.distance_to_dm);
						}
					}
#endif					
					else if (bdtype == -46)	
					{
						SetBCType_46(amesh.m_level[nlevel].m_data[dmg[i].cell],
												 amesh.m_level[nlevel].m_data[dmg[i].refcell], 
												 dmg[i].nbface,
												 amesh.m_level[nlevel].m_geom[dmg[i].cell].boxcenter);
					}
					else if (bdtype == -42)	
					{
						SetBCType_42(amesh.m_level[nlevel].m_data[dmg[i].cell],
												 amesh.m_level[nlevel].m_data[dmg[i].refcell], 
												 dmg[i].nbface,
												 amesh.m_level[nlevel].m_geom[dmg[i].refcell].boxcenter);
					}
					else if (bdtype == -45) SetBCType_45(amesh.m_level[nlevel].m_data[dmg[i].cell], 
																										 amesh.m_level[nlevel].m_data[dmg[i].refcell], dmg[i].nbface);
					else if (bdtype == -12)	SetBCType_12(amesh.m_level[nlevel].m_data[dmg[i].cell], 
																								 amesh.m_level[nlevel].m_data[dmg[i].refcell], dmg[i].nbface);
					else if (bdtype == -13)	SetBCType_13(amesh.m_level[nlevel].m_data[dmg[i].cell], 
																								 		 amesh.m_level[nlevel].m_data[dmg[i].refcell], dmg[i].nbface);
					else if (bdtype == -49) SetBCType_49(amesh.m_level[nlevel].m_data[dmg[i].cell], 
																								 		 amesh.m_level[nlevel].m_data[dmg[i].refcell], dmg[i].nbface);
					else if (bdtype == -10)	SetBCType_10(amesh.m_level[nlevel].m_data[dmg[i].cell], 
																										 amesh.m_level[nlevel].m_data[dmg[i].refcell], dmg[i].nbface);
					else if (bdtype == -125) SetBCType_125_Inputinletbdlayer(amesh.m_level[nlevel].m_data[dmg[i].cell], 
																			 amesh.m_level[nlevel].m_data[dmg[i].refcell], 
																			 dmg[i].nbface,
																			 amesh.m_level[nlevel].m_geom[dmg[i].cell].boxcenter[2],
																			 amesh.m_level[nlevel].m_box[dmg[i].cell].lowpt[2],
																			 amesh);
					else 
					{
						printf("The Boundary condition type %d has not been identified!!!"
							"Please check it\n!!! If it is 10 then two faces with index 10 can "
							"not connect with each other!!!", bdtype); 
						MPI_Abort(MPI_COMM_WORLD, 52);
					}
					// if (node == 1 && (amesh.m_level[nlevel].m_box[dmg[i].cell].ix() == 241 
					// 	             || amesh.m_level[nlevel].m_box[dmg[i].cell].ix() == 242))
					// {
					// 	printf("Bd condition for box (%d,%d,%d) it bc type is %d u is %f "
					// 		"ref is box (%d,%d,%d) u is %f\n", 
					// 		amesh.m_level[nlevel].m_box[dmg[i].cell].ix(),
					// 		amesh.m_level[nlevel].m_box[dmg[i].cell].iy(),
					// 		amesh.m_level[nlevel].m_box[dmg[i].cell].iz(),bdtype,
					// 		amesh.m_level[nlevel].m_data[dmg[i].cell].u,
					// 		amesh.m_level[nlevel].m_box[dmg[i].refcell].ix(),
					// 		amesh.m_level[nlevel].m_box[dmg[i].refcell].iy(),
					// 		amesh.m_level[nlevel].m_box[dmg[i].refcell].iz(),
					// 		amesh.m_level[nlevel].m_data[dmg[i].refcell].u);
					// }
				}
#if DIM == 2
			}
#endif							
			}
#ifdef TEMPORAL_REFINE
			}
#endif						
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}

	/*pressure inlet*/
	void SetBCType_13(FlowVariables & bcvar, FlowVariables & refvar, const int & faceindex)
	{
		double speed0 = speed(refvar);
		bcvar.u = speed0*cos(inlet_flow_angle);
		bcvar.v = speed0*sin(inlet_flow_angle);
		bcvar.w = 0.0;
		//bcvar.roe = refvar.roe;
		bcvar.T = facevar[faceindex].T-
			kineticenergy(bcvar)/Cp;
		bcvar.p = facevar[faceindex].p*pow((bcvar.T/facevar[faceindex].T), 3.5);
		Get_roe_ideal_gas(bcvar);
		Get_E(bcvar);
#ifdef TURBULENCE
		bcvar.viseddy = facevar[faceindex].viseddy;
#endif								
	}
	/*pressure outlet*/
	void SetBCType_12(FlowVariables & bcvar, FlowVariables & refvar, const int & faceindex)
	{		
		bcvar = refvar;
		bcvar.p = facevar[faceindex].p + backpressure;
		Get_roe_ideal_gas(bcvar);
		//refvar.showdata("pressure outlet reference value");
	}

	/*outflow*/
	void SetBCType_50(FlowVariables & bcvar, FlowVariables & refvar, const int & faceindex)
	{
		bcvar = refvar;
	}
	/*no-slipping wall*/
	void SetBCType_49(FlowVariables & bcvar, FlowVariables & refvar, const int & faceindex)
	{
		bcvar.u = 0.0;
		bcvar.v = 0.0;
		bcvar.w = 0.0;
		bcvar.p = refvar.p;
		bcvar.T = refvar.T;
		Get_roe_ideal_gas(bcvar);
		Get_E(bcvar);
#ifdef TURBULENCE
		bcvar.viseddy = 0.0;
#endif								
	}
	/*no-slipping wall for frame motion*/
	void SetBCType_46(FlowVariables & bcvar, FlowVariables & refvar, const int & faceindex, Pointxyz & refboxcenter)
	{
		Pointxyz boxvel = frame_rotate_vel.cross(refboxcenter - frame_rotate_origin) +frame_trans_vel;
		bcvar.u = boxvel[0];
		bcvar.v = boxvel[1];
		bcvar.w = boxvel[2];
		bcvar.p = refvar.p;
		bcvar.T = refvar.T;
		Get_roe_ideal_gas(bcvar);
		Get_E(bcvar);
#ifdef TURBULENCE
		bcvar.viseddy = 0.0;
#endif														
	}
	/*no-slipping wall for a moving wall*/
	void SetBCType_42(FlowVariables & bcvar, FlowVariables & refvar, const int & faceindex, Pointxyz & refboxcenter)
	{
		Pointxyz boxvel = face_force_rotate[faceindex].cross(refboxcenter - face_force_rot_center[faceindex]) + face_force_move[faceindex];
		bcvar.u = boxvel[0];
		bcvar.v = boxvel[1];
		bcvar.w = boxvel[2];
		bcvar.p = refvar.p;
		bcvar.T = refvar.T;
		Get_roe_ideal_gas(bcvar);
		Get_E(bcvar);
#ifdef TURBULENCE
		bcvar.viseddy = 0.0;
#endif														
	}
	/*Slipping-wall*/
	void SetBCType_51(FlowVariables & bcvar, FlowVariables & refvar, const int & faceindex, Pointxyz & facenmv)
	{
		double un = refvar.u*facenmv[0]+refvar.v*facenmv[1]+refvar.w*facenmv[2];
		bcvar.u = refvar.u - un*facenmv[0];
		bcvar.v = refvar.v - un*facenmv[1];
		bcvar.w = refvar.w - un*facenmv[2];
		bcvar.p = refvar.p;
		bcvar.T = refvar.T;
		Get_roe_ideal_gas(bcvar);
		Get_E(bcvar);
#ifdef TURBULENCE
		bcvar.viseddy = refvar.viseddy;
#endif		
	}
	/*fixed values*/
	void SetBCType_10(FlowVariables & bcvar, FlowVariables & refvar, const int & faceindex)
	{
		bcvar = facevar[faceindex];					
	}
	/*Fixed velocity others zero gradient*/
	void SetBCType_45(FlowVariables & bcvar, FlowVariables & refvar, const int & faceindex)
	{
		bcvar.u = facevar[faceindex].u;
		bcvar.v = facevar[faceindex].v;
		bcvar.w = facevar[faceindex].w;
		bcvar.T = facevar[faceindex].T;
		bcvar.p = refvar.p;
		Get_roe_ideal_gas(bcvar);	
		bcvar.e = refvar.e;
#ifdef TURBULENCE
		bcvar.viseddy = facevar[faceindex].viseddy;
#endif						
	}

	void SetBDFaceVars();

	void SetBCType_125_Inputinletbdlayer(FlowVariables & bcvar, FlowVariables & refvar, const int & faceindex, double & z0, int & zindex, Mesh & amesh);

};

#endif