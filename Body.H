#ifndef _Body_H_
#define _Body_H_

#include <iostream>
#include <fstream>
#include <cstdio>

#include "DataArray.H"
#include "Mesh.H"
#include "Cell_iterator.H"
#include "Box.H"
#include "BCValues.H"
#include "Turbulence.H"

#define BODY_SEARCH_NUM 4

#define WALLINTPFLAG_NOTINFECTED

extern const bool moveflag;

const double patch_out_distance = 0.08;

struct PatchNeib
{
	int nbnum;
#ifdef QUA_ELEMENT	
	int nb[15];
#endif
#ifdef TRI_ELEMENT
	int nb[20];
#endif	
	
	PatchNeib()
	{
		nbnum = 0;	
#ifdef QUA_ELEMENT			
		for (int i = 0; i < 15; ++i)
#endif
#ifdef TRI_ELEMENT
		for (int i = 0; i < 20; ++i)
#endif					
		{
			nb[i] = -1;
		}
	}		
};

struct PatchamongNodes
{
	vector<vector<BoxlocChange> > pat_others;
	vector<BoxlocChange> recv_pats;
	int total_exc_num;

	PatchamongNodes()
	{
		pat_others.resize(nodenum);
		total_exc_num = 0;
	}

	void initsendarray()
	{
		for (int i = 0; i < nodenum; ++i)
		{
			pat_others[i].resize(0);
		}
	}

	void initrecvarray()
	{
		recv_pats.resize(0);
	}

	void newsendmeg(int & pat0, int & tonode)
	{
		for (int i = 0; i < nodenum; ++i)
		{
			pat_others[i].push_back(BoxlocChange(pat0, tonode));
		}
	}

	void dataexchange()
	{
		BcastNewInfo_Node(pat_others, recv_pats, 1);
		int recvpatnum = recv_pats.size();
		MPI_Allreduce(&recvpatnum, &total_exc_num, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	}
};

class Body
{
	friend class Pltoutput;
	friend class AMR;
public:

	Body()
	{
		if (nrank == 0) printf("Construct a body!!!\n");
	}

	Body(const int & i0)
	{
		this->bodyindex = i0;
	}

	// ~Body()
	// {
	// 	if (nrank == 0) printf("Body %d was freed!!!\n", bodyindex);
	// }

	inline double & cx()
	{
		return bodycenter[0];
	}
	inline double & cy()
	{
		return bodycenter[1];
	}
	inline double & cz()
	{
		return bodycenter[2];
	}
	inline static double & clength()
	{
		return c_length;
	}

	inline Pointxyz & bodybc()
	{
		return bodycenter;
	}

	void SetBodyCenter(const double & ox, const double & oy, const double & oz)
	{
		this->bodycenter[0] = ox;
		this->bodycenter[1] = oy;
		this->bodycenter[2] = oz;
	}

	void SetBodyOffset(const double & ox, const double & oy, const double & oz)
	{
		this->bcoffset[0] = ox;
		this->bcoffset[1] = oy;
		this->bcoffset[2] = oz;
	}

	void SetBodyScale(const double & sx, const double & sy, const double & sz)
	{
		this->scale[0] = sx;
		this->scale[1] = sy;
		this->scale[2] = sz;
	}

	void ReadShapeFile(string filename)
	{
#ifdef WALLINTPFLAG_NOTSOLID
#ifdef WALLINTPFLAG_NOTINFECTED
		printf("Two wall interpolation flag can not be defined at the same time!!!\n");
		MPI_Abort(MPI_COMM_WORLD, 158);
#endif
#endif
#ifdef WALLINTPFLAG_NOTINFECTED
#ifdef WALLINTPFLAG_NOTSOLID
		printf("Two wall interpolation flag can not be defined at the same time!!!\n");
		MPI_Abort(MPI_COMM_WORLD, 158);
#endif
#endif
#ifndef WALLINTPFLAG_NOTINFECTED
#ifndef WALLINTPFLAG_NOTSOLID
		printf("Please define the wall interpolation flag!!!\n");
		MPI_Abort(MPI_COMM_WORLD, 158);
#endif
#endif				
		ifstream shapefile(filename.c_str());
		if (!shapefile.good())
		{
			printf("File %s does not exist!!!\n", filename.c_str());
			MPI_Abort(MPI_COMM_WORLD,81);
		}
		// shapefile.open(filename)
		int bptnum, belenum;
		shapefile >> bptnum;
		shapefile >> belenum;
		allpoint.setnum_nocopy(bptnum, 0);
		patch.setnum_nocopy(belenum, 0);
		/*Read and scale original point values*/
		if (srank == 0)
		{
			for (int i = 0; i < bptnum; ++i)
			{
				for (int j = 0; j < 3; ++j)
				{
					shapefile >> allpoint[i][j];
					allpoint[i][j] *= scale[j];
					allpoint[i][j] += bcoffset[j];
				}
			}
			/*Read patch connection*/
			for (int i = 0; i < belenum; ++i)
			{
#ifdef QUA_ELEMENT				
				for (int j = 0; j < 4; ++j)
#endif
#ifdef TRI_ELEMENT
				for (int j = 0; j < 3; ++j)
#endif									
				{
					shapefile >> patch[i].corpt[j];
					//printf("patch %d coprt %d is %d\n", i, j, patch[i].corpt[j]);
					if (patch[i].corpt[j] < 1)
					{
						printf("Patch %d corpt %d is %d!!! Please check the body mesh!!!\n",
							i,j,patch[i].corpt[j]);
						MPI_Abort(MPI_COMM_WORLD, 178);
					}
					for (int j0 = 0; j0 < j; ++j0)
					{
						if (patch[i].corpt[j] == patch[i].corpt[j0]+1)
						{
							printf("Patch %d corpt %d and %d is %d!!! Please check the body mesh!!!\n",
								i,j0,j,patch[i].corpt[j]);
							MPI_Abort(MPI_COMM_WORLD, 178);
						}
					}
					patch[i].corpt[j] -= 1;
				}
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
		shapefile.close();	
		CompressPointArray();
		/*Compute point offset with respect to the body center*/
		surfptoff.resize(allpoint.pe()-allpoint.ps());
		ComptSurfptOff();
		/*Compute patch parameters*/
		for (int i = this->patch.ps(); i < this->patch.pe(); ++i)
		{
			ComptPatchCenter(this->patch[i]);
			ComptPatchArea((this->patch[i]));
			patch[i].pc.get_theta();
		}
		MPI_Barrier(share_comm);
#if DIM == 2		
		//SortPatchCenterLocation(0, patch.size()-1, 0);
		//MPI_Barrier(share_comm);
#endif		
		for (int i = this->patch.ps(); i < this->patch.pe(); ++i)
		{
			ComptPatchNv_Init(this->patch[i], i);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		double ttf = surfacearea();
		//if (nrank == 0) printf("The body total surface area is %f\n", ttf);
		InitVelZero();		
		// Get_Space_Limit();
		CollectPatchNeib();
		if (nrank == 0) printf("Finish ReadShapeFile for body %d\n", bodyindex);
	}

	void CompressPointArray();

	void ComptSurfptOff()
	{
		for (int i = 0; i < surfptoff.size(); ++i)
		{
			int pt0 = i + allpoint.ps();
			surfptoff[i][0] = allpoint[pt0][0]-bodycenter[0];
			surfptoff[i][1] = allpoint[pt0][1]-bodycenter[1];
			surfptoff[i][2] = allpoint[pt0][2]-bodycenter[2];
		}
		MPI_Barrier(share_comm);
	}

	void CopyBody(Body & abody)
	{
		this->surfptoff.resize(abody.surfptoff.size());
		/*Copy the point offset to the new body*/
		for (int i = 0; i < abody.surfptoff.size(); ++i)
		{
			surfptoff[i] = abody.surfptoff[i];
		}
		/*Copy the topology of connection to the new body*/
		this->allpoint.setnum_nocopy(abody.allpoint.size(), 0);
		this->patch.setnum_nocopy(abody.patch.size(), 0);
		for (int i = patch.ps(); i < patch.pe(); ++i)
		{
#ifdef QUA_ELEMENT			
			for (int j = 0; j < 4; ++j)
#endif
#ifdef TRI_ELEMENT				
			for (int j = 0; j < 3; ++j)	
#endif				
			{
				this->patch[i].corpt[j] = abody.patch[i].corpt[j];
			}
		}
		/*Compute the point values of the new body with respect to the body center*/
		for (int i = this->allpoint.ps(); i < this->allpoint.pe(); ++i)
		{
			int i0 = i-this->allpoint.ps();
			for (int j = 0; j < 3; ++j)
			{
				this->allpoint[i][j] = this->bodycenter[j]+surfptoff[i-i0][j];
			}
		}
		/*Compute patch parameters*/
		for (int i = this->patch.ps(); i < this->patch.pe(); ++i)
		{
			ComptPatchCenter(this->patch[i]);
			ComptPatchArea((this->patch[i]));
			patch[i].pc.get_theta();
		}
		MPI_Barrier(share_comm);
		for (int i = this->patch.ps(); i < this->patch.pe(); ++i)
		{
			ComptPatchNv_Init(this->patch[i], i);
		}
		InitVelZero();
		// Get_Space_Limit();
		CollectPatchNeib();
	}
	void ComptPatchCenter(Surfpatch & apatch)
	{
		/*compute the center values*/
		apatch.pc[0] = 0.0;
		apatch.pc[1] = 0.0;
		apatch.pc[2] = 0.0;
		for (int j = 0; j < 3; ++j)
		{
#ifdef QUA_ELEMENT		
			for (int i = 0; i < 4; ++i)
			{
				apatch.pc[j] += allpoint[apatch.corpt[i]][j]/4.0;
			}
#endif
#ifdef TRI_ELEMENT
			for (int i = 0; i < 3; ++i)
			{
				apatch.pc[j] += allpoint[apatch.corpt[i]][j]/3.0;
			}
#endif
		}
		// double ptlen = apatch.pc.length();
		// if (ptlen < 0.000001)
		// {
		// 		printf("The point coordinate is close to zero!!! Point is (%f,%f,%f)\n",
		// 			apatch.pc[0], apatch.pc[1], apatch.pc[2]);
		// 		MPI_Abort(MPI_COMM_WORLD, 205);
		// }
	}

	void ComptPatchNv(Surfpatch & apatch);

	// inline void ComptPatchNv_Init(Surfpatch & apatch, const int & pat0)
	// {
	// 	Pointxyz n1 = (allpoint[apatch.corpt[1]] - allpoint[apatch.corpt[0]]);
	// 	Pointxyz n2 = (allpoint[apatch.corpt[2]] - allpoint[apatch.corpt[0]]);
	// 	apatch.nv = n1.cross(n2);
	// 	apatch.nv.normalize();
	// 	int patsize = patch.size();
	// 	for (int i = 0; i < patsize; ++i)
	// 	{
	// 		Pointxyz stoc = (patch[i].pc - apatch.pc);
	// 		if (stoc.length()*0.8 < stoc.dot(apatch.nv)-0.000001)
	// 		{
	// 			apatch.reverse_nv_dir = true;
	// 			apatch.nv *= -1.0;
	// 			break;
	// 		}
	// 	}
	// }

	inline void ComptPatchNv_Init(Surfpatch & apatch, const int & pat0)
	{
		Pointxyz n1 = (allpoint[apatch.corpt[1]] - allpoint[apatch.corpt[0]]);
		Pointxyz n2 = (allpoint[apatch.corpt[2]] - allpoint[apatch.corpt[0]]);
		apatch.nv = n1.cross(n2);
		apatch.nv.normalize();
#ifndef PASSAGE_ANGLE		
		if (apatch.nv.dot(apatch.pc - bodycenter) < 0.0)
		{
			apatch.reverse_nv_dir = true;
			apatch.nv *= -1.0;
		}
#endif		
	}

#ifdef QUA_ELEMENT
	inline void ComptPatchArea(Surfpatch & apatch)
	{
		Pointxyz n1 = allpoint[apatch.corpt[1]] - allpoint[apatch.corpt[0]];
		Pointxyz n2 = allpoint[apatch.corpt[3]] - allpoint[apatch.corpt[0]];
		Pointxyz n3 = allpoint[apatch.corpt[1]] - allpoint[apatch.corpt[2]];
		Pointxyz n4 = allpoint[apatch.corpt[3]] - allpoint[apatch.corpt[2]];
		Pointxyz area1 = n1.cross(n2);
		Pointxyz area2 = n3.cross(n4);
		apatch.area = (area1.length()+area2.length())*0.5;
	}
#endif
#ifdef TRI_ELEMENT
	inline void ComptPatchArea(Surfpatch & apatch)
	{
		Pointxyz n1 = allpoint[apatch.corpt[1]] - allpoint[apatch.corpt[0]];
		Pointxyz n2 = allpoint[apatch.corpt[2]] - allpoint[apatch.corpt[0]];
		Pointxyz area1 = n1.cross(n2);
		apatch.area = area1.length()*0.5;
	}
#endif		

	// void EvaluateMaxPatDis(const int & bid)
	// {
	// 	if (nrank == 0)
	// 	{
	// 		double maxdis = 0.0;
	// 		for (int i = 1; i < patch.pe(); ++i)
	// 		{
	// 			double ptd = patch[i].patdis(patch[i-1]);
	// 			if (ptd > maxdis)
	// 			{
	// 				maxdis = ptd;
	// 			}
	// 		}
	// 		printf("Max patch distance is %e\n", maxdis);
	// 		if (maxdis > dh[lastlevel][0] ||
	// 				maxdis > dh[lastlevel][1] ||
	// 				maxdis > dh[lastlevel][2])
	// 		{
	// 			printf("The max distance of body %d is coarser then the mesh!!!Please improve the body!!!\n", bid);
	// 			MPI_Abort(MPI_COMM_WORLD,242);
	// 		}
	// 	}
	// }

	void ComptPatchParams();

	void CollectPatchNeib();

	void Attachbox_Init(Mesh & amesh);

	void Attachbox_NewLevel(Mesh & amesh, DataArray<rftag> & mtag, DataArray<Boxson<int> > & mson);

	static void Renew_Attachbox(vector<Body> & abody, Mesh & amesh);

	static void Check_Attachbox_willberemoved(vector<Body> & abody, Mesh & amesh);

	void CheckAttachBoxTag(DataArray<rftag> & mtag, const int & ilevel, Mesh & amesh)
	{
#ifndef LOCALBODY		
		int pps = patch.ps();
		int ppe = patch.pe();
		for (int i = pps; i < ppe; ++i)
		{
			double ptlen = patch[i].pc.length();
			int i0 = i - pps;
			int pat0 = patchbox[i0];
			if (mtag[pat0].tag == -1 && amesh.m_level[ilevel].m_box.isnormal(pat0))
			{
				printf("Patch (%f,%f,%f) is close to level %d box (%d,%d,%d) center (%f,%f,%f) box type is %d isnormal %d signdis is %f the box should be tagged!!!\n",
					patch[i].pc[0],patch[i].pc[1],patch[i].pc[2],ilevel,
					amesh.m_level[ilevel].m_box[pat0].ix(),
					amesh.m_level[ilevel].m_box[pat0].iy(),
					amesh.m_level[ilevel].m_box[pat0].iz(),
					amesh.m_level[ilevel].m_geom[pat0].boxcenter[0],
					amesh.m_level[ilevel].m_geom[pat0].boxcenter[1],
					amesh.m_level[ilevel].m_geom[pat0].boxcenter[2],
					amesh.m_level[ilevel].m_box[pat0].type,
					amesh.m_level[ilevel].m_box.isnormal(pat0),
					amesh.m_level[ilevel].m_box[pat0].pair.signdis);
				MPI_Abort(MPI_COMM_WORLD, 497);
			}
		}
#endif		
	}

	void CheckBodyPatchCenter();

	void CycleNearCell(Mesh & amesh, int & atbox, double & atdis, const int & pat0)
	{
		int oldbox = -1;
		while (oldbox != atbox)
		{
			oldbox = atbox;
			for (Point_iterator p(0,3); p.end(); ++p)
			{
				int myneib = amesh.m_level[lastlevel].m_box[oldbox].neib[p.i][p.j][p.k];					
				//if (amesh.m_level[lastlevel].m_box.isnormal(myneib))
				if (myneib > -1)
				{						
					double dis = PatchDistoPoint(amesh.m_level[lastlevel].m_geom[myneib].boxcenter, pat0);
					if (dis < atdis && amesh.m_level[lastlevel].m_box[myneib].type != Dmghost)
					{
						atbox = myneib;
						atdis = dis;
					}
				}				
			}
		}
	}

	void ModifyBoxtoPatch(Mesh & amesh);

	void ModifyBoxtoPatch_Somepatches(Mesh & amesh, vector<Boxloc> & bdloc, const int & locnum, vector<int> & bdloc_tag);

	static void DistanceSign_boxtowall(Mesh & amesh, vector<Body> & abody, const int & ilevel);

	void FindClosePatchtoCell(BoxtoWall & abw, Pointxyz & cellcenter)
	{
		PatchDistoBox(abw, cellcenter, abw.patch);
		int oldpatch = -1;
		int ps0 = abw.patch;
		Pointxyz dxyz;
		while(oldpatch != ps0)
		{
			oldpatch = abw.patch;
			int nbnum = patneib[ps0].nbnum;
#ifdef QUA_ELEMENT			
			if (nbnum > 15)
			{
				printf("Patch %d (%f,%f,%f) has %d neibs...\n", ps0, patch[ps0].pc[0], patch[ps0].pc[1], patch[ps0].pc[2], nbnum);
			}
#endif
#ifdef TRI_ELEMENT				
			Assert(nbnum <= 20, "The patch neib number must in range!!!", 887);
#endif				
			for (int n0 = 0; n0 < nbnum; ++n0)
			{
				int pat0 = patneib[ps0].nb[n0];
				Assert(pat0 > -1, "The patch neib index error!!!", 496);
				PatchDistoBox(abw, cellcenter, pat0);
			}
			ps0 = abw.patch;
		}
#ifndef PASSAGE_ANGLE		
		dxyz = cellcenter - patch[abw.patch].pc;
		PeriodicLength(dxyz);
#else
		Pointxyz newbcxyz;
		PeriodicAnnulaLength(cellcenter, patch[abw.patch].pc, newbcxyz);
		dxyz = newbcxyz - patch[abw.patch].pc;
#endif		
		abw.signdis = patch[abw.patch].nv.dot(dxyz);
		abw.distance_to_body = abw.signdis;
	}

	void FindClosePatchtoCell_CycleAllPatches(BoxtoWall & abw, Pointxyz & cellcenter)
	{
		int patnum = patch.size();
		for (int i = 0; i < patnum; ++i)
		{
			PatchDistoBox(abw, cellcenter, i);
		}
#ifndef PASSAGE_ANGLE		
		Pointxyz dxyz = cellcenter - patch[abw.patch].pc;
		PeriodicLength(dxyz);
#else
		Pointxyz newbcxyz;
		PeriodicAnnulaLength(cellcenter, patch[abw.patch].pc, newbcxyz);
		Pointxyz dxyz = newbcxyz - patch[abw.patch].pc;
#endif
		abw.signdis = patch[abw.patch].nv.dot(dxyz);
		abw.distance_to_body = abw.signdis;
	}

	void PatchDistoBox(BoxtoWall & abw, Pointxyz & cellcenter, const int & pat0)
	{
#ifndef PASSAGE_ANGLE		
		Pointxyz dxyz = cellcenter - patch[pat0].pc;
		PeriodicLength(dxyz);
#else
		Pointxyz newbcxyz;
		PeriodicAnnulaLength(cellcenter, patch[pat0].pc, newbcxyz);
		Pointxyz dxyz = newbcxyz - patch[pat0].pc;	
#endif				
		double mydis = dxyz.length();
		if (mydis < abw.distance)
		{
			abw.distance = mydis;
			abw.body = bodyindex;
			abw.patch = pat0;
		}
	}

	double PatchDistoPoint(Pointxyz & cellcenter, const int & pat0)
	{
#ifndef PASSAGE_ANGLE		
		Pointxyz dxyz = cellcenter - patch[pat0].pc;
		PeriodicLength(dxyz);
#else
		Pointxyz newbcxyz;
		PeriodicAnnulaLength(cellcenter, patch[pat0].pc, newbcxyz);
		Pointxyz dxyz = newbcxyz - patch[pat0].pc;
#endif				
		return dxyz.length();
	}

	inline Pointxyz & bc()
	{
		return bodycenter;
	}

	void Printout(FILE * a_file)
	{
		fprintf(a_file, "zone T = Body%d\n", bodyindex);
#ifdef QUA_ELEMENT		
		fprintf(a_file, "N=%d, E=%d, F=fepoint, ET=quadrilateral\n", allpoint.size(), patch.size());
#endif
#ifdef TRI_ELEMENT		
		fprintf(a_file, "N=%d, E=%d, F=fepoint, ET=triangle\n", allpoint.size(), patch.size());
#endif
		for (int i = 0; i < allpoint.size(); ++i)
		{
			fprintf(a_file, "%16.8f%16.8f%16.8f\n", allpoint[i][0], allpoint[i][1], allpoint[i][2]);
		}
		for (int i = 0; i < patch.size(); ++i)
		{
#ifdef QUA_ELEMENT			
			fprintf(a_file, "%d %d %d %d\n", patch[i].corpt[0]+1, patch[i].corpt[1]+1, patch[i].corpt[2]+1, patch[i].corpt[3]+1);
#endif
#ifdef TRI_ELEMENT
			fprintf(a_file, "%d %d %d\n", patch[i].corpt[0]+1, patch[i].corpt[1]+1, patch[i].corpt[2]+1);
#endif					
		}
	}

	void Printout_soln(FILE * a_file, const int & varnum)
	{
		fprintf(a_file, "zone T = Body%d\n", bodyindex);
#ifdef QUA_ELEMENT		
		fprintf(a_file, "N=%d, E=%d, F=fepoint, ET=quadrilateral\n", allpoint.size(), patch.size());
#endif
#ifdef TRI_ELEMENT		
		fprintf(a_file, "N=%d, E=%d, F=fepoint, ET=triangle\n", allpoint.size(), patch.size());
#endif
		for (int i = 0; i < allpoint.size(); ++i)
		{
			for (int j = 0; j < varnum; ++j)
			{
				fprintf(a_file, " 0");
			}
			fprintf(a_file, "\n");
		}
	}

	//static void MovingBody(Mesh & amesh, vector<Body> & abody);

	// static void ImmerseWallBoundary(Mesh & amesh, vector<Body> & abody)
	// {
	// 	lastlevel = amesh.MyCurNum()-1;
	// 	DataArray<Box> & boxarray = amesh.ShowLevelBox(lastlevel);
	// 	amesh.InitWallInfect();
	// 	GiveAFlag("Finish InitWallInfect!!!", 5);
	// 	for (int i = 0; i < abody.size(); ++i)
	// 	{
	// 		abody[i].ModifyBoxtoPatch(amesh);
	// 		GiveAFlag("Finish ModifyBoxtoPatch!!!", 5);
	// 		abody[i].ModifyNearWallCellDistance(amesh);
	// 		GiveAFlag("Finish ModifyNearWallCellDistance!!!", 5);
	// 	}
	// 	DistanceSign(amesh, abody);
	// 	GiveAFlag("Finish DistanceSign!!!", 5);
	// 	CheckBdCells(amesh);
	// 	GiveAFlag("Finish CheckBdCells!!!", 5);
	// 	amesh.DataExchangeCells_Wall();
	// }

	static void ImmerseWallBoundary(Mesh & amesh, vector<Body> & abody)
	{
		Body::InfectWallBox_Init_Negative(amesh, abody);
		GiveAFlag("Finish infect box first time...!!!", 5);
		for (int i = 0; i < abody.size(); ++i)
		{					
			abody[i].InitPatchsNode(amesh);
			GiveAFlag("Finish initialize the patches node!!!", 5);
			abody[i].ModifyBoxtoPatch(amesh);
			GiveAFlag("Finish ModifyBoxtoPatch!!!", 5);
			abody[i].FindPatchHangingNode(amesh);
			GiveAFlag("Finish FindPatchHangingNode!!!", 5);
		}
		Body::FindImageCell(amesh, abody);
		GiveAFlag("Finish FindImageCell!!!", 5);
		amesh.DataExchangeCells_Wall();
		GiveAFlag("Finish finding the data exchange cells near wall...", 5);
		Body::FindImageCell_GhostIBCell(amesh, abody);
		GiveAFlag("Finish finding the interpolation node for the exchanged ib cells...", 5);
		amesh.ComptDistancetoDomain();
		GiveAFlag("Finish ComptDistancetoDomain...", 5);
	}

	void FindPatchHangingNode(Mesh & amesh);

	void FindHangingNode_Onepatch(Mesh & amesh, const int & i);

	static void SetBodyParams(vector<Body> & abody, Mesh & amesh);

	static void SetBodyParams(vector<Body> & abody, Mesh & amesh, vector<Body> & localbody);

	void InitVelZero()
	{
		force_tra_vel.zero();
		force_rot_vel.zero();
		free_tra_vel.zero();
		free_rot_vel.zero();
		inlineosc.setvalue(0.0,0.0,0.0);
		inlineosc_new.setvalue(0.0,0.0,0.0);
		rotosc.setvalue(0.0,0.0,0.0);
		rotosc_new.setvalue(0.0,0.0,0.0);
		angle_location.setvalue(0.0, 0.0, 0.0);
	}

  	inline void Surface_Point_Velocity(Pointxyz & apt, Pointxyz & avel)
  	{
			Pointxyz r0 = apt - frame_rotate_origin;
			Pointxyz r1 = apt - bodycenter;
			avel = free_tra_vel + force_tra_vel + frame_trans_vel;
#ifdef PASSAGE_ANGLE					
			avel += (force_rot_vel+frame_rotate_vel).cross(r0);
			avel += free_rot_vel.cross(r1);
#else
			avel += frame_rotate_vel.cross(r0);
			avel += (force_rot_vel+free_rot_vel).cross(r1);
#endif						
  	}

  	static void IBConditions_Image(Mesh & amesh, vector<Body> & abody);

  	static void FindImageCell(Mesh & amesh, vector<Body> & abody);

  	static void FindImageCell_GhostIBCell(Mesh & amesh, vector<Body> & abody);

  	static int HGCellVars(Mesh & amesh, 
  												const int & ibbc, 
  												Pointxyz & hgpt, 
  												FlowVariables & hgfv,
  												Pointxyz & nmv,
  												const Pointxyz & atp = Pointxyz(0.0, 0.0, 0.0),
  												const double & hl = 0.0);

  	static void LinearDist(FlowVariables & hgfv, 
  												 IBCell & ibpt, 
													 FlowVariables & ibfv_n0, 
													 Pointxyz & wallvel, 
													 Pointxyz & wallnv,
													 double & ibsigndis);

  	void ComptWallForce(Mesh & amesh);

  	static void ComptIBHglength(Mesh & amesh, const int & ilevel, const int & ibox, double & hglength, Pointxyz & patnv);

  	static void InitIBCellParams(Mesh & amesh, 
  															 vector<Body> & abody);

  	static void ComptDisFromBoxtoBody(Mesh & amesh, vector<Body> & abody, const int & ilevel);

  	static void RenewLevelDistance(Mesh & amesh, vector<Body> & abody, const int & level0,
  		vector<int> & move_from_in_to_out,
		vector<int> & move_from_out_to_in);

  	static int & Attachlevel()
  	{
  		return lastlevel;
  	}

  	static int & Bodynum()
  	{
  		return bodynum;
  	}

  	inline double & draglift(const int & i0)
  	{
  		return force[i0];
  	}

  	inline double & shearforce(const int & i0)
  	{
  		return tau_force[i0];
  	}

  	inline double & pressureforce(const int & i0)
  	{
  		return pre_force[i0];
  	}

  	inline double & bodymoment(const int & i0)
  	{
  		return moment[i0];
  	}

  	// void Get_Space_Limit();

  	inline DataArray<Surfpatch> & bodypatch()
  	{
  		return patch;
  	}

  	double surfacearea()
  	{
  		int s0 = patch.ps();
  		int se = patch.pe();
  		double ttf = 0.0;
  		for (int i = s0; i < se; ++i)
  		{
  			ttf += patch[i].area;
  		}
  		MPI_Allreduce(MPI_IN_PLACE, &ttf, 1, MPI_DOUBLE, MPI_SUM, share_comm);
  		return ttf;
  	}

  	//void IntpHgUt(Mesh & amesh, const int & lastlevel, const int & patindex, Pointxyz & patv, double & patchut);

  	void FindNearBody()
  	{
  		int bodysign[BODY_SEARCH_NUM] = {-1,1, -2, 2};
  		for (int i = 0; i < BODY_SEARCH_NUM; ++i)
  		{
  			search_near_body[i] = -1;
  			int ns0 = bodyindex + bodysign[i];
  			if (periodic[1])
  			{
  				if (ns0 < 0) ns0 += bodynum;
  				else if (ns0 > bodynum - 1) ns0 -= bodynum;
  			}
  			if (ns0 > -1 && ns0 < bodynum)
  			{
  				search_near_body[i] = ns0;
  			}
  		}
  	}

  	void ComptNewLocation();

  	void Motion_Rule(Mesh & amesh, double & dt00);

  	void Rotate_Axis_X(double & rot_angle, Pointxyz rot_bc = Pointxyz (0.0,0.0,0.0));

  	void Rotate_Axis_Z(double & rot_angle, Pointxyz rot_bc = Pointxyz (0.0,0.0,0.0));

  	inline Pointxyz & moveoffset()
  	{
  		return inlineosc;
  	}

  	inline Pointxyz & Body_Low_Limit()
  	{
  		return xyz_min;
  	}

  	inline Pointxyz & Body_High_Limit()
  	{
  		return xyz_max;
  	}

  	static void ComptDisFromBoxtoBody_Init(Mesh & amesh, vector<Body> & abody);

  	static void ComptDisFromBoxtoBody_Newlevel(Mesh & amesh, DataArray<rftag> & mtag, DataArray<Boxson<int> > & mson, vector<Body> & abody);

  	static void InfectWallBox_Init_Negative(Mesh & amesh, vector<Body> & abody);

  	static void MovingBody_Image(Mesh & amesh, vector<Body> & abody, double & dt00);

  	//static void RenewIBCellsDistance(Mesh & amesh, vector<Body> & abody, vector<double> & celldis);

  	static void RemoveIBCells(Mesh & amesh, vector<Body> & abody, vector<int> & move_from_in_to_out,
		vector<int> & move_from_out_to_in);

  	static void InfectNewIBCells(Mesh & amesh, vector<Body> & abody, vector<int> & move_from_in_to_out,
		vector<int> & move_from_out_to_in);

  	void InitPatchsNode(Mesh & amesh);

  	static void FindWallfaces(vector<Body> & abody, Mesh & amesh);

  	static void ModifyWallfaceflux(vector<Body> & abody, Mesh & amesh);

  	void set_pointnum(const int & ptn)
  	{
  		allpoint.setnum_nocopy(ptn,0);
  	}

  	void set_patchnum(const int & ptn)
  	{
  		patch.setnum_nocopy(ptn,0);
  		patneib.setnum_nocopy(ptn,0);
  	}

  	void set_bodyindex(const int & i0)
  	{
  		bodyindex = i0;
  	}

  	void SplitBodytoNode(Mesh & amesh, Body & localbody);

  	void SortPatchCenterLocation(const int & left, const int & right, const int & dir);

private:
	int bodyindex;
	DataArray<Pointxyz>  allpoint;
	DataArray<Surfpatch> patch;
	DataArray<PatchNeib> patneib;

	/*Patch parameters*/
	vector<Pointxyz>     surfptoff;
	/*The following two arraies must be modified after amr*/
	vector<int>			 patchbox;
	vector<double>   pdistobox;
	/*Body Parameters*/
	Pointxyz			 bodycenter;
	Pointxyz       		 bcoffset;
	Pointxyz			 inlineosc;
	Pointxyz			 inlineosc_new;

	Pointxyz			 rotosc;
	Pointxyz			 rotosc_new;

	Pointxyz			 scale;

	/*motion velocity*/
	Pointxyz			 force_tra_vel;
	Pointxyz			 force_rot_vel;
	Pointxyz			 free_tra_vel;
	Pointxyz			 free_rot_vel;

	Pointxyz  			 angle_location;
	Pointxyz 			 force;
	Pointxyz 			 moment;
	Pointxyz             pre_force;
	Pointxyz 			 tau_force;
	Pointxyz			 pre_moment;
	Pointxyz			 tau_moment;

	Pointxyz xyz_min;
	Pointxyz xyz_max;

	double mass;
	double damping;
	double fn;
	double stiff;
	double inertia;
	double rotdamping;
	double rotstiff;

	int search_near_body[BODY_SEARCH_NUM];

	static double        c_length;
	static double        hg_ratio;
	static double        force_hg_ratio[2];
	static int 					 lastlevel;
	static int 					 bodynum;
	/*The near wall box and its distance to the wall*/
	
};



#endif
