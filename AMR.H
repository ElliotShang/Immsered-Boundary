#ifndef _AMR_H
#define _AMR_H
#include <vector>
#include <iostream>
#include <string>


#include "AMRLevel.H"
#include "Mesh.H"
#include "Meshlevel.H"
#include "Check.H"
//#include "Pltoutput.H"
#include "AMRmpi.H"
#include "AMRSpaceTime.H"
#include "TwoLevelOp.H"
#include "MeshOp.H"
#include "Body.H"
#include "Taskdist.H"

using namespace std;

extern bool adpflag;

class AMR
{

public:

	AMR(Mesh & a_mesh):
	m_mesh(a_mesh)
	{
		m_level.reserve(m_mesh.max_level_num);
		CreateAMRLevel();
		flowdataexist = false;
		m_mesh.CheckFaceBoxIndex();
	}

	void CreateAMRLevel()
	{
		int levelindex = m_level.size();
		if (levelindex == 0)
		{
			m_level.push_back(AMRLevel(m_mesh.m_level[levelindex].m_box, 
																 m_mesh.m_level[levelindex].blockpair,
																 m_mesh.m_level[levelindex].dghost,
																 levelindex));
			m_level.back().coarselevel = NULL;
		}
		else if (levelindex > 0)
		{
			m_mesh.AddMeshlevel();
			m_level.push_back(AMRLevel(m_mesh.m_level[levelindex].m_box, 
																 m_mesh.m_level[levelindex].blockpair,
																 m_mesh.m_level[levelindex].dghost,
																 levelindex));
			m_level.back().coarselevel = &m_level[levelindex-1];
			m_level[levelindex-1].finelevel = &m_level[levelindex];
		}
		Assert(m_level.size() == levelindex+1, "Mesh and AMR level connect error!!!", 33);
		ShowAllRankData("Current Level number", levelindex, 2);
		m_level.back().twodflag = level_twod_flag[levelindex];
		m_level.back().power_ratio = level_power_ratio[levelindex];
	}

	void ADPRule(vector<Body> & abody);

	void CheckRemoveorNot(vector<Body> & abody)
	{
		if (m_level.back().m_box.size() == 0)
		{
			m_level.back().clearlevel();
			m_level.pop_back();
			m_level.back().finelevel = NULL;
			m_mesh.RemoveLastMeshlevel();
			if (abody.size() > 0)
			{
				printf("Remove the internal layer but there are bodies...\n");
				MPI_Abort(MPI_COMM_WORLD, 78);
			}
		}
	}

	void CheckLevelRefineFlag()
	{
		morelevel = false;
		if (m_level.size() != m_mesh.cur_level_num)
		{
			MPI_Abort(MPI_COMM_WORLD, 64);
		}
		if (m_level[0].toderefine)
		{
			m_level[0].toderefine = false;
		}
		if (m_level.back().torefine && m_mesh.Canberefined())
		{
			CreateAMRLevel();
			m_level.back().toderefine = true;
			morelevel = true;
		}
		else
		{
			m_level.back().torefine = false;
		}
		if (init_adjust_level > m_mesh.cur_level_num - 1)
		{
			printf("The initial adjust level is larger than the current number of level!!!\n");
			MPI_Abort(MPI_COMM_WORLD, 86);
		}
		if (m_level[init_adjust_level].toderefine)
		{
			printf("The initial adjust level %d can not be derefined!!!\n", init_adjust_level);
			MPI_Abort(MPI_COMM_WORLD, 90);
		}
	}

	void AllocateDataMem()
	{
		for (int i = 0; i < m_mesh.cur_level_num; ++i)
		{
			int bn = m_level[i].m_box.realsize();
			m_mesh.m_level[i].m_data.setnum_nocopy(bn, m_level[i].m_box.arrayghostnum());
		}		
		Flowexist();
	}

	void SetRefineFlag();

	void SetupNewAMR(vector<Body> & abody)
	{
		for (int i = 0; i < init_level_num-1; ++i)
		{
			PerformAdaptive(abody);
		}
		if (init_level_num == 1)
		{
			m_mesh.DataExchangeCells(0);
			FindModifyFace(m_mesh.mdyface);
		}
	}


	void TagtheBox(const int & ilevel, vector<Body> & abody);

	void PerformAdaptive(vector<Body> & abody);

/*Write in 2020-1-2*/
	template <class T>	
	void ComptDatainNewBox(const int & ilevel, vector<T> & coarsedata, vector<Boxson<T> > & finedata)
	{	
		int coarseboxnum = m_level[ilevel].dfnum_procs[srank];
		int fineboxsonnum = m_level[ilevel].rfnum_procs[srank];
		//printf("Level %d refine data %d derefine data %d\n", ilevel, fineboxsonnum, coarseboxnum);
		if (coarseboxnum > 0)
		{
			//vector<T>(0).swap(coarsedata);
			coarsedata.resize(coarseboxnum);
			if (!m_level[ilevel].twodflag)
			{
				for (int j = 0; j < coarseboxnum; ++j)
				{				
					ComptAdjointData(ilevel, m_level[ilevel].derfbox[j], coarsedata[j]);
				}
			}
			else
			{
				for (int j = 0; j < coarseboxnum; ++j)
				{				
					ComptAdjointData_2d(ilevel, m_level[ilevel].derfbox[j], coarsedata[j]);
				}
			}
		}
		GiveLevelFlag("Finish compute coarse data!!!", ilevel, 5);
		if (fineboxsonnum > 0)
		{
			//vector<Boxson<T> >(0).swap(finedata);
			finedata.resize(fineboxsonnum);
			if (!m_level[ilevel+1].twodflag)
			{
				for (int j = 0; j < fineboxsonnum; ++j)
				{
					Assert(m_level[ilevel].m_box.isnorg(m_level[ilevel].rfbox[j]), "The refined box is not in range!!!", 162);
					ComptAdjointData(ilevel, m_level[ilevel].rfbox[j], finedata[j]);
				}
			}
			else
			{
				for (int j = 0; j < fineboxsonnum; ++j)
				{
					Assert(m_level[ilevel].m_box.isnorg(m_level[ilevel].rfbox[j]), "The refined box is not in range!!!", 162);
					ComptAdjointData_2d(ilevel, m_level[ilevel].rfbox[j], finedata[j]);
				}
			}
		}
		GiveLevelFlag("Finish compute fine data!!!", ilevel, 5);
		MPI_Barrier(MPI_COMM_WORLD);
	}
// /*Write in 2020-1-2*/	
	inline void ComptAdjointData(const int & sonlevel, const int & sonid, FlowVariables & newcoarsedata)
	{
		Restriction(m_mesh, sonlevel, sonid, newcoarsedata, tlop);
	}	
	inline void ComptAdjointData(const int & clevel, const int & cindex, Boxson<FlowVariables> & newfinedata)
	{
		if (m_mesh.m_level[clevel].m_box[cindex].type == Normalcell) 
		{
			Prolongation(m_mesh, clevel, cindex, newfinedata, tlop);
		}
		else
		{
			for (Point_iterator p210(0,2); p210.end(); ++p210)
			{
				newfinedata.son[p210.i][p210.j][p210.k] = m_mesh.m_level[clevel].m_data[cindex];
			}
		}
	}
	inline void ComptAdjointData_2d(const int & sonlevel, const int & sonid, FlowVariables & newcoarsedata)
	{
		Restriction_2d(m_mesh, sonlevel, sonid, newcoarsedata, tlop);
	}	
	inline void ComptAdjointData_2d(const int & clevel, const int & cindex, Boxson<FlowVariables> & newfinedata)
	{
		if (m_mesh.m_level[clevel].m_box[cindex].type == Normalcell) 
		{	
			Prolongation_2d(m_mesh, clevel, cindex, newfinedata, tlop);
		}
		else
		{
			for (Point_iterator_2d p228(0,2); p228.end(); ++p228)
			{
				newfinedata.son[p228.i][p228.j][p228.k] = m_mesh.m_level[clevel].m_data[cindex];
			}
		}
	}
	/*compute the box geom parameters must consider the face parameters*/
	/*in the future work, when the mesh is not square, this part must be included*/
	void ComptAdjointData(const int & flevel, const int & fbox, BoxCellGeom & newcoarsedata);

	void ComptAdjointData(const int & rflevel, const int & rfbox0, Boxson<BoxCellGeom> & newfinedata);

	void ComptAdjointData_2d(const int & flevel, const int & fbox, BoxCellGeom & newcoarsedata);

	void ComptAdjointData_2d(const int & rflevel, const int & rfbox0, Boxson<BoxCellGeom> & newfinedata);
/*------------------------------------------------------------------------------------------------------------------------*/
/*Write in 2020-1-9*/
	void PutDataBoxintoArray(const int & i)
	{			
		GiveAFlag("Finish SwitchDataLocation!!!", 5);
		if (NULL != m_level[i].finelevel)
		{
#ifndef IMPORT_MESH			
			IncludeNewFineData(i, m_mesh.m_level[i+1].m_geom, m_level[i].finegeom);
#endif			
			if (flowdataexist)
			{
				IncludeNewFineData(i, m_mesh.m_level[i+1].m_data, m_level[i].finedata);
			}				
		}
		if (NULL != m_level[i].coarselevel)
		{
#ifndef IMPORT_MESH			
			IncludeNewCoarseData(i, m_mesh.m_level[i-1].m_geom, m_level[i].coarsegeom);
#endif			
			if (flowdataexist)
			{
				IncludeNewCoarseData(i, m_mesh.m_level[i-1].m_data, m_level[i].coarsedata);
			}					
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}

	void ComputeNewcellGeom(const int & ilevel)
	{
		int flevel = ilevel+1;
		int clevel = ilevel-1;
		if (!level_twod_flag[flevel])
		{
			for (int j = 0; j < m_level[ilevel].rfnum_procs[srank]; ++j)
			{
				int i0 = m_level[ilevel].rfbox[j];
				int tag0 = m_level[ilevel].m_tag[i0].tag;
				for (Point_iterator p(0,2); p.end(); ++p)
				{
					int sonid = m_level[ilevel].allson[tag0].son[p.i][p.j][p.k];
					ComptCellParams_Pts_AMR(flevel,sonid);
				}
			}
		}
		else
		{
			for (int j = 0; j < m_level[ilevel].rfnum_procs[srank]; ++j)
			{
				int i0 = m_level[ilevel].rfbox[j];
				int tag0 = m_level[ilevel].m_tag[i0].tag;
				for (Point_iterator_2d p(0,2); p.end(); ++p)
				{
					int sonid = m_level[ilevel].allson[tag0].son[p.i][p.j][p.k];
					ComptCellParams_Pts_AMR(flevel,sonid);
				}
			}
		}
		for (int j = 0; j < m_level[ilevel].dfnum_procs[srank]; ++j)
		{
			int dfid = m_level[ilevel].derfbox[j];
			int momid = m_level[ilevel].m_tag[dfid].detag;
			ComptCellParams_Pts_AMR(clevel,momid);
		}
		MPI_Barrier(share_comm);
	}

	void ComputeNewCellDistancetoDomian(const int & ilevel);

	void SwitchNormalGhost(const int & i)
	{
		//MPI_Win_fence(0, m_level[i].m_box.arraywin());
		m_level[i].m_box.SwitchTwoData();
		//MPI_Win_fence(0, m_level[i].m_box.arraywin());
		MPI_Barrier(share_comm);
		SwitchDataLocation(i, m_mesh.m_level[i].m_geom);
		if (flowdataexist)
		{
			SwitchDataLocation(i, m_mesh.m_level[i].m_data);
		}
		MPI_Barrier(share_comm);
#ifdef DEBUG
		for (int gi = 0; gi < m_level[i].m_box.ghost_index.size(); ++gi)
		{
			int i0 = m_level[i].m_box.ghost_index[gi];
			if (!m_level[i].m_box.isghost(i0))
			{
				PRINTFinLEVEL("The ghost_index %d is Box %d (%d,%d,%d) but it is not a ghost!!!", 
					i, gi, i0, m_level[i].m_box[i0].ix(), m_level[i].m_box[i0].iz(), m_level[i].m_box[i0].iz());
				MPI_Abort(MPI_COMM_WORLD, 246);
			}
			if (m_level[i].m_box[i0].neib[1][1][1] != i0)
			{
				PRINTFinLEVEL("Ghost box %d (%d,%d,%d) index is %d\n",
					i,i0,m_level[i].m_box[i0].ix(), m_level[i].m_box[i0].iz(), m_level[i].m_box[i0].iz(),m_level[i].m_box[i0].neib[1][1][1]);
				MPI_Abort(MPI_COMM_WORLD, 252);
			}
		}
#endif				
	}

	template <class T>
	void SwitchDataLocation(const int & level, DataArray<T> & m_array)
	{
		// if (newarrsize > oldarrsize)
		// {
		// 	m_array.setnum(newarrsize);
		// }
		//MPI_Win_fence(0, m_array.arraywin());
		for (int j = 0; j < m_level[level].m_box.newloctag.size(); ++j)
		{
			m_array.switcholdnew(m_level[level].m_box.newloctag[j].oldloc, 
				m_level[level].m_box.newloctag[j].newloc);
		}
		//PI_Win_fence(0, m_array.arraywin());
	}
	/*After the holes in the box array have been filled, the data also needs to be switched!!!*/
	/*At this moment, the new data has been computed, but has not been put into the array!!!*/
	/*The normal and ghost cells are mixed!!!*/
	/*The data switch is achieved through the following three function*/

	/*This template function switch the flow data and other data only in normal cell*/
	template <class T>
	void InitNormalDataSwitch(const int & ilevel, DataArray<T> & m_array)
	{
		int olddsize = m_array.realsize();
		int n0 = m_level[ilevel].m_box.newloctag.size();
		for (int i = 0; i < n0; ++i)
		{
			Assert(m_level[ilevel].m_box.newloctag[i].newloc < m_level[ilevel].m_box.newloctag[i].oldloc, 
				"Initial location switch error!!!", 378);
			int oldloc = m_level[ilevel].m_box.newloctag[i].oldloc;
			if (oldloc < olddsize)
			{
				int newloc = m_level[ilevel].m_box.newloctag[i].newloc;
				m_array[newloc] = m_array[oldloc];
			}			 
		}
	}
	/*This function switch the data in both the normal and ghost cells*/
	template <class T>
	void InitGeomDataSwitch(const int & ilevel, DataArray<T> & m_array)
	{
		for (int i = 0; i < m_level[ilevel].m_box.newloctag.size(); ++i)
		{
			if (m_level[ilevel].m_box.newloctag[i].oldloc < m_array.realsize())
			{
				Assert(m_level[ilevel].m_box.newloctag[i].newloc < m_level[ilevel].m_box.newloctag[i].oldloc, 
					"Initial location switch error!!!", 378);
				m_array[m_level[ilevel].m_box.newloctag[i].newloc] = 
					m_array[m_level[ilevel].m_box.newloctag[i].oldloc];
			}
		}
	}
	/*Assembly function*/
	void AllDataInitSwitch(const int & ilevel)
	{
		InitNormalDataSwitch(ilevel, m_mesh.m_level[ilevel].m_geom);
		if (flowdataexist)
		{
			InitNormalDataSwitch(ilevel, m_mesh.m_level[ilevel].m_data);
		}
		MPI_Barrier(share_comm);
		for (int i = init_adjust_level; i < m_mesh.cur_level_num; ++i)
		{
			int newarrsize = m_level[i].m_box.realsize();
			int gn0 = m_level[i].m_box.arrayghostnum();
			if (flowdataexist)
			{
				m_mesh.m_level[i].m_data.setnum_copy(newarrsize, gn0);
			}			
			m_mesh.m_level[i].m_geom.setnum_copy(newarrsize, gn0);
		}
	}
	/*---------------------------------------*/
	template <class T>	
	void IncludeNewFineData(const int & level, DataArray<T> & f_array, 
		vector<Boxson<T> > & newfinedata)
	{			
		if (newfinedata.size() > 0)
		{
			if (!level_twod_flag[level+1])
			{
				for (int j = 0; j < m_level[level].rfnum_procs[srank]; ++j)
				{
					int i0 = m_level[level].rfbox[j];
					int tag0 = m_level[level].m_tag[i0].tag;
					for (Point_iterator p(0,2); p.end(); ++p)
					{
						int sonid = m_level[level].allson[tag0].son[p.i][p.j][p.k];
						Assert(sonid < m_level[level+1].m_box.realsize(), 
							"I can not give a right location for the data of refine!!!", 181);
						f_array[sonid] = newfinedata[j].son[p.i][p.j][p.k];
					}
				}
			}
			else
			{
				for (int j = 0; j < m_level[level].rfnum_procs[srank]; ++j)
				{
					int i0 = m_level[level].rfbox[j];
					int tag0 = m_level[level].m_tag[i0].tag;
					for (Point_iterator_2d p(0,2); p.end(); ++p)
					{
						int sonid = m_level[level].allson[tag0].son[p.i][p.j][p.k];
						Assert(sonid < m_level[level+1].m_box.realsize(), 
							"I can not give a right location for the data of refine!!!", 181);
						f_array[sonid] = newfinedata[j].son[p.i][p.j][p.k];
					}
				}
			}
		}
	}

	template <class T>
	void IncludeNewCoarseData(const int & level, DataArray<T> & c_array, vector<T> & newcoarsedata)
	{	
		if (newcoarsedata.size() > 0)
		{
			for (int j = 0; j < m_level[level].dfnum_procs[srank]; ++j)
			{
				int dfid = m_level[level].derfbox[j];
				int momid = m_level[level].m_tag[dfid].detag;
				Assert(momid < m_level[level-1].m_box.realsize(), 
					"I can not give a right location for the data of refine!!!", 181);
				c_array[momid] = newcoarsedata[j];
			}
		}
	}
/*-------------------------------------------------------------------------*/
	void Flowexist()
	{
		flowdataexist = true;
	}

	void MakeNewPoints(vector<int> & newbox, 
										 const int & aimlevel, 
										 vector<newpointloc> & newpts);

	void ManageLevelPoints(const int & ilevel);

	void CheckBoxNeib(const int & level_n)
	{
		// int bs, be;
		// ArrayOrder_s(0, m_level[level_n].m_box.realsize(), bs, be, sprocs, srank);
		for (int i = m_level[level_n].m_box.ps(); i < m_level[level_n].m_box.pe(); ++i)
		{
			for (Point_iterator p(0,3); p.end(); ++p)
			{
				int myneib = m_level[level_n].m_box[i].neib[p.i][p.j][p.k];
				if (m_level[level_n].m_box.outrange(myneib))
				{
					PRINTFinLEVEL("Box %d (%d,%d,%d) neib (%d, %d, %d) is %d index not in range!!!\n",
						level_n, i, 
						m_level[level_n].m_box[i].ix(), 
						m_level[level_n].m_box[i].iy(),
						m_level[level_n].m_box[i].iz(),
						p.i, p.j, p.k, myneib);
					MPI_Abort(MPI_COMM_WORLD, 110);
				}
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
		if (nrank == 0)
		{
			printf("$$$All ranks has passed level %d Neib check!!!\n", level_n);
		}
	}

	void CheckNeib27(const int & level_n)
	{
		for (int i = m_level[level_n].m_box.ps(); i < m_level[level_n].m_box.pe(); ++i)
		{
			for (Point_iterator p(0,3); p.end(); ++p)
			{
				int myneib = m_level[level_n].m_box[i].neib[p.i][p.j][p.k];
				if (m_level[level_n].m_box.outrange(myneib))
				{
					printf("[%d]Level %d box %d (%d,%d,%d) is %d OUT OF RANGE!!!\n",
					nrank, level_n, i, p.i, p.j, p.k, myneib);
					MPI_Abort(MPI_COMM_WORLD, 897);
				}
			}
		}
		//printf("***level %d rank %d has %d fully filled box!!!\n", level_n, nrank, fullfilled);
		//printf("***level %d rank %d has %d boundary box!!!\n", level_n, nrank, bdboxnum);
	}

	void CheckPointIndex()
	{
		for (int i = 0; i < m_mesh.cur_level_num; ++i)
		{
			for (int ptn = m_mesh.m_level[i].m_point.ps(); ptn < m_mesh.m_level[i].m_point.pe(); ++ptn)
			{
				if (ptn != m_mesh.m_level[i].m_point[ptn].index)
				{
					printf("$$$Rank %d Level %d Point %d index is %d\n", nrank, i, ptn, 
						m_mesh.m_level[i].m_point[ptn].index);
					MPI_Abort(MPI_COMM_WORLD, 784);
				}
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
		if (nrank == 0)
		{
			printf("XiXi~~All processors have passed the point index check!!!\n");
		}
	}

	void CheckGhostLocation()
	{
		for (int i = 0; i < m_mesh.cur_level_num-1; ++i)
		{
			if (!level_twod_flag[i+1])
			{
				for (int g = m_level[i].f_pro_ghost.ps(); g < m_level[i].f_pro_ghost.pe(); ++g)
				{
					int ci0 = m_level[i].f_pro_ghost[g].ci;
					int cx = m_level[i].m_box[ci0].ix();
					int cy = m_level[i].m_box[ci0].iy();
					int cz = m_level[i].m_box[ci0].iz();
					if (!m_level[i].m_box.isnormal(ci0))
					{
						printf("Level %d prolongation cell %d (%d,%d,%d) is not a normal cell type is %d tag is %d detag is %d!!!\n", 
							i, ci0, cx, cy, cz, m_level[i].m_box[ci0].type, m_level[i].m_tag[ci0].tag, m_level[i].m_tag[ci0].detag);
						MPI_Abort(MPI_COMM_WORLD, 1027);
					}
					for (Point_iterator p(0,2); p.end(); ++p)
					{
						int fi0 = m_level[i].f_pro_ghost[g].fi.son[p.i][p.j][p.k];
						if (!m_level[i+1].m_box.isghost(fi0))
						{
							printf("Level %d prolongation cell %d its son (%d,%d,%d) box %d is not a ghost of level %d\n", 
								i,ci0,p.i,p.j,p.k,fi0,i+1);
							MPI_Abort(MPI_COMM_WORLD,1036);
						}
						int fix = m_level[i+1].m_box[fi0].ix();
						int fiy = m_level[i+1].m_box[fi0].iy();
						int fiz = m_level[i+1].m_box[fi0].iz();
						if (fix-2*cx != p.i || fiy-2*cy != p.j || fiz-2*cz != p.k)
						{
							printf("Level %d prolongation cell %d is at (%d,%d,%d) but its son (%d,%d,%d) is at (%d,%d,%d)\n",
								i,ci0,cx,cy,cz,p.i,p.j,p.k,fix,fiy,fiz);
							MPI_Abort(MPI_COMM_WORLD,1048);
						}
					}
				}
				for (int ig = 0; ig < ighost; ++ig)
				{
					for (int g = m_level[i].f_res_ghost[ig].ps(); g < m_level[i].f_res_ghost[ig].pe(); ++g)
					{
						int ci0 = m_level[i].f_res_ghost[ig][g].ci;
						int cx = m_level[i].m_box[ci0].ix();
						int cy = m_level[i].m_box[ci0].iy();
						int cz = m_level[i].m_box[ci0].iz();
						if (!m_level[i].m_box.isghost(ci0))
						{
							printf("Level %d restriction cell %d is not a ghost cell!!!\n", i, ci0);
							MPI_Abort(MPI_COMM_WORLD, 1049);
						}
						for (Point_iterator p(0,2); p.end(); ++p)
						{
							int fi0 = m_level[i].f_res_ghost[ig][g].fi.son[p.i][p.j][p.k];
							if (!m_level[i+1].m_box.isnormal(fi0))
							{
								printf("Level %d restriction cell %d its son (%d,%d,%d) box %d is not a normal cell of level %d\n",
									i,ci0,p.i,p.j,p.k,fi0,i+1);
							}
							int fix = m_level[i+1].m_box[fi0].ix();
							int fiy = m_level[i+1].m_box[fi0].iy();
							int fiz = m_level[i+1].m_box[fi0].iz();
							if (fix-2*cx != p.i || fiy-2*cy != p.j || fiz-2*cz != p.k)
							{
								printf("Level %d prolongation cell %d is at (%d,%d,%d) but its son (%d,%d,%d) is at (%d,%d,%d)\n",
									i,ci0,cx,cy,cz,p.i,p.j,p.k,fix,fiy,fiz);
								MPI_Abort(MPI_COMM_WORLD,1081);
							}
						}
					}
				}
			}
			else
			{
				for (int g = m_level[i].f_pro_ghost.ps(); g < m_level[i].f_pro_ghost.pe(); ++g)
				{
					int ci0 = m_level[i].f_pro_ghost[g].ci;
					int cx = m_level[i].m_box[ci0].ix();
					int cy = m_level[i].m_box[ci0].iy();
					int cz = m_level[i].m_box[ci0].iz();
					if (!m_level[i].m_box.isnormal(ci0))
					{
						printf("Level %d prolongation cell %d (%d,%d,%d) is not a normal cell type is %d tag is %d detag is %d!!!\n", 
							i, ci0, cx, cy, cz, m_level[i].m_box[ci0].type, m_level[i].m_tag[ci0].tag, m_level[i].m_tag[ci0].detag);
						MPI_Abort(MPI_COMM_WORLD, 1027);
					}
					for (Point_iterator_2d p(0,2); p.end(); ++p)
					{
						int fi0 = m_level[i].f_pro_ghost[g].fi.son[p.i][p.j][p.k];
						if (!m_level[i+1].m_box.isghost(fi0))
						{
							printf("Level %d prolongation cell %d its son (%d,%d,%d) box %d is not a ghost of level %d\n", 
								i,ci0,p.i,p.j,p.k,fi0,i+1);
							MPI_Abort(MPI_COMM_WORLD,1036);
						}
						int fix = m_level[i+1].m_box[fi0].ix();
						int fiy = m_level[i+1].m_box[fi0].iy();
						int fiz = m_level[i+1].m_box[fi0].iz();
						if (fix-2*cx != p.i || fiy-2*cy != p.j || fiz-cz != 0)
						{
							printf("Level %d prolongation cell %d is at (%d,%d,%d) but its son (%d,%d,%d) is at (%d,%d,%d)\n",
								i,ci0,cx,cy,cz,p.i,p.j,p.k,fix,fiy,fiz);
							MPI_Abort(MPI_COMM_WORLD,1048);
						}
					}
				}
				GiveLevelFlag("Finsh check pro ghost pair!!!",i,5);
				for (int ig = 0; ig < ighost; ++ig)
				{
					for (int g = m_level[i].f_res_ghost[ig].ps(); g < m_level[i].f_res_ghost[ig].pe(); ++g)
					{
						int ci0 = m_level[i].f_res_ghost[ig][g].ci;
						int cx = m_level[i].m_box[ci0].ix();
						int cy = m_level[i].m_box[ci0].iy();
						int cz = m_level[i].m_box[ci0].iz();
						if (!m_level[i].m_box.isghost(ci0))
						{
							printf("Level %d restriction cell %d is not a ghost cell!!!\n", i, ci0);
							MPI_Abort(MPI_COMM_WORLD, 1049);
						}
						for (Point_iterator_2d p(0,2); p.end(); ++p)
						{
							int fi0 = m_level[i].f_res_ghost[ig][g].fi.son[p.i][p.j][p.k];
							if (!m_level[i+1].m_box.isnormal(fi0))
							{
								printf("Level %d restriction cell %d its son (%d,%d,%d) box %d is not a normal cell of level %d\n",
									i,ci0,p.i,p.j,p.k,fi0,i+1);
							}
							int fix = m_level[i+1].m_box[fi0].ix();
							int fiy = m_level[i+1].m_box[fi0].iy();
							int fiz = m_level[i+1].m_box[fi0].iz();
							if (fix-2*cx != p.i || fiy-2*cy != p.j || fiz-cz != 0)
							{
								printf("Level %d prolongation cell %d is at (%d,%d,%d) but its son (%d,%d,%d) is at (%d,%d,%d)\n",
									i,ci0,cx,cy,cz,p.i,p.j,p.k,fix,fiy,fiz);
								MPI_Abort(MPI_COMM_WORLD,1081);
							}
						}
					}
				}
			}
		}
		for (int i = 0; i < m_mesh.cur_level_num; ++i)
		{
			int bkps = m_level[i].blockpair.ps();
			int bkpe = m_level[i].blockpair.pe();
			for (int i0 = bkps; i0 < bkpe; ++i0)
			{
				int in0 = m_level[i].blockpair[i0].innode;
				if (m_level[i].m_box[in0].type != Normalcell)
				{
					printf("Level %d Block pair %d (%d,%d,%d) type is %d\n", i, in0, 
						m_level[i].m_box[in0].ix(),m_level[i].m_box[in0].iy(),m_level[i].m_box[in0].iz(),
						m_level[i].m_box[in0].type);
					MPI_Abort(MPI_COMM_WORLD, 714);
				}
			}
		}
	}

	/*-----------------------------------------*/
	//这个函数的作用是为了清除res pair插值范围内的detag
	//从而能够保证
	void CheckDetagLocation(const int & ilevel)
	{
		m_level[ilevel].SynGhostTag();
		//MPI_Win_fence(0, m_level[ilevel].m_tag.arraywin());
		Point dp[6] = {Point(0,1,1),Point(2,1,1),Point(1,0,1),Point(1,2,1),Point(1,1,0),Point(1,1,2)};
		int bs = m_level[ilevel].f_res_ghost[0].ps();
		int be = m_level[ilevel].f_res_ghost[0].pe();
		// for (int i = bs; i < be; ++i)
		// {
		// 		int ci0 = m_level[ilevel].f_pro_ghost[i].ci;
		// 		for (Point_iterator pn(0,3); pn.end(); ++pn)
		// 		{
		// 			int an0 = m_level[ilevel].m_box[ci0].neib[pn.i][pn.j][pn.k];
		// 			if (an0 > -1)
		// 			{
		// 				m_level[ilevel].m_tag[an0].detag = -1;
		// 				//printf("Recover the tag at L%dB%d\n", ilevel, an0);
		// 			}
		// 		}
		// }
		if (!m_level[ilevel].twodflag)
		{
			for (int i = bs; i < be; ++i)
			{
				int ci0 = m_level[ilevel].f_res_ghost[0][i].ci;
				int ix0  = m_level[ilevel].m_box[ci0].ix()&1;
				int iy0  = m_level[ilevel].m_box[ci0].iy()&1;
				int iz0  = m_level[ilevel].m_box[ci0].iz()&1;
				int ss[3] = {1-ix0, 1-iy0, 1-iz0};
				int se[3] = {3-ix0, 3-iy0, 3-iz0};
				for (int nx = ss[0]; nx < se[0]; ++nx)
				{
					for (int ny = ss[1]; ny < se[1]; ++ny)
					{
						for (int nz = ss[2]; nz < se[2]; ++nz)
						{
							int anb = m_level[ilevel].m_box[ci0].neib[nx][ny][nz];
							if (anb > -1)
							{
								int anbx = m_level[ilevel].m_box[anb].ix()&1;
								int anby = m_level[ilevel].m_box[anb].iy()&1;
								int anbz = m_level[ilevel].m_box[anb].iz()&1;
								for (Point_iterator p(0,2); p.end(); ++p)
								{
									int anb2 = m_level[ilevel].m_box[anb].neib[anbx+p.i][anby+p.j][anbz+p.k];
									if (anb2 > -1)
									{
										if (m_level[ilevel].m_box[anb2].type == Blockghost)
										{
											m_level[ilevel].m_tag[anb2].detag = -5;
										}
										else
										{
											if (m_level[ilevel].m_tag[anb2].detag == 0)
											{
												printf("L%d res pair neib will be coarsened!!! Please increase the gap between level %d and its fine level!!!\n",
													ilevel, ilevel);
												MPI_Abort(MPI_COMM_WORLD, 772);
											}
											m_level[ilevel].m_tag[anb2].detag = -1;
										}
									}
								}
							}
						}
					}
				}
				// m_level[ilevel].m_tag[ci0].detag = -1;
				for (int di = 0; di < 6; ++di)
				{
					int an1 = ci0;
					for (int i00 = 0; i00 < 4; ++i00)
					{
						if (an1 > -1)
						{
							an1 = m_level[ilevel].m_box[an1].neib[dp[di][0]][dp[di][1]][dp[di][2]];
						}
						else break;
					}
					if (an1 > -1)
					{
						if (m_level[ilevel].m_box[an1].type == Blockghost)
						{
							m_level[ilevel].m_tag[an1].detag = -5;
						}
						else
						{

							if (m_level[ilevel].m_tag[an1].detag == 0)
							{
								printf("L%d res pair neib will be coarsened!!! Please increase the gap between level %d and its fine level!!!\n",
										ilevel, ilevel);
								MPI_Abort(MPI_COMM_WORLD, 772);
							}		
							m_level[ilevel].m_tag[an1].detag = -1;
						}
					}
				}
			}
		}
		else
		{
			int dnz_num[6] = {4,4,4,4,2,2};
			for (int i = bs; i < be; ++i)
			{
				int ci0 = m_level[ilevel].f_res_ghost[0][i].ci;
				int ix0  = m_level[ilevel].m_box[ci0].ix()&1;
				int iy0  = m_level[ilevel].m_box[ci0].iy()&1;
				int ss[2] = {1-ix0, 1-iy0};
				int se[2] = {3-ix0, 3-iy0};
				for (int nx = ss[0]; nx < se[0]; ++nx)
				{
					for (int ny = ss[1]; ny < se[1]; ++ny)
					{
						int anb = m_level[ilevel].m_box[ci0].neib[nx][ny][1];
						if (anb > -1)
						{
							int anbx = m_level[ilevel].m_box[anb].ix()&1;
							int anby = m_level[ilevel].m_box[anb].iy()&1;
							for (Point_iterator_2d p(0,2); p.end(); ++p)
							{
								for (int nz0 = 0; nz0 < 3; ++nz0)
								{
									int anb2 = m_level[ilevel].m_box[anb].neib[anbx+p.i][anby+p.j][nz0];
									if (anb2 > -1)
									{
										if (m_level[ilevel].m_box[anb2].type == Blockghost)
										{
											m_level[ilevel].m_tag[anb2].detag = -5;
										}
										else
										{
											if (m_level[ilevel].m_tag[anb2].detag == 0)
											{
												printf("L%d res pair neib will be coarsened!!! Please increase the gap between level %d and its fine level!!!\n",
												ilevel, ilevel);
												MPI_Abort(MPI_COMM_WORLD, 772);
											}
											m_level[ilevel].m_tag[anb2].detag = -1;
										}
									}
								}
							}
						}
					}
				}
				// m_level[ilevel].m_tag[ci0].detag = -1;
				for (int di = 0; di < 6; ++di)
				{
					int an1 = ci0;
					for (int i00 = 0; i00 < dnz_num[di]; ++i00)
					{
						if (an1 > -1)
						{
							an1 = m_level[ilevel].m_box[an1].neib[dp[di][0]][dp[di][1]][dp[di][2]];
						}
						else break;
					}
					if (an1 > -1)
					{
						if (m_level[ilevel].m_box[an1].type == Blockghost)
						{
							m_level[ilevel].m_tag[an1].detag = -5;
						}
						else
						{
							if (m_level[ilevel].m_tag[an1].detag == 0)
							{
								printf("L%d res pair neib will be coarsened!!! Please increase the gap between level %d and its fine level!!!\n",
										ilevel, ilevel);
								MPI_Abort(MPI_COMM_WORLD, 772);
							}
							m_level[ilevel].m_tag[an1].detag = -1;
						}
					}
				}
			}
		}
		MPI_Barrier(share_comm);
		be = m_level[ilevel].ghost_target.size();
		vector<vector<int> > informbkp(nodenum);
		vector<int> recvbkp;
		for (int i = 0; i < be; ++i)
		{
			int i0 = m_level[ilevel].ghost_target[i];
			Assert(m_level[ilevel].m_box[i0].type == Blockghost, "The ghost target is not a block ghost!!!", 793);
			if (m_level[ilevel].m_tag[i0].detag == -5)
			{
				informbkp[m_level[ilevel].originbp[i].node].push_back(m_level[ilevel].originbp[i].index);
				m_level[ilevel].m_tag[i0].detag = -1;
			}
		}
		MPI_Barrier(share_comm);
		BcastNewInfo_Node(informbkp, recvbkp, 0);
		be = recvbkp.size();
		for (int i = 0; i < be; ++i)
		{
			int i0 = recvbkp[i];
			int g0 = m_level[ilevel].blockpair[i0].innode;
			if (m_level[ilevel].m_tag[g0].detag == 0)
			{
				printf("L%d res pair neib will be coarsened!!! Please increase the gap between level %d and its fine level!!!\n",
										ilevel, ilevel);
				MPI_Abort(MPI_COMM_WORLD, 772);
			}
			m_level[ilevel].m_tag[g0].detag = -1;
		}
		MPI_Barrier(MPI_COMM_WORLD);
		m_level[ilevel].TrimDetag();
#ifdef DEBUG
	if (!m_level[ilevel].twodflag)
	{
		for (int i0 = m_level[ilevel].m_tag.ps(); i0 < m_level[ilevel].m_tag.pe(); ++i0)
		{
			if (m_level[ilevel].m_tag[i0].detag != -1)
			{
				int sx00 = m_level[ilevel].m_box[i0].ix() & 1;
				int sy00 = m_level[ilevel].m_box[i0].iy() & 1;
				int sz00 = m_level[ilevel].m_box[i0].iz() & 1;
				for (int nx = 1-sx00; nx < 3-sx00; ++nx)
				{
					for (int ny = 1-sy00; ny < 3-sy00; ++ny)
					{
						for (int nz = 1-sz00; nz < 3-sz00; ++nz)
						{
							int sneib = m_level[ilevel].m_box[i0].neib[nx][ny][nz];
							if (sneib < 0)
							{
								PRINTFinLEVEL("Box odd %d %d %d box is %d neib %d %d %d is %d\n", ilevel, sx00, sy00, sz00, i0, nx, ny, nz, sneib);
							}
							Assert(sneib > -1, "The four box must be non-negative 79!!!", 79);
							if (m_level[ilevel].m_tag[sneib].detag != m_level[ilevel].m_tag[i0].detag)
							{
								printf("Box (%d,%d,%d) detag is 0 but its neib (%d,%d,%d) detag is %d\n",
									m_level[ilevel].m_box[i0].ix(), m_level[ilevel].m_box[i0].iy(), m_level[ilevel].m_box[i0].iz(), 
									nx, ny, nz, m_level[ilevel].m_tag[sneib].detag);
								MPI_Abort(MPI_COMM_WORLD, 84);
							}
						}
					}
				}	
			}
		}
	}
	else
	{
		for (int i0 = m_level[ilevel].m_tag.ps(); i0 < m_level[ilevel].m_tag.pe(); ++i0)
		{
			if (m_level[ilevel].m_tag[i0].detag != -1)
			{
				int sx00 = m_level[ilevel].m_box[i0].ix() & 1;
				int sy00 = m_level[ilevel].m_box[i0].iy() & 1;
				for (int nx = 1-sx00; nx < 3-sx00; ++nx)
				{
					for (int ny = 1-sy00; ny < 3-sy00; ++ny)
					{
						int sneib = m_level[ilevel].m_box[i0].neib[nx][ny][1];
						if (sneib < 0)
						{
							PRINTFinLEVEL("2d level Box odd %d %d box is %d neib %d %d %d is %d\n", 
								ilevel, sx00, sy00, i0, nx, ny, 1, sneib);
						}
						Assert(sneib > -1, "The four box must be non-negative 79!!!", 79);
						if (m_level[ilevel].m_tag[sneib].detag != m_level[ilevel].m_tag[i0].detag)
						{
							printf("Box (%d,%d,%d) detag is 0 but its neib (%d,%d,%d) detag is %d\n",
								m_level[ilevel].m_box[i0].ix(), m_level[ilevel].m_box[i0].iy(), m_level[ilevel].m_box[i0].iz(), 
								nx, ny, 1, m_level[ilevel].m_tag[sneib].detag);
							MPI_Abort(MPI_COMM_WORLD, 84);
						}
					}
				}	
			}
		}
	}
#endif
}

	//void FaceAdjust();

	// void RenewFaceIndex(const int & ilevel);

	//void MakeNewFace(const int & i, vector<Face> & newfineface, vector<Face> & newcface);

	//void FaceArrayHole(const int & ilevel);

	//void GetRemoveFaceInfoRefine();

	//void GetRemoveFaceInfoDerefine();

	// void PrintAllOnesideFace();

	void FaceSideBox(const int & ilevel, vector<Face> & newfineface, 
		const int & fn, const int & bn);

	void TagGhostCell();

	void InterfaceExchange();

	void InterfaceProlongation(const int & nlevel);

	void InterfaceRestriction(const int & nlevel);

	void FindModifyFace(vector<DataArray<Facepair> > & mdyface);

	void ConstructGhostFace();

	

	// void RecoverSlidingPartPoint(const int & nlevel)
	// {
	// 	for (int i = 0; i < m_level[nlevel].drf_sliding_res_to_pro.size(); ++i)
	// 	{
	// 		int ci0 = m_level[nlevel].drf_sliding_res_to_pro[i].ci;
	// 		for (Point_iterator p(0,2); p.end(); ++p)
	// 		{
	// 			int ptn = m_level[nlevel].m_box[ci0].pts[p.i][p.j][p.k];
	// 			if (m_mesh.m_level[nlevel].m_point[ptn].index == -1)
	// 			{
	// 				printf("Recover box (%d,%d,%d) point (%d,%d,%d) index is %d\n", 
	// 					m_level[nlevel].m_box[ci0].ix(),
	// 					m_level[nlevel].m_box[ci0].iy(),
	// 					m_level[nlevel].m_box[ci0].iz(),
	// 					p.i,p.j,p.k,ptn);
	// 				m_mesh.m_level[nlevel].m_point[ptn].index = ptn;
	// 			}
	// 		}
	// 	}
	// }

// 	void BalanceGhostPair()
// 	{
// 		for (int i = init_adjust_level; i < m_mesh.cur_level_num; ++i)
// 		{
// 			taskbalance.LoadBalance(m_level[i].f_pro_ghost);
// 			taskbalance.LoadBalance(m_level[i].f_res_ghost[0]);
// 			taskbalance.LoadBalance(m_level[i].f_res_ghost[1]);
// #ifdef DEBUG
// 			printf("[%d]Level %d Prolongation pair %d\n", nrank, i, (int)m_level[i].f_pro_ghost.size());
// 			printf("[%d]Level %d Restriction 0 pair %d\n", nrank, i, (int)m_level[i].f_res_ghost[0].size());
// 			printf("[%d]Level %d Restriction 0 pair %d\n", nrank, i, (int)m_level[i].f_res_ghost[1].size());
// #endif			
// 		}	
// 	}

private:

	Mesh &              m_mesh;
	vector<AMRLevel>    m_level;
	bool                flowdataexist;
	int init_adjust_level;
	TwoLevelOp tlop;
	bool morelevel;
	// Taskdist 	taskbalance;
};


#endif
