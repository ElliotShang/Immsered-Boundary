#ifndef _Mesh_H_
#define _Mesh_H_

#include <iostream>
#include <vector>
//#include <cmath>
#include <unistd.h>

#include "Box.H"
#include "Domain.H"
#include "Check.H"
#include "AMRmpi.H"
#include "AMRDatatype.H"
#include "Meshlevel.H"
#include "AMRSpaceTime.H"
#include "MeshOp.H"
#include "Face.H"
#include "GhostPair.H"
#include "Wallface.H"

using namespace std;

extern string meshfile;
extern Pointxyz meshscale;

extern vector<bool> level_twod_flag;
extern vector<int> level_power_ratio;

extern bool gridisnew;

extern Pointxyz frame_rotate_vel;
extern Pointxyz frame_rotate_origin;
extern Pointxyz frame_trans_vel;

void SetTwodFlag();

void SetPower_Ratio();

struct newfaceinfo
{
	int box;
	int facedir;
	int face0;
	newfaceinfo()
	{
		box = -1;
		facedir = -1;
		face0 = -1;
	}
};

class Mesh
{

protected:
	int                  max_level_num;
	int                  cur_level_num;
	vector<Meshlevel>    m_level;
	Domain &             m_dm;
	DataArray<IBCell> 				 m_dis;
	DataArray<wallface> m_wallface;
	DataArray<int>             infectbox;
	vector<DataArray<Facepair> > mdyface;
	Nodedata nd_mdis;
	Nodedata rev_mdis;
	Nodedata meshnd;
	vector<vector<Dmcellface> > facecenter[3][2];

public:

	friend class Pltoutput;
	friend class AMR;
	friend class NS_Solver;
	friend class BCValues;
	friend class TwoLevelOp;
	friend class Body;

	Mesh(Domain & a_dm, int n):
	m_dm(a_dm)
	{
		Mesh_Init(n);
		m_level.reserve(n);
		GiveAFlag("add a Meshlevel...", 5);
		m_level.push_back(Meshlevel());
		GiveAFlag("Finish add a Meshlevel!!!", 5);
		BuiltBackMesh(a_dm);
		GiveAFlag("Finish build back mesh!!!", 5);
	}

	void Mesh_Init(int n)
	{
		max_level_num = n;
		cur_level_num = 1;
	}

	void AddMeshlevel()
	{
		m_level.push_back(Meshlevel());
		cur_level_num += 1;			
	}

	inline bool Canberefined()
	{
		if (cur_level_num < max_level_num)
		{
			return true;
		}
		else
		{
			return false;
		}
	}

	inline bool Canbederefined()
	{
		if (cur_level_num > 1)
		{
			return true;
		}
		else
		{
			return false;
		}
	}

	void RemoveLastMeshlevel()
	{
		m_level.back().clearlevel();
		m_level.pop_back();
		cur_level_num -= 1;		
	}

	inline DataArray<Box> & ShowLevelBox(const int & level_n)
	{
		return m_level[level_n].m_box;
	}

	inline DataArray<Box> & LastLayerBox()
	{
		return m_level[cur_level_num-1].m_box;
	}

	inline DataArray<Box> * ShowLevelBoxPtr(const int & level_n)
	{
		return &m_level[level_n].m_box;
	}

	inline Meshlevel * ShowLevel(const int & level_n)
	{
		return &m_level[level_n];
	}

	inline const int & BoxNeib(const int & aboxlevel, const int & aboxindex, 
		const int & ni, const int & nj, const int & nk)
	{
		return m_level[aboxlevel].m_box[aboxindex].neib[ni][nj][nk];
	}

	inline DataArray<FlowVariables> * ShowLevelDataPtr(const int & level_n)
	{
		return &m_level[level_n].m_data;
	}

	inline const int & MyCurNum()
	{
		return cur_level_num;
	}

	inline const int & MyMaxNum()
	{
		return max_level_num;
	}

	inline vector<Point> & FaceNormDir()
	{
		return m_dm.facedir;
	}

	void BuiltBackMesh(Domain & a_dm);

	inline const int & LevelBoxStart(const int & i)
	{
		return m_level[i].m_box.ps();
	}

	inline const int & LevelBoxEnd(const int & i)
	{
		return m_level[i].m_box.pe();
	}

	inline Box & LevelBox(const int & i, const int & bj)
	{
		return m_level[i].m_box[bj];
	}

	inline const Box * LevelBoxPtr(const int & i, const int & bj)
	{
		return &m_level[i].m_box[bj];
	}

	inline FlowVariables & LevelBoxData(const int & i, const int & bj)
	{
		return m_level[i].m_data[bj];
	}

	inline DataArray<FlowVariables> & LevelDataArray(const int & i)
	{
		return m_level[i].m_data;
	}

	inline FlowVariables * LevelBoxDataPtr(const int & i, const int & bj)
	{
		return &m_level[i].m_data[bj];
	}

	inline const bool & PeriodicX()
	{
		return m_dm.PeriodicX();
	}

	inline const bool & PeriodicY()
	{
		return m_dm.PeriodicY();
	}

	inline const bool & PeriodicZ()
	{
		return m_dm.PeriodicZ();
	}

	inline Point & Boxnum()
	{
		return m_dm.DmBoxNum();
	}
	inline Domain & Meshdomain()
	{
		return m_dm;
	}
	inline DataArray<Domainghost> & MeshGhostArray(const int & i0)
	{
		return m_level[i0].dghost;
	}
	/*------------------------------------------------------*/
	inline const int & Facestart(const int & i)
	{
		return m_level[i].m_face.ps();
	}

	inline const int & Faceend(const int & i)
	{
		return m_level[i].m_face.pe();
	}

	inline Face & Levelface(const int & ilevel, const int & ibox)
	{
		return m_level[ilevel].m_face[ibox];
	}

	inline DataArray<Face> & Meshface(const int & ilevel)
	{
		return m_level[ilevel].m_face;
	}

	inline double & Facearea(const int & ilevel, const int & iface)
	{
		return m_level[ilevel].m_face[iface].area;
	}
	/*-----------------------------------------------------*/
	inline int & GridPoint(const int & ilevel, const int & ibox, 
		const int & pi, const int & pj, const int & pk)
	{
		return m_level[ilevel].m_box[ibox].pts[pi][pj][pk];
	}

	inline double distoboxpt(const int & alevel, const int & abox, 
		const int & pi, const int & pj, const int & pk, const Pointxyz & apt)
	{
		int pti = m_level[alevel].m_box[abox].pts[pi][pj][pk];
		Pointxyz bpt = m_level[alevel].m_point[pti].xyz - apt;
		return bpt.length();
	}

	inline Pointxyz & GridPointXYZ(const int & ilevel, const int & ibox, 
		const int & pi, const int & pj, const int & pk)
	{
		return m_level[ilevel].m_point[m_level[ilevel].m_box[ibox].pts[pi][pj][pk]].xyz;
	}

  inline Pointxyz & bc(const int & ilevel, const int & bn)
	{
		return m_level[ilevel].m_geom[bn].boxcenter;
	}

  inline CellMCoef & boxmc(const int & ilevel, const int & ibox)
	{
		return m_level[ilevel].m_geom[ibox].keisa;
	}

	inline double & boxv(const int & ilevel, const int & ibox)
	{
		return m_level[ilevel].m_geom[ibox].v;
	}

	inline int & boxface(const int & ilevel, const int & ibox, const int & dir, const int & f0)
	{
		return m_level[ilevel].m_box[ibox].faces[dir][f0];
	}

	inline Face & Boxsideface(const int & ilevel, const int & ibox, const int & dir, const int & f0)
	{
		int fi = boxface(ilevel,ibox,dir,f0);
		return m_level[ilevel].m_face[fi];
	}

	inline double distobc(const int & ilevel, const int & bn, 
		const double & x0, const double & y0, const double & z0)
  {
  	Pointxyz inputxyz(x0, y0, z0);
#ifndef PASSAGE_ANGLE  	
  	Pointxyz dxyz = m_level[ilevel].m_geom[bn].boxcenter - inputxyz;
  	PeriodicLength(dxyz);
#else
  	Pointxyz newbcxyz;
  	PeriodicAnnulaLength(m_level[ilevel].m_geom[bn].boxcenter, inputxyz, newbcxyz);
  	Pointxyz dxyz = newbcxyz - inputxyz;
#endif  	
    return dxyz.length();
  }

  inline double distobc(const int & ilevel, const int & bn, 
		Pointxyz & ap)
  {
#ifndef PASSAGE_ANGLE  	
  	Pointxyz dxyz = m_level[ilevel].m_geom[bn].boxcenter - ap;
  	PeriodicLength(dxyz);
#else
  	Pointxyz newbcxyz;
  	PeriodicAnnulaLength(m_level[ilevel].m_geom[bn].boxcenter, ap, newbcxyz);
  	Pointxyz dxyz = newbcxyz - ap;
#endif
    return dxyz.length();
  }

  void ComptMassFlow(double & massflow);

  void ComptMassFlow_Outlet(double & massflow);

  void ComptAveInletVel(double & invel);

  void ComptSectionAvePre(const int & level, double * p0, double * T0, const int & array_size);

  void ComptSectionAvePre(const int & level, double & p0, double & T0, int z0);

  void ComptSectionAvePre(const int & level, double * p0, double * T0, int * z0, int ptnum);

  void GetPointVar(const int & ilevel, Pointxyz & mypt, const int & varindex, double & avar);

  void GetPointVar(const int & ilevel, Point & mypt, const int & varindex, double & avar);

  void ComptDistancetoDomain();

	void CheckBoxCenter();

	void CheckFaceBoxIndex();

	void CheckBoxOnlyOne();

	void CheckBox(const int & level_n);

	void CheckPointxyz();

	void CheckBoxPointIndex(const int & i);

	void ReverseIndexCheck(const int & level_n);

	void ComputeLevelGeom(const int & ilevel);

	void CheckHgIntpcells(const int & ilevel, HGCell & hginfo);

	inline int Infectboxnum()
  {
  	return infectbox.size();
  }

  inline int & BoxInfectIndex(const int & i0)
  {
  	return infectbox[i0];
  }

  inline int & ibcellifc(const int & i0)
  {
  	return infectbox[m_dis[i0].ci];
  }

  inline void InfectABox(const int & i0)
  {
  	infectbox[i0] = -2;
  }

  inline double & InfectedBoxDis(const int & abox)
	{
		return m_level[cur_level_num-1].m_box[infectbox[abox]].pair.signdis;
	}

	inline DataArray<int> & IfcArray()
	{
		return infectbox;
	}

  void InitWallInfect()
  {
  	infectbox.setnum_nocopy(m_level[cur_level_num-1].m_box.realsize(), 0);
    for (int i = infectbox.ps(); i < infectbox.pe(); ++i)
    {
      infectbox[i] = -1;
    }
    //m_dis.setnum(0);
    MPI_Barrier(share_comm);
  }

  void CountInfectedBox();

  void RenewInfectIndex()
  {
  	int ifs = infectbox.ps();
  	int ife = infectbox.pe();
  	for (int i = ifs; i < ife; ++i)
  	{
  		infectbox[i] = -1;
  	}
  	MPI_Barrier(share_comm);
  	int mds = m_dis.ps();
  	int mde = m_dis.pe();
  	for (int i = mds; i < mde; ++i)
  	{
  		infectbox[m_dis[i].ci] = i;
  	}
  	MPI_Barrier(share_comm);
  }

  void AdjustInfectedArray()
  {
  	if (infectbox.size() > 0)
  	{
  		infectbox.setnum_nocopy(m_level[cur_level_num-1].m_box.realsize(), 0);
  		int ifs = infectbox.ps();
  		int ife = infectbox.pe();
  		for (int i = ifs; i < ife; ++i)
  		{
  			infectbox[i] = -1;
  		}
  		MPI_Barrier(share_comm);
  		int mds = m_dis.ps();
  		int mde = m_dis.pe();
  		for (int i = mds; i < mde; ++i)
  		{
  			int ci0 = m_dis[i].ci;
  			m_dis[i].ci = m_level[cur_level_num-1].m_box[ci0].neib[1][1][1];
  			Assert(m_dis[i].ci > -1, "The ib array main cell should be non-negative!!!", 485);
  			infectbox[m_dis[i].ci] = i;

  			ci0 = m_dis[i].hgc.closecell;
  			Assert(ci0 > -1, "The ib cell hg close cell should be non-negative!!!", 489);
  			m_dis[i].hgc.closecell =  m_level[cur_level_num-1].m_box[ci0].neib[1][1][1];
  			for (Point_iterator p(0,2); p.end(); ++p)
  			{
  				int an0 = m_dis[i].hgc.intpcell[p.i][p.j][p.k];
  				//Assert(an0 > -1, "The ib cell hg intp cell should be non-negative!!!", 494);
  				if (an0 > -1) m_dis[i].hgc.intpcell[p.i][p.j][p.k] = m_level[cur_level_num-1].m_box[an0].neib[1][1][1];
  			}
  		}
  	}
  }

  inline IBCell & ibcell(const int & i0)
  {
  	return m_dis[i0];
  }

  inline IBCell & ifccell(const int & b0)
  {
  	return m_dis[infectbox[b0]];
  }

  inline DataArray<IBCell> & ibarray()
  {
  	return m_dis;
  }

  inline bool isghostcell(const int & ilevel, const int & ibox)
  {
  	return m_level[ilevel].m_box.isghost(ibox);
  }

  inline bool isnormalcell(const int & ilevel, const int & ibox)
  {
  	return m_level[ilevel].m_box.isnormal(ibox);
  }

  void ConstructBoxFace(const int & ilevel);

  void ComputeBlockPairTheta(const int & ilevel);

  void ConstructBoxPoint(DataArray<Box> & abox, 
												 DataArray<mPoint> & apt, 
												 const int & ilevel,
												 const int & startbox,
												 const int & endbox);

  void DataExchangeCells(const int & nlevel);

  void DataExchangeCells_Wall();

  void DataExchange(Nodedata & nd, const int & ilevel);

  void DataExchange_alllevels(const int & time_step);

  void ConstructGhostPoint(vector<vector<vector<Pointxyz> > > & meshpts, Point & ptnums_g);

  bool InTransferRange(const int & level0, const int & box0)
  {
  	int r0 = pow(2,level0);
  	int r0z = level_power_ratio[level0];
  	int ix0 = m_level[level0].m_box[box0].ix();
  	int iy0 = m_level[level0].m_box[box0].iy();
  	int iz0 = m_level[level0].m_box[box0].iz();
  	if (ix0 < m_dm.dmbound[0][1]*r0+dtln && ix0 > m_dm.dmbound[0][0]*r0-1-dtln &&
  		iy0 < m_dm.dmbound[1][1]*r0+dtln && iy0 > m_dm.dmbound[1][0]*r0-1-dtln &&
  		iz0 < m_dm.dmbound[2][1]*r0z+dtln && iz0 > m_dm.dmbound[2][0]*r0z-1-dtln)
  	{
  		return true;
  	}
  	else
  	{
  		printf("N%dL%dB(%d,%d,%d) domain range:x(%d,%d) y(%d,%d) z(%d,%d)\n",
  			node, level0, ix0, iy0, iz0,
  			m_dm.dmbound[0][0],m_dm.dmbound[0][1],
  			m_dm.dmbound[1][0],m_dm.dmbound[1][1],
  			m_dm.dmbound[2][0],m_dm.dmbound[2][1]);
  		return false;
  	}
  }

  bool IsNormalcellButNotBlockpair(const int & level0, const int & box0)
  {
  	int r0 = pow(2,level0);
  	int r0z = level_power_ratio[level0];
  	int ix0 = m_level[level0].m_box[box0].ix();
  	int iy0 = m_level[level0].m_box[box0].iy();
  	int iz0 = m_level[level0].m_box[box0].iz();
  	if (ix0 < m_dm.dmbound[0][1]*r0-dtln && ix0 > m_dm.dmbound[0][0]*r0-1+dtln &&
  			iy0 < m_dm.dmbound[1][1]*r0-dtln && iy0 > m_dm.dmbound[1][0]*r0-1+dtln &&
  			iz0 < m_dm.dmbound[2][1]*r0z-dtln && iz0 > m_dm.dmbound[2][0]*r0z-1+dtln)
  	{
  		return true;
  		printf("N%dL%dB(%d,%d,%d) IsNormalcellButNotBlockpair domain range:x(%d,%d) y(%d,%d) z(%d,%d)\n",
  			node, level0, ix0, iy0, iz0,
  			m_dm.dmbound[0][0],m_dm.dmbound[0][1],
  			m_dm.dmbound[1][0],m_dm.dmbound[1][1],
  			m_dm.dmbound[2][0],m_dm.dmbound[2][1]);
  	}
  	else
  	{
  		return false;
  	}
  }

  void SynNormal_BlockDistance(const int & ilevel, vector<int> & move_from_in_to_out, vector<int> & move_from_out_to_in);

  void CounterSliceCellNumber(const int & dir, vector<int> & nodestart, vector<int> & nodeend);

  void ConstructBoxPoint_FromImportMesh();

  void ReadPlot3dFile(vector<vector<vector<Pointxyz> > > & meshpts, Point & ptnums_g, double * dh0);

  void ReadAlphaFlowGrid(vector<vector<vector<Pointxyz> > > & meshpts, Point & ptnums_g, double * dh0);

  void ReadUnstructGrid(vector<vector<vector<Pointxyz> > > & meshpts, Point & ptnums_g, double * dh0);

  void ComputeDmgFaceAngle(const int & ilevel);

  inline Nodedata & LevelDataWin(const int & ilevel)
  {
  	return m_level[ilevel].nd;
  }

  void ComptHgIntpCoef(HGCell & hginfo, Pointxyz & dkeisa, const int & ilevel);

  void LocateHGCell(const int & ilevel, HGCell & hginfo, Pointxyz & patc, int & hgnode, Pointxyz & dkeisa, int & rev_dir);

  void FindHGCellIntpCell(HGCell & hginfo, const int & ilevel, Pointxyz & dkeisa, int & hgnode);

  void LocatePointinComptDomain_AllDomain(Pointxyz & apt, const int & ilevel, int & init_box, Pointxyz & dkeisa);

  void ReversePeriodicSide(const int & ilevel, int & init_box, int & rev_dir);

  //void LocatePointinComptDomain_NormalRange(Pointxyz & apt, const int & ilevel, int & init_box, int & apt_node);

  void IntpHGgcellVar(HGCell & myhg, const int & lastlevel);

  void IntpcellNotSolid(bool & checkresult, const int & mylevel, HGCell & hgc);

  void IntpcellNotSolid_NotInfect(bool & checkresult, const int & mylevel, HGCell & hgc);

  void CopyIBCellDistance(vector<double> & celldis);

  void ComptCelldistancetoDomain(const int & ilevel, const int & ibox);

  void ComptCelldistancetoSection(const int & ilevel, const int & ibox, const int & sec_d1, const int & sec_d2);

  void ConstructDomainCellFace(vector<vector<vector<Pointxyz> > > & meshpts);

  void GetCellDomainNormalVector(const int & ilevel, const int & ibox, Pointxyz & cellnmv);

  int FindNeibBlock(Box & lastbox, const int & ilevel);

  void TagIBCellOtherNode();

  void DataExchange_Pairinfo(Nodedata & nd, const int & ilevel);

  void MarkSpecialBox();

};





#endif
