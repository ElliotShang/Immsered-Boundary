#ifndef _AMRLevel_H_
#define _AMRLevel_H_

#include <vector>
#include <iostream>
#include <mpi.h>

#include "Box.H"
#include "Check.H"
#include "Point_iterator.H"
#include "AMRmpi.H"
#include "DataArray.H"
#include "FlowVariables.H"
#include "GhostPair.H"
#include "AMRDatatype.H"


using namespace std;
struct newpointloc
{
	int boxindex;
	int vectindex;
	int boxpointx;
	int boxpointy;
	int boxpointz;
	int arrindex;

	newpointloc()
	{}

	newpointloc(const int & i1, const int & i2, const int & i3, const int & i4, const int & i5)
	{
		boxindex = i1;
		boxpointx = i2;
		boxpointy = i3;
		boxpointz = i4;
		arrindex = -1;
		vectindex = i5;
	}
};

class AMRLevel
{
	friend class AMR;
	friend class Pltoutput;

public:

	AMRLevel(DataArray<Box> & a_box, DataArray<NodePair> & andp, DataArray<Domainghost> & admg, int & n):
	m_box(a_box), blockpair(andp), dmghost(admg), cur_level(n)
	{
		rfnum_procs.resize(sprocs);
		rfnum_start.resize(sprocs);
		dfnum_procs.resize(sprocs);
		dfnum_start.resize(sprocs);
		torefine = false;
		toderefine = false;
		finelevel = NULL;
		coarselevel = NULL;
		twodflag = false;
	}

	AMRLevel & operator=(const AMRLevel & a_level)
	{
		this->m_box = a_level.m_box;
		this->cur_level = a_level.cur_level;
		this->rfnum_procs.resize(sprocs);
		this->rfnum_start.resize(sprocs);
		this->dfnum_procs.resize(sprocs);
		this->dfnum_start.resize(sprocs);
		this->torefine = a_level.torefine;
		this->toderefine = a_level.toderefine;
		return *this;
	}
	void Inittag();

	void AdjustMesh();

	void RefineLevel();

	void Assignnewloc(const bool & holeflag);

	void DerefineLevel();

	void MakeSmallBox();

	void ConnectRefinedBox();

	void RenewGhostIndex();

	void EndMeshHoleConnection();

	void MergeBox();

	void ConnectCoarseBox();

	void clearlevel();

	inline const int & LevelID();

	inline void BoxOddeven(const int & i0, int & sx, int & sy, int & sz)
	{
		sx = m_box[i0].ix() & 1;
		sy = m_box[i0].iy() & 1;
		sz = m_box[i0].iz() & 1;
	}

	inline void BoxOddeven(const int & i0, int & sx, int & sy)
	{
		sx = m_box[i0].ix() & 1;
		sy = m_box[i0].iy() & 1;
	}

	void TagRefineGhost();

	void TagDerefineGhost();

	void TagDerefineGhost_2d();

	void BuildNewPair_Refine();

	void BuildNewPair_Derefine();

	void FillGhostPair_Refine();

	void FillGhostPair_Derefine();

	void RemoveGhostandNormal_Refine();

	void RemoveGhostandNormal_Derefine();

	void RemoveArrayGhost(DataArray<Box> & mybox);

	void IncludeNewboxinGhost_Refine();

	void IncludeNewboxinGhost_Derefine();

	void TrimDetag();

	void RenewProandResIndex();

	void SynGhostTag();

	void RenewLevelProGhostIndex();

	void RenewBlockPairIndex();

	void ReorderBlockpair();

	void AdjustOldPairInFluid_Refine();

	void AdjustOldPairInGhost_Refine();

	void CheckExitingProPair();

	void RemoveExpiredRestriction_Refine();

	void RemoveRefinedBlockCell_Refine();

	void ConstructBlockProPair_Refine();

	void MarkProPairBeforeSyn();

	void AdjustOldPair_Coarsen();

	void RemoveExpiredProlongation_Coarsen();

	void CheckExitingResPair();

	void RemoveExpiredLevelProGhost_Coarsen();

	//void RemoveExpiredLevelProFluid_Coarsen();

	void RemoveFineBlockPair_Coarsen();

	void ConstructNewDomainGhost(vector<Point> & facenmv);

	void RenewGhostRefDmghost(vector<Point> & facenmv);

	void MarkDomainGhgost_Coarsen();

	void MarkResPairBeforeSyn();

	void MarkallPair();

	void ClearExpiredDomainGhost();

	void RenewDomainGhost();

	void AdjustProPair();

	void CheckResProPair_Debug();

	void TellBlockPairNewTag();

	void TellBlockPairNewDetag();

	void ConstructBlockProPair_Coarsen();

	void CheckNewProResPair(CtoFPair & ap);

	void PrintoutData();

	void SynNewGhostTag();

	void CheckNumerofNormalghost();

	void CheckLevelProGhostType();

	void CheckRefineTagLocation();

	void CreatGhostProPair();

	void InitPair();

	void AssociateCelltoBlockPair(const int & mylevel);

	void MarkAllProResPair();

	void MoveBox_Refine(const int & dir, const int & ddir, int & tagnum);

	void MoveBox_Refine_Rotor(const int & dir, const int & ddir, int & tagnum, const int & rs_sep);

	void MoveBox_Derefine(const int & dir, const int & ddir, int & detagnum);

	void MoveBox_Derefine_Rotor(const int & dir, const int & ddir, int & detagnum, const int & rs_sep);

	void cleartemparray()
	{
		newrefinebox.resize(0);
		newcoarsebox.resize(0);
		newrefinebox_ghosttag.resize(0);
		newcoarsebox_ghosttag.resize(0);
		coarsedata.resize(0);
		finedata.resize(0);
		coarsegeom.resize(0);
		finegeom.resize(0);
		rfbox.resize(0);
		derfbox.resize(0);
		taggedbox.resize(0);

		allson.setnum_nocopy(0,0);

		rf_pro_switch_ghost_layer[0].resize(0);
		rf_pro_switch_ghost_layer[1].resize(0);
		rf_pro_keep_f_remove_c.resize(0);
		rf_res_remove_ghost_ccell[0].resize(0);
		rf_res_remove_ghost_ccell[1].resize(0);
		rf_res_move2layer.resize(0);
		rf_recover_pair.resize(0);

		rf_recover_block.resize(0);
		rf_new_fine_bp_ghost.resize(0);
		rf_new_level_pro_ghost.resize(0);
		rf_remove_blockghost.resize(0);
		rf_remove_levelproghost.resize(0);
		rf_oldres_to_pro.resize(0);

		drf_res_switch_ghost_layer[0].resize(0);
		drf_res_switch_ghost_layer[1].resize(0);
		drf_res_keep_c_remove_f[0].resize(0);
		drf_res_keep_c_remove_f[1].resize(0);
		drf_remove_pro_ghost_fcell.resize(0);
		drf_res_move2layer.resize(0);
		drf_recover_pair[0].resize(0);
		drf_recover_pair[1].resize(0);
		drf_remove_respair.resize(0);

		drf_remove_blockghost.resize(0);
		drf_remove_level_pro_ghost.resize(0);
		drf_recover_level_pro_ghost.resize(0);
		drf_new_coarse_bpghost.resize(0);
		drf_new_level_pro_ghost.resize(0);
		
		newfinept.resize(0);
		newcoarsept.resize(0);
		tocoarse.resize(0);
		tofine.resize(0);

		rfnum_procs[srank] = 0;
		dfnum_procs[srank] = 0;
	}

	void MoveBox_Y(const int & ydir);

private:

	DataArray<Box> &  					m_box;
	DataArray<rftag>					m_tag;
	bool 								torefine;
	bool                    			toderefine;
	vector<int>             			rfbox;
	DataArray<Boxson<int> >   	allson;
	vector<int> 								derfbox;
	AMRLevel *      					finelevel;
	AMRLevel *      					coarselevel;
	int             					cur_level;
	vector<int>							rfnum_procs;
	vector<int>							rfnum_start;
	vector<int>							dfnum_start;
	vector<int>							dfnum_procs;
	int 								totrfnum;
	int                 totdfnum;

	vector<int> newrefinebox;
	vector<int> newcoarsebox;
	vector<int> newrefinebox_ghosttag;
	vector<int> newcoarsebox_ghosttag;
	// vector<newpointloc> newcoarsept;
	// vector<newpointloc> newfinept;
	/*New data box for flow variables*/
	vector<FlowVariables>		coarsedata;
	vector<Boxson<FlowVariables> > finedata;
	/*New data box for cell geometry*/
	vector<BoxCellGeom>		coarsegeom;
	vector<Boxson<BoxCellGeom> > finegeom;

	vector<int> taggedbox;

	DataArray<CtoFPair> f_res_ghost[ighost];
	DataArray<CtoFPair> f_pro_ghost;
	DataArray<int> g_pro;

	vector<Boxloc> originbp;
	vector<int> ghost_target;

	DataArray<NodePair> & blockpair;
	DataArray<Domainghost> & dmghost;
	DataArray<CtoFPair> level_pro_ghost;

	vector<CtoFPair> new_f_res_ghost[ighost];
	vector<CtoFPair> new_f_pro_ghost;

	vector<int> rf_pro_switch_ghost_layer[2];
	vector<int> rf_pro_keep_f_remove_c;
	vector<int> rf_res_remove_ghost_ccell[2];
	vector<int> rf_res_move2layer;
	vector<int> rf_recover_pair;
	
	vector<int> rf_new_fine_bp_ghost;
	vector<int> rf_new_level_pro_ghost;
	vector<int> rf_recover_block;
	vector<int> rf_remove_levelproghost;
	vector<int> rf_remove_blockghost;
	vector<int> rf_oldres_to_pro;

	vector<int> drf_res_switch_ghost_layer[2];
	vector<int> drf_res_keep_c_remove_f[2];
	vector<int> drf_remove_pro_ghost_fcell;
	vector<int> drf_res_move2layer;
	vector<int> drf_recover_pair[2];
	vector<CtoFPair> drf_remove_respair;

	vector<int> drf_remove_level_pro_ghost;
	vector<int> drf_recover_level_pro_ghost;
	vector<int> drf_remove_blockghost;
	vector<int> drf_new_coarse_bpghost;
	vector<int> drf_new_level_pro_ghost;
	
	/*new points*/
	vector<mPoint> tocoarse;
	vector<mPoint> tofine;

	vector<newpointloc> newcoarsept;
	vector<newpointloc> newfinept;

	bool twodflag;
	int power_ratio;

};

#endif
