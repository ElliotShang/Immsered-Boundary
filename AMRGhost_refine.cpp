#include "AMRLevel.H"
/*----------------------Notation---------------------------*/
/*----------------------Step of refine---------------------*/
//1. Transfer information about tag of the adjacent blocks;
//-----------function SynGhostTag
//2. Extend the local tagged cells;
//-----------function TagRefineGhost
//3. Modify tags based on previous pairs;
//----->4.1 In the fluid domain:
//--------->4.1.1 some pro pairs will refined and the pair will be removed;
//--------->4.1.2 some pro pairs will be refined and be transfered to res pair;
//--------->4.1.3 some res pairs will be removed. -->RemoveExpiredRestriction_Refine
//--------------->1 Each block needs to remove the local part;
//--------------->2 No need to consider other blocks.
//--------->4.1.4 some level_pro_fluid pairs will also be recovered
//--------->4.1.5 cells to be cleared:
//---------------> 1. which are marked as tag = 0 and not occupied by block pairs; -->BuildNewPair_Refine MarkRefinedBlockPair_Refine
//---------------> 2. mother cell of (4.1.3) and not occupied by block pairs; -->RemoveExpiredRestriction_Refine;
//----->4.2 In the patch ghost domain:
//--------->4.2.1 recover marked cells which are occupied by level_pro_ghost to -1; -->AdjustOldPairInGhost_Refine
//--------->4.2.2 make new level_pro_ghost pair; -->ConstructNewfineBP_Refine
//--------->4.2.3 send the fine cells of the level_pro_ghost pairs to the adjacent blocks to make block pairs; ConstructNewfineBP_Refine
//--------->4.2.4 DONT REMOVE ANY cells of level i in the patch ghost domain.
//4. Construct new pairs;
//5. Make new fine cells;
//6. Fill the new pairs;
//7. Include the sub-cells in the modified cells.
//8. Connect fine cells;

void AMRLevel::TagRefineGhost()
{
	MarkResPairBeforeSyn();
	TellBlockPairNewTag();
	/*---tag = -2 prolongation*/
	/*---tag = 1  first restriction layer*/
	/*---tag = 2  second restriction layer*/
	Point dpt[6] = {Point(0,1,1),Point(2,1,1),Point(1,0,1),Point(1,2,1),Point(1,1,0),Point(1,1,2)};
	MPI_Win_fence(0, m_tag.arraywin());
	for (int i0 = m_tag.ps(); i0 < m_tag.pe(); ++i0)
	{
		if (m_tag[i0].tag == 0)
		{
			//PRINTFinLEVEL("Box %d (%d,%d,%d) has been tagged!!!",cur_level,i0,m_box[i0].ix(),m_box[i0].iy(),m_box[i0].iz());
			for (Point_iterator p(0,3); p.end(); ++p)
			{
				int aneib = m_box[i0].neib[p.i][p.j][p.k];
				if (m_box.isnorg(aneib))
				{
					if (m_tag[aneib].tag == -1 || m_tag[aneib].tag == -2)
					{												
						m_tag[aneib].tag = -2;
						//printf("a -2 tag\n");
						if (m_box[aneib].type != Dmghost)
						{
							m_tag[i0].tag = 1;
						}
					}
				}
				else if (m_box.isnew(aneib))
				{
					PRINTFinLEVEL("Tag 0 Box (%d,%d,%d) neib (%d,%d,%d) is new (%d,%d,%d)", cur_level,
						m_box[i0].ix(), m_box[i0].iy(), m_box[i0].iz(), p.i, p.j, p.k,
						m_box[aneib].ix(), m_box[aneib].iy(), m_box[aneib].iz());
					MPI_Abort(MPI_COMM_WORLD,64);
				}
			}
			if (m_tag[i0].tag == 0)
			{
				for (int di = 0; di < 6; ++di)
				{
					int aneib = m_box[i0].neib[dpt[di][0]][dpt[di][1]][dpt[di][2]];
					if (aneib > -1)
					{
						int aneib2 = m_box[aneib].neib[dpt[di][0]][dpt[di][1]][dpt[di][2]];
						if (m_box.isnorg(aneib2))
						{
							if (m_tag[aneib2].tag == -1 || m_tag[aneib2].tag == -2)
							{
								if (m_box[aneib2].type != Dmghost)
								{
									m_tag[i0].tag = 2;
								}
							}
						}
						else if (m_box.isnew(aneib2))
						{
							Assert(m_box[aneib2].type != Dmghost, "Box type error when expanding the refine tag!!!", 80);
							m_tag[i0].tag = 2;
						}
					}
				}
			}
		}
	}
	MPI_Win_fence(0, m_tag.arraywin());	
	MPI_Barrier(share_comm);
	TellBlockPairNewTag();
}

void AMRLevel::AdjustOldPairInFluid_Refine()
{
	RemoveExpiredRestriction_Refine();
	GiveLevelFlag("Finish RemoveExpiredRestriction_Refine!!!", cur_level, 5);
	RemoveRefinedBlockCell_Refine(); //tag at the block pair locations are set to be 50
	GiveLevelFlag("Finish MarkRefinedBlockCell_Refine!!!", cur_level, 5);
	AdjustOldPairInGhost_Refine();
	CheckExitingProPair(); // tag = 50 is utilized to identify a block pair;
	GiveAFlag("Finish CheckExitingProPair!!!", 5);

}

void AMRLevel::CheckExitingProPair()
{
	for (int i = f_pro_ghost.ps(); i < f_pro_ghost.pe(); ++i)
	{
		int ci0 = f_pro_ghost[i].ci;
		if (m_box.isnorg(ci0))
		{
			if (m_tag[ci0].tag == -2)
			{
				m_tag[ci0].tag = -1;
				rf_recover_pair.push_back(i);
				//PRINTFinLEVEL("rf_recover_pair box (%d,%d,%d)",cur_level,m_box[ci0].ix(),m_box[ci0].iy(),m_box[ci0].iz());
			}
			else if (m_tag[ci0].tag > -1)
			{	
				if (m_tag[ci0].tag == 0)
				{
						rf_pro_keep_f_remove_c.push_back(i);			
				}
				else if (m_tag[ci0].tag == 1 || m_tag[ci0].tag == 2)
				{
					rf_pro_switch_ghost_layer[m_tag[ci0].tag-1].push_back(i);
				}
				m_tag[ci0].tag = -1;
			}
		}
	}
	MPI_Barrier(share_comm);
	int se0 = drf_remove_respair.size();
	if (!finelevel->twodflag)
	{
		for (int i = 0; i < se0; ++i)
		{
			int ci0 = drf_remove_respair[i].ci;
			Assert(m_box.isghost(ci0), "The removed res pair main box should be in the ghost range!!!", 199);
			if (m_tag[ci0].tag == -2)
			{
				m_tag[ci0].tag = -1;
				for (Point_iterator p(0,2); p.end(); ++p)
				{
					finelevel->m_box.ghost_index.push_back(drf_remove_respair[i].fi.son[p.i][p.j][p.k]);
				}
				rf_oldres_to_pro.push_back(i);
			}
			else
			{
				Assert(m_tag[ci0].tag == -1, "A removed res pair tag can only be -1!!!", 210);
				for (Point_iterator p(0,2); p.end(); ++p)
				{
					finelevel->m_box.givehole(drf_remove_respair[i].fi.son[p.i][p.j][p.k]);
				}
			}
		}
	}
	else
	{
		for (int i = 0; i < se0; ++i)
		{
			int ci0 = drf_remove_respair[i].ci;
			Assert(m_box.isghost(ci0), "The removed res pair main box should be in the ghost range!!!", 199);
			if (m_tag[ci0].tag == -2)
			{
				m_tag[ci0].tag = -1;
				for (Point_iterator_2d p(0,2); p.end(); ++p)
				{
					finelevel->m_box.ghost_index.push_back(drf_remove_respair[i].fi.son[p.i][p.j][p.k]);
				}
				rf_oldres_to_pro.push_back(i);
			}
			else
			{
				Assert(m_tag[ci0].tag == -1, "A removed res pair tag can only be -1!!!", 210);
				for (Point_iterator_2d p(0,2); p.end(); ++p)
				{
					finelevel->m_box.givehole(drf_remove_respair[i].fi.son[p.i][p.j][p.k]);
				}
			}
		}
	}
	MPI_Barrier(share_comm);
}


void AMRLevel::RemoveExpiredRestriction_Refine()
{
	Point ap[6] = {Point(0,1,1),Point(2,1,1),Point(1,0,1),Point(1,2,1),Point(1,1,0),Point(1,1,2)};
	//int keepnum = 0; int resholenum = 0; int move2num = 0;
	for (int i = f_res_ghost[0].ps(); i < f_res_ghost[0].pe(); ++i)
	{
		int ci0 = f_res_ghost[0][i].ci;
		if (m_box.isghost(ci0))
		{
			Assert(m_tag[ci0].tag == -20, "Error when removing res pair", 165);
			// bool removeflag = true;
			for (Point_iterator pn(0,3); pn.end(); ++pn)
			{
				int aneib = m_box[ci0].neib[pn.i][pn.j][pn.k];
				if (m_box.isnorg(aneib)) //not new cell from coarsening 
				{
					if ((m_tag[aneib].tag == -1 || m_tag[aneib].tag == -2) && m_box[aneib].type != Dmghost)
					{
						// removeflag = false;
						//keepnum++;
						// PRINTFinLEVEL("Res (%d,%d,%d) neib (%d,%d,%d)",cur_level,
						// 			m_box[ci0].ix(),m_box[ci0].iy(),m_box[ci0].iz(),
						// 			m_box[aneib].ix(),m_box[aneib].iy(),m_box[aneib].iz());
						// break;
						goto NEXTRES;
					}
				}
				else if (m_box.isnew(aneib))
				{
					// removeflag = false;
					goto NEXTRES;
				}
			}
			// if (removeflag)
			{
				bool copyflag = false;
				for (int di = 0; di < 6; ++di)
				{
					int aneib = m_box[ci0].neib[ap[di][0]][ap[di][1]][ap[di][2]];
					if (aneib > -1)
					{
						int an1 = m_box[aneib].neib[ap[di][0]][ap[di][1]][ap[di][2]];
						//Assert(an1 > -1, "A negative index for res 0 neib!!!", 182);
						if (m_box.isnorg(an1))
						{
							if ((m_tag[an1].tag == -1 || m_tag[an1].tag == -2) && m_box[an1].type != Dmghost)
							{
								copyflag = true;
								// PRINTFinLEVEL("Res copy (%d,%d,%d) neib (%d,%d,%d)",cur_level,
								// 	m_box[ci0].ix(),m_box[ci0].iy(),m_box[ci0].iz(),
								// 	m_box[an1].ix(),m_box[an1].iy(),m_box[an1].iz());
							}
						}
						else if (m_box.isnew(an1))
						{
							copyflag = true;
						}
					}
				}
				if (copyflag)
				{
					rf_res_move2layer.push_back(i);
					//move2num++;
				}
				else
				{
					m_tag[ci0].tag = -21;
					rf_res_remove_ghost_ccell[0].push_back(i);
					//resholenum++;
				}
			}
			NEXTRES:;
		}
	}
	MPI_Barrier(share_comm);
	for (int i = f_res_ghost[1].ps(); i < f_res_ghost[1].pe(); ++i)
	{
		int ci0 = f_res_ghost[1][i].ci;
		if (m_box.isghost(ci0))
		{
			Assert(m_tag[ci0].tag = -20, "Error when removing res pair", 216);
			bool removeflag = true;
			for (int di = 0; di < 6; ++di)
			{
				int aneib = m_box[ci0].neib[ap[di][0]][ap[di][1]][ap[di][2]];
				if (aneib > -1)
				{
					int aneib2 = m_box[aneib].neib[ap[di][0]][ap[di][1]][ap[di][2]];
					if (m_box.isnorg(aneib2))
					{
						if ((m_tag[aneib2].tag == -1 || m_tag[aneib2].tag == -2) && m_box[aneib2].type != Dmghost)
						{
							removeflag = false;
							break;
						}
					}
					else if (m_box.isnew(aneib2))
					{
						removeflag = false;
						break;
					}
				}
			}
			if (removeflag)
			{
				//DONTREMOVEBLOCKPOINT
				m_tag[ci0].tag = -21;
				rf_res_remove_ghost_ccell[1].push_back(i);
				//printf("Remove res 1 is (%d,%d,%d)\n", m_box[ci0].ix(), m_box[ci0].iy(), m_box[ci0].iz());
			}
		}
	}
	GiveAFlag("Finish second layer!!!",5);
	MPI_Barrier(share_comm);
	TellBlockPairNewTag();
	// for (int i = m_box.gps(); i < m_box.gpe(); ++i)
	// {
	// 	int ix = m_box[i].ix();
	// 	int iy = m_box[i].iy();
	// 	int iz = m_box[i].iz();
	// 	if (ix == 34 && iy == 0 && iz == 0)
	// 	{
	// 		printf("My tag is %d\n", m_tag[i].tag);
	// 	}
	// }
#ifdef DEBUG
	for (int ig = 0; ig < ighost; ++ig)
	{
		for (int i = f_res_ghost[ig].ps(); i < f_res_ghost[ig].pe(); ++i)
		{
			int ci0 = f_res_ghost[ig][i].ci;
			if (m_box.isghost(ci0))
			{
				if (m_tag[ci0].tag != -20 && m_tag[ci0].tag != -21)
				{
					PRINTFinLEVEL("The res tag at box (%d,%d,%d) has been changed to %d during refine!!!",
						cur_level,m_box[ci0].ix(),m_box[ci0].iy(),m_box[ci0].iz(),m_tag[ci0].tag);
					MPI_Abort(MPI_COMM_WORLD,267);
				}
			}
		}
	}
#endif			
}


void AMRLevel::AdjustOldPairInGhost_Refine()
{
	Point ap[6] = {Point(0,1,1),Point(2,1,1),Point(1,0,1),Point(1,2,1),Point(1,1,0),Point(1,1,2)};
	for (int i = level_pro_ghost.ps(); i < level_pro_ghost.pe(); ++i)
	{
		int i0 = level_pro_ghost[i].ci;
		if (m_box.isghost(i0))
		{
			if (m_tag[i0].tag == -1 || m_tag[i0].tag == -20)
			{
				goto NEXTPAIR;
			}
			if (m_tag[i0].tag == -21)
			{
				rf_remove_levelproghost.push_back(i);
			}
			else if (m_tag[i0].tag != -1 && m_tag[i0].tag != -20)
			{
				if (m_tag[i0].tag == 0)
				{
					rf_remove_levelproghost.push_back(i);
				}
				rf_recover_block.push_back(i);
				m_tag[i0].tag = -1;
			}
			NEXTPAIR:;
		}
	}
	MPI_Barrier(share_comm);
}
/*-----This function executed after tag extension---------*/
void AMRLevel::RemoveRefinedBlockCell_Refine()
{
	for (int n0 = blockpair.ps(); n0 < blockpair.pe(); ++n0)
	{
		NodePair * nptr = &blockpair[n0];
		int a0 = nptr->innode;
		if (m_box.isnorg(a0))
		{
			if (m_tag[a0].tag == 0 || m_tag[a0].tag == -21)
			{
				blockpair.givehole(n0);
			}
		}
	}
}

void AMRLevel::BuildNewPair_Refine()
{
	//MPI_Barrier(MPI_COMM_WORLD);
	MPI_Barrier(share_comm);
	int mbps = m_box.ps();
	int mbpe = m_box.pe();
	for (int i = mbps; i < mbpe; ++i)
	{
		if (m_tag[i].tag == -2)
		{
			new_f_pro_ghost.push_back(CtoFPair(i));
			taggedbox.push_back(i);
		}
		else if (m_tag[i].tag == 1 || m_tag[i].tag == 2)
		{
			//PRINTFinLEVEL("box %d ", cur_level, i);
			int r0 = m_tag[i].tag-1;
			new_f_res_ghost[r0].push_back(CtoFPair(i));
			m_tag[i].tag = 0;
			m_box.ghost_index.push_back(i);
			//PRINTFinLEVEL("a box %d ", cur_level, i);
		}
		else if (m_tag[i].tag == 0)
		{
			m_box.givehole(i);
		}
	}
	MPI_Barrier(share_comm);
	GiveLevelFlag("Finish check refine tag!!!", cur_level, 5);
	int mbgps = m_box.gps();
	int mbgpe = m_box.gpe();
	for (int i = mbgps; i < mbgpe; ++i)
	{
		if (m_tag[i].tag == -1)
		{
			goto NEXTCELL;
		}
		else
		{
			if (m_tag[i].tag == 0)
			{
				Assert(m_box[i].type == Blockghost, "A ghost with tag = 0 can not be a normal cell!!!", 576);
				rf_remove_blockghost.push_back(i);
				m_tag[i].tag = -2;
				taggedbox.push_back(i);
			}
			else if (m_tag[i].tag != -20 && m_tag[i].tag != -21)
			{
				if (m_box[i].type == Blockghost || m_box[i].type == Dmghost)
				{
					//PRINTFinLEVEL("new level_pro_ghost %d tag size %d", cur_level, i, m_tag[i].tag);
					rf_new_level_pro_ghost.push_back(i);
				}
				/*Another part of cells are son of old prolongation*/
				else
				{
					// PRINTFinLEVEL("This situation is not possible tag is %d!!!", cur_level, m_tag[i].tag);
					// MPI_Abort(MPI_COMM_WORLD, 634);
					Assert(m_tag[i].tag==-2 && m_box[i].type == Normalcell, "A tagged normal ghost should be -2!!!", 591);
					new_f_pro_ghost.push_back(CtoFPair(i));
				}
				taggedbox.push_back(i);	
				m_tag[i].tag = -2;
			}
		}
		NEXTCELL:;
	}
	MPI_Barrier(MPI_COMM_WORLD);
	for (int i = 0; i < ghost_target.size(); ++i)
	{
		int i0 = ghost_target[i];
		if (m_box.isnorg(i0))
		{
			//if (m_tag[i0].tag != -1 && m_tag[i0].tag != -20)
			if (m_tag[i0].tag == -2)
			{
				rf_new_fine_bp_ghost.push_back(i);
				// PRINTFinLEVEL("new block pair in ghost (%d,%d,%d)",
				// 	cur_level, 
				// 	m_box[i0].ix(),
				// 	m_box[i0].iy(),
				// 	m_box[i0].iz());
			}
		}
	}
	
	MPI_Barrier(share_comm);
	PRINTFinLEVEL("Refine new pro pair is %d",cur_level,(int)new_f_pro_ghost.size());
	PRINTFinLEVEL("Refine new res 1 pair is %d",cur_level, (int)new_f_res_ghost[0].size());
	PRINTFinLEVEL("Refine new res 2 pair is %d",cur_level, (int)new_f_res_ghost[1].size());
	PRINTFinLEVEL("Refine ghost number is %d before add",cur_level, (int)m_box.ghost_index.size());
	PRINTFinLEVEL("Refine new block pair in ghost %d",cur_level,(int)rf_new_fine_bp_ghost.size());
	PRINTFinLEVEL("Refine new level pro ghost %d",cur_level,(int)rf_new_level_pro_ghost.size());
}

void AMRLevel::ConstructBlockProPair_Refine()
{
	vector<vector<S_newfinebp0> > senddata(nodenum);
	vector<S_newfinebp0> recvdata;
	vector<int> recvnode;
	//printf("Node %d start new_level_pro_ghost\n", node);
	/*------------------------------------------------------------*/
	int new_level_pro_size = rf_new_level_pro_ghost.size();
	vector<CtoFPair> new_level_pro_ghost(new_level_pro_size);
	for (int i = 0; i < new_level_pro_size; ++i)
	{
		int i0 = rf_new_level_pro_ghost[i];
		new_level_pro_ghost[i].ci = i0;
		Assert(m_tag[i0].tag > -1, "A level pro ghost has a invalid refine tag!!! ", 444);
		new_level_pro_ghost[i].fi = allson[m_tag[i0].tag];
	}
	/*--------------From the newly refined normal cells-----------*/
	int fine_bp_size = rf_new_fine_bp_ghost.size();
	if (!finelevel->twodflag)
	{
		for (int i = 0; i < fine_bp_size; ++i)
		{
			int i0 = rf_new_fine_bp_ghost[i];
			int g0 = ghost_target[i0];
			int tag0 = m_tag[g0].tag;
			Assert(tag0 > -1, "The new block pair was not tagged!!!", 620);
			int node_dest = originbp[i0].node;
			senddata[node_dest].push_back(S_newfinebp0(originbp[i0].index,allson[tag0]));
			for (Point_iterator p(0,2); p.end(); ++p)
			{
				int fi0 = allson[tag0].son[p.i][p.j][p.k];
				finelevel->m_box[fi0].type = Blockghost;
#ifdef DEBUG
				senddata[node_dest].back().sonpt.son[p.i][p.j][p.k] = finelevel->m_box[fi0].lowpt;
#endif						
			}
		}
	}
	else
	{
		for (int i = 0; i < fine_bp_size; ++i)
		{
			int i0 = rf_new_fine_bp_ghost[i];
			int g0 = ghost_target[i0];
			int tag0 = m_tag[g0].tag;
			Assert(tag0 > -1, "The new block pair was not tagged!!!", 620);
			int node_dest = originbp[i0].node;
			senddata[node_dest].push_back(S_newfinebp0(originbp[i0].index,allson[tag0]));
			for (Point_iterator_2d p(0,2); p.end(); ++p)
			{
				int fi0 = allson[tag0].son[p.i][p.j][p.k];
				finelevel->m_box[fi0].type = Blockghost;
#ifdef DEBUG
				senddata[node_dest].back().sonpt.son[p.i][p.j][p.k] = finelevel->m_box[fi0].lowpt;
#endif						
			}
		}
	}
	PRINTFinLEVEL("new level pro ghost number %d fine block pair size %d",cur_level,new_level_pro_size,fine_bp_size);
	/*------------------------------------------------------------*/
	BcastNewInfo_NodeRank(senddata, recvdata, recvnode, 5);
	int rd = recvdata.size();
	if (!finelevel->twodflag)
	{
		vector<NodePair> new_bkp(8*rd);
		for (int i = 0; i < rd; ++i)
		{
			int i0 = recvdata[i].aimbkp;
			int c0 = blockpair[i0].innode;
			int t0 = m_tag[c0].tag;
			Assert(t0 > -1, "error in build new block pair!!!", 75);
			int s = 0;
			for (Point_iterator p(0,2); p.end(); ++p)
			{
				int fbkp = i*8+s;
				new_bkp[fbkp].innode = allson[t0].son[p.i][p.j][p.k];
#ifdef DEBUG			
				int rf0 = pow(2, cur_level+1);
				Point s1 = finelevel->m_box[allson[t0].son[p.i][p.j][p.k]].lowpt;
				Point s2 = recvdata[i].sonpt.son[p.i][p.j][p.k];
#if DIM == 3			
				Assert((s1[0]==s2[0] || (abs(s1[0]-s2[0]) == highpt.xy[0]*rf0 && periodic[0])) && 
							 (s1[1]==s2[1] || (abs(s1[1]-s2[1]) == highpt.xy[1]*rf0 && periodic[1])) &&
							 (s1[2]==s2[2] || (abs(s1[2]-s2[2]) == highpt.xy[2]*finelevel->power_ratio && periodic[2])),
							 "The new block pair during refining does not match!!!", 657);
#elif DIM == 2
				Assert((s1[0]==s2[0] || (abs(s1[0]-s2[0]) == highpt.xy[0]*rf0 && periodic[0])) && 
							 (s1[1]==s2[1] || (abs(s1[1]-s2[1]) == highpt.xy[1]*rf0 && periodic[1])),
							 "The new block pair during refining does not match!!!", 657);
#endif						 			
#endif						
				new_bkp[fbkp].outnode.node = recvnode[i];
				new_bkp[fbkp].outnode.index = recvdata[i].son0.son[p.i][p.j][p.k];
#ifdef PASSAGE_ANGLE				
				new_bkp[fbkp].theta = blockpair[i0].theta;
#endif				
				++s;
			}
		}
		finelevel->blockpair.Addnew(new_bkp);
	}
	else
	{
		vector<NodePair> new_bkp(4*rd);
		for (int i = 0; i < rd; ++i)
		{
			int i0 = recvdata[i].aimbkp;
			int c0 = blockpair[i0].innode;
			int t0 = m_tag[c0].tag;
			Assert(t0 > -1, "error in build new block pair!!!", 75);
			int s = 0;
			for (Point_iterator_2d p(0,2); p.end(); ++p)
			{
				int fbkp = i*4+s;
				new_bkp[fbkp].innode = allson[t0].son[p.i][p.j][p.k];
#ifdef DEBUG			
				int rf0 = pow(2, cur_level+1);
				Point s1 = finelevel->m_box[allson[t0].son[p.i][p.j][p.k]].lowpt;
				Point s2 = recvdata[i].sonpt.son[p.i][p.j][p.k];
#if DIM == 3			
				Assert((s1[0]==s2[0] || (abs(s1[0]-s2[0]) == highpt.xy[0]*rf0 && periodic[0])) && 
							 (s1[1]==s2[1] || (abs(s1[1]-s2[1]) == highpt.xy[1]*rf0 && periodic[1])) &&
							 (s1[2]==s2[2] || (abs(s1[2]-s2[2]) == highpt.xy[2]*finelevel->power_ratio && periodic[2])),
							 "The new block pair during refining does not match!!!", 657);
#elif DIM == 2
				Assert((s1[0]==s2[0] || (abs(s1[0]-s2[0]) == highpt.xy[0]*rf0 && periodic[0])) && 
							 (s1[1]==s2[1] || (abs(s1[1]-s2[1]) == highpt.xy[1]*rf0 && periodic[1])),
							 "The new block pair during refining does not match!!!", 657);
#endif						 			
#endif						
				new_bkp[fbkp].outnode.node = recvnode[i];
				new_bkp[fbkp].outnode.index = recvdata[i].son0.son[p.i][p.j][p.k];
#ifdef PASSAGE_ANGLE				
				new_bkp[fbkp].theta = blockpair[i0].theta;
#endif				
				++s;
			}
		}
		finelevel->blockpair.Addnew(new_bkp);
	}
	finelevel->SynNewGhostTag();
	finelevel->blockpair.DirectlyReduceNew();
	level_pro_ghost.Addnew(new_level_pro_ghost);
	level_pro_ghost.DirectlyReduceNew();
	GiveLevelFlag("finish build ghost pro ghost", cur_level, 5);
}


