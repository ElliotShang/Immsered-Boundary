#include "AMR.H"

void AMRLevel::SynGhostTag()
{
	vector<vector<S_tag> > local_tag(nodenum);
	vector<S_tag> recvtag;
	vector<int> recvnode;
	// if (option == 0)
	// {
		for (int n0 = blockpair.ps(); n0 < blockpair.pe(); ++n0)
		{
			NodePair * nptr = &blockpair[n0];
			local_tag[nptr->outnode.node].push_back(S_tag(nptr->outnode.index, n0));
			if (nptr->outnode.index == 6916839)
			{
				printf("N%dR%d the %d blockpair: outnode %d remotecell %d local cell %d (%d,%d,%d)\n",
					node, srank, n0, nptr->outnode.node, nptr->outnode.index, nptr->innode,
					m_box[nptr->innode].ix(),
					m_box[nptr->innode].iy(),
					m_box[nptr->innode].iz());
			}
#ifdef DEBUG
			local_tag[nptr->outnode.node].back().bxyz = m_box[nptr->innode].lowpt;
#endif						
			// local_tag[nptr->outnode.node].back().aimcell = nptr->outnode.index;
			// local_tag[nptr->outnode.node].back().localbkp = n0;
			// local_tag[nptr->outnode.node].back().bkptag[0] = m_tag[nptr->innode].tag;
			// local_tag[nptr->outnode.node].back().bkptag[1] = m_tag[nptr->innode].detag;
			// if (m_tag[nptr->innode].tag != -1)
			// {
			// 	printf("N%d Block pair box %d (%d,%d,%d) SEND tag %d\n",
			// 		node,nptr->innode,
			// 		m_box[nptr->innode].ix(),
			// 		m_box[nptr->innode].iy(),
			// 		m_box[nptr->innode].iz(),
			// 		m_tag[nptr->innode].tag);
			// }
		}
		//printf("Start to BcastNewInfo_Node!!!\n");
		BcastNewInfo_NodeRank(local_tag, recvtag, recvnode, 3);
		//printf("Finish BcastNewInfo_Node!!!\n");
		int s = 0;
		int recvsizetag = recvtag.size();
		// while (s < recvsizetag)
		// {
		// 	m_tag[recvtag[s].aimcell].tag = recvtag[s].bkptag[0];
		// 	m_tag[recvtag[s].aimcell].detag = recvtag[s].bkptag[1];
		// 	// if (m_tag[recvtag[s].aimcell].tag != -1 && ts == 1)
		// 	// {
		// 	// 	PRINTFinLEVEL("Block pair box(%d,%d,%d) receive tag %d\n",
		// 	// 		cur_level,
		// 	// 		m_box[recvtag[s].aimcell].ix(),
		// 	// 		m_box[recvtag[s].aimcell].iy(),
		// 	// 		m_box[recvtag[s].aimcell].iz(),
		// 	// 		m_tag[recvtag[s].aimcell].tag);
		// 	// }
		// 	++s;
		// }	
	// }
	// else if (option == 1)
	// {
	// 	for (int n0 = bkpair.ps(); n0 < bkpair.pe(); ++n0)
	// 	{
	// 		NodePair * nptr = &bkpair[n0];
	// 		local_tag[nptr->outnode.node].push_back(S_tag());
	// 		local_tag[nptr->outnode.node].back().aimcell = nptr->outnode.index;
	// 		local_tag[nptr->outnode.node].back().localbkp = n0;
	// 		local_tag[nptr->outnode.node].back().bkptag = m_tag[nptr->innode].tag;
	// 		BcastNewInfo_Node(local_tag, recvtag, nodecomm, recvnode_tag);
	// 		int s = 0;
	// 		recvsizetag = recvtag.size();
	// 		while (s < recvsizetag)
	// 		{
	// 			m_tag[recvtag[s].aimcell].detag = recvtag[s].bkptag; 
	// 			++s;
	// 		}	
	// 	}
	// }
	originbp.resize(recvsizetag);
	ghost_target.resize(recvsizetag);
	// s = 0;
	while(s < recvsizetag)
	{
		ghost_target[s] = recvtag[s].aimcell;
		// if (!m_box.isghost(ghost_target[s]))
		// {
		// 	PRINTFinLEVEL("The ghost target is not ghost box is %d (%d,%d,%d) gps %d gpe %d",
		// 		cur_level,ghost_target[s],
		// 		m_box[ghost_target[s]].ix(),
		// 		m_box[ghost_target[s]].iy(),
		// 		m_box[ghost_target[s]].iz(),
		// 		m_box.gps(), m_box.gpe());
		// }		
		Assert(m_box.isghost(ghost_target[s]) || m_box.isnew(ghost_target[s]), "A ghost target is not a ghost!!!", 74);
		Assert(m_box[ghost_target[s]].type == Blockghost, "A ghost target is not a block ghost!!!", 74);
		// if (!((m_box[ghost_target[s]].lowpt[0] == recvtag[s].bxyz[0] ||
		// 			 (abs(m_box[ghost_target[s]].lowpt[0]-recvtag[s].bxyz[0]) == highpt.xy[0]*pow(2,cur_level) && periodic[0])) &&
		// 			 (m_box[ghost_target[s]].lowpt[1] == recvtag[s].bxyz[1] ||
		// 			 (abs(m_box[ghost_target[s]].lowpt[1]-recvtag[s].bxyz[1]) == highpt.xy[1]*pow(2,cur_level) && periodic[1])) &&
		// 			 (m_box[ghost_target[s]].lowpt[2] == recvtag[s].bxyz[2] ||
		// 			 (abs(m_box[ghost_target[s]].lowpt[2]-recvtag[s].bxyz[2]) == highpt.xy[2]*pow(2,cur_level) && periodic[2]))))
		// {
		// 	printf("local ghost is (%d,%d,%d) and received ghost is (%d,%d,%d)\n",
		// 		m_box[ghost_target[s]].lowpt[0],m_box[ghost_target[s]].lowpt[1],m_box[ghost_target[s]].lowpt[2],
		// 		recvtag[s].bxyz[0],recvtag[s].bxyz[1],recvtag[s].bxyz[2]);
		// }
#if DIM == 3
#ifdef DEBUG		
		if (!((m_box[ghost_target[s]].lowpt[0] == recvtag[s].bxyz[0] ||
			(abs(m_box[ghost_target[s]].lowpt[0]-recvtag[s].bxyz[0]) == highpt.xy[0]*pow(2,cur_level) && periodic[0])) &&
			(m_box[ghost_target[s]].lowpt[1] == recvtag[s].bxyz[1] ||
			(abs(m_box[ghost_target[s]].lowpt[1]-recvtag[s].bxyz[1]) == highpt.xy[1]*pow(2,cur_level) && periodic[1])) &&
			(m_box[ghost_target[s]].lowpt[2] == recvtag[s].bxyz[2] ||
			(abs(m_box[ghost_target[s]].lowpt[2]-recvtag[s].bxyz[2]) == highpt.xy[2]*power_ratio && periodic[2]))))
			{
				printf("Level%d Ghost target is Box(%d,%d,%d) but the recv ghost xyz is (%d,%d,%d)!!!\n",
					power_ratio, m_box[ghost_target[s]].lowpt[0], m_box[ghost_target[s]].lowpt[1], m_box[ghost_target[s]].lowpt[2],
					recvtag[s].bxyz[0],recvtag[s].bxyz[1],recvtag[s].bxyz[2]);
				MPI_Abort(MPI_COMM_WORLD, 128);
			}
#endif						
		Assert((m_box[ghost_target[s]].lowpt[0] == recvtag[s].bxyz[0] ||
					 (abs(m_box[ghost_target[s]].lowpt[0]-recvtag[s].bxyz[0]) == highpt.xy[0]*pow(2,cur_level) && periodic[0])) &&
					 (m_box[ghost_target[s]].lowpt[1] == recvtag[s].bxyz[1] ||
					 (abs(m_box[ghost_target[s]].lowpt[1]-recvtag[s].bxyz[1]) == highpt.xy[1]*pow(2,cur_level) && periodic[1])) &&
					 (m_box[ghost_target[s]].lowpt[2] == recvtag[s].bxyz[2] ||
					 (abs(m_box[ghost_target[s]].lowpt[2]-recvtag[s].bxyz[2]) == highpt.xy[2]*power_ratio && periodic[2])),
					 "A ghost target xyz does not match the original box!!!", 104);
#elif DIM == 2
		Assert((m_box[ghost_target[s]].lowpt[0] == recvtag[s].bxyz[0] ||
					 (abs(m_box[ghost_target[s]].lowpt[0]-recvtag[s].bxyz[0]) == highpt.xy[0]*pow(2,cur_level) && periodic[0])) &&
					 (m_box[ghost_target[s]].lowpt[1] == recvtag[s].bxyz[1] ||
					 (abs(m_box[ghost_target[s]].lowpt[1]-recvtag[s].bxyz[1]) == highpt.xy[1]*pow(2,cur_level) && periodic[1])),
					 "A ghost target xyz does not match the original box!!!", 104);
#endif		
		// PRINTFinLEVEL("Received ghost_target %d is box %d (%d,%d,%d)",
		// 	cur_level,s,ghost_target[s],m_box[ghost_target[s]].ix(),m_box[ghost_target[s]].iy(),m_box[ghost_target[s]].iz());
		originbp[s] = Boxloc(recvnode[s], recvtag[s].localbkp);
		++s;
		// if (node == 0)
		// {
		// 	if (m_box[ghost_target[s-1]].ix() < 20)
		// 	{
		// 		printf("syn node 0 ghost_target (%d,%d,%d) is not a ghost!!!\n",
		// 			m_box[ghost_target[s-1]].ix(),
		// 			m_box[ghost_target[s-1]].iy(),
		// 			m_box[ghost_target[s-1]].iz());
		// 	}
		// 	Assert(m_box[ghost_target[s-1]].ix() > 19, "node 0 error", 758);
		// }
		// else if (node == 1)
		// {
		// 	Assert(m_box[ghost_target[s-1]].ix() < 20, "syn node 1 ghost_target is not a ghost!!!", 762);
		// }
	}
}

void AMRLevel::SynNewGhostTag()
{
	vector<vector<S_tag> > local_tag(nodenum);
	vector<S_tag> recvtag;
	vector<int> recvnode;
	int bs,be;
	ArrayOrder_s(blockpair.newstart(), blockpair.newend(), bs, be, sprocs, srank);
	for (int n0 = bs; n0 < be; ++n0)
	{
		NodePair * nptr = &blockpair[n0];
		local_tag[nptr->outnode.node].push_back(S_tag(nptr->outnode.index, n0));
#ifdef DEBUG
		local_tag[nptr->outnode.node].back().bxyz = m_box[nptr->innode].lowpt;
#endif						
	}
	//printf("Start to BcastNewInfo_Node!!!\n");
	BcastNewInfo_NodeRank(local_tag, recvtag, recvnode, 3);
	//printf("Finish BcastNewInfo_Node!!!\n");
	int s = 0;
	int recvsizetag = recvtag.size();
	// s = 0;
	int gsnum = ghost_target.size();
	while(s < recvsizetag)
	{
		ghost_target.push_back(recvtag[s].aimcell);
#ifdef DEBUG
		Point s1 = m_box[ghost_target[s+gsnum]].lowpt;
		Point s2 = recvtag[s].bxyz;
		int rf0 = pow(2, cur_level);
#if DIM == 3		
		Assert((s1[0]==s2[0] || (abs(s1[0]-s2[0]) == highpt.xy[0]*rf0 && periodic[0])) && 
					 (s1[1]==s2[1] || (abs(s1[1]-s2[1]) == highpt.xy[1]*rf0 && periodic[1])) &&
					 (s1[2]==s2[2] || (abs(s1[2]-s2[2]) == highpt.xy[2]*power_ratio && periodic[2])), 
					 "A ghost target xyz does not match the original box!!!", 104);
#elif DIM == 2
		Assert((s1[0]==s2[0] || (abs(s1[0]-s2[0]) == highpt.xy[0]*rf0 && periodic[0])) && 
					 (s1[1]==s2[1] || (abs(s1[1]-s2[1]) == highpt.xy[1]*rf0 && periodic[1])), 
					 "A ghost target xyz does not match the original box!!!", 104);
#endif					 	
#endif		
		// PRINTFinLEVEL("Received ghost_target %d is box %d (%d,%d,%d)",
		// 	cur_level,s,ghost_target[s],m_box[ghost_target[s]].ix(),m_box[ghost_target[s]].iy(),m_box[ghost_target[s]].iz());
		originbp.push_back(Boxloc(recvnode[s], recvtag[s].localbkp));
		++s;
	}
}


void AMRLevel::FillGhostPair_Refine()
{
	for (int i = 0; i < new_f_pro_ghost.size(); ++i)
	{
		int ci0 = new_f_pro_ghost[i].ci;
		int mytag = m_tag[ci0].tag;
		Assert(mytag > -1, "The ghost cell has not been tagged!!!", 85);
		new_f_pro_ghost[i].fi = allson[mytag];
	}
	for (int ig = 0; ig < ighost; ++ig)
	{
		for (int i = 0; i < new_f_res_ghost[ig].size(); ++i)
		{
			int ci0 = new_f_res_ghost[ig][i].ci;
			int mytag = m_tag[ci0].tag;
			Assert(mytag > -1, "The ghost cell has not been tagged!!!", 85);
			new_f_res_ghost[ig][i].fi = allson[mytag];
		}
	}
	int se225 = rf_oldres_to_pro.size();
	for (int i = 0; i < se225; ++i)
	{
		new_f_pro_ghost.push_back(drf_remove_respair[rf_oldres_to_pro[i]]);
	}
	// printf("Level %d prolongation pair is %d\n", cur_level, (int)f_pro_ghost.size());
	// printf("Level %d restriction pair 1 is %d\n", cur_level, (int)f_res_ghost[0].size());
	// printf("Level %d restriction pair 2 is %d\n", cur_level, (int)f_res_ghost[1].size());
}

void AMRLevel::FillGhostPair_Derefine()
{
	for (int i = 0; i < new_f_pro_ghost.size(); ++i)
	{
		int fi0 = new_f_pro_ghost[i].fi.son[0][0][0];
		new_f_pro_ghost[i].ci = finelevel->m_tag[fi0].detag;
		//printf("new pro derf %d ci is %d\n", i, f_pro_ghost[i0].ci);
		Assert(finelevel->m_tag[fi0].detag > -1, "Negative index for prolongation ghost of derefine!!!", 366);
#ifdef DEBUG
		CheckNewProResPair(new_f_pro_ghost[i]);
#endif		
	}
	GiveLevelFlag("Finish Pro part in FillGhostPair_Derefine!!!",cur_level,5);
	for (int ig = 0; ig < ighost; ++ig)
	{
		for (int i = 0; i < new_f_res_ghost[ig].size(); ++i)
		{
			int fi0 = new_f_res_ghost[ig][i].fi.son[0][0][0];
			new_f_res_ghost[ig][i].ci = finelevel->m_tag[fi0].detag;
			int ci0 = new_f_res_ghost[ig][i].ci;
			Assert(finelevel->m_tag[fi0].detag > -1, "Negative index for restriction ghost of derefine!!!", 366);
#ifdef DEBUG
			CheckNewProResPair(new_f_res_ghost[ig][i]);
#endif						
		}
	}
	GiveLevelFlag("Finish Res part in FillGhostPair_Derefine!!!",cur_level,5);
	ConstructBlockProPair_Coarsen();
}

void AMRLevel::RemoveGhostandNormal_Refine()
{
	if (NULL != finelevel)
	{
		//MPI_Win_fence(0, finelevel->m_box.arraywin());
		for (int ig = 0; ig < ighost; ++ig)
		{
			for (int i = 0; i < rf_res_remove_ghost_ccell[ig].size(); ++i)
			{
				int i0 = rf_res_remove_ghost_ccell[ig][i];
				m_box.givehole(f_res_ghost[ig][i0].ci);
				m_box[f_res_ghost[ig][i0].ci].neib[1][1][1] = -1;
				f_res_ghost[ig].givehole(i0);
				/*---The holes have been given in determination---*/
			}
		}
		for (int i = 0; i < rf_res_move2layer.size(); ++i)
		{
			int i0 = rf_res_move2layer[i];
			new_f_res_ghost[1].push_back(f_res_ghost[0][i0]);
			f_res_ghost[0].givehole(i0);
		}
		for (int i = 0; i < rf_remove_blockghost.size(); ++i)
		{
			int i0 = rf_remove_blockghost[i];
			m_box.givehole(i0);
			m_box[i0].neib[1][1][1] = -1;	
		}
		for (int i = 0; i < rf_remove_levelproghost.size(); ++i)
		{
			int i0 = rf_remove_levelproghost[i];
			level_pro_ghost.givehole(i0);
			m_box.givehole(level_pro_ghost[i0].ci);
			m_box[level_pro_ghost[i0].ci].neib[1][1][1] = -1;
		}
		/*-----------------------------------------------------------------------*/
		if (!finelevel->twodflag)
		{
			for (int ig = 0; ig < ighost; ++ig)
			{
				for (int i = 0; i < rf_pro_switch_ghost_layer[ig].size(); ++i)
				{
					int i0 = rf_pro_switch_ghost_layer[ig][i];
					int ci0 = f_pro_ghost[i0].ci;
					m_box.ghost_index.push_back(ci0);
					for (Point_iterator p(0,2); p.end(); ++p)
					{
						int fi0 = f_pro_ghost[i0].fi.son[p.i][p.j][p.k];
						finelevel->m_box[fi0].neib[1][1][1] = -1;
					}
					f_pro_ghost.givehole(i0);
					new_f_res_ghost[ig].push_back(f_pro_ghost[i0]);
				}
			}
			for (int i = 0; i < rf_pro_keep_f_remove_c.size(); ++i)
			{
				int i0 = rf_pro_keep_f_remove_c[i];
				int ci0 = f_pro_ghost[i0].ci;
				m_box.givehole(ci0);
				for (Point_iterator p(0,2); p.end(); ++p)
				{
					int fi0 = f_pro_ghost[i0].fi.son[p.i][p.j][p.k];
					finelevel->m_box[fi0].neib[1][1][1] = -1;
				}
				f_pro_ghost.givehole(i0);
			}
		}
		else
		{
			for (int ig = 0; ig < ighost; ++ig)
			{
				for (int i = 0; i < rf_pro_switch_ghost_layer[ig].size(); ++i)
				{
					int i0 = rf_pro_switch_ghost_layer[ig][i];
					int ci0 = f_pro_ghost[i0].ci;
					m_box.ghost_index.push_back(ci0);
					for (Point_iterator_2d p(0,2); p.end(); ++p)
					{
						int fi0 = f_pro_ghost[i0].fi.son[p.i][p.j][p.k];
						finelevel->m_box[fi0].neib[1][1][1] = -1;
					}
					f_pro_ghost.givehole(i0);
					new_f_res_ghost[ig].push_back(f_pro_ghost[i0]);
				}
			}
			for (int i = 0; i < rf_pro_keep_f_remove_c.size(); ++i)
			{
				int i0 = rf_pro_keep_f_remove_c[i];
				int ci0 = f_pro_ghost[i0].ci;
				m_box.givehole(ci0);
				for (Point_iterator_2d p(0,2); p.end(); ++p)
				{
					int fi0 = f_pro_ghost[i0].fi.son[p.i][p.j][p.k];
					finelevel->m_box[fi0].neib[1][1][1] = -1;
				}
				f_pro_ghost.givehole(i0);
			}
		}
		// PRINTFinLEVEL("Refine res 0 array hole number: %d", cur_level, (int)f_res_ghost[0].Holenum());
		// PRINTFinLEVEL("Refine res 1 array hole number: %d", cur_level, (int)f_res_ghost[1].Holenum());
		// PRINTFinLEVEL("Refine pro 0 array hole number: %d", cur_level, (int)f_pro_ghost.Holenum());
		AdjustProPair();
		// RemoveGhostPair();
	}
}

void AMRLevel::RemoveGhostandNormal_Derefine()
{
	if (!finelevel->twodflag)
	{
		for (int ig = 0; ig < ighost; ++ig)
		{
			for (int i = 0; i < drf_res_switch_ghost_layer[ig].size(); ++i)
			{
				int i0 = drf_res_switch_ghost_layer[ig][i];
				int ci0 = f_res_ghost[ig][i0].ci;
				m_box[ci0].neib[1][1][1] = -1;
				for (Point_iterator p(0,2); p.end(); ++p)
				{
					int fi0 = f_res_ghost[ig][i0].fi.son[p.i][p.j][p.k];
					finelevel->m_box.ghost_index.push_back(fi0);
				}
				new_f_pro_ghost.push_back(f_res_ghost[ig][i0]);
				f_res_ghost[ig].givehole(i0);
			}
		}
		for (int i = 0; i < drf_remove_pro_ghost_fcell.size(); ++i)
		{
			int i0 = drf_remove_pro_ghost_fcell[i];
			for (Point_iterator p(0,2); p.end(); ++p)
			{
				int fi0 = f_pro_ghost[i0].fi.son[p.i][p.j][p.k];
				finelevel->m_box.givehole(fi0);
				finelevel->m_box[fi0].neib[1][1][1] = -1;
			}
			f_pro_ghost.givehole(i0);	
		}
		for (int i = 0; i < drf_remove_level_pro_ghost.size(); ++i)
		{
			int i0 = drf_remove_level_pro_ghost[i];
			level_pro_ghost.givehole(i0);
			for (Point_iterator p(0,2); p.end(); ++p)
			{
				int fi0 = level_pro_ghost[i0].fi.son[p.i][p.j][p.k];
				finelevel->m_box.givehole(fi0);
				finelevel->m_box[fi0].neib[1][1][1] = -1;
			}
		}
	}
	else
	{
		for (int ig = 0; ig < ighost; ++ig)
		{
			for (int i = 0; i < drf_res_switch_ghost_layer[ig].size(); ++i)
			{
				int i0 = drf_res_switch_ghost_layer[ig][i];
				int ci0 = f_res_ghost[ig][i0].ci;
				m_box[ci0].neib[1][1][1] = -1;
				for (Point_iterator_2d p(0,2); p.end(); ++p)
				{
					int fi0 = f_res_ghost[ig][i0].fi.son[p.i][p.j][p.k];
					finelevel->m_box.ghost_index.push_back(fi0);
				}
				new_f_pro_ghost.push_back(f_res_ghost[ig][i0]);
				f_res_ghost[ig].givehole(i0);
			}
		}
		for (int i = 0; i < drf_remove_pro_ghost_fcell.size(); ++i)
		{
			int i0 = drf_remove_pro_ghost_fcell[i];
			for (Point_iterator_2d p(0,2); p.end(); ++p)
			{
				int fi0 = f_pro_ghost[i0].fi.son[p.i][p.j][p.k];
				finelevel->m_box.givehole(fi0);
				finelevel->m_box[fi0].neib[1][1][1] = -1;
			}
			f_pro_ghost.givehole(i0);	
		}
		for (int i = 0; i < drf_remove_level_pro_ghost.size(); ++i)
		{
			int i0 = drf_remove_level_pro_ghost[i];
			level_pro_ghost.givehole(i0);
			for (Point_iterator_2d p(0,2); p.end(); ++p)
			{
				int fi0 = level_pro_ghost[i0].fi.son[p.i][p.j][p.k];
				finelevel->m_box.givehole(fi0);
				finelevel->m_box[fi0].neib[1][1][1] = -1;
			}
		}
	}
	/*-----------------------------------------------------------------------*/
	for (int ig = 0; ig < ighost; ++ig)
	{
		for (int i = 0; i < drf_res_keep_c_remove_f[ig].size(); ++i)
		{
			int i0 = drf_res_keep_c_remove_f[ig][i];
			int ci0 = f_res_ghost[ig][i0].ci;
			m_box[ci0].neib[1][1][1] = -1;
			f_res_ghost[ig].givehole(i0);
		}
	}
	for (int i = 0; i < drf_res_move2layer.size(); ++i)
	{
		int i0 = drf_res_move2layer[i];
		new_f_res_ghost[0].push_back(f_res_ghost[1][i0]);
		f_res_ghost[1].givehole(i0);
	}
	
	int ds379 = drf_remove_blockghost.size();
	for (int i = 0; i < ds379; ++i)
	{
		int i0 = drf_remove_blockghost[i];
		finelevel->m_box[i0].neib[1][1][1] = -1;
	}
	MPI_Barrier(share_comm);
	// PRINTFinLEVEL("Coarsen res 0 array hole number: %d", cur_level, (int)f_res_ghost[0].Holenum());
	// PRINTFinLEVEL("Coarsen res 1 array hole number: %d", cur_level, (int)f_res_ghost[1].Holenum());
	// PRINTFinLEVEL("Coarsen pro 0 array hole number: %d", cur_level, (int)f_pro_ghost.Holenum());
	// for (int i = 0; i < new_f_res_ghost[1].size(); ++i)
	// {
	// 	int ci0 = new_f_res_ghost[1][i].ci;
	// 	PRINTFinLEVEL("new_f_res_ghost 1 is box (%d,%d,%d)",cur_level,
	// 		m_box[ci0].ix(),
	// 		m_box[ci0].iy(),
	// 		m_box[ci0].iz())
	// }
	AdjustProPair();
	// RemoveGhostPair();
}

void AMRLevel::RemoveArrayGhost(DataArray<Box> & mybox)
{
	//MPI_Win_lock_all(0, mybox.arraywin());
	int gn = mybox.ghost_index.size();
	for (int i = 0; i < gn; ++i)
	{
		int fc0 = mybox.ghost_index[i];
		// printf("Level %d ghost %d point is (%d,%d,%d) index is %d\n", cur_level, i, 
		// 	mybox[fc0].ix(), mybox[fc0].iy(), mybox[fc0].iz(), mybox[fc0].neib[1][1][1].index);
		if (mybox[fc0].neib[1][1][1] == -1)
		{
			mybox.ghost_index[i] = mybox.ghost_index.back();
			mybox.ghost_index.pop_back();
			--i;
			--gn;
			mybox[fc0].neib[1][1][1] = fc0;
			//printf("Level %d ghost %d point is (%d,%d,%d)\n", cur_level, i, mybox[fc0].ix(), mybox[fc0].iy(), mybox[fc0].iz());
		}
	}
	//MPI_Win_unlock_all(mybox.arraywin());
	MPI_Barrier(share_comm);
	int bs, be;
	ArrayOrder_s(0, mybox.realsize(), bs, be, sprocs, srank);
#ifdef DEBUG	
	for (int i = bs; i < be; ++i)
	{
		if (mybox[i].neib[1][1][1] == -1)
		{
			PRINTFinLEVEL("Box %d (%d,%d,%d) main index is -1 after clear the ghost index!!!", 
				cur_level, i, mybox[i].ix(), mybox[i].iy(), mybox[i].iz());
			MPI_Abort(MPI_COMM_WORLD, 524);
		}
	}
#endif	
}

void AMRLevel::IncludeNewboxinGhost_Refine()
{
	vector<Boxson<int> > ghostson;
	int origin_rfboxnum = rfbox.size();
	// printf("Level %d Pro switch ghost pair is %d %d rf_pro_keep_f_remove_c is %d\n", 
	// 	cur_level, (int)rf_pro_switch_ghost_layer[0].size(), (int)rf_pro_switch_ghost_layer[1].size(),
	// 	(int)rf_pro_keep_f_remove_c.size());
	PRINTFinLEVEL("rf_pro_switch_ghost_layer[0] size %d", cur_level, (int)rf_pro_switch_ghost_layer[0].size());
	PRINTFinLEVEL("rf_pro_switch_ghost_layer[1] size %d", cur_level, (int)rf_pro_switch_ghost_layer[1].size());
	PRINTFinLEVEL("rf_pro_keep_f_remove_c size %d", cur_level, (int)rf_pro_keep_f_remove_c.size());
	PRINTFinLEVEL("rf_new_level_pro_ghost size %d", cur_level, (int)rf_new_level_pro_ghost.size());
	PRINTFinLEVEL("rf_res_remove_ghost_ccell[0] size %d", cur_level, (int)rf_res_remove_ghost_ccell[0].size());
	PRINTFinLEVEL("rf_res_remove_ghost_ccell[1] size %d", cur_level, (int)rf_res_remove_ghost_ccell[1].size());
	PRINTFinLEVEL("rf_res_move2layer size %d",cur_level,(int)rf_res_move2layer.size());
	PRINTFinLEVEL("rf_recover_pair size %d", cur_level, (int)rf_recover_pair.size());
	PRINTFinLEVEL("rf_recover_block size %d", cur_level, (int)rf_recover_block.size());
	PRINTFinLEVEL("rf_remove_levelproghost size %d", cur_level, (int)rf_remove_levelproghost.size());
	if (!finelevel->twodflag)
	{
		/*----------cells used to be a pro pair and marked as a res pair---------------*/
		for (int ig = 0; ig < ighost; ++ig)
		{
			for (int i = 0; i < rf_pro_switch_ghost_layer[ig].size(); ++i)
			{
				int i0 = rf_pro_switch_ghost_layer[ig][i];
				rfbox.push_back(f_pro_ghost[i0].ci);
				Assert(m_tag[f_pro_ghost[i0].ci].tag == -1, "non zero tag in including new ghost box!!!", 370);
				ghostson.push_back(f_pro_ghost[i0].fi);
				for (Point_iterator p(0,2); p.end(); ++p)
				{
					newrefinebox.push_back(f_pro_ghost[i0].fi.son[p.i][p.j][p.k]);
				}
			}
		}
		/*---------cells used to be pro pairs and marked as refined cells--------------*/
		for (int i = 0; i < rf_pro_keep_f_remove_c.size(); ++i)
		{
			int i0 = rf_pro_keep_f_remove_c[i];
			rfbox.push_back(f_pro_ghost[i0].ci);
			Assert(m_tag[f_pro_ghost[i0].ci].tag == -1, "non zero tag in including new ghost box!!!", 386);
			ghostson.push_back(f_pro_ghost[i0].fi);
			for (Point_iterator p(0,2); p.end(); ++p)
			{
				newrefinebox.push_back(f_pro_ghost[i0].fi.son[p.i][p.j][p.k]);
			}
		}
		/*---------cells used to be pro pairs and still marked as pro pairs----------*/
		for (int i = 0; i < rf_recover_pair.size(); ++i)
		{
			int i0 = rf_recover_pair[i];
			rfbox.push_back(f_pro_ghost[i0].ci);
			Assert(m_tag[f_pro_ghost[i0].ci].tag == -1, "non zero tag in including new ghost box!!!", 414);
			ghostson.push_back(f_pro_ghost[i0].fi);
			for (Point_iterator p(0,2); p.end(); ++p)
			{
				newrefinebox.push_back(f_pro_ghost[i0].fi.son[p.i][p.j][p.k]);
			}
		}
		/*----------cells in the block conjunction-----------*/
		for (int i = 0; i < rf_recover_block.size(); ++i)
		{
			int i0 = rf_recover_block[i];
			rfbox.push_back(level_pro_ghost[i0].ci);
			Assert(m_tag[level_pro_ghost[i0].ci].tag == -1, "non zero tag in including new ghost box!!!", 427);
			ghostson.push_back(level_pro_ghost[i0].fi);
			for (Point_iterator p(0,2); p.end(); ++p)
			{
				newrefinebox.push_back(level_pro_ghost[i0].fi.son[p.i][p.j][p.k]);
			}
		}
		int se540 = rf_oldres_to_pro.size();
		for (int i = 0; i < se540; ++i)
		{
			int i0 = rf_oldres_to_pro[i];
			rfbox.push_back(drf_remove_respair[i0].ci);
			if (m_tag[drf_remove_respair[i0].ci].tag != -1)
			{
				PRINTFinLEVEL("The removed res pair main box tag is %d!!!", cur_level, m_tag[drf_remove_respair[i0].ci].tag);
			}
			Assert(m_tag[drf_remove_respair[i0].ci].tag == -1, "non zero tag in including new ghost box!!!", 545);
			ghostson.push_back(drf_remove_respair[i0].fi);
			for (Point_iterator p(0,2); p.end(); ++p)
			{
				newrefinebox.push_back(drf_remove_respair[i0].fi.son[p.i][p.j][p.k]);
			}
		}
	}
	else
	{
		/*----------cells used to be a pro pair and marked as a res pair---------------*/
		for (int ig = 0; ig < ighost; ++ig)
		{
			for (int i = 0; i < rf_pro_switch_ghost_layer[ig].size(); ++i)
			{
				int i0 = rf_pro_switch_ghost_layer[ig][i];
				rfbox.push_back(f_pro_ghost[i0].ci);
				Assert(m_tag[f_pro_ghost[i0].ci].tag == -1, "non zero tag in including new ghost box!!!", 370);
				ghostson.push_back(f_pro_ghost[i0].fi);
				for (Point_iterator_2d p(0,2); p.end(); ++p)
				{
					newrefinebox.push_back(f_pro_ghost[i0].fi.son[p.i][p.j][p.k]);
				}
			}
		}
		/*---------cells used to be pro pairs and marked as refined cells--------------*/
		for (int i = 0; i < rf_pro_keep_f_remove_c.size(); ++i)
		{
			int i0 = rf_pro_keep_f_remove_c[i];
			rfbox.push_back(f_pro_ghost[i0].ci);
			Assert(m_tag[f_pro_ghost[i0].ci].tag == -1, "non zero tag in including new ghost box!!!", 386);
			ghostson.push_back(f_pro_ghost[i0].fi);
			for (Point_iterator_2d p(0,2); p.end(); ++p)
			{
				newrefinebox.push_back(f_pro_ghost[i0].fi.son[p.i][p.j][p.k]);
			}
		}
		/*---------cells used to be pro pairs and still marked as pro pairs----------*/
		for (int i = 0; i < rf_recover_pair.size(); ++i)
		{
			int i0 = rf_recover_pair[i];
			rfbox.push_back(f_pro_ghost[i0].ci);
			Assert(m_tag[f_pro_ghost[i0].ci].tag == -1, "non zero tag in including new ghost box!!!", 414);
			ghostson.push_back(f_pro_ghost[i0].fi);
			for (Point_iterator_2d p(0,2); p.end(); ++p)
			{
				newrefinebox.push_back(f_pro_ghost[i0].fi.son[p.i][p.j][p.k]);
			}
		}
		/*----------cells in the block conjunction-----------*/
		for (int i = 0; i < rf_recover_block.size(); ++i)
		{
			int i0 = rf_recover_block[i];
			rfbox.push_back(level_pro_ghost[i0].ci);
			Assert(m_tag[level_pro_ghost[i0].ci].tag == -1, "non zero tag in including new ghost box!!!", 427);
			ghostson.push_back(level_pro_ghost[i0].fi);
			for (Point_iterator_2d p(0,2); p.end(); ++p)
			{
				newrefinebox.push_back(level_pro_ghost[i0].fi.son[p.i][p.j][p.k]);
			}
		}
		int se540 = rf_oldres_to_pro.size();
		for (int i = 0; i < se540; ++i)
		{
			int i0 = rf_oldres_to_pro[i];
			rfbox.push_back(drf_remove_respair[i0].ci);
			if (m_tag[drf_remove_respair[i0].ci].tag != -1)
			{
				PRINTFinLEVEL("The removed res pair main box tag is %d!!!", cur_level, m_tag[drf_remove_respair[i0].ci].tag);
			}
			Assert(m_tag[drf_remove_respair[i0].ci].tag == -1, "non zero tag in including new ghost box!!!", 545);
			ghostson.push_back(drf_remove_respair[i0].fi);
			for (Point_iterator_2d p(0,2); p.end(); ++p)
			{
				newrefinebox.push_back(drf_remove_respair[i0].fi.son[p.i][p.j][p.k]);
			}
		}
	}
	int existsonnum = allson.size();
	int ghostboxnum = rfbox.size() - rfnum_procs[srank];
	//PRINTFinLEVEL("^^^^^^^^^^^^Include more new fine box number %d",cur_level,ghostboxnum);
	vector<int> rank_ghostbox_num(sprocs);
	vector<int> rank_ghostbox_start(sprocs,0);
	MPI_Allgather(&ghostboxnum, 1, MPI_INT, &rank_ghostbox_num[0], 1, MPI_INT, share_comm);
	ArrayProcsStart(rank_ghostbox_num, rank_ghostbox_start);
	for (int i = 0; i < ghostboxnum; ++i)
	{
		Assert(m_tag[rfbox[i+origin_rfboxnum]].tag == -1, "non zero tag in including new ghost box!!!", 444);
		m_tag[rfbox[i+origin_rfboxnum]].tag = existsonnum + rank_ghostbox_start[srank] + i;
	}
	allson.Addnew(ghostson);
	allson.DirectlyReduceNew();
	MPI_Barrier(share_comm);
	ConstructBlockProPair_Refine();
}

void AMRLevel::IncludeNewboxinGhost_Derefine()
{
	Boxson<int> * ghostson;
	PRINTFinLEVEL("drf_res_switch_ghost_layer [0]: %d; [1]: %d", cur_level, 
		(int)drf_res_switch_ghost_layer[0].size(),(int)drf_res_switch_ghost_layer[1].size());
	PRINTFinLEVEL("drf_res_keep_c_remove_f [0]: %d; [1]: %d", cur_level, 
		(int)drf_res_keep_c_remove_f[0].size(),(int)drf_res_keep_c_remove_f[1].size());
	PRINTFinLEVEL("drf_recover_pair [0]: %d; [1]: %d", cur_level, 
		(int)drf_recover_pair[0].size(),(int)drf_recover_pair[1].size());
	PRINTFinLEVEL("drf_recover_level_pro_ghost: %d",cur_level,(int)drf_recover_level_pro_ghost.size());
	PRINTFinLEVEL("drf_remove_level_pro_ghost: %d",cur_level,(int)drf_remove_level_pro_ghost.size());
	PRINTFinLEVEL("drf_res_move2layer: %d",cur_level,(int)drf_res_move2layer.size());
	/*---------------Cells used to be in res pair and fine cells will be new ghosts------------*/
	if (!finelevel->twodflag)
	{
		for (int ig = 0; ig < ighost; ++ig)
		{
			for (int i = 0; i < drf_res_switch_ghost_layer[ig].size(); ++i)
			{
				int i0 = drf_res_switch_ghost_layer[ig][i];
				finelevel->derfbox.push_back(f_res_ghost[ig][i0].fi.son[0][0][0]);
				ghostson = &f_res_ghost[ig][i0].fi;
				finelevel->newcoarsebox.push_back(f_res_ghost[ig][i0].ci);
				for (Point_iterator p(0,2); p.end(); ++p)
				{
					finelevel->m_tag[ghostson->son[p.i][p.j][p.k]].detag = f_res_ghost[ig][i0].ci;
				}
			}
			/*---------------Cells used to be in res pair and fine cells will be deleted------------*/
			for (int i = 0; i < drf_res_keep_c_remove_f[ig].size(); ++i)
			{
				int i0 = drf_res_keep_c_remove_f[ig][i];
				finelevel->derfbox.push_back(f_res_ghost[ig][i0].fi.son[0][0][0]);
				ghostson = &f_res_ghost[ig][i0].fi;
				finelevel->newcoarsebox.push_back(f_res_ghost[ig][i0].ci);
				for (Point_iterator p(0,2); p.end(); ++p)
				{
					finelevel->m_tag[ghostson->son[p.i][p.j][p.k]].detag = f_res_ghost[ig][i0].ci;
				}
			}
		}
		/*-----cells used to be in res pair and still marked as res pair*/
		for (int ig = 0; ig < ighost; ++ig)
		{
			for (int i = 0; i < drf_recover_pair[ig].size(); ++i)
			{
				int i0 = drf_recover_pair[ig][i];
				finelevel->derfbox.push_back(f_res_ghost[ig][i0].fi.son[0][0][0]);
				finelevel->newcoarsebox.push_back(f_res_ghost[ig][i0].ci);
				for (Point_iterator p(0,2); p.end(); ++p)
				{
					finelevel->m_tag[f_res_ghost[ig][i0].fi.son[p.i][p.j][p.k]].detag = f_res_ghost[ig][i0].ci;
				}
			}
		}
		/*-----cells in the block conjunction-------------------------*/
		for (int i = 0; i < drf_recover_level_pro_ghost.size(); ++i)
		{
			int i0 = drf_recover_level_pro_ghost[i];
			finelevel->derfbox.push_back(level_pro_ghost[i0].fi.son[0][0][0]);
			finelevel->newcoarsebox.push_back(level_pro_ghost[i0].ci);
			for (Point_iterator p(0,2); p.end(); ++p)
			{
				finelevel->m_tag[level_pro_ghost[i0].fi.son[p.i][p.j][p.k]].detag = level_pro_ghost[i0].ci;
			}
		}
	}
	else
	{
		for (int ig = 0; ig < ighost; ++ig)
		{
			for (int i = 0; i < drf_res_switch_ghost_layer[ig].size(); ++i)
			{
				int i0 = drf_res_switch_ghost_layer[ig][i];
				finelevel->derfbox.push_back(f_res_ghost[ig][i0].fi.son[0][0][0]);
				ghostson = &f_res_ghost[ig][i0].fi;
				finelevel->newcoarsebox.push_back(f_res_ghost[ig][i0].ci);
				for (Point_iterator_2d p(0,2); p.end(); ++p)
				{
					finelevel->m_tag[ghostson->son[p.i][p.j][p.k]].detag = f_res_ghost[ig][i0].ci;
				}
			}
			/*---------------Cells used to be in res pair and fine cells will be deleted------------*/
			for (int i = 0; i < drf_res_keep_c_remove_f[ig].size(); ++i)
			{
				int i0 = drf_res_keep_c_remove_f[ig][i];
				finelevel->derfbox.push_back(f_res_ghost[ig][i0].fi.son[0][0][0]);
				ghostson = &f_res_ghost[ig][i0].fi;
				finelevel->newcoarsebox.push_back(f_res_ghost[ig][i0].ci);
				for (Point_iterator_2d p(0,2); p.end(); ++p)
				{
					finelevel->m_tag[ghostson->son[p.i][p.j][p.k]].detag = f_res_ghost[ig][i0].ci;
				}
			}
		}
		/*-----cells used to be in res pair and still marked as res pair*/
		for (int ig = 0; ig < ighost; ++ig)
		{
			for (int i = 0; i < drf_recover_pair[ig].size(); ++i)
			{
				int i0 = drf_recover_pair[ig][i];
				finelevel->derfbox.push_back(f_res_ghost[ig][i0].fi.son[0][0][0]);
				finelevel->newcoarsebox.push_back(f_res_ghost[ig][i0].ci);
				for (Point_iterator_2d p(0,2); p.end(); ++p)
				{
					finelevel->m_tag[f_res_ghost[ig][i0].fi.son[p.i][p.j][p.k]].detag = f_res_ghost[ig][i0].ci;
				}
			}
		}
		/*-----cells in the block conjunction-------------------------*/
		for (int i = 0; i < drf_recover_level_pro_ghost.size(); ++i)
		{
			int i0 = drf_recover_level_pro_ghost[i];
			finelevel->derfbox.push_back(level_pro_ghost[i0].fi.son[0][0][0]);
			finelevel->newcoarsebox.push_back(level_pro_ghost[i0].ci);
			for (Point_iterator_2d p(0,2); p.end(); ++p)
			{
				finelevel->m_tag[level_pro_ghost[i0].fi.son[p.i][p.j][p.k]].detag = level_pro_ghost[i0].ci;
			}
		}
	}
	MPI_Barrier(share_comm);
}

void AMRLevel::RenewLevelProGhostIndex()
{
	for (int i = level_pro_ghost.ps(); i < level_pro_ghost.pe(); ++i)
	{
		int i0 = level_pro_ghost[i].ci;
		level_pro_ghost[i].ci = m_box[i0].neib[1][1][1];
	}
	// for (int i = level_pro_fluid.ps(); i < level_pro_fluid.pe(); ++i)
	// {
	// 	int i0 = level_pro_fluid[i].ci;
	// 	level_pro_fluid[i].ci = m_box[i0].neib[1][1][1];
	// }
	if (NULL != coarselevel)
	{
		if (!twodflag)
		{
			for (int i = coarselevel->level_pro_ghost.ps(); i < coarselevel->level_pro_ghost.pe(); ++i)
			{
				for (Point_iterator p(0,2); p.end(); ++p)
				{
					int fi0 = coarselevel->level_pro_ghost[i].fi.son[p.i][p.j][p.k];
					coarselevel->level_pro_ghost[i].fi.son[p.i][p.j][p.k] = m_box[fi0].neib[1][1][1];
				}
			}
		}
		else
		{
			for (int i = coarselevel->level_pro_ghost.ps(); i < coarselevel->level_pro_ghost.pe(); ++i)
			{
				for (Point_iterator_2d p(0,2); p.end(); ++p)
				{
					int fi0 = coarselevel->level_pro_ghost[i].fi.son[p.i][p.j][p.k];
					coarselevel->level_pro_ghost[i].fi.son[p.i][p.j][p.k] = m_box[fi0].neib[1][1][1];
				}
			}
		}
		// for (int i = coarselevel->level_pro_fluid.ps(); i < coarselevel->level_pro_fluid.pe(); ++i)
		// {
		// 	for (Point_iterator p(0,2); p.end(); ++p)
		// 	{
		// 		int fi0 = coarselevel->level_pro_fluid[i].fi.son[p.i][p.j][p.k];
		// 		coarselevel->level_pro_fluid[i].fi.son[p.i][p.j][p.k] = fi0;
		// 	}
		// }
	}
	MPI_Barrier(share_comm);
}

void AMRLevel::RenewBlockPairIndex()
{
	vector<vector<S_RemoteBKPInfo> > senddata(nodenum);
	//PRINTFinLEVEL("ghost_target size %d",cur_level,(int)ghost_target.size());
	for (int i = 0; i < ghost_target.size(); ++i)
	{
		int g0 = ghost_target[i];
		if (g0 > -1)
		{
			int i0 = m_box[g0].neib[1][1][1];
			// if (cur_level == 1 && i0 > 129000)
			// {
			// PRINTFinLEVEL("ghost target old is %d new is %d", cur_level, g0, i0);}
			if (g0 != i0)
			{
				senddata[originbp[i].node].push_back(S_RemoteBKPInfo(originbp[i].index, i0));
#ifdef DEBUG				
				senddata[originbp[i].node].back().bxyz = m_box[g0].lowpt;
#endif				
			}
			ghost_target[i] = i0;
			//PRINTFinLEVEL("No %d send to node %d is its block pair %d out cell is box %d", 
			//	cur_level,(int)senddata[originbp[i].node].size()-1,originbp[i].node,originbp[i].index,m_box[g0].neib[1][1][1]);
		}
	}
	vector<S_RemoteBKPInfo> recvdata;
	BcastNewInfo_Node(senddata, recvdata, 4);
	for (int i = 0; i < recvdata.size(); ++i)
	{
		int bkp0 = recvdata[i].bkpindex;
		blockpair[bkp0].outnode.index = recvdata[i].localcell;
#ifdef DEBUG
		Point s1 = recvdata[i].bxyz;
		Point s2 = m_box[blockpair[bkp0].innode].lowpt;
		int rf0 = pow(2, cur_level);
#if DIM == 3		
		Assert((s1[0]==s2[0] || (abs(s1[0]-s2[0]) == highpt.xy[0]*rf0 && periodic[0])) && 
					 (s1[1]==s2[1] || (abs(s1[1]-s2[1]) == highpt.xy[1]*rf0 && periodic[1])) &&
					 (s1[2]==s2[2] || (abs(s1[2]-s2[2]) == highpt.xy[2]*power_ratio && periodic[2])), 
					 "A ghost target xyz does not match the original box!!!", 104);
#elif DIM == 2
		Assert((s1[0]==s2[0] || (abs(s1[0]-s2[0]) == highpt.xy[0]*rf0 && periodic[0])) && 
					 (s1[1]==s2[1] || (abs(s1[1]-s2[1]) == highpt.xy[1]*rf0 && periodic[1])), 
					 "A ghost target xyz does not match the original box!!!", 104);	
#endif					 	
#endif
		//Assert(recvdata[i].bxyz==m_box[blockpair[bkp0].innode].lowpt, "Error in renewing ghost index!!!",878);	
		//PRINTFinLEVEL("Block pair %d out node is %d", cur_level, bkp0, recvdata[i].localcell);
	}
	MPI_Barrier(share_comm);
	// for (int i = 0; i < blockpair.hole.size(); ++i)
	// {
	// 	int i0 = blockpair.hole[i];
	// 	int inn0 = blockpair[i0].innode;
	// 	if (inn0 > -1)
	// 	{
	// 		if (m_tag[inn0].tag == 0)
	// 		{
	// 		printf("N%d Level %d Block pair hole %d (%d,%d,%d) new neib 1 %d tag is %d total box %d gps %d gpe %d\n", 
	// 			node,cur_level,inn0,m_box[inn0].ix(),m_box[inn0].iy(),m_box[inn0].iz(),m_box[inn0].neib[1][1][1],m_tag[inn0].tag,
	// 			m_box.realsize(), m_box.gps(), m_box.gpe());}
	// 	}
	// }	
	for (int i = blockpair.ps(); i < blockpair.pe(); ++i)
	{
		int in0 = blockpair[i].innode;
		// if (blockpair[i].outnode.index == 27744)
		// {
		// 	PRINTFinLEVEL("#######the box 27744 is the pair for (%d,%d,%d)",cur_level,
		// 		m_box[blockpair[i].innode].ix(),
		// 		m_box[blockpair[i].innode].iy(),
		// 		m_box[blockpair[i].innode].iz());
		// 	// MPI_Abort(MPI_COMM_WORLD,871);
		// }
		if (in0 > -1)
		{
			blockpair[i].innode = m_box[in0].neib[1][1][1];
		}
#ifdef DEBUG		
		if (blockpair[i].innode == -1)
		{
			if (blockpair[i].outnode.index != -1)
			{
				PRINTFinLEVEL("The in cell %d (%d,%d,%d) detag is %d type is %d will be removed but the out node is not !!!",cur_level,
					in0, m_box[in0].ix(), m_box[in0].iy(), m_box[in0].iz(), m_tag[in0].detag, m_box[in0].type);
				MPI_Abort(MPI_COMM_WORLD,893);
			}
		}
#endif		
	}
	MPI_Barrier(share_comm);
// #ifdef	DEBUG
// 	for (int i = 0; i < blockpair.hole.size(); ++i)
// 	{
// 		int i0 = blockpair.hole[i];
// 		if (blockpair[i0].innode != -1)
// 		{
// 			printf("N%d A hole bkp innode is %d is box (%d,%d,%d) tag is %d\n",
// 				node, 
// 				blockpair[i0].innode, 
// 				m_box[blockpair[i0].innode].ix(), 
// 				m_box[blockpair[i0].innode].iy(),
// 				m_box[blockpair[i0].innode].iz(),
// 				m_tag[blockpair[i0].innode].tag);
// 			printf("Error in check the block pair innode index!!!");
// 			MPI_Abort(MPI_COMM_WORLD,892);
// 		}
// 	}
// 	MPI_Barrier(share_comm);
// #endif
}
struct S_state
{
	int remotebkp;
	int state;

	S_state()
	{}

	S_state(const int & rc0, const int & s0)
	{
		remotebkp = rc0;
		state = s0;
	}
};
void AMRLevel::ReorderBlockpair()
{
	blockpair.holeplan(true);
	blockpair.CompressArray();
}

void AMRLevel::RenewGhostIndex()
{
	RenewProandResIndex();
	GiveLevelFlag("Finish renew ghost target", cur_level, 5);
	RenewBlockPairIndex();
	RenewLevelProGhostIndex();
	RenewDomainGhost();
}
/*2020-7-6 end*/
void AMRLevel::ConstructNewDomainGhost(vector<Point> & facenmv)
{
	if (NULL != finelevel)
	{
		//printf("Node %d\n", node);
		vector<Domainghost> finedmg(1000*cur_level);
		finedmg.resize(0);
		if (!finelevel->twodflag)
		{
			for (int i = dmghost.ps(); i < dmghost.pe(); ++i)
			{
				int & d0 = dmghost[i].cell;
				Assert(d0 > -1 && d0 < m_box.realsize(), "The domain ghost not in range!!!", 963);
				if (m_box.isnorg(d0))
				{
					if (m_tag[d0].tag > -1 && m_tag[d0].tag < totrfnum)
					{
						int tag0 = m_tag[d0].tag;
						Assert(dmghost[i].refcell > -1 && dmghost[i].refcell < m_box.realsize(),
							"The refined domain ghost refcell not in range", 1184);
						Assert(tag0 < allson.size(), "Error of tag0 in ConstructNewDomainGhost!!!", 979);
						Assert((dmghost[i].nbface > -1 && dmghost[i].nbface < facenmv.size()), "Face error", 984);
						Point & dir0 = facenmv[dmghost[i].nbface];
						for (Point_iterator p(0,2); p.end(); ++p)
						{
							int s0 = allson[tag0].son[p.i][p.j][p.k];
							Assert((s0 > -1 && s0 < finelevel->m_box.realsize()), 
								"son index not in range in ConstructNewDomainGhost!!!", 984);
							finelevel->m_box[s0].type = Dmghost;
							int r0 = finelevel->m_box[s0].neib[1+dir0[0]][1+dir0[1]][1+dir0[2]];
							Assert(r0 > -1, "New DomainGhost has a negative ref cell!!!", 904);
							finedmg.push_back(Domainghost(s0, dmghost[i].nbface, r0));
						}
					}
				}
			}
		}
		else
		{
			for (int i = dmghost.ps(); i < dmghost.pe(); ++i)
			{
				int & d0 = dmghost[i].cell;
				Assert(d0 > -1 && d0 < m_box.realsize(), "The domain ghost not in range!!!", 963);
				if (m_box.isnorg(d0))
				{
					if (m_tag[d0].tag > -1 && m_tag[d0].tag < totrfnum)
					{
						int tag0 = m_tag[d0].tag;
						Assert(dmghost[i].refcell > -1 && dmghost[i].refcell < m_box.realsize(),
							"The refined domain ghost refcell not in range", 1184);
						Assert(tag0 < allson.size(), "Error of tag0 in ConstructNewDomainGhost!!!", 979);
						Assert((dmghost[i].nbface > -1 && dmghost[i].nbface < facenmv.size()), "Face error", 984);
						Point & dir0 = facenmv[dmghost[i].nbface];
						for (Point_iterator_2d p(0,2); p.end(); ++p)
						{
							int s0 = allson[tag0].son[p.i][p.j][p.k];
							Assert((s0 > -1 && s0 < finelevel->m_box.realsize()), 
								"son index not in range in ConstructNewDomainGhost!!!", 984);
							finelevel->m_box[s0].type = Dmghost;
							int r0 = finelevel->m_box[s0].neib[1+dir0[0]][1+dir0[1]][1+dir0[2]];
							Assert(r0 > -1, "New DomainGhost has a negative ref cell!!!", 904);
							finedmg.push_back(Domainghost(s0, dmghost[i].nbface, r0));
						}
					}
				}
			}
		}
		MPI_Barrier(share_comm);
		int be117 = finedmg.size();
		for (int i = 0; i < be117; ++i)
		{
			if (finelevel->m_box[finedmg[i].refcell].type == Dmghost)
			{
				Point & dir0 = facenmv[finedmg[i].nbface];
				finedmg[i].refcell = finelevel->m_box[finedmg[i].refcell].neib[1+dir0[0]][1+dir0[1]][1+dir0[2]];
				Assert(finedmg[i].refcell > -1, "The new ref cell for a dmghost must be non-negative!!!", 1177);
			}
		}
		MPI_Barrier(share_comm);
		GiveLevelFlag("Start renew negative refcell!!!", cur_level, 5);
		for (int i = dmghost.ps(); i < dmghost.pe(); ++i)
		{
			if (dmghost[i].refcell < 0)
			{
				Point & dir0 = facenmv[dmghost[i].nbface];
				dmghost[i].refcell = m_box[dmghost[i].cell].neib[1+dir0[0]][1+dir0[1]][1+dir0[2]];
				if (dmghost[i].refcell > -1)
				{
					if (m_box[dmghost[i].refcell].type == Dmghost)
					{
						dmghost[i].refcell = m_box[dmghost[i].refcell].neib[1+dir0[0]][1+dir0[1]][1+dir0[2]];
					}
				}
			}
		}
		GiveLevelFlag("Start add new dmghost array!!!", cur_level, 5);
		finelevel->dmghost.Addnew(finedmg);
		//finelevel->dmghost.DirectlyReduceNew();
		finelevel->dmghost.holeplan(true);
		finelevel->dmghost.CompressArray();	
		PRINTFinLEVELRANK0("!!!Domain ghost number: %d new number: %d", 
			cur_level+1, finelevel->dmghost.size(), (int)finedmg.size());
	}
}

void AMRLevel::ClearExpiredDomainGhost()
{
	for (int i = dmghost.ps(); i < dmghost.pe(); ++i)
	{
		int i0 = dmghost[i].cell;
		Assert(i0 > -1 && i0 < m_tag.size(), "The main cell index of the domain ghost out of range!!!", 5);
		if (m_tag[i0].detag == 0)
		{
			dmghost.givehole(i);
		}
		Assert(m_tag[i0].detag != 1, "The domain ghost detag can not be 1!!!", 1225);
	}
}

void AMRLevel::RenewDomainGhost()
{
	for (int i = dmghost.ps(); i < dmghost.pe(); ++i)
	{
		int c0 = dmghost[i].cell;
		int rc0 = dmghost[i].refcell;
		//Assert(c0 > -1 && c0 < m_box.realsize(), "The main cell of the domain ghost not in range!!!", 1255);
		if (c0 > -1)
		{
			dmghost[i].cell = m_box[c0].neib[1][1][1];
		}
		//Assert(dmghost[i].cell > -1 && dmghost[i].cell < m_box.realsize(), 
		//	"The new location of the domain ghost main cell must be non-negative!!!", 1257);
		// if (dmghost[i].cell < 0)
		// {
		// 	printf("The domain ghost is (%d,%d,%d)\n", m_box[c0].ix(), m_box[c0].iy(),m_box[c0].iz());
		// }
		if (rc0 > -1)
		{
			dmghost[i].refcell = m_box[rc0].neib[1][1][1];
		}
	}
#ifdef DEBUG
	int hole_sum = 0;
	for (int i = dmghost.ps(); i < dmghost.pe(); ++i)
	{
		if (dmghost[i].cell == -1)
		{
			++hole_sum;
		}
	}
	MPI_Allreduce(MPI_IN_PLACE, &hole_sum, 1, MPI_INT, MPI_SUM, share_comm);
	int dahn = dmghost.Holenum();
	MPI_Allreduce(MPI_IN_PLACE, &dahn, 1, MPI_INT, MPI_SUM, share_comm);
	if (dahn != hole_sum && srank == 0)
	{
		PRINTFinLEVEL("total -1 main cell number %d total hole number %d", cur_level, hole_sum, dahn);
		MPI_Abort(MPI_COMM_WORLD, 1246);
	}
#endif	
}

void AMRLevel::RenewGhostRefDmghost(vector<Point> & facenmv)
{
	for (int i = dmghost.ps(); i < dmghost.pe(); ++i)
	{
		int r0 = dmghost[i].cell;
		Assert(r0 > -1 && r0 < m_box.realsize(), "The main cell of a domain ghost not in range when renew the ref cell!!!", 1274);
		if (m_box.isghost(r0))
		{
			Point & dir0 = facenmv[dmghost[i].nbface];
			r0 = m_box[r0].neib[1+dir0[0]][1+dir0[1]][1+dir0[2]];
			//Assert(m_box.isnormal(dmghost[i].refcell), "The new ghost ref cell is not normal yet!!!", 1156);
			if (m_box.isnormal(r0))
			{
				dmghost[i].refcell = r0;
			}
		}
	}
}

void AMRLevel::AdjustProPair()
{
	level_pro_ghost.holeplan(true);
	// level_pro_fluid.holeplan(true);
	level_pro_ghost.CompressArray();
	// level_pro_fluid.CompressArray();
	f_pro_ghost.Addnew(new_f_pro_ghost);
	f_pro_ghost.DirectlyReduceNew();
	f_pro_ghost.holeplan(true);
	f_pro_ghost.CompressArray();
	new_f_pro_ghost.resize(0);
	f_res_ghost[0].Addnew(new_f_res_ghost[0]);
	f_res_ghost[0].DirectlyReduceNew();
	f_res_ghost[0].holeplan(true);
	f_res_ghost[0].CompressArray();
	new_f_res_ghost[0].resize(0);
	f_res_ghost[1].Addnew(new_f_res_ghost[1]);
	f_res_ghost[1].DirectlyReduceNew();
	f_res_ghost[1].holeplan(true);
	f_res_ghost[1].CompressArray();
	new_f_res_ghost[1].resize(0);		
}

void AMRLevel::CheckResProPair_Debug()
{
	Point dpt[6] = {Point(0,1,1),Point(2,1,1),Point(1,0,1),Point(1,2,1),Point(1,1,0),Point(1,1,2)};
	for (int i = f_res_ghost[0].ps(); i < f_res_ghost[0].pe(); ++i)
	{
		int ci0 = f_res_ghost[0][i].ci;
		for (int d0 = 0 ; d0 < 6; ++d0)
		{
			int an1 = m_box[ci0].neib[dpt[d0][0]][dpt[d0][1]][dpt[d0][2]];
			if (an1 < 0)
			{
				int an2 = m_box[ci0].neib[2-dpt[d0][0]][2-dpt[d0][1]][2-dpt[d0][2]];
				if (m_box.isnormal(an2))
				{
					PRINTFinLEVEL("CheckResProPair_Debug failed!!! Res 0 main box (%d,%d,%d) has negative neib %d AT (%d,%d,%d)",
						cur_level,m_box[ci0].ix(),m_box[ci0].iy(),m_box[ci0].iz(),an1,dpt[d0][0],dpt[d0][1],dpt[d0][2]);
					MPI_Abort(MPI_COMM_WORLD,990);
				}
			}
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if (nrank == 0)
	{
		printf("###All processors passed res and pro pair check!!!\n");
	}
}

void AMRLevel::MarkallPair()
{
	//Inittag();
	InitPair();
	///-----------------------------------
	int bs = f_pro_ghost.ps();
	int be = f_pro_ghost.pe();
	for (int i = bs; i < be; ++i)
	{
		int ci0 = f_pro_ghost[i].ci;
	// 	m_tag[ci0].tag = 0;
		m_box[ci0].ptype = f_pro;
	}
	bs = f_res_ghost[0].ps();
	be = f_res_ghost[0].pe();
	for (int i = bs; i < be; ++i)
	{
		int ci0 = f_res_ghost[0][i].ci;
		//m_tag[ci0].tag = 1;
		m_box[ci0].ptype = f_res1;
	}
	bs = f_res_ghost[1].ps();
	be = f_res_ghost[1].pe();
	for (int i = bs; i < be; ++i)
	{
		int ci0 = f_res_ghost[1][i].ci;
		//m_tag[ci0].tag = 1;
		m_box[ci0].ptype = f_res2;
	}
	//TellBlockPairNewTag();
	bs = level_pro_ghost.ps();
	be = level_pro_ghost.pe();
	for (int i = bs; i < be; ++i)
	{
		int ci0 = level_pro_ghost[i].ci;
		if (m_tag[ci0].tag != -1)
		{
			if (m_tag[ci0].tag > 0) m_box[ci0].ptype = ghost_pro;
			else if (m_tag[ci0].tag == -2) m_box[ci0].ptype = ghost_res1;
			else if (m_tag[ci0].tag == -3) m_box[ci0].ptype = ghost_res2;
			else
			{
				printf("The box pair type error!!! tag is %d\n", m_tag[ci0].tag);
				MPI_Abort(MPI_COMM_WORLD, 1322);
			}
		}
	}
	///-----------------------------------------------
	// if (NULL != coarselevel)
	// {
	// 	int bs = coarselevel->f_res_ghost[0].ps();
	// 	int be = coarselevel->f_res_ghost[0].pe();
	// 	for (int i = bs; i < be; ++i)
	// 	{
	// 		Boxson<int> & fi0 = coarselevel->f_res_ghost[0][i].fi;
	// 		if (twodflag)
	// 		{
	// 			for (Point_iterator_2d p(0,2); p.end(); ++p)
	// 			{
	// 				int cfi0 = fi0.son[p.i][p.j][p.k];
	// 				Assert((cfi0 > -1 && cfi0 < m_box.size()), "The box to be marked res1 is not in range!!!", 1358);
	// 				m_box[cfi0].ptype = f_res1;
	// 				m_box[cfi0].pindex = i;
	// 			}
	// 		}
	// 		else
	// 		{
	// 			for (Point_iterator p(0,2); p.end(); ++p)
	// 			{
	// 				int cfi0 = fi0.son[p.i][p.j][p.k];
	// 				Assert((cfi0 > -1 && cfi0 < m_box.size()), "The box to be marked res1 is not in range!!!", 1358);
	// 				m_box[cfi0].ptype = f_res1;
	// 				m_box[cfi0].pindex = i;
	// 			}
	// 		}
	// 	}
	// 	bs = coarselevel->f_res_ghost[1].ps();
	// 	be = coarselevel->f_res_ghost[1].pe();
	// 	for (int i = bs; i < be; ++i)
	// 	{
	// 		Boxson<int> & fi0 = coarselevel->f_res_ghost[1][i].fi;
	// 		if (twodflag)
	// 		{
	// 			for (Point_iterator_2d p(0,2); p.end(); ++p)
	// 			{
	// 				int cfi0 = fi0.son[p.i][p.j][p.k];
	// 				Assert((cfi0 > -1 && cfi0 < m_box.size()), "The box to be marked res2 is not in range!!!", 1358);
	// 				m_box[cfi0].ptype = f_res2;
	// 				m_box[cfi0].pindex = i;
	// 			}
	// 		}
	// 		else
	// 		{
	// 			for (Point_iterator p(0,2); p.end(); ++p)
	// 			{
	// 				int cfi0 = fi0.son[p.i][p.j][p.k];
	// 				Assert((cfi0 > -1 && cfi0 < m_box.size()), "The box to be marked res2 is not in range!!!", 1358);
	// 				m_box[cfi0].ptype = f_res2;
	// 				m_box[cfi0].pindex = i;
	// 			}
	// 		}
	// 	}
	// }
	// MPI_Barrier(share_comm);
}

void AMRLevel::MarkResPairBeforeSyn()
{
	for (int ig = 0; ig < ighost; ++ig)
	{
		for (int i = f_res_ghost[ig].ps(); i < f_res_ghost[ig].pe(); ++i)
		{
			int ci0 = f_res_ghost[ig][i].ci;
			Assert(ci0 > -1 && ci0 < m_box.realsize(), "The res main cell not in range when marking!!!",1354);
			if (m_box.isghost(ci0))
			{
				m_tag[ci0].tag = -20;
			}
		}
	}
	MPI_Barrier(share_comm);
}

void AMRLevel::TellBlockPairNewTag()
{
	vector<vector<S_newtag> > senddata(nodenum);
	for (int i = blockpair.ps(); i < blockpair.pe(); ++i)
	{
		int i0 = blockpair[i].innode;
		if (m_box.isnorg(i0))
		{
			if (m_tag[i0].tag != -1)
			{
				senddata[blockpair[i].outnode.node].push_back(S_newtag(blockpair[i].outnode.index,m_tag[i0].tag));
			}
		}
	}
	vector<S_newtag> recvdata;
	BcastNewInfo_Node(senddata, recvdata, 2);
	for (int i = 0; i < recvdata.size(); ++i)
	{
		int i0 = recvdata[i].remotebkp;
		m_tag[i0].tag = recvdata[i].newtag;
	}
	MPI_Barrier(share_comm);
}

void AMRLevel::TellBlockPairNewDetag()
{
	vector<vector<S_newtag> > senddata(nodenum);
	for (int i = blockpair.ps(); i < blockpair.pe(); ++i)
	{
		int i0 = blockpair[i].innode;
		Assert(m_box.isnorg(i0), "The new deta location error for telling other blocks!!!", 1417);
		if (m_tag[i0].detag != -1)
		{
			senddata[blockpair[i].outnode.node].push_back(S_newtag(blockpair[i].outnode.index, m_tag[i0].detag));
		}
	}
	vector<S_newtag> recvdata;
	BcastNewInfo_Node(senddata, recvdata, 2);
	for (int i = 0; i < recvdata.size(); ++i)
	{
		int i0 = recvdata[i].remotebkp;
		m_tag[i0].detag = recvdata[i].newtag;
	}
	MPI_Barrier(share_comm);
}

void AMRLevel::CheckNewProResPair(CtoFPair & ap)
{
	int cxyz[3] = {m_box[ap.ci].ix(),m_box[ap.ci].iy(),m_box[ap.ci].iz()};
	if (!finelevel->twodflag)
	{
		for (Point_iterator p(0,2); p.end(); ++p)
		{
			int fi0 = ap.fi.son[p.i][p.j][p.k];
			int fxyz[3] = {finelevel->m_box[fi0].ix(),finelevel->m_box[fi0].iy(),finelevel->m_box[fi0].iz()};
			if (fxyz[0]!=cxyz[0]*2+p.i)
			{
				PRINTFinLEVEL("Error for a new CtoFPair x!!! Main box (%d,%d,%d) son (%d,%d,%d) is box (%d,%d,%d)",
					cur_level,cxyz[0],cxyz[1],cxyz[2],p.i,p.j,p.k,fxyz[0],fxyz[1],fxyz[2]);
				MPI_Abort(MPI_COMM_WORLD,1161);
			}
			if (fxyz[1]!=cxyz[1]*2+p.j)
			{
				PRINTFinLEVEL("Error for a new CtoFPair y!!! Main box (%d,%d,%d) son (%d,%d,%d) is box (%d,%d,%d)",
					cur_level,cxyz[0],cxyz[1],cxyz[2],p.i,p.j,p.k,fxyz[0],fxyz[1],fxyz[2]);
				MPI_Abort(MPI_COMM_WORLD,1167);
			}
			if (fxyz[2]!=cxyz[2]*2+p.k)
			{
				PRINTFinLEVEL("Error for a new CtoFPair z!!! Main box (%d,%d,%d) son (%d,%d,%d) is box (%d,%d,%d)",
					cur_level,cxyz[0],cxyz[1],cxyz[2],p.i,p.j,p.k,fxyz[0],fxyz[1],fxyz[2]);
				MPI_Abort(MPI_COMM_WORLD,1173);
			}
		}
	}
	else
	{
		for (Point_iterator_2d p(0,2); p.end(); ++p)
		{
			int fi0 = ap.fi.son[p.i][p.j][p.k];
			int fxyz[3] = {finelevel->m_box[fi0].ix(),finelevel->m_box[fi0].iy(),finelevel->m_box[fi0].iz()};
			if (fxyz[0]!=cxyz[0]*2+p.i)
			{
				PRINTFinLEVEL("Error for a new CtoFPair x!!! Main box (%d,%d,%d) son (%d,%d,%d) is box (%d,%d,%d)",
					cur_level,cxyz[0],cxyz[1],cxyz[2],p.i,p.j,p.k,fxyz[0],fxyz[1],fxyz[2]);
				MPI_Abort(MPI_COMM_WORLD,1161);
			}
			if (fxyz[1]!=cxyz[1]*2+p.j)
			{
				PRINTFinLEVEL("Error for a new CtoFPair y!!! Main box (%d,%d,%d) son (%d,%d,%d) is box (%d,%d,%d)",
					cur_level,cxyz[0],cxyz[1],cxyz[2],p.i,p.j,p.k,fxyz[0],fxyz[1],fxyz[2]);
				MPI_Abort(MPI_COMM_WORLD,1167);
			}
			if (fxyz[2]!=cxyz[2])
			{
				PRINTFinLEVEL("Error for a new CtoFPair z!!! Main box (%d,%d,%d) son (%d,%d,%d) is box (%d,%d,%d)",
					cur_level,cxyz[0],cxyz[1],cxyz[2],p.i,p.j,p.k,fxyz[0],fxyz[1],fxyz[2]);
				MPI_Abort(MPI_COMM_WORLD,1173);
			}
		}
	}
}

void AMRLevel::PrintoutData()
{
	if (NULL != coarselevel)
	{
		for (int i = coarselevel->level_pro_ghost.ps(); i < coarselevel->level_pro_ghost.pe() && node == 0; ++i)
		{
			int i0 = coarselevel->level_pro_ghost[i].ci;
			int ix0 = coarselevel->m_box[i0].ix();
			int iy0 = coarselevel->m_box[i0].iy();
			int iz0 = coarselevel->m_box[i0].iz();
			int fi0 = coarselevel->level_pro_ghost[i].fi.son[0][0][0];
			int fix = m_box[fi0].ix();
			int fiy = m_box[fi0].iy();
			int fiz = m_box[fi0].iz();
			PRINTFinLEVEL("level_pro_ghost Pair %d Box %d (%d,%d,%d) fine box 0 is %d,%d,%d",
				cur_level,i,i0,ix0,iy0,iz0,fix,fiy,fiz);
		}
	}
}

void AMRLevel::MarkAllProResPair()
{
	Inittag();
	//MPI_Barrier(share_comm);
	int bs = f_pro_ghost.ps();
	int be = f_pro_ghost.pe();
	for (int i = bs; i < be; ++i)
	{
		int ci0 = f_pro_ghost[i].ci;
		m_tag[ci0].tag = i+1;
	}
	for (int ig = 0; ig < ighost; ++ig)
	{
		bs = f_res_ghost[ig].ps();
		be = f_res_ghost[ig].pe();
		for (int i = bs; i < be; ++i)
		{
			int ci0 = f_res_ghost[ig][i].ci;
			m_tag[ci0].tag = -(ig+2);
		}
	}
	MPI_Barrier(share_comm);
	TellBlockPairNewTag();
}

void AMRLevel::CreatGhostProPair()
{
	g_pro.setnum_nocopy(0,0);	
	vector<int> new_g_pro;
	int bs = level_pro_ghost.ps();
	int be = level_pro_ghost.pe();
	for (int i = bs; i < be; ++i)
	{
		int ci0 = level_pro_ghost[i].ci;
		if (m_tag[ci0].tag > 0)
		{
			// PRINTFinLEVEL("level_pro_ghost for g_pro box (%d,%d,%d) type is %d!!!", cur_level,
			// 	m_box[ci0].ix(), m_box[ci0].iy(), m_box[ci0].iz(), m_box[ci0].type);
			for (Point_iterator p(0,3); p.end(); ++p)
			{
				int an0 = m_box[ci0].neib[p.i][p.j][p.k];
				if (an0 > -1)
				{
					if (m_box[an0].type == Normalcell)
					{
						new_g_pro.push_back(i);
						break;
					}
				}
			}
		}
	}
	g_pro.Addnew(new_g_pro);
	g_pro.DirectlyReduceNew();
}
