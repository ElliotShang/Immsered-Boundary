#include "AMRLevel.H"
#include "AMRSpaceTime.H"

#include <cmath>

void AMRLevel::Inittag()
{
	int mbnum = m_box.realsize();
	m_tag.setnum_nocopy(mbnum, 0);
	int bs = m_tag.ps();
	int be = m_tag.pe();
	for (int i = bs; i < be; ++i)
	{
		m_tag[i].tag = -1;
		m_tag[i].detag = -1;
	}
	MPI_Barrier(share_comm);
}

void AMRLevel::InitPair()
{
	int bs = m_box.ps();
	int be = m_box.pe();
	for (int i = bs; i < be; ++i)
	{
		m_box[i].ptype = none;
	}
	bs = m_box.gps();
	be = m_box.gpe();
	for (int i = bs; i < be; ++i)
	{
		m_box[i].ptype = none;
	}
	MPI_Barrier(share_comm);
}

void AMRLevel::AssociateCelltoBlockPair(const int & mylevel)
{		
	int bs, be;
	m_box.GlobalOrder(bs,be);
	for (int i = bs; i < be; ++i)
	{
		m_box[i].bkpid = -1;
	}
	MPI_Barrier(share_comm);
	bs = blockpair.ps();
	be = blockpair.pe();
	for (int i = bs; i < be; ++i)
	{
		NodePair & mybkp = blockpair[i];
		// for (int j = 0; j < blockpair.size(); ++j)
		// {
		// 	if (mybkp.outnode.node == blockpair[j].outnode.node &&
		// 		mybkp.outnode.index == blockpair[j].outnode.index &&
		// 		mybkp.innode == blockpair[j].innode && i!=j)
		// 	{
		// 		printf("N%dL%d blockpair %d and %d are the same innode %d (%d,%d,%d)"
		// 			"outnode is N%dB%d (%d,%d,%d)!!!\n",
		// 			node, mylevel, i,j,
		// 			mybkp.innode, m_box[mybkp.innode].ix(),
		// 			m_box[mybkp.innode].iy(),
		// 			m_box[mybkp.innode].iz(),
		// 			mybkp.outnode.node, mybkp.outnode.index,
		// 			m_box[mybkp.outnode.index].ix(),
		// 			m_box[mybkp.outnode.index].iy(),
		// 			m_box[mybkp.outnode.index].iz());
		// 		MPI_Abort(MPI_COMM_WORLD, 67);
		// 	}
		// }
		if (mybkp.outnode.node == node)
		{
			// if (mybkp.outnode.index == mybkp.innode)
			// {
			// 	printf("N%d The same-node blockpair outnode and innode both is %d (%d,%d,%d)\n",
			// 		node, mybkp.innode, m_box[mybkp.innode].ix(), m_box[mybkp.innode].iy(), m_box[mybkp.innode].iz());
			// 	MPI_Abort(MPI_COMM_WORLD, 57);
			// }
			if (m_box[mybkp.outnode.index].bkpid != -1)
			{
				printf("N%dL%d The blockpair is in the same node. Innode %d(%d,%d,%d) outnode %d(%d,%d,%d)"
					" but the outnode is occupied by box %d(%d,%d,%d)\n",
					node, mylevel,mybkp.innode,
					m_box[mybkp.innode].ix(),
					m_box[mybkp.innode].iy(),
					m_box[mybkp.innode].iz(),
					mybkp.outnode.index,
					m_box[mybkp.outnode.index].ix(),
					m_box[mybkp.outnode.index].iy(),
					m_box[mybkp.outnode.index].iz(),
					m_box[mybkp.outnode.index].bkpid,
					m_box[m_box[mybkp.outnode.index].bkpid].ix(),
					m_box[m_box[mybkp.outnode.index].bkpid].iy(),
					m_box[m_box[mybkp.outnode.index].bkpid].iz());
				MPI_Abort(MPI_COMM_WORLD, 67);
			}
			m_box[mybkp.outnode.index].bkpid = mybkp.innode;
		}
	}	
}

void AMRLevel::RefineLevel()
{
	// const int ptn0 = 2;
	// int icx0[ptn0] = {135,135};
	// int icy0[ptn0] = {32,33};
	// int icz0[ptn0] = {0,0};
	// for (int i = m_tag.ps(); i < m_tag.pe(); ++i)
	// 	{
	// 		int i0 = i;
	// 		int ix = m_box[i0].ix();
	// 		int iy = m_box[i0].iy();
	// 		int iz = m_box[i0].iz();
	// 		for (int p0s = 0; p0s < ptn0; ++p0s)
	// 		{
	// 			if (ix == icx0[p0s] && iy == icy0[p0s] && iz == icz0[p0s])
	// 			{
	// 				printf("Level %d before TagRefineGhost Box %d (%d,%d,%d) tag is %d\n", cur_level, i0, ix, iy,iz, m_tag[i0].tag);
	// 			}
	// 		}
	// 	}
	TagRefineGhost();
	GiveAFlag("Finish tag refine ghost!!!",5);
	AdjustOldPairInFluid_Refine();
	GiveAFlag("Finish AdjustOldPairInFluid_Refine", 5);
	BuildNewPair_Refine();
	GiveAFlag("Finish BuildNewPair_Refine!!!", 5);
	MakeSmallBox();
	GiveAFlag("Finish MakeSmallBox!!!", 5);
	FillGhostPair_Refine();
	GiveLevelFlag("Finish FillGhostPair_Refine!!!", cur_level, 5);
	IncludeNewboxinGhost_Refine();
	ConnectRefinedBox();
	GiveAFlag("Finish connect refined box!!!", 5);
	RemoveGhostandNormal_Refine();
	RemoveArrayGhost(m_box);
	RemoveArrayGhost(finelevel->m_box);
}

void AMRLevel::DerefineLevel()
{
	// const int ptn = 3;
	// int icx0[ptn] = {204, 204, 204};
	// int icy0[ptn] = {106, 108, 104};
	// int icz0[ptn] = {4, 4, 4};
	// for (int i = finelevel->m_tag.ps(); i < finelevel->m_tag.pe(); ++i)
	// 	{
	// 		int i0 = i;
	// 		int ix = finelevel->m_box[i0].ix();
	// 		int iy = finelevel->m_box[i0].iy();
	// 		int iz = finelevel->m_box[i0].iz();
	// 		for (int p0 = 0; p0 < ptn; ++p0)
	// 		{
	// 			if (ix == icx0[p0] && iy == icy0[p0] && iz == icz0[p0])
	// 		{
	// 			printf("N%d Level %d Before tag derefine ghost Box %d (%d,%d,%d) detag is %d type is %d\n", node, cur_level, i0, ix, iy,iz, finelevel->m_tag[i0].detag,finelevel->m_box[i0].type);
	// 		}
	// 		}	
	// 	}
	if (!finelevel->twodflag) TagDerefineGhost();
	else TagDerefineGhost_2d();
		// int bs, be;
		// finelevel->m_box.GlobalOrder(bs, be);
		// for (int i = bs; i < be; ++i)
		// {
		// 		if ((finelevel->m_box[i].ix() >= 314 && finelevel->m_box[i].ix() <= 315) &&
		// 			(finelevel->m_box[i].iy() >= 158 && finelevel->m_box[i].iy() <= 159) &&
		// 			(finelevel->m_box[i].iz() >= 26 && finelevel->m_box[i].iz() <= 29))
		// 		{
		// 			printf("Finelevel derefine Box (%d,%d,%d) detag is %d\n",
		// 				finelevel->m_box[i].ix(),
		// 				finelevel->m_box[i].iy(),
		// 				finelevel->m_box[i].iz(),
		// 				finelevel->m_tag[i].detag);
		// 		}
		// }
#ifdef DEBUG
	for (int ig = 0; ig < ighost; ++ig)
	{
		for (int i = finelevel->f_res_ghost[ig].ps(); i < finelevel->f_res_ghost[ig].pe(); ++i)
		{
			int ci0 = finelevel->f_res_ghost[ig][i].ci;
			if (finelevel->m_tag[ci0].detag !=-1)
			{
				PRINTFinLEVEL("The second finer level res %d main box (%d,%d,%d) detag is %d during coarsen", 
					cur_level, ig,
					finelevel->m_box[ci0].ix(),
					finelevel->m_box[ci0].iy(),
					finelevel->m_box[ci0].iz(),
					finelevel->m_tag[ci0].detag);
				MPI_Abort(MPI_COMM_WORLD,82);
			}
		}
	}
#endif		
	GiveAFlag("TagDerefineGhost!!!", 5);
	AdjustOldPair_Coarsen();
	BuildNewPair_Derefine();
	GiveAFlag("Start to merge box", 5);
	finelevel->MergeBox();
	GiveAFlag("Finish merge box", 5);
	FillGhostPair_Derefine();
	GiveLevelFlag("Finish FillGhostPair_Derefine!!!",cur_level,5);
	IncludeNewboxinGhost_Derefine();
	GiveLevelFlag("Finish IncludeNewboxinGhost_Derefine!!!",cur_level,5);
	finelevel->ConnectCoarseBox();
	GiveAFlag("Finish ConnectCoarseBox!!!", 5);
	RemoveGhostandNormal_Derefine();
	//finelevel->EndMeshHoleConnection();
	RemoveArrayGhost(m_box);
	RemoveArrayGhost(finelevel->m_box);
}

	void AMRLevel::AdjustMesh()
	{
		if (NULL != finelevel)
		{
			if (finelevel->toderefine)
			{
				DerefineLevel();		
			}
			finelevel->toderefine = false;
		}
		GiveAFlag("Finish de-refine level!!!", 5);
		if (torefine)
		{
			CheckRefineTagLocation();
			RefineLevel();
			torefine = false;
		}
		GiveAFlag("Finish refine level!!!", 5);
		MPI_Barrier(MPI_COMM_WORLD);
	}

void AMRLevel::RenewProandResIndex()
{
	/*2020-2-26*/
	/*Renew the ghost location in the m_box array*/
	//MPI_Win_fence(0, m_box.arraywin());
	vector<int> & boxghost = m_box.arrayghost();
	for (int i = 0; i < boxghost.size(); ++i)
	{
		int gi0 = boxghost[i];
		boxghost[i] = m_box[gi0].neib[1][1][1];
	}
	//GiveAFlag("Renew box ghost!!!", 5);
	/*Renew the ghost array tag*/
	for (int ig = 0; ig < ighost; ++ig)
	{
		for (int i = f_res_ghost[ig].ps(); i < f_res_ghost[ig].pe(); ++i)
		{
			int ci0 = f_res_ghost[ig][i].ci;
			//printf("Level res layer %d pair %d cbox %d new index %d\n", ig, i, ci0, m_box[ci0].neib[1][1][1].index);
			f_res_ghost[ig][i].ci = m_box[ci0].neib[1][1][1];
			Assert(f_res_ghost[ig][i].ci > -1, "the new index for the res main should be non-negative!!!", 107);
		}
	}
	GiveAFlag("Renew res ghost!!!", 5);
	for (int i = f_pro_ghost.ps(); i < f_pro_ghost.pe(); ++i)
	{
		int ci0 = f_pro_ghost[i].ci;
		//printf("Level pro pair %d cbox %d new index %d size %d\n", i, ci0, m_box[ci0].neib[1][1][1].index, m_box.realsize());
		f_pro_ghost[i].ci = m_box[ci0].neib[1][1][1];
		Assert(f_pro_ghost[i].ci > -1, "the new index for the pro main should be non-negative!!!", 116);
	}
	// for (int i = 0; i < drf_sliding_res_to_pro.size(); ++i)
	// {
	// 	int ci0 = drf_sliding_res_to_pro[i].ci;
	// 	drf_sliding_res_to_pro[i].ci = m_box[ci0].neib[1][1][1];
	// }
	GiveAFlag("Renew my ghost relation!!!", 5);
	if (NULL != coarselevel)
	{
		/*renew the refine box index*/
		if (!twodflag)
		{
			for (int i = 0; i < coarselevel->rfbox.size(); ++i)
			{
				int i0 = coarselevel->rfbox[i];
				int tag0 = coarselevel->m_tag[i0].tag;
				for (Point_iterator p(0,2); p.end(); ++p)
				{
					int sonid = coarselevel->allson[tag0].son[p.i][p.j][p.k];
					// printf("coarselevel son old index is %d\n", coarselevel->myson[i].son[p.i][p.j][p.k]);
					coarselevel->allson[tag0].son[p.i][p.j][p.k] = m_box[sonid].neib[1][1][1];
					Assert(m_box[sonid].neib[1][1][1] > -1, "the new index for the refined son should be non-negative!!!", 136);
					// printf("coarselevel son new index is %d\n", coarselevel->myson[i].son[p.i][p.j][p.k]);
				}
			}
			/*renew the ghost index stored in the coarser level*/
			for (int ig = 0; ig < ighost; ++ig)
			{
				for (int i = coarselevel->f_res_ghost[ig].ps(); 
								 i < coarselevel->f_res_ghost[ig].pe(); ++i)
				{
					for (Point_iterator p(0,2); p.end(); ++p)
					{
						int fi0 = coarselevel->f_res_ghost[ig][i].fi.son[p.i][p.j][p.k];
						coarselevel->f_res_ghost[ig][i].fi.son[p.i][p.j][p.k] = m_box[fi0].neib[1][1][1];
#ifdef DEBUG					
						if (!(m_box[fi0].neib[1][1][1] > -1))
						{
							int ci0 = coarselevel->f_res_ghost[ig][i].ci;
							PRINTFinLEVEL("RES %d negative index son main box (%d,%d,%d)",cur_level,ig,
								coarselevel->m_box[ci0].ix(),
								coarselevel->m_box[ci0].iy(),
								coarselevel->m_box[ci0].iz());
						}
#endif					
						Assert(m_box[fi0].neib[1][1][1] > -1, "the new index for the res son should be non-negative!!!", 150);
					}
				}
			}
			for (int i = coarselevel->f_pro_ghost.ps(); 
							 i < coarselevel->f_pro_ghost.pe(); ++i)
			{
				for (Point_iterator p(0,2); p.end(); ++p)
				{
					int fi0 = coarselevel->f_pro_ghost[i].fi.son[p.i][p.j][p.k];
					coarselevel->f_pro_ghost[i].fi.son[p.i][p.j][p.k] = m_box[fi0].neib[1][1][1];
#ifdef DEBUG				
					if (m_box[fi0].neib[1][1][1] < 0)
					{
						int ci0 = coarselevel->f_pro_ghost[i].ci;
						PRINTFinLEVEL("Pro son %d (%d,%d,%d) detag %d", cur_level,
							fi0, m_box[fi0].ix(), m_box[fi0].iy(), m_box[fi0].iz(), m_tag[fi0].detag);
					}
#endif				
					Assert(m_box[fi0].neib[1][1][1] > -1, "the new index for the pro son should be non-negative!!!", 161);
				}
			}
		}
		else //the finer level is 2d-refined
		{
			for (int i = 0; i < coarselevel->rfbox.size(); ++i)
			{
				int i0 = coarselevel->rfbox[i];
				int tag0 = coarselevel->m_tag[i0].tag;
				for (int nx = 0; nx < 2; ++nx)
				{
					for (int ny = 0; ny < 2; ++ny)
					{
						int sonid = coarselevel->allson[tag0].son[nx][ny][0];
						// printf("coarselevel son old index is %d\n", coarselevel->myson[i].son[p.i][p.j][p.k]);
						coarselevel->allson[tag0].son[nx][ny][0] = m_box[sonid].neib[1][1][1];
						Assert(m_box[sonid].neib[1][1][1] > -1, "the new index for the refined son should be non-negative!!!", 136);
						// printf("coarselevel son new index is %d\n", coarselevel->myson[i].son[p.i][p.j][p.k]);
					}
				}
			}
			/*renew the ghost index stored in the coarser level*/
			for (int ig = 0; ig < ighost; ++ig)
			{
				for (int i = coarselevel->f_res_ghost[ig].ps(); 
								 i < coarselevel->f_res_ghost[ig].pe(); ++i)
				{
					for (Point_iterator_2d p(0,2); p.end(); ++p)
					{
						int fi0 = coarselevel->f_res_ghost[ig][i].fi.son[p.i][p.j][p.k];
						coarselevel->f_res_ghost[ig][i].fi.son[p.i][p.j][p.k] = m_box[fi0].neib[1][1][1];
#ifdef DEBUG					
						if (!(m_box[fi0].neib[1][1][1] > -1))
						{
							int ci0 = coarselevel->f_res_ghost[ig][i].ci;
							PRINTFinLEVEL("RES %d negative index son main box (%d,%d,%d)",cur_level,ig,
								coarselevel->m_box[ci0].ix(),
								coarselevel->m_box[ci0].iy(),
								coarselevel->m_box[ci0].iz());
						}
#endif					
						Assert(m_box[fi0].neib[1][1][1] > -1, "the new index for the res son should be non-negative!!!", 150);
					}
				}
			}
			for (int i = coarselevel->f_pro_ghost.ps(); 
							 i < coarselevel->f_pro_ghost.pe(); ++i)
			{
				for (Point_iterator_2d p(0,2); p.end(); ++p)
				{
					int fi0 = coarselevel->f_pro_ghost[i].fi.son[p.i][p.j][p.k];
					coarselevel->f_pro_ghost[i].fi.son[p.i][p.j][p.k] = m_box[fi0].neib[1][1][1];
#ifdef DEBUG				
					if (m_box[fi0].neib[1][1][1] < 0)
					{
						int ci0 = coarselevel->f_pro_ghost[i].ci;
						PRINTFinLEVEL("Pro son %d (%d,%d,%d) detag %d", cur_level,
							fi0, m_box[fi0].ix(), m_box[fi0].iy(), m_box[fi0].iz(), m_tag[fi0].detag);
					}
#endif				
					Assert(m_box[fi0].neib[1][1][1] > -1, "the new index for the pro son should be non-negative!!!", 161);
				}
			}
		}
	}
	
	if (NULL != finelevel)
	{
		for (int i = 0; i < finelevel->dfnum_procs[srank]; ++i)
		{
			int dfboxid = finelevel->derfbox[i];
			int momboxid = finelevel->m_tag[dfboxid].detag;
			finelevel->m_tag[dfboxid].detag = m_box[momboxid].neib[1][1][1];
			Assert(m_box[momboxid].neib[1][1][1] > -1, "the new index for the coarsened detag should be non-negative!!!", 181);
		}
	}
	//MPI_Win_fence(0, m_box.arraywin());
	MPI_Barrier(share_comm);
/*End for the add in 2020-2-26*/
}

/*--------------------------------------------------*/
/*This function makes small boxes at each refined box.*/
/*--------------------------------------------------*/
void AMRLevel::MakeSmallBox()
	{
		vector<Box> newbox;
		//DEFINE_MPI_BOXSON();
		vector<Boxson<int> > 				myson(0);
		int tagboxsize = taggedbox.size();
		if (!finelevel->twodflag)
		{
			for (int i0 = 0; i0 < tagboxsize; ++i0)
			{
				int i = taggedbox[i0];
				if (m_tag[i].tag != -1)
				{
					rfbox.push_back(i);
					myson.push_back(Boxson<int>(-1));
					int lptx = m_box[i].ix();
					int lpty = m_box[i].iy();
					int lptz = m_box[i].iz();
					int mylpx, mylpy, mylpz;
					for (Point_iterator p(0,2); p.end(); ++p)
					{
						mylpx = 2*lptx+p.i;
						mylpy = 2*lpty+p.j;
						mylpz = 2*lptz+p.k;
						newbox.push_back(Box(mylpx, mylpy, mylpz));
						Box & anewbox = newbox.back();
						anewbox.type = Normalcell;
						anewbox.pair = m_box[i].pair;
						anewbox.solid = m_box[i].solid;
						// anewbox.pair.body = m_box[i].pair.body;
						// anewbox.pair.patch = m_box[i].pair.patch;
					}
					/*for the cell to be refined*/
					if (m_tag[i].tag == 0)
					{
						for (Point_iterator p(0,2); p.end(); ++p)
						{
							/*-1 means the new box is not a ghost*/
							newrefinebox_ghosttag.push_back(-1);
						}
					}
					/*For the cell to be the ghost*/
					else
					{
						// if (m_tag[i].tag != -2)
						// {
						// 	PRINTFinLEVEL("Box (%d,%d,%d) tag %d",cur_level,m_box[i].ix(),m_box[i].iy(),m_box[i].iz(),m_tag[i].tag);
						// }
						Assert(m_tag[i].tag == -2, "the refine tag is not 0 or 2!!!", 234);
						for (Point_iterator p(0,2); p.end(); ++p)
						{
							/*0 means the new box is a ghost*/
							newrefinebox_ghosttag.push_back(0);
						}
					}
					m_tag[i].tag = rfbox.size() - 1;
				}
			}
		}
		else
		{
			for (int i0 = 0; i0 < tagboxsize; ++i0)
			{
				int i = taggedbox[i0];
				if (m_tag[i].tag != -1)
				{
					rfbox.push_back(i);
					myson.push_back(Boxson<int>(-1));
					int lptx = m_box[i].ix();
					int lpty = m_box[i].iy();
					int lptz = m_box[i].iz();
					int mylpx, mylpy, mylpz;
					for (Point_iterator_2d p(0,2); p.end(); ++p)
					{
						mylpx = 2*lptx+p.i;
						mylpy = 2*lpty+p.j;
						mylpz = lptz;
						newbox.push_back(Box(mylpx, mylpy, mylpz));
						Box & anewbox = newbox.back();
						anewbox.type = Normalcell;
						anewbox.pair = m_box[i].pair;
						anewbox.solid = m_box[i].solid;
						// anewbox.pair.body = m_box[i].pair.body;
						// anewbox.pair.patch = m_box[i].pair.patch;
					}
					/*for the cell to be refined*/
					if (m_tag[i].tag == 0)
					{
						for (Point_iterator_2d p(0,2); p.end(); ++p)
						{
							/*-1 means the new box is not a ghost*/
							newrefinebox_ghosttag.push_back(-1);
						}
					}
					/*For the cell to be the ghost*/
					else
					{
						// if (m_tag[i].tag != -2)
						// {
						// 	PRINTFinLEVEL("Box (%d,%d,%d) tag %d",cur_level,m_box[i].ix(),m_box[i].iy(),m_box[i].iz(),m_tag[i].tag);
						// }
						Assert(m_tag[i].tag == -2, "the refine tag is not 0 or 2!!!", 234);
						for (Point_iterator_2d p(0,2); p.end(); ++p)
						{
							/*0 means the new box is a ghost*/
							newrefinebox_ghosttag.push_back(0);
						}
					}
					m_tag[i].tag = rfbox.size() - 1;
				}
			}
		}
		int refine_box_num = rfbox.size();
		ShowAllRankData("num_refined", refine_box_num, 4);
		MPI_Allgather(&refine_box_num, 1, MPI_INT,
			&rfnum_procs[0], 1, MPI_INT, share_comm);
		CountTotalNum(rfnum_procs, totrfnum);
		MPI_Barrier(share_comm);
		if (totrfnum > 0)
		{
			ArrayProcsStart(rfnum_procs, rfnum_start);
			ShowAllRankData("rfnum_start", rfnum_start[srank], 4);
			int s0;
			if (!finelevel->twodflag)
			{
				for (int i = 0; i < refine_box_num; ++i)
				{
					int i0 = rfbox[i];
					m_tag[i0].tag += rfnum_start[srank];		
					s0 = 0;
					for (Point_iterator p(0,2); p.end(); ++p)
					{
						myson[i].son[p.i][p.j][p.k] = finelevel->m_box.realsize()+(rfnum_start[srank]+i)*8+s0;
						++s0;
						newrefinebox.push_back(myson[i].son[p.i][p.j][p.k]);
					}
				}
			}
			else
			{
				for (int i = 0; i < refine_box_num; ++i)
				{
					int i0 = rfbox[i];
					m_tag[i0].tag += rfnum_start[srank];		
					s0 = 0;
					for (Point_iterator_2d p(0,2); p.end(); ++p)
					{
						myson[i].son[p.i][p.j][p.k] = finelevel->m_box.realsize()+(rfnum_start[srank]+i)*4+s0;
						++s0;
						newrefinebox.push_back(myson[i].son[p.i][p.j][p.k]);
					}
				}
			}
		}
		allson.Addnew(myson);
		allson.DirectlyReduceNew();
		finelevel->m_box.Addnew(newbox, newrefinebox_ghosttag);
		MPI_Barrier(share_comm);
	}

void AMRLevel::ConnectRefinedBox()
	{
		int tss[3] = {0, 0, 1};
		int tse[3] = {0, 1, 1};

		int sns[4] = {0, 1, 0, 2};
		int sne[4] = {0, 2, 1, 2};

		int ms[4] = {1, 0, 1, 0};
		//int judflag;
		//MPI_Win_fence(0, allson.arraywin());
		if (!finelevel->twodflag)
		{
			for (int i = 0; i < rfbox.size(); ++i)
			{
				int box0 = rfbox[i];
				//PRINTFinLEVEL("refine box %d is (%d,%d,%d)", cur_level, i, m_box[box0].ix(), m_box[box0].iy(), m_box[box0].iz());
				Assert(m_tag[box0].tag > -1, "A refined box tag error!!!", 442);
				for (Point_iterator m(0,3); m.end(); ++m)
				{
					int momneibloc = m_box[box0].neib[m.i][m.j][m.k];
					if (m_box.isnorg(momneibloc))
					{
						//Assert(m_box.isnorg(momneibloc), "Connect refined box error!!!", 418);
						if (m_tag[momneibloc].tag > -1)
						{	
							int tag0 = m_tag[momneibloc].tag;
							for (int tsi = tss[m.i]; tsi <= tse[m.i]; ++tsi)
							{
								for (int tsj = tss[m.j]; tsj <= tse[m.j]; ++tsj)
								{
									for (int tsk = tss[m.k]; tsk <= tse[m.k] ; ++tsk)
									{
										int sonid = allson[m_tag[box0].tag].son[tsi][tsj][tsk];
										int subx = m.i+tsi;
										int suby = m.j+tsj;
										int subz = m.k+tsk;							
										for (int sni = sns[subx]; sni <= sne[subx]; ++sni)
										{
											for (int snj = sns[suby]; snj <= sne[suby]; ++snj)
											{
												for (int snk = sns[subz]; snk <= sne[subz]; ++snk)
												{
													int sonneibloc = allson[tag0].
														son[ms[tsi+sni]][ms[tsj+snj]][ms[tsk+snk]]; //something wrong!!!	
													Assert(sonneibloc > -1, "A negative son in connect refined box!!!", 514);
													finelevel->m_box[sonid].setneib(sni, snj, snk, sonneibloc);
													finelevel->m_box[sonneibloc].setneib(2-sni, 2-snj, 2-snk, sonid);
												}
											}
										}
									}
								}
							}
						}
					}		
				}
			}
		}
		else
		{
			for (int i = 0; i < rfbox.size(); ++i)
			{
				int box0 = rfbox[i];
				//PRINTFinLEVEL("refine box %d is (%d,%d,%d)", cur_level, i, m_box[box0].ix(), m_box[box0].iy(), m_box[box0].iz());
				Assert(m_tag[box0].tag > -1, "A refined box tag error!!!", 442);
				for (Point_iterator m(0,3); m.end(); ++m)
				{
					int momneibloc = m_box[box0].neib[m.i][m.j][m.k];
					if (m_box.isnorg(momneibloc))
					{
						//Assert(m_box.isnorg(momneibloc), "Connect refined box error!!!", 418);
						if (m_tag[momneibloc].tag > -1)
						{	
							int tag0 = m_tag[momneibloc].tag;
							for (int tsi = tss[m.i]; tsi <= tse[m.i]; ++tsi)
							{
								for (int tsj = tss[m.j]; tsj <= tse[m.j]; ++tsj)
								{
									int sonid = allson[m_tag[box0].tag].son[tsi][tsj][0];
									int subx = m.i+tsi;
									int suby = m.j+tsj;
									for (int sni = sns[subx]; sni <= sne[subx]; ++sni)
									{
										for (int snj = sns[suby]; snj <= sne[suby]; ++snj)
										{
											int sonneibloc = allson[tag0].
												son[ms[tsi+sni]][ms[tsj+snj]][0]; //something wrong!!!	
											Assert(sonneibloc > -1, "A negative son in connect refined box!!!", 514);
											finelevel->m_box[sonid].setneib(sni, snj, m.k, sonneibloc);
											finelevel->m_box[sonneibloc].setneib(2-sni, 2-snj, 2-m.k, sonid);
										}
									}
								}
							}
						}
					}		
				}
			}
		}
		MPI_Barrier(share_comm);
		GiveAFlag("Connect refined part!!!",5);
	}

int psns[3] = {1, 1, 1};
int psne[3] = {1, 2, 1};

int psums[3] = {3, 2, 1};
int psume[3] = {3, 3, 1};

void AMRLevel::TrimDetag()
{
	int bs = m_tag.ps();
	int be = m_tag.pe();
	if (!twodflag)
	{
		for (int i = bs; i < be; ++i)
		{	
			Assert(m_tag[i].detag != -5, "Detag -5 in trim !!!", 608);
			if (m_tag[i].detag > -1)
			{
				bool derf_this_cell = true;
				int sx = m_box[i].ix() & 1;
				/*In 2D condition, sx must be 0*/
				int sy = m_box[i].iy() & 1;
				int sz = m_box[i].iz() & 1;
				int ss[3] = {1-sx, 1-sy, 1-sz};
				int se[3] = {3-sx, 3-sy, 3-sz};
				for (int nx = ss[0]; nx < se[0]; ++nx)
				{
					for (int ny = ss[1]; ny < se[1]; ++ny)
					{
						for (int nz = ss[2]; nz < se[2]; ++nz)
						{
							int an0 = m_box[i].neib[nx][ny][nz];
							Assert(an0>-1, "The four cells must be non-negative when TrimDetag!!!", 526);
							if (m_tag[an0].detag == -1)
							{
								derf_this_cell = false;
								break;
							}
						}
					}
				}
				if (!derf_this_cell)
				{
					for (int nx = ss[0]; nx < se[0]; ++nx)
					{
						for (int ny = ss[1]; ny < se[1]; ++ny)
						{
							for (int nz = ss[2]; nz < se[2]; ++nz)
							{
								int an0 = m_box[i].neib[nx][ny][nz];
								m_tag[an0].detag = -1;
							}
						}
					}
				}
			}
		}
	}
	else
	{
		for (int i = bs; i < be; ++i)
		{	
			Assert(m_tag[i].detag != -5, "Detag -5 in trim !!!", 608);
			if (m_tag[i].detag > -1)
			{
				bool derf_this_cell = true;
				int sx = m_box[i].ix() & 1;
				/*In 2D condition, sx must be 0*/
				int sy = m_box[i].iy() & 1;
				int ss[2] = {1-sx, 1-sy};
				int se[2] = {3-sx, 3-sy};
				for (int nx = ss[0]; nx < se[0]; ++nx)
				{
					for (int ny = ss[1]; ny < se[1]; ++ny)
					{
						int an0 = m_box[i].neib[nx][ny][1];
						Assert(an0>-1, "The four cells must be non-negative when TrimDetag!!!", 526);
						if (m_tag[an0].detag == -1)
						{
							derf_this_cell = false;
							break;
						}
					}
				}
				if (!derf_this_cell)
				{
					for (int nx = ss[0]; nx < se[0]; ++nx)
					{
						for (int ny = ss[1]; ny < se[1]; ++ny)
						{
							int an0 = m_box[i].neib[nx][ny][1];
							m_tag[an0].detag = -1;
						}
					}
				}
			}
		}
	}
	MPI_Barrier(share_comm);
}

	void AMRLevel::MergeBox()	
	{	
		vector<Box> newbox;
		int bs = m_tag.ps();
		int be = m_tag.pe();
		if (!twodflag)
		{
			for (int i = bs; i < be; ++i)
			{	
				Assert(m_tag[i].detag < 1, "Invalid detag when merging boxes!!!", 533);
				if (m_tag[i].detag == 0 || m_tag[i].detag == -2)
				{
					int sx = m_box[i].ix() & 1;
					/*In 2D condition, sx must be 0*/
					int sy = m_box[i].iy() & 1;
					int sz = m_box[i].iz() & 1;
					if (sx + sy + sz == 0)
					{
						derfbox.push_back(i);					
						int bx = m_box[i].ix()/2; 
						int by = m_box[i].iy()/2; 
						int bz = m_box[i].iz()/2;
						newbox.push_back(Box(bx, by, bz));
						Box & thebackbox = newbox.back();
						thebackbox.type = Normalcell;
						// thebackbox.pair.body = m_box[i].pair.body;
						// thebackbox.pair.patch = m_box[i].pair.patch;
						thebackbox.pair = m_box[i].pair;
						thebackbox.solid = m_box[i].solid;
						if (m_tag[i].detag == 0)
						{
							newcoarsebox_ghosttag.push_back(-1);
						}
						else
						{					
							newcoarsebox_ghosttag.push_back(0);
						}
						m_tag[i].detag = derfbox.size()-1;
					}		
				}
			}
		}
		else
		{
			for (int i = bs; i < be; ++i)
			{	
				Assert(m_tag[i].detag < 1, "Invalid detag when merging boxes!!!", 533);
				if (m_tag[i].detag == 0 || m_tag[i].detag == -2)
				{
					int sx = m_box[i].ix() & 1;
					/*In 2D condition, sx must be 0*/
					int sy = m_box[i].iy() & 1;
					if (sx + sy == 0)
					{
						derfbox.push_back(i);					
						int bx = m_box[i].ix()/2; 
						int by = m_box[i].iy()/2; 
						newbox.push_back(Box(bx, by, m_box[i].iz()));
						Box & thebackbox = newbox.back();
						thebackbox.type = Normalcell;
						// thebackbox.pair.body = m_box[i].pair.body;
						// thebackbox.pair.patch = m_box[i].pair.patch;
						thebackbox.pair = m_box[i].pair;
						thebackbox.solid = m_box[i].solid;
						if (m_tag[i].detag == 0)
						{
							newcoarsebox_ghosttag.push_back(-1);
						}
						else
						{					
							newcoarsebox_ghosttag.push_back(0);
						}
						m_tag[i].detag = derfbox.size()-1;
					}		
				}
			}
		}
		/*detag in other processors has been modified and the data should be synchronized before futher calculation*/
		MPI_Barrier(MPI_COMM_WORLD);
		int derefi_box_num = derfbox.size();
		//ShowAllRankData("num_derefined", derefi_box_num, 5);
		MPI_Allgather(&derefi_box_num, 1, MPI_INT,
			&dfnum_procs[0], 1, MPI_INT, share_comm);
		CountTotalNum(dfnum_procs, totdfnum);
		// ShowAllRankData("dfnum_procs", dfnum_procs[srank], 5);
		if (totdfnum > 0)
		{
			ArrayProcsStart(dfnum_procs, dfnum_start);
			// ShowAllRankData("dfnum_start", dfnum_start[srank], -1);
			int existcboxnum = coarselevel->m_box.realsize();
			if (!twodflag)
			{
				for (int i = 0; i < dfnum_procs[srank]; ++i)
				{
					int i0 = derfbox[i];
					int mydetag = m_tag[i0].detag + dfnum_start[srank];
					m_tag[i0].detag	= existcboxnum+mydetag;
					newcoarsebox.push_back(m_tag[i0].detag);
					newbox[i].neib[1][1][1] = m_tag[i0].detag;			
					for (Point_iterator p(1,3); p.end(); ++p)
					{
						int sonneib = m_box[i0].neib[p.i][p.j][p.k];
						m_tag[sonneib].detag = m_tag[i0].detag;
					}
				}
			}
			else
			{
				for (int i = 0; i < dfnum_procs[srank]; ++i)
				{
					int i0 = derfbox[i];
					int mydetag = m_tag[i0].detag + dfnum_start[srank];
					m_tag[i0].detag	= existcboxnum+mydetag;
					newcoarsebox.push_back(m_tag[i0].detag);
					newbox[i].neib[1][1][1] = m_tag[i0].detag;			
					for (Point_iterator_2d p(1,3); p.end(); ++p)
					{
						int sonneib = m_box[i0].neib[p.i][p.j][1];
						m_tag[sonneib].detag = m_tag[i0].detag;
					}
				}
			}
		}
		coarselevel->m_box.Addnew(newbox, newcoarsebox_ghosttag);
	}

	void AMRLevel::ConnectCoarseBox()
	{

		int fnb[3] = {1, 1, 2};
		//printf("derfbox size is %d newcoarsebox size is %d\n", (int)derfbox.size(),(int)newcoarsebox.size());
		if (!twodflag)
		{
			for (int i = 0; i < derfbox.size(); ++i)
			{
				int i0 = derfbox[i];
				int cbox0 = newcoarsebox[i];
				for (Point_iterator p(0,3); p.end(); ++p)
				{					
					int rvny = 2-p.j; int rvnz = 2-p.k; int rvnx = 2-p.i;	
					int fineneib = m_box[i0].neib[fnb[p.i]][fnb[p.j]][fnb[p.k]];
					int fnn = m_box[fineneib].neib[p.i][p.j][p.k];
					if (fnn > -1)
					{
						if (m_tag[fnn].detag > -1)
						{
							coarselevel->m_box[cbox0].setneib(p.i, p.j, p.k, m_tag[fnn].detag);
							coarselevel->m_box[m_tag[fnn].detag].setneib(rvnx, rvny, rvnz, cbox0);							
						}
					}
				}
			}
		}
		else
		{
			for (int i = 0; i < derfbox.size(); ++i)
			{
				int i0 = derfbox[i];
				int cbox0 = newcoarsebox[i];
				for (Point_iterator p(0,3); p.end(); ++p)
				{					
					int rvny = 2-p.j; int rvnz = 2-p.k; int rvnx = 2-p.i;
					int fineneib = m_box[i0].neib[fnb[p.i]][fnb[p.j]][1];
					int fnn = m_box[fineneib].neib[p.i][p.j][p.k];
					if (fnn > -1)
					{
						if (m_tag[fnn].detag > -1)
						{
							coarselevel->m_box[cbox0].setneib(p.i, p.j, p.k, m_tag[fnn].detag);
							coarselevel->m_box[m_tag[fnn].detag].setneib(rvnx, rvny, rvnz, cbox0);							
						}
					}
				}
			}
		}
		MPI_Barrier(share_comm);
	}

inline const int & AMRLevel::LevelID()
{
	return cur_level;
}

// bool AMRLevel::notnewbox(const int & a_box) const
// {
// 	if (a_box < m_tag.size() && a_box > -1)
// 	{
// 		return true;
// 	}
// 	else
// 	{
// 		return false;
// 	}
// }

// bool AMRLevel::notnewbox(const Boxloc * a_box) const
// {
// 	if (a_box->index < m_tag.size() && a_box->index > -1)
// 	{
// 		return true;
// 	}
// 	else
// 	{
// 		return false;
// 	}
// }

void AMRLevel::clearlevel()
{
	m_tag.cleararray();
}

void AMRLevel::CheckRefineTagLocation()
	{
		if (NULL != coarselevel)
		{
			if (!twodflag)
			{
				for (int ig = 0; ig < ighost; ++ig)
				{
					int bs = coarselevel->f_res_ghost[ig].ps();
					int be = coarselevel->f_res_ghost[ig].pe();
					for (int i = bs; i < be; ++i)
					{
						for (Point_iterator p(0,2); p.end(); ++p)
						{
							int fi0 = coarselevel->f_res_ghost[ig][i].fi.son[p.i][p.j][p.k];
							if (m_box.isnorg(fi0))
							{
								if (m_tag[fi0].tag == 0)
								{
									printf("Level %d box (%d,%d,%d) signdis %f is a son of coarselevel res pair but it is going to be refined!!! Please enlarge this level or shrink the fine level!!!!!!\n",
										cur_level,
										m_box[fi0].ix(),m_box[fi0].iy(),m_box[fi0].iz(),m_box[fi0].pair.signdis);
									MPI_Abort(MPI_COMM_WORLD, 1009);
								}
								m_tag[fi0].tag = -1;
							}
						}
					}
				}
			}
			else
			{
				for (int ig = 0; ig < ighost; ++ig)
				{
					int bs = coarselevel->f_res_ghost[ig].ps();
					int be = coarselevel->f_res_ghost[ig].pe();
					for (int i = bs; i < be; ++i)
					{
						for (Point_iterator_2d p(0,2); p.end(); ++p)
						{
							int fi0 = coarselevel->f_res_ghost[ig][i].fi.son[p.i][p.j][p.k];
							if (m_box.isnorg(fi0))
							{
								if (m_tag[fi0].tag == 0)
								{
									printf("Level %d box (%d,%d,%d) signdis %f is a son of coarselevel res pair but it is going to be refined!!! Please enlarge this level or shrink the fine level!!!!!!\n",
										cur_level,
										m_box[fi0].ix(),m_box[fi0].iy(),m_box[fi0].iz(),m_box[fi0].pair.signdis);
									// for (Point_iterator q(0,3); q.end(); ++q)
									// {
									// 	int anb = m_box[fi0].neib[q.i][q.j][q.k];
									// 	if (m_box.isnorg(anb))
									// 	{
									// 		printf("Box (%d,%d,%d) tag is %d\n", m_box[anb].ix(),
									// 			m_box[anb].iy(),m_box[anb].iz(),m_tag[anb].tag);
									// 	}
									// }
									MPI_Abort(MPI_COMM_WORLD, 1009);
								}
								m_tag[fi0].tag = -1;
							}
						}
					}
				}
			}
		}
		//MPI_Win_fence(0, m_level[ilevel].m_tag.arraywin());
		MPI_Barrier(share_comm);
		int tag0 = taggedbox.size();
		for (int i = 0; i < tag0; ++i)
		{
			int i0 = taggedbox[i];
			if (m_tag[i0].tag == -1)
			{
				taggedbox[i] = taggedbox.back();
				taggedbox.pop_back();
				--i;
				--tag0;
			}
		}
		MPI_Barrier(share_comm);
	}

void AMRLevel::MoveBox_Refine(const int & dir, const int & ddir, int & tagnum)
{
	Point nbbox(1,1,1);
	*(&nbbox[0]+dir) = 1-ddir;
	int bs = f_pro_ghost.ps();
	int be = f_pro_ghost.pe();
	for (int i = bs; i < be; ++i)
	{
		int ci0 = f_pro_ghost[i].ci;
		int nb_neg = m_box[ci0].neib[nbbox[0]][nbbox[1]][nbbox[2]];
		if (m_box[nb_neg].ptype == f_res1 || m_box[nb_neg].ptype == ghost_res1)
		{
			m_tag[ci0].givetag();
			taggedbox.push_back(ci0);
			++tagnum;
		}
	}
	MPI_Barrier(share_comm);
}

void AMRLevel::MoveBox_Refine_Rotor(const int & dir, const int & ddir, int & tagnum, const int & rs_sep)
{
	Point nbbox(1,1,1);
	*(&nbbox[0]+dir) = 1-ddir;
	int bs = f_pro_ghost.ps();
	int be = f_pro_ghost.pe();
	for (int i = bs; i < be; ++i)
	{
		int ci0 = f_pro_ghost[i].ci;
		if (m_box[ci0].ix() <= rs_sep)
		{
			int nb_neg = m_box[ci0].neib[nbbox[0]][nbbox[1]][nbbox[2]];
			if (m_box[nb_neg].ptype == f_res1 || m_box[nb_neg].ptype == ghost_res1)
			{
				m_tag[ci0].givetag();
				taggedbox.push_back(ci0);
				++tagnum;
			}
		}
	}
	MPI_Barrier(share_comm);
}

void AMRLevel::MoveBox_Derefine(const int & dir, const int & ddir, int & detagnum)
{
	Point nbbox(1,1,1);
	*(&nbbox[0]+dir) = 1-ddir;
	int bs = coarselevel->f_res_ghost[0].ps();
	int be = coarselevel->f_res_ghost[0].pe();
	if (twodflag)
	{
		for (int i = bs; i < be; ++i)
		{
			int ci0 = coarselevel->f_res_ghost[0][i].ci;
			int nb_neg = coarselevel->m_box[ci0].neib[nbbox[0]][nbbox[1]][nbbox[2]];
			if (nb_neg > -1)
			{
				if (coarselevel->m_box[nb_neg].ptype == f_pro || coarselevel->m_box[nb_neg].ptype == ghost_pro)
				{
					for (Point_iterator_2d p(0,2); p.end(); ++p)
					{
						int fi0 = coarselevel->f_res_ghost[0][i].fi.son[p.i][p.j][p.k];
						m_tag[fi0].givedetag();
						++detagnum;
					}
				}
			}
		}
	}
	else
	{
		for (int i = bs; i < be; ++i)
		{
			int ci0 = coarselevel->f_res_ghost[0][i].ci;
			int nb_neg = coarselevel->m_box[ci0].neib[nbbox[0]][nbbox[1]][nbbox[2]];
			if (nb_neg > -1)
			{
				if (coarselevel->m_box[nb_neg].ptype == f_pro || coarselevel->m_box[nb_neg].ptype == ghost_pro)
				{
					for (Point_iterator p(0,2); p.end(); ++p)
					{
						int fi0 = coarselevel->f_res_ghost[0][i].fi.son[p.i][p.j][p.k];
						m_tag[fi0].givedetag();
						++detagnum;
					}
				}
			}
		}
	}
	MPI_Barrier(share_comm);
}

void AMRLevel::MoveBox_Derefine_Rotor(const int & dir, const int & ddir, int & detagnum, const int & rs_sep)
{
	Point nbbox(1,1,1);
	*(&nbbox[0]+dir) = 1-ddir;
	int bs = coarselevel->f_res_ghost[0].ps();
	int be = coarselevel->f_res_ghost[0].pe();
	if (twodflag)
	{
		for (int i = bs; i < be; ++i)
		{
			int ci0 = coarselevel->f_res_ghost[0][i].ci;
			if (coarselevel->m_box[ci0].ix() <= rs_sep)
			{
				int nb_neg = coarselevel->m_box[ci0].neib[nbbox[0]][nbbox[1]][nbbox[2]];
				if (nb_neg > -1)
				{
					if (coarselevel->m_box[nb_neg].ptype == f_pro || coarselevel->m_box[nb_neg].ptype == ghost_pro)
					{
						for (Point_iterator_2d p(0,2); p.end(); ++p)
						{
							int fi0 = coarselevel->f_res_ghost[0][i].fi.son[p.i][p.j][p.k];
							m_tag[fi0].givedetag();
							++detagnum;
						}
					}
				}
			}
		}
	}
	else
	{
		for (int i = bs; i < be; ++i)
		{
			int ci0 = coarselevel->f_res_ghost[0][i].ci;
			if (coarselevel->m_box[ci0].ix() <= rs_sep)
			{
				int nb_neg = coarselevel->m_box[ci0].neib[nbbox[0]][nbbox[1]][nbbox[2]];
				if (nb_neg > -1)
				{
					if (coarselevel->m_box[nb_neg].ptype == f_pro || coarselevel->m_box[nb_neg].ptype == ghost_pro)
					{
						for (Point_iterator p(0,2); p.end(); ++p)
						{
							int fi0 = coarselevel->f_res_ghost[0][i].fi.son[p.i][p.j][p.k];
							m_tag[fi0].givedetag();
							++detagnum;
						}
					}
				}
			}
		}
	}
	MPI_Barrier(share_comm);
}


