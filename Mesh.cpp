#include "Mesh.H"

vector<bool> level_twod_flag;
vector<int> level_power_ratio;

#ifndef PASSAGE_ANGLE
#define NewFlowMesg(ilevel, bkp)\
  {\
    m_level[ilevel].nd.cell_extrac[bkp.outnode.node].push_back(ExtractMesg(bkp.outnode.index,bkp.innode));\
  }
#else 
#define NewFlowMesg(ilevel, bkp)\
  {\
    m_level[ilevel].nd.cell_extrac[bkp.outnode.node].push_back(ExtractMesg(bkp.outnode.index,bkp.innode));\
    m_level[ilevel].nd.cell_extrac[bkp.outnode.node].back().rotangle = bkp.theta;\
  }
#endif    
void Mesh::ConstructBoxFace(const int & ilevel)
{
	//MPI_Win_fence(0,abox.arraywin());
	//MPI_Win_fence(0,aface.arraywin());
	Point nb[3][2] = {{Point(0,1,1),Point(2,1,1)},
										{Point(1,0,1),Point(1,2,1)},
										{Point(1,1,0),Point(1,1,2)}};
	vector<Face> nface;
	//MPI_Win_fence(0,abox.arraywin());
	int bs = m_level[ilevel].m_box.ps();
	int be = m_level[ilevel].m_box.pe();
	for (int f0 = 0; f0 < 2; ++f0)
	{
		for (int fd = 0; fd < 3; ++fd)
		{	
			for (int i = bs; i < be; ++i)
			{	
				int aneib = m_level[ilevel].m_box[i].neib[nb[fd][f0][0]][nb[fd][f0][1]][nb[fd][f0][2]];	
				if (m_level[ilevel].m_box[i].faces[fd][f0] < 0)
				{
					Assert(aneib > -1, "The face side of a normal box must be non-negative!!!", 30);
					if (aneib > -1)
					{
						if (m_level[ilevel].m_box[aneib].faces[fd][1-f0] < 0)
						{
							m_level[ilevel].m_box[aneib].faces[fd][1-f0] = 0;
							m_level[ilevel].m_box[i].faces[fd][f0] = 0;

							nface.push_back(Face());
							Face & lastface = nface.back();
							lastface[1-f0] = i;
							lastface[f0] = aneib;
							lastface.fnv = fd;
							m_level[ilevel].NewFaceArea(lastface, ilevel, fd, i, f0);
							// printf("N%dL%dB%d (%d,%d,%d) Dir %d Face %d area is %f keisa (%f,%f,%f)\n", 
							// 	node, ilevel,
							// 	i,
							// 	m_level[ilevel].m_box[i].ix(),
							// 	m_level[ilevel].m_box[i].iy(),
							// 	m_level[ilevel].m_box[i].iz(), 
							// 	fd, f0, nface.back().area,nface.back().keisa[0],nface.back().keisa[1],nface.back().keisa[2]);													
						}
						else
						{
							int s0 = m_level[ilevel].m_box[aneib].faces[fd][1-f0];
							m_level[ilevel].m_box[i].faces[fd][f0] = s0;
							m_level[ilevel].m_face[s0][1-f0] = i;
						}
					}
				}
				else
				{
					if (m_level[ilevel].m_box[aneib].faces[fd][1-f0] < 0)
					{
						int & fid = m_level[ilevel].m_box[aneib].faces[fd][1-f0];
						int s0 = m_level[ilevel].m_box[i].faces[fd][f0];
						fid = s0;
						m_level[ilevel].m_face[fid][f0] = aneib;
					}
				}
			}
		}
		MPI_Barrier(share_comm);
	}
	m_level[ilevel].IncludeNewFace(nface);
}

struct mptloc
{
	int b0;
	int px[3];

	mptloc()
	{}

	mptloc(const int & i1, const int & i2, const int & i3, const int & i4)
	{
		b0 = i1;
		px[0] = i2;
		px[1] = i3;
		px[2] = i4;
	}
};
void Mesh::ConstructBoxPoint(DataArray<Box> & abox, 
														 DataArray<mPoint> & apt, 
														 const int & ilevel,
														 const int & startbox,
														 const int & endbox)
	{
		vector<mptloc> npt(0);
		int neibpts[4] = {1,0,1,0};
		for (int i = startbox; i < endbox; ++i)
		{
			for (Point_iterator q(0,2); q.end(); ++q)
			{
				int ptindex = abox[i].pts[q.i][q.j][q.k];
				if (ptindex == -1)
				{
					bool ptflag = true;
					int nxe = q.i+2; int nye = q.j+2; int nze = q.k+2;
					for (int ni = q.i; ni < nxe; ++ni)
					{
						for (int nj = q.j; nj < nye; ++nj)
						{
#if DIM == 3							
							for (int nk = q.k; nk < nze; ++nk)
							{
								int ptneib = abox[i].neib[ni][nj][nk];
								if (ptneib > -1)
								{									
									if (abs(abox[ptneib].ix()-abox[i].ix()) > 1 ||
											abs(abox[ptneib].iy()-abox[i].iy()) > 1 ||
											abs(abox[ptneib].iz()-abox[i].iz()) > 1)
									{}
									else
									{
										int neib_pti = abox[ptneib].pts[neibpts[q.i+ni]][neibpts[q.j+nj]][neibpts[q.k+nk]];
										if (neib_pti > -1)
										{
											abox[i].pts[q.i][q.j][q.k] = neib_pti;
											ptflag = false;
											goto FOUNDAPOINT;
										}
									}
								}
							}
#else 
							int nk = 1;
							int ptneib = abox[i].neib[ni][nj][nk];
							if (ptneib > -1)
							{									
								if (abs(abox[ptneib].ix()-abox[i].ix()) > 1 ||
										abs(abox[ptneib].iy()-abox[i].iy()) > 1)
								{}
								else
								{
									int neib_pti = abox[ptneib].pts[neibpts[q.i+ni]][neibpts[q.j+nj]][neibpts[q.k+nk]];
									if (neib_pti > -1)
									{
										abox[i].pts[q.i][q.j][q.k] = neib_pti;
										ptflag = false;
										goto FOUNDAPOINT;
									}
								}
							}
#endif													
						}
					}
					FOUNDAPOINT:;
					if (ptflag)
					{
						abox[i].pts[q.i][q.j][q.k] = 0;
						npt.push_back(mptloc(i, q.i, q.j, q.k));
					}
				}
			}
		}
		int ptnum_procs = npt.size();
		vector<int> newptsnum(sprocs,0);
		vector<int> ptnum_arrstart(sprocs,0);
		int totnewptnum;
		int existptnum = apt.realsize();
		MPI_Barrier(share_comm);
		//PRINTFinLEVEL("Back ground mesh local new point number is %d", 0, ptnum_procs);
		MPI_Allgather(&ptnum_procs, 1, MPI_INT, &newptsnum[0], 1, MPI_INT, share_comm);
		CountTotalNum(newptsnum, totnewptnum);
		ArrayProcsStart(newptsnum, ptnum_arrstart);
		npt.resize(ptnum_procs);
		vector<mPoint> makepoint(ptnum_procs);
		for (int ptn = 0; ptn < ptnum_procs; ++ptn)
		{
			int pt0 = existptnum+ptn+ptnum_arrstart[srank];
			int b0 = npt[ptn].b0;
			makepoint[ptn].index = pt0;
			makepoint[ptn][0] = double(abox[b0].ix()+npt[ptn].px[0])*dh[ilevel][0];
			makepoint[ptn][1] = double(abox[b0].iy()+npt[ptn].px[1])*dh[ilevel][1];
			makepoint[ptn][2] = double(abox[b0].iz()+npt[ptn].px[2])*dh[ilevel][2];
			//PRINTFinLEVEL("Back Mesh new point %d is (%f,%f,%f)", 0, ptn, makepoint[ptn][0], makepoint[ptn][1], makepoint[ptn][2]);
			int nxe = npt[ptn].px[0]+2; int nye = npt[ptn].px[1]+2;		 
			int nze = npt[ptn].px[2]+2;			
			for (int ni = npt[ptn].px[0]; ni < nxe; ++ni)
			{
				for (int nj = npt[ptn].px[1]; nj < nye; ++nj)
				{
#if DIM == 3					
					for (int nk = npt[ptn].px[2]; nk < nze; ++nk)
					{
						int ptneib = abox[b0].neib[ni][nj][nk];
						if (ptneib > -1)
						{
							if (abs(abox[ptneib].ix()-abox[b0].ix()) > 1 ||
									abs(abox[ptneib].iy()-abox[b0].iy()) > 1 ||
									abs(abox[ptneib].iz()-abox[b0].iz()) > 1)
							{}
							else
							{
								abox[ptneib].
									pts[neibpts[npt[ptn].px[0]+ni]]
								   	[neibpts[npt[ptn].px[1]+nj]]
								   	[neibpts[npt[ptn].px[2]+nk]] = pt0;
							}
						}					
					}
#else
					int nk = 1;
					int ptneib = abox[b0].neib[ni][nj][nk];
					if (ptneib > -1)
					{
						if (abs(abox[ptneib].ix()-abox[b0].ix()) > 1 ||
								abs(abox[ptneib].iy()-abox[b0].iy()) > 1)
						{}
						else
						{
							abox[ptneib].
								pts[neibpts[npt[ptn].px[0]+ni]]
								 	 [neibpts[npt[ptn].px[1]+nj]]
								   [neibpts[npt[ptn].px[2]+nk]] = pt0;
						}
					}
#endif					
				}
			}
		}
		apt.Addnew(makepoint);
		apt.DirectlyReduceNew();
		MPI_Barrier(share_comm);
	}

void Mesh::BuiltBackMesh(Domain & a_dm)
{
	a_dm.splitdomain_tonode(m_level[0].m_box, m_level[0].blockpair);
	CheckBox(0);
	GiveAFlag("finish splitting domain to node!", 5);
#ifdef IMPORT_MESH
	ConstructBoxPoint_FromImportMesh();
	GiveAFlag("Finish ConstructBoxPoint_FromImportMesh!!!", 5);
#else
	int newboxstart = 0;
	int newboxend = 0;
	m_level[0].m_box.GlobalOrder(newboxstart, newboxend);
	ConstructBoxPoint(m_level[0].m_box, m_level[0].m_point, 0, newboxstart, newboxend);
#endif	
	GiveAFlag("finish ConstructBoxPoint!!!", 5);
	ConstructBoxFace(0);
	GiveAFlag("finish ConstructBoxFace!!!", 5);
   	m_level[0].RenewBoxFace(0);
  	ComputeLevelGeom(0);
  	GiveAFlag("Finish ComputeLevelGeom for the back mesh!!!", 5);
#ifdef PASSAGE_ANGLE
  	ComputeBlockPairTheta(0);
#endif  	
  	m_dm.FaceProperty(m_level[0].m_box, m_level[0].m_geom, m_level[0].dghost);
	ComputeDmgFaceAngle(0);
  	CheckFaceBoxIndex();
#ifndef IMPORT_MESH 
  	CheckBoxCenter();
  	CheckPointxyz();
#endif
  	CheckBoxPointIndex(0);
  	GiveAFlag("Finish CheckBoxPointIndex for the back mesh!!!", 5);
  	ReverseIndexCheck(0);
  	GiveAFlag("Finish ReverseIndexCheck for the back mesh!!!", 5);
  	m_level[0].CheckGhost(0);
  	GiveAFlag("Finish back ground mesh check!!!", 5);
	/*----Sort the box array on order to move the ghost box to the tail of the array*/
}

void Mesh::ComputeLevelGeom(const int & ilevel)
{
	m_level[ilevel].m_geom.setnum_nocopy(m_level[ilevel].m_box.realsize(),
																			 m_level[ilevel].m_box.realsize()-m_level[ilevel].m_box.size());
#ifndef IMPORT_MESH	
	for (int i = m_level[ilevel].m_geom.ps(); i < m_level[ilevel].m_geom.pe(); ++i)
	{
		ComptCellParams(m_level[ilevel].m_geom[i], ilevel, i);
	}
	for (int i = m_level[ilevel].m_geom.gps(); i < m_level[ilevel].m_geom.gpe(); ++i)
	{
		ComptCellParams(m_level[ilevel].m_geom[i], ilevel, i);
  }
#endif

#ifdef IMPORT_MESH
  PRINTFinLEVEL("level geom start %d end %d", 0, m_level[ilevel].m_geom.ps(), m_level[ilevel].m_geom.pe());
	for (int i = m_level[ilevel].m_geom.ps(); i < m_level[ilevel].m_geom.pe(); ++i)
	{
		ComptCellParams_Pts_Mesh(ilevel, i);
	}
	for (int i = m_level[ilevel].m_geom.gps(); i < m_level[ilevel].m_geom.gpe(); ++i)
	{
		ComptCellParams_Pts_Mesh(ilevel, i);
  }
	GiveAFlag("Finish compute the normal cell m_geom for the imported mesh!!!", 5);
	MPI_Barrier(MPI_COMM_WORLD);
#endif    
}

void Mesh::DataExchangeCells(const int & nlevel)
  {
    m_level[nlevel].nd.cell_extrac.resize(nodenum);
    vector<Direct_Mesg> local_direct_extrac;
    int rcn = pow(2, nlevel);
    int rcnz = level_power_ratio[nlevel];
    int bprange[3][2] = {{m_dm.dmbound[0][0]*rcn+dtln,m_dm.dmbound[0][1]*rcn-dtln-1},
                         {m_dm.dmbound[1][0]*rcn+dtln,m_dm.dmbound[1][1]*rcn-dtln-1},
                         {m_dm.dmbound[2][0]*rcnz+dtln,m_dm.dmbound[2][1]*rcnz-dtln-1}};
    for (int i = 0; i < nodenum; ++i)
    {
      m_level[nlevel].nd.cell_extrac[i].resize(0);
    }
    for (int i = m_level[nlevel].blockpair.ps(); i < m_level[nlevel].blockpair.pe(); ++i)
    {
    	int in0 = m_level[nlevel].blockpair[i].innode;
      Point & inxyz = m_level[nlevel].m_box[in0].lowpt; 
#if DIM == 3      			
      if (inxyz[0] < bprange[0][0] || inxyz[0] > bprange[0][1] ||
          inxyz[1] < bprange[1][0] || inxyz[1] > bprange[1][1] ||
          inxyz[2] < bprange[2][0] || inxyz[2] > bprange[2][1])
#elif DIM == 2
			if (inxyz[0] < bprange[0][0] || inxyz[0] > bprange[0][1] ||
          inxyz[1] < bprange[1][0] || inxyz[1] > bprange[1][1])
#endif    	
    	{	
    		if (m_level[nlevel].blockpair[i].outnode.node != node)
    		{
        	NewFlowMesg(nlevel, m_level[nlevel].blockpair[i]);
    		}
    		else
    		{
    			local_direct_extrac.push_back(Direct_Mesg(m_level[nlevel].blockpair[i].innode, 
    																												m_level[nlevel].blockpair[i].outnode.index));
#ifdef PASSAGE_ANGLE
    			local_direct_extrac.back().rotangle = m_level[nlevel].blockpair[i].theta; 			
#ifdef DEBUG
					int mytoid = m_level[nlevel].blockpair[i].outnode.index;
					int myfromid = m_level[nlevel].blockpair[i].innode;    	
					if (abs(m_level[nlevel].m_box[mytoid].iy()-m_level[nlevel].m_box[myfromid].iy()) != highpt.xy[1]*level_grid_ratio[nlevel][1])
					{
						printf("N%d need to transfer L%dB(%d,%d,%d) to Box(%d,%d,%d)!!!Error because they are not periodic pair!!!\n",
							node,nlevel,
							m_level[nlevel].m_box[myfromid].ix(),
							m_level[nlevel].m_box[myfromid].iy(),
							m_level[nlevel].m_box[myfromid].iz(),
							m_level[nlevel].m_box[mytoid].ix(),
							m_level[nlevel].m_box[mytoid].iy(),
							m_level[nlevel].m_box[mytoid].iz());
						MPI_Abort(MPI_COMM_WORLD, 662);
					}
#endif	
					// adm0.rotangle = ComptPointAngle_Rotate_X(m_level[nlevel].m_geom[adm0.toid].boxcenter, m_level[nlevel].m_geom[adm0.fromid].boxcenter);
					// if (adm0.rotangle < 0.0) adm0.rotangle = -PASSAGE_ANGLE;
					// else adm0.rotangle = PASSAGE_ANGLE;  
#endif					  			
    		}
    	}
    }
    //printf("N%dR%dL%d Direct_Mesg number is %d!!!\n",node, srank, nlevel, (int)m_level[nlevel].nd.direct_extrac.size());
    MPI_Barrier(MPI_COMM_WORLD);
    GiveAFlag("Finish create data exchange win for a layer...",5);
    m_level[nlevel].nd.CreateWin();
    m_level[nlevel].nd.Creat_direct_extract_array(local_direct_extrac);
  }

  void Mesh::DataExchangeCells_Wall()
  {
  	vector<Direct_Mesg> local_direct_extrac;
    nd_mdis.cell_extrac.resize(nodenum);
    rev_mdis.cell_extrac.resize(nodenum);
    rev_mdis.ibothernodes.resize(nodenum);
    int rcn = pow(2, cur_level_num-1);
    int rcnz = level_power_ratio[cur_level_num-1];
    int bprange[3][2] = {{m_dm.dmbound[0][0]*rcn+dtln,m_dm.dmbound[0][1]*rcn-dtln-1},
                         {m_dm.dmbound[1][0]*rcn+dtln,m_dm.dmbound[1][1]*rcn-dtln-1},
                         {m_dm.dmbound[2][0]*rcnz+dtln,m_dm.dmbound[2][1]*rcnz-dtln-1}};
    
    for (int i = 0; i < nodenum; ++i)
    {
      nd_mdis.cell_extrac[i].resize(0);
      rev_mdis.cell_extrac[i].resize(0);
      rev_mdis.ibothernodes[i].resize(0);
    }
    int lastlevel = cur_level_num-1;
    int bps = m_level[lastlevel].blockpair.ps();
    int bpe = m_level[lastlevel].blockpair.pe();
    for (int i = bps; i < bpe; ++i)
    {
    	int in0 = m_level[lastlevel].blockpair[i].innode;
    	Point & inxyz = m_level[lastlevel].m_box[in0].lowpt;
      NodePair & bkp = m_level[lastlevel].blockpair[i];
#if DIM == 3      			
      if (inxyz[0] < bprange[0][0] || inxyz[0] > bprange[0][1] ||
          inxyz[1] < bprange[1][0] || inxyz[1] > bprange[1][1] ||
          inxyz[2] < bprange[2][0] || inxyz[2] > bprange[2][1])
#elif DIM == 2
			if (inxyz[0] < bprange[0][0] || inxyz[0] > bprange[0][1] ||
          inxyz[1] < bprange[1][0] || inxyz[1] > bprange[1][1])
#endif          			        			
      { 
      	if (infectbox[in0] > -1)
      	{ 	
    			int n0 = m_dis[infectbox[in0]].mynode;
    			if (n0 != node && bkp.outnode.node != node)
      		{
#ifdef DEBUG      					
						rev_mdis.ibothernodes[n0].push_back(ExtractMesg_xyz(bkp.outnode.index, in0, inxyz, m_level[lastlevel].m_box[in0].pair.signdis,
									m_level[lastlevel].m_geom[in0].boxcenter));
#ifdef PASSAGE_ANGLE						
						rev_mdis.ibothernodes[n0].back().rotangle = bkp.theta;
#endif						
						printf("N%dR%d send to node %d the %d ib cell is: bkp %d remotecell %d localcell %d xyz (%d,%d,%d)!!!\n",
							node, srank, n0, (int)rev_mdis.ibothernodes[n0].size(), i, bkp.outnode.index, in0, inxyz[0], inxyz[1], inxyz[2]);
#else
						rev_mdis.ibothernodes[n0].push_back(ExtractMesg_xyz(bkp.outnode.index, in0, inxyz));
#ifdef PASSAGE_ANGLE						
						rev_mdis.ibothernodes[n0].back().rotangle = bkp.theta;
#endif						
#endif																			      					
      		}
      		//if (n0 == -2) the ib cell will be transferred to the ghost cell in the same node...  		
    			if (m_level[lastlevel].blockpair[i].outnode.node != node)
    			{	
      			nd_mdis.cell_extrac[bkp.outnode.node].push_back(ExtractMesg(bkp.outnode.index,bkp.innode));
#ifdef PASSAGE_ANGLE
						nd_mdis.cell_extrac[bkp.outnode.node].back().rotangle = m_level[lastlevel].blockpair[i].theta;
#endif						    			     					
      		}
      		else
      		{
      			local_direct_extrac.push_back(Direct_Mesg(m_level[lastlevel].blockpair[i].innode,
      																						 m_level[lastlevel].blockpair[i].outnode.index));
#ifdef PASSAGE_ANGLE  
						local_direct_extrac.back().rotangle = m_level[lastlevel].blockpair[i].theta;  			
#ifdef DEBUG
					int mytoid = m_level[lastlevel].blockpair[i].outnode.index;
					int myfromid = m_level[lastlevel].blockpair[i].innode;     	
						if (abs(m_level[lastlevel].m_box[mytoid].iy()-m_level[lastlevel].m_box[myfromid].iy()) != highpt.xy[1]*level_grid_ratio[lastlevel][1])
						{
							printf("N%d need to transfer L%dB(%d,%d,%d) to Box(%d,%d,%d)!!!Error because they are not periodic pair!!!\n",
								node,lastlevel,
								m_level[lastlevel].m_box[myfromid].ix(),
								m_level[lastlevel].m_box[myfromid].iy(),
								m_level[lastlevel].m_box[myfromid].iz(),
								m_level[lastlevel].m_box[mytoid].ix(),
								m_level[lastlevel].m_box[mytoid].iy(),
								m_level[lastlevel].m_box[mytoid].iz());
							MPI_Abort(MPI_COMM_WORLD, 662);
						}
#endif							
						// adm0.rotangle = ComptPointAngle_Rotate_X(m_level[lastlevel].m_geom[adm0.toid].boxcenter, m_level[lastlevel].m_geom[adm0.fromid].boxcenter); 
						// if (adm0.rotangle < 0.0) adm0.rotangle = -PASSAGE_ANGLE;
						// else adm0.rotangle = PASSAGE_ANGLE;
#endif	      			
      		}
      	}
    	}     
    }
    //printf("N%dR%d Wall Direct_Mesg number is %d!!!\n",node, srank, (int)m_level[lastlevel].nd.direct_extrac.size());
    MPI_Barrier(MPI_COMM_WORLD);
    GiveAFlag("Start TagIBCellOtherNode...",5);
    TagIBCellOtherNode();
    GiveAFlag("Start create data win for the wall layer...",5);
    nd_mdis.CreateWin();
    nd_mdis.Creat_direct_extract_array(local_direct_extrac);
    rev_mdis.CreateWin();
    rev_mdis.direct_extrac.setnum_nocopy(0,0);
    GiveAFlag("Finish create data win for the wall layer...",5);
  }
  void Mesh::TagIBCellOtherNode()
  {        
    rev_mdis.recvibcells.resize(0);
    rev_mdis.recvib_src_node.resize(0);
    BcastNewInfo_NodeRank(rev_mdis.ibothernodes, rev_mdis.recvibcells, rev_mdis.recvib_src_node, 7);
    GiveAFlag("Finish transfer the reverse IB cells!!!", 5);
    int recvnum = rev_mdis.recvibcells.size();
#ifdef DEBUG    
    if (recvnum > 0) printf("N%dR%d receive %d ib cells from other node!!!\n", node, srank,recvnum);
#endif    
    for (int i = 0; i < recvnum; ++i)
    {
      	ExtractMesg_xyz & i0 = rev_mdis.recvibcells[i];
      	int c0 = i0.remotecell;
      	/*----------Reason for this if--------------*/
      	/*-----The might be block pair due to periodicity-------*/
      	/*-----But this is not what we need---------------------*/
      	if (!(c0 > -1 && c0 < m_level[cur_level_num-1].m_box.realsize()))
      	{
      		printf("N%dR%d Recv %d ib cell innode is %d (%d,%d,%d) level box realsize is %d!!!\n",
      			node, srank, i, c0, i0.xyz[0], i0.xyz[1], i0.xyz[2],
      			m_level[cur_level_num-1].m_box.realsize());
      		MPI_Abort(MPI_COMM_WORLD, 481);
      	}
      	Assert((c0 > -1 && c0 < m_level[cur_level_num-1].m_box.realsize()), "The recv ib cell index error!!!", 476);
      	if (i0.xyz == m_level[cur_level_num-1].m_box[i0.remotecell].lowpt)
      	{
      		int if0 = infectbox[c0];
      		if (if0 < 0)
      		{
      			printf("N%d The recv ib cell (%d,%d,%d) signdis %f from N%d is not infected!!! Decteded in TagIBCellOtherNode!!!\n",
      				node, i0.xyz[0], i0.xyz[1], i0.xyz[2], 
      				rev_mdis.recvib_src_node[i]);
      			MPI_Abort(MPI_COMM_WORLD, 384);
      		}
      	  int n0 = rev_mdis.recvib_src_node[i];
      	  Assert((n0 > -1 && n0 < nodenum), "The recv src node error!!!", 490);  	
      		rev_mdis.cell_extrac[n0].push_back(ExtractMesg(i0.localcell, c0));
#ifdef PASSAGE_ANGLE      		
      		rev_mdis.cell_extrac[n0].back().rotangle = -i0.rotangle;
#endif      		
      		Assert((if0 > -1 && if0 < m_dis.size()), "The recv ib m_dis error!!!", 491);
      		m_dis[if0].mynode = rev_mdis.recvib_src_node[i];
      	}   
    }
    GiveAFlag("Finish tag the reverse IB cells!!!", 5);
#ifdef DEBUG    
    for (int i = 0; i < nodenum; ++i)
    {    	
    	if (rev_mdis.cell_extrac[i].size() > 0)
    	{
    		printf("N%dR%d Marked ghost IB cell number to Node %d: %d\n", node, srank, i, (int)rev_mdis.cell_extrac[i].size());
    	}
    }
#endif    
  }

void Mesh::DataExchange_alllevels(const int & time_step)
{
#ifdef SHOWTIME
	double starttime = MPI_Wtime();
#endif	
	for (int i = 0; i < cur_level_num; ++i)
	{
#ifdef TEMPORAL_REFINE
		if (time_step%marching_step[i] == marching_left_step[i])
		{
#endif
		DataExchange(m_level[i].nd, i);
#ifdef TEMPORAL_REFINE
		}
#endif		
	}
#ifdef SHOWTIME
	double endtime = MPI_Wtime();
	step_exchange_time += endtime - starttime;
#endif 	
}

  void Mesh::ComptMassFlow(double & massflow)
  {
  	massflow = 0.0;
  	double area0 = 0.0;
  	int bs = m_level[0].m_box.ps();
  	int be = m_level[0].m_box.pe();
  	for (int i = bs; i < be; ++i)
  	{
  		if (m_level[0].m_box[i].ix() == 0)
  		{
  			int f0 = m_level[0].m_box[i].faces[0][0];
  			int rnb = m_level[0].m_box[i].neib[0][1][1];
  			double v0 = 0.5*((m_level[0].m_data[i].u+m_level[0].m_data[rnb].u)*m_level[0].m_face[f0].keisa[0]+
  											 (m_level[0].m_data[i].v+m_level[0].m_data[rnb].v)*m_level[0].m_face[f0].keisa[1]+
  											 (m_level[0].m_data[i].w+m_level[0].m_data[rnb].w)*m_level[0].m_face[f0].keisa[2]);
  			massflow += m_level[0].m_data[i].roe*v0*m_level[0].m_face[f0].area;
  			area0+= m_level[0].m_face[f0].area;
  		}
  	}
  	MPI_Allreduce(MPI_IN_PLACE, &massflow, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  	MPI_Allreduce(MPI_IN_PLACE, &area0, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  	massflow *= roe_ref*u_ref*area_ref;
  	//if (ts < 100) printf("total mass flow is %f total area is %f!!!\n", massflow, area0);
  }

  void Mesh::ComptMassFlow_Outlet(double & massflow)
  {
  	massflow = 0.0;
  	double area0 = 0.0;
  	int bs = m_level[0].m_box.ps();
  	int be = m_level[0].m_box.pe();
  	for (int i = bs; i < be; ++i)
  	{
  		if (m_level[0].m_box[i].ix() == highpt.xy[0]-1)
  		{
  			int f0 = m_level[0].m_box[i].faces[0][1];
  			int rnb = m_level[0].m_box[i].neib[2][1][1];
  			double v0 = 0.5*((m_level[0].m_data[i].u+m_level[0].m_data[rnb].u)*m_level[0].m_face[f0].keisa[0]+
  											 (m_level[0].m_data[i].v+m_level[0].m_data[rnb].v)*m_level[0].m_face[f0].keisa[1]+
  											 (m_level[0].m_data[i].w+m_level[0].m_data[rnb].w)*m_level[0].m_face[f0].keisa[2]);
  			massflow += m_level[0].m_data[i].roe*v0*m_level[0].m_face[f0].area;
  			area0+= m_level[0].m_face[f0].area;
  		}
  	}
  	MPI_Allreduce(MPI_IN_PLACE, &massflow, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  	MPI_Allreduce(MPI_IN_PLACE, &area0, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  	massflow *= roe_ref*u_ref*area_ref;
  	//if (ts < 100) printf("total mass flow is %f total area is %f!!!\n", massflow, area0);
  }

  void Mesh::ComptAveInletVel(double & invel)
  {
  	invel = 0.0;
  	int cellnum = 0;
  	int bs = m_level[0].m_box.ps();
  	int be = m_level[0].m_box.pe();
  	for (int i = bs; i < be; ++i)
  	{
  		if (m_level[0].m_box[i].ix() == 0)
  		{
  			int f0 = m_level[0].m_box[i].faces[0][0];
  			double v0 = m_level[0].m_data[i].u*m_level[0].m_face[f0].keisa[0]+
  						m_level[0].m_data[i].v*m_level[0].m_face[f0].keisa[1]+
  						m_level[0].m_data[i].w*m_level[0].m_face[f0].keisa[2];
  			invel += v0;
  			++cellnum;
  		}
  	}
  	MPI_Allreduce(MPI_IN_PLACE, &invel, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  	MPI_Allreduce(MPI_IN_PLACE, &cellnum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  	invel /= double(cellnum);
  }

  void Mesh::GetPointVar(const int & ilevel, Pointxyz & mypt, const int & varindex, double & avar)
  {
  	int bs = m_level[ilevel].m_box.ps();
  	int be = m_level[ilevel].m_box.pe();
  	int var_num = 0;
  	avar = 0.0;
  	for (int i = bs; i < be; ++i)
  	{
  		if (abs(m_level[ilevel].m_geom[i].boxcenter[0]-mypt[0]) < dh[ilevel][0]*0.5 &&
  			abs(m_level[ilevel].m_geom[i].boxcenter[1]-mypt[1]) < dh[ilevel][1]*0.5 &&
  			abs(m_level[ilevel].m_geom[i].boxcenter[2]-mypt[2]) < dh[ilevel][2]*0.5)
  		{
  			avar = m_level[ilevel].m_data[i][varindex];
  			++var_num;
  		}
  	}
  	MPI_Allreduce(MPI_IN_PLACE, &var_num, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  	MPI_Allreduce(MPI_IN_PLACE, &avar, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  	if (var_num > 1)
  	{
  		printf("Fail to get point (%f,%f,%f) var %d because %d cells at the point!!!\n",
  			mypt[0], mypt[1], mypt[2], varindex, var_num);
  		MPI_Abort(MPI_COMM_WORLD, 496);
  	}
  }

  void Mesh::GetPointVar(const int & ilevel, Point & mypt, const int & varindex, double & avar)
  {
  	int bs = m_level[ilevel].m_box.ps();
  	int be = m_level[ilevel].m_box.pe();
  	int var_num = 0;
  	avar = 0.0;
  	for (int i = bs; i < be; ++i)
  	{
  		if (m_level[ilevel].m_box[i].lowpt == mypt) 
  		{
  			avar = m_level[ilevel].m_data[i][varindex];
  			++var_num;
  		}
  	}
  	MPI_Allreduce(MPI_IN_PLACE, &var_num, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  	MPI_Allreduce(MPI_IN_PLACE, &avar, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  	if (var_num > 1)
  	{
  		printf("Fail to get point (%d,%d,%d) var %d because %d cells at the point!!!\n",
  			mypt[0], mypt[1], mypt[2], varindex, var_num);
  		MPI_Abort(MPI_COMM_WORLD, 496);
  	}
  }

  void Mesh::DataExchange(Nodedata & nd, const int & ilevel)
  {
    Assert(nd.recvdestcell.size()==nd.recvfv.size(), "Error in data exchange!!!", 5);
    for (int n0 = 0; n0 < nodenum; ++n0)
    {
      	//int nfvsize = nd.num2ranks[n0];
      	int rps = nd.global_cell_extrac[n0].ps();
      	int rpe = nd.global_cell_extrac[n0].pe();
      	//printf("Node %d level %d exchange cells with %d is %d\n", node, ilevel, n0, nd.global_cell_extrac[n0].size());
      	for (int i = rps; i < rpe; ++i)
     		{
     			int i0 = i - rps;
     			int lc0 = nd.global_cell_extrac[n0][i].localcell;
       		nd.flowmesg_s[n0][i0] = m_level[ilevel].m_data[lc0];
#ifdef PASSAGE_ANGLE
					if (nd.global_cell_extrac[n0][i].rotangle != 0.0)
					{
						Pointxyz fromvel = Pointxyz(m_level[ilevel].m_data[lc0].u, m_level[ilevel].m_data[lc0].v, m_level[ilevel].m_data[lc0].w);
						fromvel.rotate_x(-nd.global_cell_extrac[n0][i].rotangle);
						nd.flowmesg_s[n0][i0].u = fromvel[0];
						nd.flowmesg_s[n0][i0].v = fromvel[1];
						nd.flowmesg_s[n0][i0].w = fromvel[2];
					}
#endif					     		
      	}
    }
    MPI_Barrier(nodecomm);
    MPI_Request send_req[nodenum];
    MPI_Status send_sta[nodenum];
    int tot_recv_num = nd.recvfv.size();
    //ShowAllRankData("tot_recv_num", tot_recv_num, 5);
    int have_recv_num = 0; 
    for (int i = 0; i < nodenum; ++i)
    {
      if (nd.numrecv[i] > 0)
      {
        //printf("<<<<<My rank is %d level %d recvnum from rank %d is %d\n", nrank, ilevel, i, nd.numrecv[i]);
        MPI_Irecv(&nd.recvfv[have_recv_num], nd.numrecv[i], MPI_FV, i, i, nodecomm, &send_req[i]);
        have_recv_num += nd.numrecv[i];
      }
    }
//     for (int i = 0; i < nodenum; ++i)
//     {
// //      printf("My rank is %d my number to rank %d is %d\n",  nrank, i, num2ranks[i]);
//       if (nd.num2ranks[i] > 0)
//       {
//         MPI_Isend(&nd.flowmesg_s[i][0], nd.num2ranks[i], MPI_FV, i, node, nodecomm, &send_req);
//         MPI_Wait(&send_req, &recv_stats);
//       }
//     }
    for (int i = 0; i < nodenum; ++i)
    {
      
      if (nd.num2ranks[i] > 0)
      {
      	//printf("My rank is %d level %d my number to rank %d is %d\n",  nrank, ilevel, i, nd.num2ranks[i]);
        MPI_Send(&nd.flowmesg_s[i][0], nd.num2ranks[i], MPI_FV, i, node, nodecomm);
      }
    }
    // int recvhand = 0;
    // MPI_Waitany(nodenum, send_req, &recvhand, send_sta);
    for (int i = 0; i < nodenum; ++i)
    {
    	if (nd.numrecv[i] > 0)
    	{
    		MPI_Wait(&send_req[i], &send_sta[i]);
    	}
    }
    int recvsize = nd.recvfv.size();
    for (int i = 0; i < recvsize; ++i)
    {
      int i0 = nd.recvdestcell[i];
      m_level[ilevel].m_data[i0] = nd.recvfv[i];
      // if (ilevel == max_mesh_level -1 && node == 1) printf("RECVFLAG<<N%dR%d Recv %d box (%d,%d,%d) receive a data!!!\n",node, srank, i, 
      // 	m_level[ilevel].m_box[i0].ix(),
      // 	m_level[ilevel].m_box[i0].iy(),
      // 	m_level[ilevel].m_box[i0].iz());
    }
    MPI_Barrier(nodecomm);
    int deps = nd.direct_extrac.ps();
    int depe = nd.direct_extrac.pe();
    for (int i = deps; i < depe; ++i)
    {
    	int fromid = nd.direct_extrac[i].fromid;
    	int toid = nd.direct_extrac[i].toid;
    	m_level[ilevel].m_data[toid] = m_level[ilevel].m_data[fromid];
#ifdef PASSAGE_ANGLE
			Pointxyz fromvel = Pointxyz(m_level[ilevel].m_data[fromid].u, m_level[ilevel].m_data[fromid].v, m_level[ilevel].m_data[fromid].w);
			fromvel.rotate_x(-nd.direct_extrac[i].rotangle);
			m_level[ilevel].m_data[toid].u = fromvel[0];
			m_level[ilevel].m_data[toid].v = fromvel[1];
			m_level[ilevel].m_data[toid].w = fromvel[2];
#endif			    	
    }
    MPI_Barrier(share_comm);
  }

  void Mesh::DataExchange_Pairinfo(Nodedata & nd, const int & ilevel)
  {
    Assert(nd.recvdestcell.size()==nd.recvfv.size(), "Error in data exchange!!!", 5);
    for (int n0 = 0; n0 < nodenum; ++n0)
    {
      	int rps = nd.global_cell_extrac[n0].ps();
      	int rpe = nd.global_cell_extrac[n0].pe();
      	for (int i = rps; i < rpe; ++i)
     	{
     		int i0 = i - rps;
       		nd.boxpair[n0][i0].body = m_level[ilevel].m_box[nd.global_cell_extrac[n0][i].localcell].pair.body;
       		nd.boxpair[n0][i0].patch = m_level[ilevel].m_box[nd.global_cell_extrac[n0][i].localcell].pair.patch;
      	}
    }
    MPI_Barrier(nodecomm);
    MPI_Request send_req[nodenum];
    MPI_Status send_sta[nodenum];
    int have_recv_num = 0; 
    for (int i = 0; i < nodenum; ++i)
    {
      if (nd.numrecv[i] > 0)
      {
        MPI_Irecv(&nd.recvpair[have_recv_num], nd.numrecv[i], MPI_PAIRINFO, i, i, nodecomm, &send_req[i]);
        have_recv_num += nd.numrecv[i];
      }
    }
    for (int i = 0; i < nodenum; ++i)
    {
      
      if (nd.num2ranks[i] > 0)
      {
        MPI_Send(&nd.boxpair[i][0], nd.num2ranks[i], MPI_PAIRINFO, i, node, nodecomm);
      }
    }
    for (int i = 0; i < nodenum; ++i)
    {
    	if (nd.numrecv[i] > 0)
    	{
    		MPI_Wait(&send_req[i], &send_sta[i]);
    	}
    }
  }

  void Mesh::ComputeDmgFaceAngle(const int & ilevel)
  {
    int dgs = m_level[ilevel].dghost.ps();
    int dge = m_level[ilevel].dghost.pe();
    for (int i = dgs; i < dge; ++i)
    {
    	Domainghost & adg = m_level[ilevel].dghost[i];
    	Boxloc & fo = m_dm.refboxface[adg.nbface];
    	if (fo.node > -1)
    	{ 		
      	if (!adg.nmvflag)
      	{
      	  int fi0 = m_level[ilevel].m_box[adg.refcell].faces[fo.node][fo.index];
      	  if (fi0 > -1)
      	  {
      	  	adg.nmv = m_level[ilevel].m_face[fi0].keisa;
      	  }
      	  else
      	  {
      	  	Assert(fi0 == -1, "The domain ghost ref face index error!!!", 503);
      	  	if (fo.node == 0)
      	  	{
      	  		ComptFaceNmv_ThreePts(adg.nmv,
      	  													m_level[ilevel].m_point[m_level[ilevel].m_box[adg.refcell].pts[fo.index][0][0]].xyz,
																		m_level[ilevel].m_point[m_level[ilevel].m_box[adg.refcell].pts[fo.index][1][0]].xyz,
																		m_level[ilevel].m_point[m_level[ilevel].m_box[adg.refcell].pts[fo.index][0][1]].xyz);
      	  	}
      	  	else if (fo.node == 1)
      	  	{
      	  		ComptFaceNmv_ThreePts(adg.nmv,
      	  													m_level[ilevel].m_point[m_level[ilevel].m_box[adg.refcell].pts[0][fo.index][0]].xyz,
																		m_level[ilevel].m_point[m_level[ilevel].m_box[adg.refcell].pts[0][fo.index][1]].xyz,
																		m_level[ilevel].m_point[m_level[ilevel].m_box[adg.refcell].pts[1][fo.index][0]].xyz);
      	  	}
      	  	else if (fo.node == 2)
      	  	{
      	  		ComptFaceNmv_ThreePts(adg.nmv,
      	  													m_level[ilevel].m_point[m_level[ilevel].m_box[adg.refcell].pts[0][0][fo.index]].xyz,
																		m_level[ilevel].m_point[m_level[ilevel].m_box[adg.refcell].pts[1][0][fo.index]].xyz,
																		m_level[ilevel].m_point[m_level[ilevel].m_box[adg.refcell].pts[0][1][fo.index]].xyz);
      	  	}
      	  }
      	  adg.nmvflag = true;
      	}
      }
    }
  }
#if DIM == 4
void Mesh::LocateHGCell(const int & ilevel, HGCell & hginfo, Pointxyz & patc)
{
	Pointxyz & hgxyz = hginfo.pt;
	int ibox = -1;
	hginfo.closedis = 1000.0;
	double distonearcell[3][3][3];
	int nearcell[3][3][3];
	double maxdis;
	while (ibox != hginfo.closecell)
	{
		ibox = hginfo.closecell;
		maxdis = -1.0;
		for (Point_iterator p(0,3); p.end(); ++p)
		{
			int an0 = m_level[ilevel].m_box[ibox].neib[p.i][p.j][p.k];
			if (an0 > -1)
			{
// #ifdef DEBUG
// 			if (an0 < 0)
// 			{
// 				printf("The hg point is (%f,%f,%f) atp is (%f,%f,%f) hgdis %f close distance is %f close box (%d,%d,%d)\n",
// 					hgxyz[0], hgxyz[1], hgxyz[2], 
// 					hginfo.attach_pt[0], hginfo.attach_pt[1], hginfo.attach_pt[2],
// 					hginfo.hgdis,
// 					hginfo.closedis,
// 					m_level[ilevel].m_box[ibox].ix(),
// 					m_level[ilevel].m_box[ibox].iy(),
// 					m_level[ilevel].m_box[ibox].iz());
// 				printf("Error 473: The hg interpolation point must be non-negative!!!\n");
// 				MPI_Abort(MPI_COMM_WORLD, 479);
// 			}			
// #endif			
				double adis = distobc(ilevel, an0, hgxyz);
				nearcell[p.i][p.j][p.k] = an0;
				distonearcell[p.i][p.j][p.k] = adis;
				maxdis = max(maxdis, adis);
				if (adis < hginfo.closedis)
				{
					hginfo.closecell = an0;
					hginfo.closedis = adis;
				}
			}
			else
			{
				nearcell[p.i][p.j][p.k] = -1;
			}
		}
	}
	for (int i = 0; i < 4; ++i)
	{
		hginfo.distointp[i] = maxdis;
		hginfo.intpcell[i] = -1;
	}
	for (Point_iterator p(0,3); p.end(); ++p)
	{
		if (nearcell[p.i][p.j][p.k] > -1)
		{
		for (int i = 0; i < 4; ++i)
		{
			if (distonearcell[p.i][p.j][p.k] < hginfo.distointp[i] && nearcell[p.i][p.j][p.k] != hginfo.intpcell[i])
			{
				for (int i0 = 3; i0 > i; --i0)
				{
					int s0 = i0-1;
					hginfo.distointp[i0] = hginfo.distointp[s0];
					hginfo.intpcell[i0] = hginfo.intpcell[s0];
				}
				hginfo.distointp[i] = distonearcell[p.i][p.j][p.k];
				hginfo.intpcell[i] = nearcell[p.i][p.j][p.k];
				goto NEXTNEIB;
			}
		}
		}
		NEXTNEIB:;
	}
#ifdef DEBUG
	// int flownum = 0;
	// int solidnum = 0;
	for (int i = 0; i < 4; ++i)
	{
		if (hginfo.intpcell[i] == -1)
		{
			printf("Hg point (%f,%f,%f) does not have the %d interpolation point!!!\n", hginfo.pt[0], hginfo.pt[1], hginfo.pt[2], i);
			MPI_Abort(MPI_COMM_WORLD, 1045);
		}
		if (!InTransferRange(ilevel, hginfo.intpcell[i]))
		{
			printf("N%d Hg point (%f,%f,%f) the %d interpolation box (%d,%d,%d)(%f,%f,%f) is not InTransferRange!!!\n", node,
				hginfo.pt[0], hginfo.pt[1], hginfo.pt[2], i,
				m_level[ilevel].m_box[hginfo.intpcell[i]].ix(),
				m_level[ilevel].m_box[hginfo.intpcell[i]].iy(),
				m_level[ilevel].m_box[hginfo.intpcell[i]].iz(),
				m_level[ilevel].m_geom[hginfo.intpcell[i]].boxcenter[0],
				m_level[ilevel].m_geom[hginfo.intpcell[i]].boxcenter[1],
				m_level[ilevel].m_geom[hginfo.intpcell[i]].boxcenter[2]);
			MPI_Abort(MPI_COMM_WORLD, 1045);
		}
		// if (m_level[ilevel].m_box[hginfo.intpcell[i]].solid) ++solidnum;
		// else ++flownum;
		for (int i0 = i+1; i0 < 4; ++i0)
		{
			if (hginfo.intpcell[i0] == hginfo.intpcell[i])
			{
				printf("Hg point (%f,%f,%f) %d %d interpolation points are (%f,%f,%f)\n",
					i, i0,
					hginfo.pt[0], hginfo.pt[1], hginfo.pt[2],
					m_level[ilevel].m_geom[hginfo.intpcell[i]].boxcenter[0],
					m_level[ilevel].m_geom[hginfo.intpcell[i]].boxcenter[1],
					m_level[ilevel].m_geom[hginfo.intpcell[i]].boxcenter[2]);
				MPI_Abort(MPI_COMM_WORLD, 537);
			}
		}
	}
	for (int i = 0; i < 4; ++i)
	{
		int i0 = hginfo.intpcell[i];
		if (!m_level[ilevel].m_box.isnormal(i0))
		{
			if (m_level[ilevel].m_box[i0].type != Blockghost)
			{
					printf("The hanging node interpolation point should be a normal cell or block ghost!!! the ghost point is (%f,%f,%f)\n", 
						m_level[ilevel].m_geom[i0].boxcenter[0],
						m_level[ilevel].m_geom[i0].boxcenter[1],
						m_level[ilevel].m_geom[i0].boxcenter[2]);
					MPI_Abort(MPI_COMM_WORLD, 569);
			}
		}
	}	
#endif
	hginfo.Comptintpcoef();
}

void Mesh::IntpHGgcellVar(HGCell & myhg, const int & lastlevel)
{
	myhg.fv.zero();
	for (int p0 = 0; p0 < 4; ++p0)
	{
		int an0 = myhg.intpcell[p0];
		Assert(an0 > -1, "Error 580 The interpolation point must be non-negative!!!", 580);
		FlowVariables & anfv = m_level[lastlevel].m_data[an0];
		myhg.fv.T += myhg.intpcoef[p0]*anfv.T;
		myhg.fv.p += myhg.intpcoef[p0]*anfv.p;
		myhg.fv.u += myhg.intpcoef[p0]*anfv.u;
		myhg.fv.v += myhg.intpcoef[p0]*anfv.v;
		myhg.fv.w += myhg.intpcoef[p0]*anfv.w;
		Get_roe_ideal_gas(myhg.fv);
	}
}

void Mesh::IntpcellNotSolid(bool & checkresult, const int & mylevel, HGCell & hgc)
{
	checkresult = true;
	for (int i = 0; i < 4; ++i)
	{		
		int an0 = hgc.intpcell[i];
		if (m_level[mylevel].m_box[an0].solid)
		{
			checkresult = false;
			break;
		}
	}
}

void Mesh::IntpcellNotSolid_NotInfect(bool & checkresult, const int & mylevel, HGCell & hgc)
{
	checkresult = true;
	for (int i = 0; i < 4; ++i)
	{		
		int an0 = hgc.intpcell[i];
		if (m_level[mylevel].m_box[an0].solid || infectbox[an0] > -1)
		{
			checkresult = false;
			break;
		}
	}
}
#endif
#if DIM < 4
void Mesh::LocateHGCell(const int & ilevel, HGCell & hginfo, Pointxyz & patc, int & hgnode, Pointxyz & dkeisa,
	int & rev_dir)
{
	hginfo.InitIntparray();
	rev_dir = 0;
	Pointxyz & hgxyz = hginfo.pt;
	int ibox = -1;
	hginfo.closedis = 1000.0;
	Pointxyz dxyz,newbcxyz;
	int cycle_num = 0;
	int init_box = hginfo.closecell;
	hgnode = -1;
	while (ibox != hginfo.closecell)
	{
		ibox = hginfo.closecell;
#ifndef PASSAGE_ANGLE	
		dxyz = hgxyz - m_level[ilevel].m_geom[ibox].boxcenter;
		PeriodicLength(dxyz);
#else
  	PeriodicAnnulaLength(hgxyz, m_level[ilevel].m_geom[ibox].boxcenter, newbcxyz);
  	dxyz = newbcxyz - m_level[ilevel].m_geom[ibox].boxcenter;
#endif		
		dkeisa[0] = dxyz.dot(m_level[ilevel].m_geom[ibox].keisa[0]);
		dkeisa[1] = dxyz.dot(m_level[ilevel].m_geom[ibox].keisa[1]);
		dkeisa[2] = dxyz.dot(m_level[ilevel].m_geom[ibox].keisa[2]);
		int dnxyz[3];
		for (int di = 0; di < 3; ++di)
		{
			if (dkeisa[di] < -1.00001)
			{
				dnxyz[di] = 0;
			}
			else if (dkeisa[di] > 1.00001)
			{
				dnxyz[di] = 2;
			}
			else
			{
				dnxyz[di] = 1;
			}
		}
		hginfo.closecell = m_level[ilevel].m_box[ibox].neib[dnxyz[0]][dnxyz[1]][dnxyz[2]];
		if (hginfo.closecell == -1)
		{
			printf("Hg close cell is -1! Hg point (%f,%f,%f) last close box (%d,%d,%d)(%f,%f,%f) hg dis %f!!!"
				   "last dkeisa is (%f,%f,%f) initial close box is (%d,%d,%d)!!!\n",
				hgxyz[0], hgxyz[1], hgxyz[2],
				m_level[ilevel].m_box[ibox].ix(),
				m_level[ilevel].m_box[ibox].iy(),
				m_level[ilevel].m_box[ibox].iz(),
				m_level[ilevel].m_geom[ibox].boxcenter[0],
				m_level[ilevel].m_geom[ibox].boxcenter[1],
				m_level[ilevel].m_geom[ibox].boxcenter[2],
				hginfo.hgdis,
				dkeisa[0], dkeisa[1], dkeisa[2],
				m_level[ilevel].m_box[init_box].ix(),
				m_level[ilevel].m_box[init_box].iy(),
				m_level[ilevel].m_box[init_box].iz());
			hginfo.closecell = ibox;
			if (!m_level[ilevel].m_box.isghost(ibox))
			{
				printf("The hg close box neib will be -1 but the box is not a ghost!!!\n");
				MPI_Abort(MPI_COMM_WORLD, 1054);
			}
			break;
		}
		if (m_level[ilevel].m_box[hginfo.closecell].type == Dmghost)
		{
#ifdef DEBUG			
			printf("Hg close cell is Dmghost! Hg point (%f,%f,%f) last close box (%d,%d,%d)(%f,%f,%f) hg dis %f!!!"
				   	"last dkeisa is (%f,%f,%f) initial close box is (%d,%d,%d)!!! cycle_num is %d\n",
					hgxyz[0], hgxyz[1], hgxyz[2],
					m_level[ilevel].m_box[ibox].ix(),
					m_level[ilevel].m_box[ibox].iy(),
					m_level[ilevel].m_box[ibox].iz(),
					m_level[ilevel].m_geom[ibox].boxcenter[0],
					m_level[ilevel].m_geom[ibox].boxcenter[1],
					m_level[ilevel].m_geom[ibox].boxcenter[2],
					hginfo.hgdis,
					dkeisa[0], dkeisa[1], dkeisa[2],
					m_level[ilevel].m_box[init_box].ix(),
					m_level[ilevel].m_box[init_box].iy(),
					m_level[ilevel].m_box[init_box].iz(),
					cycle_num);
#endif			
			return;
		}
		++cycle_num;
		if (cycle_num > 30)
		{
			// printf("Cycle 31 times to locate the hg cell!!! The last values: dkeisa(%f,%f,%f) Please check!!!\n",
			// 	dkeisa[0],dkeisa[1],dkeisa[2]);
			printf("Cycle %d times to locate the hg cell!!!\n"
				   "Hg point (%f,%f,%f) last close box (%d,%d,%d)(%f,%f,%f) hg dis %f!!!\n"
				   "last dkeisa is (%f,%f,%f) initial close box is (%d,%d,%d)(%f,%f,%f)!!!\n", cycle_num,
				hgxyz[0], hgxyz[1], hgxyz[2],
				m_level[ilevel].m_box[ibox].ix(),
				m_level[ilevel].m_box[ibox].iy(),
				m_level[ilevel].m_box[ibox].iz(),
				m_level[ilevel].m_geom[ibox].boxcenter[0],
				m_level[ilevel].m_geom[ibox].boxcenter[1],
				m_level[ilevel].m_geom[ibox].boxcenter[2],
				hginfo.hgdis,
				dkeisa[0], dkeisa[1], dkeisa[2],
				m_level[ilevel].m_box[init_box].ix(),
				m_level[ilevel].m_box[init_box].iy(),
				m_level[ilevel].m_box[init_box].iz(),
				m_level[ilevel].m_geom[init_box].boxcenter[0],
				m_level[ilevel].m_geom[init_box].boxcenter[1],
				m_level[ilevel].m_geom[init_box].boxcenter[2]);
			//printf("last dkeisa is (%f,%f,%f) initial close box is (%d,%d,%d)!!!\n", );
		}
		if (cycle_num > 50)
		{
			MPI_Abort(MPI_COMM_WORLD, 994);
		}
	}
	bool now_box_isnormal = m_level[ilevel].m_box.isnormal(hginfo.closecell);
	if (!now_box_isnormal)
	{		
		if (m_level[ilevel].m_box[hginfo.closecell].type != Blockghost)
		{
			printf("N%d The hg pt (%f,%f,%f) is close to a ghost which is not a Blockghost!!! Init close cell (%d,%d,%d) signdis %f hgdis %f final (%d,%d,%d)\n",
				node,
				hginfo.pt[0], hginfo.pt[1], hginfo.pt[2],
				m_level[ilevel].m_box[init_box].ix(),
				m_level[ilevel].m_box[init_box].iy(),
				m_level[ilevel].m_box[init_box].iz(),
				m_level[ilevel].m_box[init_box].pair.signdis,
				hginfo.hgdis,
				m_level[ilevel].m_box[hginfo.closecell].ix(),
				m_level[ilevel].m_box[hginfo.closecell].iy(),
				m_level[ilevel].m_box[hginfo.closecell].iz());
			MPI_Abort(MPI_COMM_WORLD, 1108);
		}	
		//Point neibdir;
		hgnode = FindNeibBlock(m_level[ilevel].m_box[hginfo.closecell], ilevel);
		//hgnode = m_dm.neibblocks[neibdir[0]][neibdir[1]][neibdir[2]];
		if (hgnode == node)
		{
			int remotecell = m_level[ilevel].m_box[hginfo.closecell].bkpid;
			if (!m_level[ilevel].m_box.isnormal(remotecell))
			{
				printf("N%d The hg closecell will be move from (%d,%d,%d) to the other side %d (%d,%d,%d)!!!"
							"But the new box is not a normal box!!! Init close cell (%d,%d,%d) signdis %f hgdis %f\n", node,
							m_level[ilevel].m_box[hginfo.closecell].ix(),
							m_level[ilevel].m_box[hginfo.closecell].iy(),
							m_level[ilevel].m_box[hginfo.closecell].iz(),
							remotecell,
							m_level[ilevel].m_box[remotecell].ix(),
							m_level[ilevel].m_box[remotecell].iy(),
							m_level[ilevel].m_box[remotecell].iz(),
							m_level[ilevel].m_box[init_box].ix(),
							m_level[ilevel].m_box[init_box].iy(),
							m_level[ilevel].m_box[init_box].iz(),
							m_level[ilevel].m_box[init_box].pair.signdis,
							hginfo.hgdis);
				MPI_Abort(MPI_COMM_WORLD, 1091);
			}
			if (m_level[ilevel].m_box[remotecell].iy() < m_level[ilevel].m_box[hginfo.closecell].iy())
			{
				rev_dir = 1;
			}
			else
			{
				rev_dir = -1;
			}
			hginfo.closecell = remotecell;
		}
		else if (hgnode < 0 || hgnode >= nodenum)
		{
			printf("When transfer the hg point to another node, old node is %d but neib is node %d!!!\n",
				node, hgnode);
			MPI_Abort(MPI_COMM_WORLD, 1042);
		}
		else
		{
			return;
		}
	}
	if (abs(dkeisa[0]) > 1.001 || abs(dkeisa[1]) > 1.001 ||abs(dkeisa[2]) > 1.001)
	{
		printf("Patch center is (%f,%f,%f) hgdis %f hgpt (%f,%f,%f) close box (%f,%f,%f) dkeisa is (%f,%f,%f)\n", 
			patc[0], patc[1], patc[2], hginfo.hgdis,
			hgxyz[0], hgxyz[1], hgxyz[2],
			m_level[ilevel].m_geom[hginfo.closecell].boxcenter[0],
			m_level[ilevel].m_geom[hginfo.closecell].boxcenter[1],
			m_level[ilevel].m_geom[hginfo.closecell].boxcenter[2],
			dkeisa[0], dkeisa[1], dkeisa[2]);
		for (Point_iterator p(0,3); p.end(); ++p)
		{
			int an0 = m_level[ilevel].m_box[ibox].neib[p.i][p.j][p.k];
			if (an0 > -1)
			{
				double adis = distobc(ilevel, an0, hgxyz);
				printf("Hg pt (%f,%f,%f) to close box neib (%d,%d,%d) (%f,%f,%f) distance is %f\n",
					hgxyz[0],hgxyz[1],hgxyz[2],p.i,p.j,p.k,
					m_level[ilevel].m_geom[an0].boxcenter[0],
					m_level[ilevel].m_geom[an0].boxcenter[1],
					m_level[ilevel].m_geom[an0].boxcenter[2],
					adis);
			}
		}
		MPI_Abort(MPI_COMM_WORLD, 1119);
	}
	Point hgdir;
	for (int i = 0; i < 3; ++i)
	{
		if (dkeisa[i] > 0.0)
		{
			hgdir[i] = 1;
		}
		else
		{
			hgdir[i] = 0;
		}
	}	
	for (Point_iterator p(0,2); p.end(); ++p)
	{
		int nx = p.i+hgdir[0];
		int ny = p.j+hgdir[1];
		int nz = p.k+hgdir[2];
		int an1 = m_level[ilevel].m_box[hginfo.closecell].neib[nx][ny][nz];
		hginfo.intpcell[p.i][p.j][p.k] = an1;
		// printf("N%dR%d The hanging node intp(%d,%d,%d) is box(%d,%d,%d) signdis %f infecttag %d!!!Hg point (%f,%f,%f) close cell (%d,%d,%d)patch (%f,%f,%f) hgdis %f !!!\n",
		// 		node,srank,
		// 		p.i,p.j,p.k,
		// 		m_level[ilevel].m_box[an1].ix(),
		// 		m_level[ilevel].m_box[an1].iy(),
		// 		m_level[ilevel].m_box[an1].iz(),
		// 		m_level[ilevel].m_box[an1].pair.signdis,
		// 		infectbox[an1],
		// 		hgxyz[0],hgxyz[1],hgxyz[2],
		// 		m_level[ilevel].m_box[hginfo.closecell].ix(),
		// 		m_level[ilevel].m_box[hginfo.closecell].iy(),
		// 		m_level[ilevel].m_box[hginfo.closecell].iz(),
		// 		patc[0],patc[1],patc[2],hginfo.hgdis);
		if (an1 < 0)
		{
			printf("WARNING::N%dR%d The hanging node is out of domain!!!Hg point (%f,%f,%f) close cell (%d,%d,%d)(%f,%f,%f) patch (%f,%f,%f) hgdis %f!!!\n",
				node,srank,
				hgxyz[0],hgxyz[1],hgxyz[2],
				m_level[ilevel].m_box[hginfo.closecell].ix(),
				m_level[ilevel].m_box[hginfo.closecell].iy(),
				m_level[ilevel].m_box[hginfo.closecell].iz(),
				m_level[ilevel].m_geom[hginfo.closecell].boxcenter[0],
				m_level[ilevel].m_geom[hginfo.closecell].boxcenter[1],
				m_level[ilevel].m_geom[hginfo.closecell].boxcenter[2],
				patc[0],patc[1],patc[2],hginfo.hgdis);
			printf("dkeisa (%f,%f,%f)!!!\n", dkeisa[0], dkeisa[1], dkeisa[2]);
			printf("The hg cell initial close box is (%d,%d,%d)(%f,%f,%f)\n",
				m_level[ilevel].m_box[init_box].ix(),
				m_level[ilevel].m_box[init_box].iy(),
				m_level[ilevel].m_box[init_box].iz(),
				m_level[ilevel].m_geom[init_box].boxcenter[0],
				m_level[ilevel].m_geom[init_box].boxcenter[1],
				m_level[ilevel].m_geom[init_box].boxcenter[2]);
			MPI_Abort(MPI_COMM_WORLD, 1182);
		}
		else if (!InTransferRange(ilevel, an1))
		{
			printf("N%d Hanging node interpolation box L%dBox(%d,%d,%d) not InTransferRange.  The initial close box is (%d,%d,%d) final (%d,%d,%d)!!!\n",
				node, ilevel,
				m_level[ilevel].m_box[an1].ix(),
				m_level[ilevel].m_box[an1].iy(),
				m_level[ilevel].m_box[an1].iz(),
				m_level[ilevel].m_box[init_box].ix(),
				m_level[ilevel].m_box[init_box].iy(),
				m_level[ilevel].m_box[init_box].iz(),
				m_level[ilevel].m_box[hginfo.closecell].ix(),
				m_level[ilevel].m_box[hginfo.closecell].iy(),
				m_level[ilevel].m_box[hginfo.closecell].iz());
			MPI_Abort(MPI_COMM_WORLD, 1197);
		}
	}
}

void Mesh::FindHGCellIntpCell(HGCell & hginfo, const int & ilevel, Pointxyz & dkeisa, int & hgnode)
{
	hgnode = -1;
	Point hgdir;
	for (int i = 0; i < 3; ++i)
	{
		if (dkeisa[i] > 0.0)
		{
			hgdir[i] = 1;
		}
		else
		{
			hgdir[i] = 0;
		}
	}	
	for (Point_iterator p(0,2); p.end(); ++p)
	{
		int nx = p.i+hgdir[0];
		int ny = p.j+hgdir[1];
		int nz = p.k+hgdir[2];
		int an1 = m_level[ilevel].m_box[hginfo.closecell].neib[nx][ny][nz];
		hginfo.intpcell[p.i][p.j][p.k] = an1;
		if (an1 < 0)
		{
			// printf("WARNING::N%dR%d The hanging node is out of domain!!!Hg point (%f,%f,%f) close cell (%d,%d,%d)(%f,%f,%f) hgdis %f!!!\n",
			// 	node,srank,
			// 	hginfo.pt[0],hginfo.pt[1],hginfo.pt[2],
			// 	m_level[ilevel].m_box[hginfo.closecell].ix(),
			// 	m_level[ilevel].m_box[hginfo.closecell].iy(),
			// 	m_level[ilevel].m_box[hginfo.closecell].iz(),
			// 	m_level[ilevel].m_geom[hginfo.closecell].boxcenter[0],
			// 	m_level[ilevel].m_geom[hginfo.closecell].boxcenter[1],
			// 	m_level[ilevel].m_geom[hginfo.closecell].boxcenter[2],
			// 	hginfo.hgdis);
			// printf("dkeisa (%f,%f,%f)!!!\n", dkeisa[0], dkeisa[1], dkeisa[2]);
			// MPI_Abort(MPI_COMM_WORLD, 1182);
			//Point neibdir;
			hgnode = FindNeibBlock(m_level[ilevel].m_box[hginfo.closecell], ilevel);
			//hgnode = m_dm.neibblocks[neibdir[0]][neibdir[1]][neibdir[2]];
			if (hgnode == node)
			{
				printf("WARNING::N%dR%d The hanging node is out of domain!!!Hg point (%f,%f,%f) close cell (%d,%d,%d)(%f,%f,%f) hgdis %f!!!\n",
					node,srank,
					hginfo.pt[0],hginfo.pt[1],hginfo.pt[2],
					m_level[ilevel].m_box[hginfo.closecell].ix(),
					m_level[ilevel].m_box[hginfo.closecell].iy(),
					m_level[ilevel].m_box[hginfo.closecell].iz(),
					m_level[ilevel].m_geom[hginfo.closecell].boxcenter[0],
					m_level[ilevel].m_geom[hginfo.closecell].boxcenter[1],
					m_level[ilevel].m_geom[hginfo.closecell].boxcenter[2],
					hginfo.hgdis);
				printf("dkeisa (%f,%f,%f)!!!\n", dkeisa[0], dkeisa[1], dkeisa[2]);
				printf("The out-domain intp cell can not be located in the same node\n");
				MPI_Abort(MPI_COMM_WORLD, 1263);
			}
			return;
		}
		// else if (!InTransferRange(ilevel, an1))
		// {
		// 	printf("N%d Hanging node interpolation box L%dBox(%d,%d,%d) not InTransferRange.  The hg close box is (%d,%d,%d)!!!\n",
		// 		node, ilevel,
		// 		m_level[ilevel].m_box[an1].ix(),
		// 		m_level[ilevel].m_box[an1].iy(),
		// 		m_level[ilevel].m_box[an1].iz(),
		// 		m_level[ilevel].m_box[hginfo.closecell].ix(),
		// 		m_level[ilevel].m_box[hginfo.closecell].iy(),
		// 		m_level[ilevel].m_box[hginfo.closecell].iz());
		// 	MPI_Abort(MPI_COMM_WORLD, 1197);
		// }
	}
}


void Mesh::ComptHgIntpCoef(HGCell & hginfo, Pointxyz & dkeisa, const int & ilevel)
{
	double ca[3] = {0.0, 1.0, 0.0};
	double cb[3] = {1.0, -1.0, 1.0};
	Point hgdir;
	Pointxyz dxyz, newbcxyz;
	for (int i = 0; i < 3; ++i)
	{
		if (dkeisa[i] > 0.0)
		{
			hgdir[i] = 1;
		}
		else
		{
			hgdir[i] = 0;
		}
	}	
	for (Point_iterator p(0,2); p.end(); ++p)
	{
		int nx = p.i+hgdir[0];
		int ny = p.j+hgdir[1];
		int nz = p.k+hgdir[2];
		int an1 = hginfo.intpcell[p.i][p.j][p.k];
#ifdef IDW		
#ifndef PASSAGE_ANGLE	
		dxyz = hginfo.pt - m_level[ilevel].m_geom[an1].boxcenter;
		PeriodicLength(dxyz);
#else
  	PeriodicAnnulaLength(hginfo.pt, m_level[ilevel].m_geom[an1].boxcenter, newbcxyz);
  	dxyz = newbcxyz - m_level[ilevel].m_geom[an1].boxcenter;
#endif		
		hginfo.distointp[p.i][p.j][p.k] = dxyz.length();
#else		
		hginfo.intpcoef[p.i][p.j][p.k] = (ca[nx]+cb[nx]*abs(dkeisa[0]))*(ca[ny]+cb[ny]*abs(dkeisa[1]))*(ca[nz]+cb[nz]*abs(dkeisa[2]));
#endif		
	}
#ifdef IDW	
	hginfo.Comptintpcoef();
#endif
}

// void Mesh::LocateHGCell(const int & ilevel, HGCell & hginfo, Pointxyz & patc)
// {
// 	Pointxyz & hgxyz = hginfo.pt;
// 	int ibox = -1;
// 	hginfo.closedis = 1000.0;
// 	while (ibox != hginfo.closecell)
// 	{
// 		ibox = hginfo.closecell;
// 		for (Point_iterator p(0,3); p.end(); ++p)
// 		{
// 			int an0 = m_level[ilevel].m_box[ibox].neib[p.i][p.j][p.k];
// 			if (an0 > -1)
// 			{
// 				double adis = distobc(ilevel, an0, hgxyz);
// 				if (adis < hginfo.closedis)
// 				{
// 					hginfo.closecell = an0;
// 					hginfo.closedis = adis;
// 				}
// 			}
// 		}
// 	}
// 	Pointxyz dxyz,newbcxyz;
// #ifndef PASSAGE_ANGLE	
// 	dxyz = hgxyz - m_level[ilevel].m_geom[hginfo.closecell].boxcenter;
// 	PeriodicLength(dxyz);
// #else
//   	PeriodicAnnulaLength(hgxyz, m_level[ilevel].m_geom[hginfo.closecell].boxcenter, newbcxyz);
//   	dxyz = newbcxyz - m_level[ilevel].m_geom[hginfo.closecell].boxcenter;
// #endif	
// 	Pointxyz dkeisa;
// 	dkeisa[0] = dxyz.dot(m_level[ilevel].m_geom[hginfo.closecell].keisa[0]);
// 	dkeisa[1] = dxyz.dot(m_level[ilevel].m_geom[hginfo.closecell].keisa[1]);
// 	dkeisa[2] = dxyz.dot(m_level[ilevel].m_geom[hginfo.closecell].keisa[2]);
// 	if (abs(dkeisa[0]) > 1.0 || abs(dkeisa[1]) > 1.0 ||abs(dkeisa[2]) > 1.0)
// 	{
// 		printf("Patch center is (%f,%f,%f) hgdis %f hgpt (%f,%f,%f) close box (%f,%f,%f) dkeisa is (%f,%f,%f)\n", 
// 			patc[0], patc[1], patc[2], hginfo.hgdis,
// 			hgxyz[0], hgxyz[1], hgxyz[2],
// 			m_level[ilevel].m_geom[hginfo.closecell].boxcenter[0],
// 			m_level[ilevel].m_geom[hginfo.closecell].boxcenter[1],
// 			m_level[ilevel].m_geom[hginfo.closecell].boxcenter[2],
// 			dkeisa[0], dkeisa[1], dkeisa[2]);
// 		for (Point_iterator p(0,3); p.end(); ++p)
// 		{
// 			int an0 = m_level[ilevel].m_box[ibox].neib[p.i][p.j][p.k];
// 			if (an0 > -1)
// 			{
// 				double adis = distobc(ilevel, an0, hgxyz);
// 				printf("Hg pt (%f,%f,%f) to close box neib (%d,%d,%d) (%f,%f,%f) distance is %f\n",
// 					hgxyz[0],hgxyz[1],hgxyz[2],p.i,p.j,p.k,
// 					m_level[ilevel].m_geom[an0].boxcenter[0],
// 					m_level[ilevel].m_geom[an0].boxcenter[1],
// 					m_level[ilevel].m_geom[an0].boxcenter[2],
// 					adis);
// 			}
// 		}
// 	}
// 	Point hgdir;
// 	for (int i = 0; i < 3; ++i)
// 	{
// 		if (dkeisa[i] > 0.0)
// 		{
// 			hgdir[i] = 1;
// 		}
// 		else
// 		{
// 			hgdir[i] = 0;
// 		}
// 	}
// 	double ca[3] = {1.0, 0.0, 1.0};
// 	double cb[3] = {-1.0, 1.0, -1.0};
// 	for (Point_iterator p(0,2); p.end(); ++p)
// 	{
// 		int nx = p.i+hgdir[0];
// 		int ny = p.j+hgdir[1];
// 		int nz = p.k+hgdir[2];
// 		int an1 = m_level[ilevel].m_box[hginfo.closecell].neib[nx][ny][nz];
// 		hginfo.intpcell[p.i][p.j][p.k] = an1;
// 		hginfo.intpcoef[p.i][p.j][p.k] = (ca[nx]+cb[nx]*abs(dkeisa[0]))*(ca[ny]+cb[ny]*abs(dkeisa[1]))*(ca[nz]+cb[nz]*abs(dkeisa[2]));
// 		if (hginfo.intpcell[p.i][p.j][p.k] < 0)
// 		{
// 			printf("The interpolation error 774!!!\n");
// 			MPI_Abort(MPI_COMM_WORLD, 775);
// 		}
// 	}
// }

void Mesh::IntpHGgcellVar(HGCell & myhg, const int & lastlevel)
{
	myhg.fv.zero();
	for (Point_iterator p(0,2); p.end(); ++p)
	{
		int an0 = myhg.intpcell[p.i][p.j][p.k];
		Assert(an0 > -1, "Error 580 The interpolation point must be non-negative!!!", 580);
		FlowVariables & anfv = m_level[lastlevel].m_data[an0];
		myhg.fv.T += myhg.intpcoef[p.i][p.j][p.k]*anfv.T;
		myhg.fv.p += myhg.intpcoef[p.i][p.j][p.k]*anfv.p;
		myhg.fv.u += myhg.intpcoef[p.i][p.j][p.k]*anfv.u;
		myhg.fv.v += myhg.intpcoef[p.i][p.j][p.k]*anfv.v;
		myhg.fv.w += myhg.intpcoef[p.i][p.j][p.k]*anfv.w;
		Get_roe_ideal_gas(myhg.fv);
	}
}

void Mesh::IntpcellNotSolid(bool & checkresult, const int & mylevel, HGCell & hgc)
{
	checkresult = true;
	for (Point_iterator p(0,2); p.end(); ++p)
	{		
		int an0 = hgc.intpcell[p.i][p.j][p.k];
		Assert(an0 > -1, "Check hg intp cell non-negative fail!!!", 1373);
		if (m_level[mylevel].m_box[an0].solid)
		{
			checkresult = false;
			break;
		}
	}
}

void Mesh::IntpcellNotSolid_NotInfect(bool & checkresult, const int & mylevel, HGCell & hgc)
{
	checkresult = true;
	for (Point_iterator p(0,2); p.end(); ++p)
	{		
		int an0 = hgc.intpcell[p.i][p.j][p.k];
		if (m_level[mylevel].m_box[an0].solid || infectbox[an0] > -1)
		{
			checkresult = false;
			break;
		}
	}
}
#endif

void Mesh::CheckHgIntpcells(const int & ilevel, HGCell & hginfo)
{
#if DIM == 4	
	for (int i = 0; i < 4; ++i)
	{
		int an0 = hginfo.intpcell[i];
#elif DIM < 4
	for (Point_iterator p(0,2);p.end(); ++p)
	{
		int an0 = hginfo.intpcell[p.i][p.j][p.k];
#endif		
		if (m_level[ilevel].m_box[an0].solid)
		{
			printf("The interpolation cell is a solid cell!!!\n");
			//printf("hg xyz after periodic (%f,%f,%f)\n", newbcxyz[0], newbcxyz[1], newbcxyz[2]);
			printf("Intp cell index(%d,%d,%d) close box index(%d,%d,%d) close bc(%f,%f,%f)close box signdis %f\n", 
				m_level[ilevel].m_box[an0].ix(),
				m_level[ilevel].m_box[an0].iy(),
				m_level[ilevel].m_box[an0].iz(),
				m_level[ilevel].m_box[hginfo.closecell].ix(),
				m_level[ilevel].m_box[hginfo.closecell].iy(),
				m_level[ilevel].m_box[hginfo.closecell].iz(),
				m_level[ilevel].m_geom[hginfo.closecell].boxcenter[0],
				m_level[ilevel].m_geom[hginfo.closecell].boxcenter[1],
				m_level[ilevel].m_geom[hginfo.closecell].boxcenter[2],
				m_level[ilevel].m_box[hginfo.closecell].pair.signdis);
			MPI_Abort(MPI_COMM_WORLD, 769);
		}
  	if (!InTransferRange(ilevel, an0))
		{
			printf("N%d Hanging node interpolation box L%dBox(%d,%d,%d) not InTransferRange.  The hg close box is (%d,%d,%d)!!!\n",
				node, ilevel,
				m_level[ilevel].m_box[an0].ix(),
				m_level[ilevel].m_box[an0].iy(),
				m_level[ilevel].m_box[an0].iz(),
				m_level[ilevel].m_box[hginfo.closecell].ix(),
				m_level[ilevel].m_box[hginfo.closecell].iy(),
				m_level[ilevel].m_box[hginfo.closecell].iz());
			MPI_Abort(MPI_COMM_WORLD, 1197);
		}
	}
}



void Mesh::CopyIBCellDistance(vector<double> & celldis)
{
	int mps = m_dis.ps();
	int mpe = m_dis.pe();
	int mdnum = mpe-mps;
	celldis.resize(mdnum);
	for (int i = mps; i < mpe; ++i)
	{
		int i0= i - mps;
		celldis[i0] = m_level[cur_level_num-1].m_box[m_dis[i].ci].pair.signdis;
	}
	MPI_Barrier(share_comm);
}

void Mesh::ComptDistancetoDomain()
{
	int bs, be;
	for (int i = 0; i < cur_level_num-1; ++i)
	{
		m_level[i].m_box.GlobalOrder(bs, be);
		for (int b0 = bs; b0 < be; ++b0)
		{
			if (m_level[i].m_box[b0].type == Normalcell)
			{
				ComptCelldistancetoDomain(i, b0);
				m_level[i].m_box[b0].pair.signdis = min(m_level[i].m_box[b0].pair.signdis, 
														m_level[i].m_box[b0].pair.distance_to_dm);
			}			
		}
	}
	int level0 = cur_level_num-1;
	if (infectbox.size() > 0)
	{	
		m_level[level0].m_box.GlobalOrder(bs, be);
		for (int i = bs; i < be; ++i)
		{
			if (infectbox[i] == -1)
			{
				if (m_level[level0].m_box[i].type == Normalcell)
				{
					ComptCelldistancetoDomain(level0, i);
					m_level[level0].m_box[i].pair.signdis = min(m_level[level0].m_box[i].pair.signdis,
																m_level[level0].m_box[i].pair.distance_to_dm);
				}
				
			}
		}
	}
	else
	{
		m_level[level0].m_box.GlobalOrder(bs, be);
		for (int i = bs; i < be; ++i)
		{
			if (m_level[level0].m_box[i].type == Normalcell)
			{
				ComptCelldistancetoDomain(level0, i);
				m_level[level0].m_box[i].pair.signdis = min(m_level[level0].m_box[i].pair.signdis,
															m_level[level0].m_box[i].pair.distance_to_dm);
			}		
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
}
int dircs[3][2] = {{1,2},{0,2},{0,1}};
void Mesh::ComptCelldistancetoSection(const int & ilevel, const int & ibox, const int & sec_d1, const int & sec_d2)
{
	double distocell;
	Point ixyz(m_level[ilevel].m_box[ibox].ix()/level_refine_ratio[ilevel],
			   m_level[ilevel].m_box[ibox].iy()/level_refine_ratio[ilevel],
			   m_level[ilevel].m_box[ibox].iz()/level_power_ratio[ilevel]);
	Assert(ixyz[0] > -1 && ixyz[0] < highpt.xy[0], "Error in ComptCelldistancetoSection direction x!!!", 1360);
	Assert(ixyz[1] > -1 && ixyz[1] < highpt.xy[1], "Error in ComptCelldistancetoSection direction y!!!", 1361);
	Assert(ixyz[2] > -1 && ixyz[2] < highpt.xy[2], "Error in ComptCelldistancetoSection direction z!!!", 1362);
	Pointxyz & cellbc = m_level[ilevel].m_geom[ibox].boxcenter;
	int cnx = ixyz[dircs[sec_d1][0]];
	int cny = ixyz[dircs[sec_d1][1]];
	distocell = distobc(ilevel, ibox, facecenter[sec_d1][sec_d2][cnx][cny].fc);
	bool to_cycle = true;
	while (to_cycle)
	{
		to_cycle = false;
		int cnx0 = cnx; int cny0 = cny;
		for (int i = -1; i < 2; ++i)
		{
			for (int j = -1; j < 2; ++j)
			{
				int nx = cnx0+i;
				int ny = cny0+j;
				if (nx > -1 && nx < highpt.xy[dircs[sec_d1][0]])
				{
					if (ny > -1 && ny < highpt.xy[dircs[sec_d1][1]])
					{
						double newdis = distobc(ilevel, ibox, facecenter[sec_d1][sec_d2][nx][ny].fc);
						if (newdis < distocell)
						{
							to_cycle = true;
							distocell = newdis;
							cnx = nx;
							cny = ny;
						}
					}
				}
			}
		}
	}
#ifndef PASSAGE_ANGLE  	
  	Pointxyz dxyz = cellbc - facecenter[sec_d1][sec_d2][cnx][cny].fc;
  	PeriodicLength(dxyz);
#else
  	Pointxyz newbcxyz;
  	PeriodicAnnulaLength(cellbc, facecenter[sec_d1][sec_d2][cnx][cny].fc, newbcxyz);
  	Pointxyz dxyz = newbcxyz - facecenter[sec_d1][sec_d2][cnx][cny].fc;
#endif
	double cell_dis_to_dm = dxyz.dot(facecenter[sec_d1][sec_d2][cnx][cny].fcnmv);	
	if (cell_dis_to_dm < 0.0)
	{
		// printf("L%dB%d(%d,%d,%d)(%f,%f,%f) distance to domain face(%d,%d)Element(%d,%d)center(%f,%f,%f)nmv(%f,%f,%f) is negative distance %f signdis %f!!!\n",
		// 	ilevel,ibox,m_level[ilevel].m_box[ibox].ix(),m_level[ilevel].m_box[ibox].iy(),m_level[ilevel].m_box[ibox].iz(),
		// 	cellbc[0], cellbc[1], cellbc[2],
		// 	sec_d1,sec_d2,cnx,cny,
		// 	facecenter[sec_d1][sec_d2][cnx][cny].fc[0],
		// 	facecenter[sec_d1][sec_d2][cnx][cny].fc[1],
		// 	facecenter[sec_d1][sec_d2][cnx][cny].fc[2],
		// 	facecenter[sec_d1][sec_d2][cnx][cny].fcnmv[0],
		// 	facecenter[sec_d1][sec_d2][cnx][cny].fcnmv[1],
		// 	facecenter[sec_d1][sec_d2][cnx][cny].fcnmv[2],
		// 	distocell, cell_dis_to_dm);
		//MPI_Abort(MPI_COMM_WORLD, 1407);
		cell_dis_to_dm = abs(cell_dis_to_dm);
	}
	m_level[ilevel].m_box[ibox].pair.distance_to_dm = min(m_level[ilevel].m_box[ibox].pair.distance_to_dm,
														  cell_dis_to_dm); 
}

void Mesh::ComptCelldistancetoDomain(const int & ilevel, const int & ibox)
{
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 2; ++j)
		{
			if (m_dm.boundaryiswall[i][j])
			{
				ComptCelldistancetoSection(ilevel, ibox, i, j);
			}	
		}
	}
}

void Mesh::LocatePointinComptDomain_AllDomain(Pointxyz & apt, const int & ilevel, int & init_box, Pointxyz & dkeisa)
{
	Pointxyz dxyz,newbcxyz;
	Point dnxyz;
	int dnse[3][2];
	int oldbox = -1;
	int init0 = init_box;
	int cycle_num = 0;	
	while (oldbox != init_box)
	{
		oldbox = init_box;
#ifndef PASSAGE_ANGLE  	
  	dxyz = apt - m_level[ilevel].m_geom[oldbox].boxcenter;
  	PeriodicLength(dxyz);
#else
  	PeriodicAnnulaLength(apt, m_level[ilevel].m_geom[oldbox].boxcenter, newbcxyz);
  	dxyz = newbcxyz - m_level[ilevel].m_geom[oldbox].boxcenter;
#endif			
		dkeisa[0] = dxyz.dot(m_level[ilevel].m_geom[oldbox].keisa[0]);
		dkeisa[1] = dxyz.dot(m_level[ilevel].m_geom[oldbox].keisa[1]);
		dkeisa[2] = dxyz.dot(m_level[ilevel].m_geom[oldbox].keisa[2]);
		for (int i = 0; i < 3; ++i)
		{
			// if (dkeisa[i] > 1.0) dnxyz[i] = 2;
			// else if (dkeisa[i] < -1.0) 	dnxyz[i] = 0;
			// else dnxyz[i] = 1;
			if (dkeisa[i] > 0)
			{
				dnse[i][0] = 1;
				dnse[i][1] = 3;
			}
			else
			{
				dnse[i][0] = 0;
				dnse[i][1] = 2;
			}
		}
		double mydis = 999.0;
		for (int dx = dnse[0][0]; dx < dnse[0][1]; ++dx)
		{
			for (int dy = dnse[1][0]; dy < dnse[1][1]; ++dy)
			{
				for (int dz = dnse[2][0]; dz < dnse[2][1]; ++dz)
				{
					int anb = m_level[ilevel].m_box[oldbox].neib[dx][dy][dz];
					if (anb > -1)
					{
						// if (m_level[ilevel].m_box[anb].type != Dmghost)
						// {
#ifndef PASSAGE_ANGLE  	
  						dxyz = apt - m_level[ilevel].m_geom[anb].boxcenter;
  						PeriodicLength(dxyz);
#else
  						PeriodicAnnulaLength(apt, m_level[ilevel].m_geom[anb].boxcenter, newbcxyz);
  						dxyz = newbcxyz - m_level[ilevel].m_geom[anb].boxcenter;
#endif
  						double anb_dis = dxyz.length();
  						if (anb_dis < mydis)
  						{
  							init_box = anb;
  							mydis = anb_dis;
  						}
  					// }
					}
				}
			}
		}
	}
	if (m_level[ilevel].m_box[init_box].type == Dmghost)
	{
		printf("The patch (%f,%f,%f) is close to a Dmghost (%f,%f,%f)!!!Please check!!!\n",
			apt[0], apt[1], apt[2],
			m_level[ilevel].m_geom[init_box].boxcenter[0],
			m_level[ilevel].m_geom[init_box].boxcenter[1],
			m_level[ilevel].m_geom[init_box].boxcenter[2]);
		MPI_Abort(MPI_COMM_WORLD, 1661);
	}
}

void Mesh::ReversePeriodicSide(const int & ilevel, int & init_box, int & rev_dir)
{
	rev_dir = 0;
	if (m_level[ilevel].m_box.isghost(init_box))
	{
		if (m_level[ilevel].m_box[init_box].bkpid > -1)
		{
			int remotecell = m_level[ilevel].m_box[init_box].bkpid;
			//Assert(m_level[ilevel].m_box.isnormal(remotecell), "The point close box has been switched to the other side but the box is not a normal cell!!!", 1689);
			if (m_level[ilevel].m_box[remotecell].iy() < m_level[ilevel].m_box[init_box].iy())
			{
				rev_dir = 1;
			}
			else
			{
				rev_dir = -1;
			}
			init_box = remotecell;
		}
		else
		{
			for (Point_iterator p(0,3); p.end(); ++p)
			{
				int an1 = m_level[ilevel].m_box[init_box].neib[p.i][p.j][p.k];
				if (an1 > -1)
				{
					if (m_level[ilevel].m_box[an1].bkpid > -1)
					{
						int abkpid = m_level[ilevel].m_box[an1].bkpid;
						int bkpnb = m_level[ilevel].m_box[abkpid].neib[2-p.i][2-p.j][2-p.k];
						if (bkpnb > -1)
						{
							if (m_level[ilevel].m_box[bkpnb].iy() < m_level[ilevel].m_box[init_box].iy())
							{
								rev_dir = 1;
							}
							else
							{
								rev_dir = -1;
							}
							init_box = bkpnb;
							break;
						}
					}
				}
			}
		}
	}
}

// void Mesh::LocatePointinComptDomain_NormalRange(Pointxyz & apt, const int & ilevel, int & init_box, int & apt_node)
// {
// 	Pointxyz dxyz;
// 	Pointxyz dkeisa;
// 	Point dnxyz;
// 	int oldbox = -1;
// 	int init0 = init_box;
// 	apt_node = node;
// 	while (oldbox != init_box)
// 	{
// 		oldbox = init_box;
// #ifndef PASSAGE_ANGLE  	
//   		Pointxyz dxyz = apt - m_level[ilevel].m_geom[oldbox].boxcenter;
//   		PeriodicLength(dxyz);
// #else
//   		Pointxyz newbcxyz;
//   		PeriodicAnnulaLength(apt, m_level[ilevel].m_geom[oldbox].boxcenter, newbcxyz);
//   		Pointxyz dxyz = newbcxyz - m_level[ilevel].m_geom[oldbox].boxcenter;
// #endif		
// 		dkeisa[0] = dxyz.dot(m_level[ilevel].m_geom[oldbox].keisa[0]);
// 		dkeisa[1] = dxyz.dot(m_level[ilevel].m_geom[oldbox].keisa[1]);
// 		dkeisa[2] = dxyz.dot(m_level[ilevel].m_geom[oldbox].keisa[2]);
// 		for (int i = 0; i < 3; ++i)
// 		{
// 			if (dkeisa[i] > 1.0) dnxyz[i] = 2;
// 			else if (dkeisa[i] < -1.0) 	dnxyz[i] = 0;
// 			else dnxyz[i] = 1;
// 		}
// 		init_box = m_level[ilevel].m_box[oldbox].neib[dnxyz[0]][dnxyz[1]][dnxyz[2]];
// 		if (m_level[ilevel].m_box.isghost(init_box))
// 		{
// 			if (m_level[ilevel].m_box[init_box].type != Blockghost)
// 			{
// 				printf("The point is close to a ghost cell which is not a blockghost. Maybe the boundary of the mesh level!!!\n");
// 				MPI_Abort(MPI_COMM_WORLD, 1568);
// 			}
// 			Point neibdir;
// 			FindNeibBlock(neibdir, m_level[ilevel].m_box[init_box], ilevel);
// 			apt_node = m_dm.neibblocks[neibdir[0]][neibdir[1]][neibdir[2]];
// 			if (apt_node < 0 || apt_node > nodenum - 1)
// 			{
// 				printf("Error in LocatePointinComptDomain!!! Level %d Input point (%f,%f,%f) last close box (%d,%d,%d)(%f,%f,%f) last dkeisa (%f,%f,%f)!!! "
// 					"init box is (%d,%d,%d)(%f,%f,%f) this node is %d and the new node is %d!!!\n",
// 				ilevel, apt[0],apt[1],apt[2],
// 				m_level[ilevel].m_box[oldbox].ix(),
// 				m_level[ilevel].m_box[oldbox].iy(),
// 				m_level[ilevel].m_box[oldbox].iz(),
// 				m_level[ilevel].m_geom[oldbox].boxcenter[0],
// 				m_level[ilevel].m_geom[oldbox].boxcenter[1],
// 				m_level[ilevel].m_geom[oldbox].boxcenter[2],
// 				dkeisa[0],dkeisa[1],dkeisa[2],
// 				m_level[ilevel].m_box[init0].ix(),
// 				m_level[ilevel].m_box[init0].iy(),
// 				m_level[ilevel].m_box[init0].iz(),
// 				m_level[ilevel].m_geom[init0].boxcenter[0],
// 				m_level[ilevel].m_geom[init0].boxcenter[1],
// 				m_level[ilevel].m_geom[init0].boxcenter[2], node, apt_node);
// 				MPI_Abort(MPI_COMM_WORLD, 1534);
// 			}
// 			else if (apt_node == node)
// 			{
// 				int remotecell = m_level[ilevel].m_box[init_box].bkpid;
// 				if (!m_level[ilevel].m_box.isnormal(remotecell))
// 				{
// 					printf("The close cell for the hg point has been reversed to the other side of the domain through bkp!!!"
// 						"But the remote cell index is %d, not a normal cell!!!\n"
// 						"the last close box is N%dL%d(%d,%d,%d)\n",
// 						node, ilevel,
// 						m_level[ilevel].m_box[init_box].ix(),
// 						m_level[ilevel].m_box[init_box].iy(),
// 						m_level[ilevel].m_box[init_box].iz());
// 					MPI_Abort(MPI_COMM_WORLD, 1670);
// 				}
// 				init_box = remotecell;
// 			}
// 			else
// 			{
// 				break;
// 			}		
// 		}
// 	}
// }

// void Mesh::GetCellDomainNormalVector(const int & ilevel, const int & ibox, Pointxyz & cellnmv)
// {
// 	if (m_level[ilevel].m_box[ibox].type != Dmghost)
// 	{
// 		printf("Try to get the domain normal vector of a box which is not domain ghost!!!\n");
// 		MPI_Abort(MPI_COMM_WORLD, 1493);
// 	}
// 	int rcdir[2] = {1,-1};
// 	int max_ratio[3] = {level_refine_ratio[ilevel], level_refine_ratio[ilevel], level_power_ratio[ilevel]};
// 	int capture = 0;
// 	Point dnxyz(1,1,1);
// 	for (int i = 0; i < 3; ++i)
// 	{
// 		for (int j = 0; j < 2; ++j)
// 		{
// 			bool assignface = false;
// 			if (m_dm.boundaryiswall[i][j])
// 			{
// 				if (j == 0)
// 				{
// 					if (m_level[ilevel].m_box[ibox].lowpt.xy[i] < lowpt.xy[i]*max_ratio[i])
// 					{
// 						++capture;
// 						assignface = true;
// 					}
// 				}
// 				else if (j == 1)
// 				{
// 					if (m_level[ilevel].m_box[ibox].lowpt.xy[i] > highpt.xy[i]*max_ratio[i]-1)
// 					{
// 						++capture;
// 						assignface = true;
// 					}
// 				}
// 				if (assignface)
// 				{
// 					// dnxyz[i] = rcdir[j] + 1;
// 					int an0 = m_level[ilevel].m_box[ibox].neib[dnxyz[0]][dnxyz[1]][dnxyz[2]];
// 					// int f0 = m_level[ilevel].m_box[an0].faces[i][j];
// 					int f0 = m_level[ilevel].m_box[ibox].faces[i][1-j];
// 					Assert((f0 > -1 && f0 < m_level[ilevel].m_face.realsize()), "The domain cell face index error!!!", 1507);
// #ifdef DEBUG
// 					if (f0 < 0 || f0 > m_level[ilevel].m_face.realsize()-1)
// 					{
// 						printf("Domain box (%d,%d,%d) neib (%d,%d,%d) is Box(%d,%d,%d) dir %d side %d face is %d box high %d low %d!!!\n",
// 							m_level[ilevel].m_box[ibox].ix(),
// 							m_level[ilevel].m_box[ibox].iy(),
// 							m_level[ilevel].m_box[ibox].iz(),
// 							dnxyz[0],dnxyz[1],dnxyz[2],
// 							m_level[ilevel].m_box[an0].ix(),
// 							m_level[ilevel].m_box[an0].iy(),
// 							m_level[ilevel].m_box[an0].iz(),
// 							i,j,f0,
// 							highpt.xy[i]*max_ratio[i]-1,
// 							lowpt.xy[i]*max_ratio[i]);
// 						MPI_Abort(MPI_COMM_WORLD, 1521);						
// 					}
// #endif					
// 					cellnmv = m_level[ilevel].m_face[f0].keisa*double(rcdir[j]);
// 				}				
// 			}
// 		}
// 	}
// 	if (capture > 1 || capture == 0)
// 	{
// 		printf("A domain cell at level %d (%d,%d,%d) captured time %d error!!!\n", 
// 			ilevel,
// 			m_level[ilevel].m_box[ibox].ix(),
// 			m_level[ilevel].m_box[ibox].iy(),
// 			m_level[ilevel].m_box[ibox].iz(),
// 			capture);
// 		MPI_Abort(MPI_COMM_WORLD, 1517);
// 	}
// }

void Mesh::GetCellDomainNormalVector(const int & ilevel, const int & ibox, Pointxyz & cellnmv)
{
	int rcdir[2] = {1,-1};
	int max_ratio[3] = {level_refine_ratio[ilevel], level_refine_ratio[ilevel], level_power_ratio[ilevel]};
	int capture = 0;
	Point dnxyz(1,1,1);
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 2; ++j)
		{
			bool assignface = false;
			if (m_dm.boundaryiswall[i][j])
			{
				if (j == 0)
				{
					if (m_level[ilevel].m_box[ibox].lowpt.xy[i] < (lowpt.xy[i]+highpt.xy[i])*max_ratio[i]/2)
					{
						++capture;
						assignface = true;
					}
				}
				else if (j == 1)
				{
					if (m_level[ilevel].m_box[ibox].lowpt.xy[i] > (lowpt.xy[i]+highpt.xy[i])*max_ratio[i]/2)
					{
						++capture;
						assignface = true;
					}
				}	
				if (assignface)
				{
					int f0 = m_level[ilevel].m_box[ibox].faces[i][j];
					Assert((f0 > -1 && f0 < m_level[ilevel].m_face.realsize()), "The domain cell face index error!!!", 1507);
#ifdef DEBUG
					if (f0 < 0 || f0 > m_level[ilevel].m_face.realsize()-1)
					{
						printf("IB (%d,%d,%d) dir %d side %d face is %d box high %d low %d!!!\n",
							m_level[ilevel].m_box[ibox].ix(),
							m_level[ilevel].m_box[ibox].iy(),
							m_level[ilevel].m_box[ibox].iz(),
							i,j,f0,
							highpt.xy[i]*max_ratio[i]-1,
							lowpt.xy[i]*max_ratio[i]);
						MPI_Abort(MPI_COMM_WORLD, 1521);						
					}
#endif					
					cellnmv = m_level[ilevel].m_face[f0].keisa*double(rcdir[j]);
				}
			}
		}
	}
	if (capture > 1 || capture == 0)
	{
		printf("An IB cell at level %d (%d,%d,%d) captured time %d error!!!\n", 
			ilevel,
			m_level[ilevel].m_box[ibox].ix(),
			m_level[ilevel].m_box[ibox].iy(),
			m_level[ilevel].m_box[ibox].iz(),
			capture);
		MPI_Abort(MPI_COMM_WORLD, 1517);
	}
}

int Mesh::FindNeibBlock(Box & lastbox, const int & ilevel)
{
	Point neibdir = Point(1,1,1);
	Point & bpt = lastbox.lowpt;
	for (int i = 0; i < 3; ++i)
	{
		//printf("Dir %d dmbound is (%d,%d) level_grid_ratio is %d box point is %d\n", i, m_dm.dmbound[i][0],m_dm.dmbound[i][1],level_grid_ratio[ilevel][i],bpt[i]);
		if (bpt[i] < m_dm.dmbound[i][0]*level_grid_ratio[ilevel][i])
		{
			neibdir[i] = 0;
		}
		else if (bpt[i] >= m_dm.dmbound[i][1]*level_grid_ratio[ilevel][i])
		{
			neibdir[i] = 2;
		}
	}
	int bn = m_dm.neibblocks[neibdir[0]][neibdir[1]][neibdir[2]];
	Assert((bn > -1 && bn < nodenum), "the neib block not in node range!!!", 1935);
	return m_dm.neibblocks[neibdir[0]][neibdir[1]][neibdir[2]];
}

void Mesh::MarkSpecialBox()
{	
    // for (int level0 = 0; level0 < cur_level_num; ++level0)
    // {
    // 	m_level[level0].InitPair();
    // 	for (int i = 0; i < nodenum; ++i)
    // 	{
    //   		if (i != node)
    //   		{
    //   			int rps = m_level[level0].nd.global_cell_extrac[i].ps();
    //     		int rpe = m_level[level0].nd.global_cell_extrac[i].pe();
    //     		for (int r0 = rps; r0 < rpe; ++r0)
    //     		{
    //      			int c0 = m_level[level0].nd.global_cell_extrac[i][r0].localcell;
    //      			if (m_level[level0].m_box[c0].ptype != none)
    //      			{
    //      				printf("The send box should be not marked!!!\n");
    //      				MPI_Abort(MPI_COMM_WORLD, 1658);
    //      			}
    //      			m_level[level0].m_box[c0].ptype = f_pro;
    //      		}
    //     	}
    //   }
    //   int recvnum = m_level[level0].nd.recvdestcell.size();
    //   for (int i = 0; i < recvnum; ++i)
    //   {
    //   	int r0 = m_level[level0].nd.recvdestcell[i];
    //   	if (m_level[level0].m_box[r0].ptype != none)
    //     {
    //      	printf("The recv box should be not marked!!!\n");
    //      	MPI_Abort(MPI_COMM_WORLD, 1658);
    //     }
    //     m_level[level0].m_box[r0].ptype = f_res2;
    //   }
    // }
	// if (infectbox.size() > 0)
	// {
	// 	m_level[cur_level_num-1].InitPair();
	// for (int i = 0; i < nodenum; ++i)
 //  {
 //    if (i != node)
 //    {
 //      int rps = rev_mdis.global_cell_extrac[i].ps();
 //      int rpe = rev_mdis.global_cell_extrac[i].pe();
 //      for (int r0 = rps; r0 < rpe; ++r0)
 //      {
 //        int c0 = rev_mdis.global_cell_extrac[i][r0].localcell;
 //        if (m_level[cur_level_num-1].m_box[c0].ptype != none)
 //        {
 //         	printf("The send box should be not marked!!!\n");
 //         	MPI_Abort(MPI_COMM_WORLD, 1658);
 //        }
 //        m_level[cur_level_num-1].m_box[c0].ptype = f_pro;
 //      }
 //    }
 //  }
 //  int recvnum = rev_mdis.recvdestcell.size();
 //  for (int i = 0; i < recvnum; ++i)
 //  {
 //    int r0 = rev_mdis.recvdestcell[i];
 //    if (m_level[cur_level_num-1].m_box[r0].ptype != none)
 //    {
 //      printf("The recv box should be not marked!!!\n");
 //      MPI_Abort(MPI_COMM_WORLD, 1658);
 //    }
 //    m_level[cur_level_num-1].m_box[r0].ptype = f_res2;
 //  }
	// }
 //  MPI_Barrier(share_comm);
}

		void Mesh::CounterSliceCellNumber(const int & dir, vector<int> & nodestart, vector<int> & nodeend)
		{
				vector<int> slicecellnum(highpt.xy[dir], 0);
				for (int i = max_level_num-1; i < max_level_num; ++i)
				{
						int bs = m_level[i].m_box.ps();
						int be = m_level[i].m_box.pe();
						for (int b0 = bs; b0 < be; ++b0)
						{
								int z0 = *(&m_level[i].m_box[b0].lowpt[0]+dir)/level_grid_ratio[i][2];
								++slicecellnum[z0];
						}
				}
				MPI_Allreduce(MPI_IN_PLACE, &slicecellnum[0], highpt.xy[dir], MPI_INT, MPI_SUM, MPI_COMM_WORLD);
				int totalcellnum;
				CountTotalNum(slicecellnum, totalcellnum);
				int num_every_node = int(double(totalcellnum)/double(nodenum))+1;
				int slice0 = 0;
				int current_node = 0;
				nodestart.resize(nodenum, 0);
				nodeend.resize(nodenum, 0);
				vector<int> nodecellnum(nodenum, 0);
				nodeend[nodenum-1] = highpt.xy[dir]-1;
				while (slice0 < highpt.xy[dir])
				{
						nodecellnum[current_node] += slicecellnum[slice0];
						if (nodecellnum[current_node] >= num_every_node)
						{
							nodeend[current_node] = slice0;
							nodestart[current_node+1] = slice0+1;
							++current_node; 
						}
						++slice0;
				}
				if (current_node != nodenum-1)
				{
						printf("The last node was not detected in the slice distribution!!!\n");
						MPI_Abort(MPI_COMM_WORLD, 1794);
				}
				if (nrank == 0)
				{
						printf("Domain split suggestion:\n");
						for (int i = 0; i < nodenum; ++i)
						{
								printf("Node %d Domain start: %d end %d slice number %d cell num %d\n",
									i, nodestart[i], nodeend[i], nodeend[i]-nodestart[i]+1, nodecellnum[i]);
						}
				}
				MPI_Barrier(MPI_COMM_WORLD);
		}

	void Mesh::CountInfectedBox()
  {
    int infnum = 0;
    int totinfnum = 0;
    int infectednum = m_dis.size();
    vector<int> infnum_procs(sprocs);
    vector<int> infstart_procs(sprocs);
    for (int i = infectbox.ps(); i < infectbox.pe(); ++i)
    {
      if (infectbox[i] == -2)
      {
        infnum += 1;
      }
    }
    MPI_Allgather(&infnum, 1, MPI_INT,
      &infnum_procs[0], 1, MPI_INT, share_comm);
    CountTotalNum(infnum_procs, totinfnum);
    ArrayProcsStart(infnum_procs, infstart_procs);
    infnum = 0;
    PRINTFinLEVELRANK0("Total new infected box number is %d...",cur_level_num-1,totinfnum);
#ifdef DEBUG    
    if (srank == 0) printf("N%d Newly infected box number is %d!!!\n", node, totinfnum);
#endif    
    for (int i = infectbox.ps(); i < infectbox.pe(); ++i)
    {
      if (infectbox[i] == -2)
      {
        infectbox[i] = infnum + infstart_procs[srank]+infectednum;
        ++infnum;
      }
    }
    vector<IBCell> newibcell(infnum);
    m_dis.Addnew(newibcell);
    m_dis.DirectlyReduceNew();
    MPI_Barrier(share_comm);
  }


  void Mesh::SynNormal_BlockDistance(const int & ilevel,
  	vector<int> & move_from_in_to_out,
		vector<int> & move_from_out_to_in)
  {
  	int bps, bpe;
  	m_level[ilevel].m_box.GlobalOrder(bps, bpe);
  	for (int i = bps; i < bpe; ++i)
		{
			if (m_level[ilevel].m_box[i].ix() == 564 &&
				m_level[ilevel].m_box[i].iy() == -2 &&
				m_level[ilevel].m_box[i].iz() == 85)
			{
				printf("N%d Box (403,-3,79) before syn close body %d patch %d signdis %f\n",node,
					m_level[ilevel].m_box[i].pair.body,
					m_level[ilevel].m_box[i].pair.patch,
					m_level[ilevel].m_box[i].pair.signdis);
			}
		}
  	vector<vector<Block_dis> > senddata(nodenum);
		for (int i = m_level[ilevel].blockpair.ps(); i < m_level[ilevel].blockpair.pe(); ++i)
		{
			int i0 = m_level[ilevel].blockpair[i].innode;
			BoxtoWall & btw = m_level[ilevel].m_box[i0].pair;
			senddata[m_level[ilevel].blockpair[i].outnode.node].push_back(
				Block_dis(m_level[ilevel].blockpair[i].outnode.index,
									btw.body,
									btw.patch,
									btw.distance,
									btw.signdis));
		}
		vector<Block_dis> recvdata;
		BcastNewInfo_Node(senddata, recvdata, 8);
		for (int i = 0; i < recvdata.size(); ++i)
		{
			int i0 = recvdata[i].remotebkp;
			BoxtoWall & btw = m_level[ilevel].m_box[i0].pair;
			btw.body = recvdata[i].b0;
			btw.patch = recvdata[i].p0;
			btw.distance = recvdata[i].distance;
			btw.signdis = recvdata[i].signdis;
			if (btw.signdis < 0.0)
			{
				if (!m_level[ilevel].m_box[i0].solid && m_level[ilevel].m_box[i0].type != Dmghost)
				{
					move_from_out_to_in.push_back(i0);
				}
				m_level[ilevel].m_box[i0].solid = true;
			}
			else
			{
				if (m_level[ilevel].m_box[i0].solid && m_level[ilevel].m_box[i0].type != Dmghost)
				{
					move_from_in_to_out.push_back(i0);
				}
				m_level[ilevel].m_box[i0].solid = false;
			}
		}
		MPI_Barrier(share_comm);
		for (int i = bps; i < bpe; ++i)
		{
			if (m_level[ilevel].m_box[i].ix() == 564 &&
				m_level[ilevel].m_box[i].iy() == -2 &&
				m_level[ilevel].m_box[i].iz() == 85)
			{
				printf("N%d Box (403,-3,79) after syn close body %d patch %d signdis %f\n",
					node,
					m_level[ilevel].m_box[i].pair.body,
					m_level[ilevel].m_box[i].pair.patch,
					m_level[ilevel].m_box[i].pair.signdis);
			}
		}
  }

  void Mesh::ComputeBlockPairTheta(const int & ilevel)
  {
#ifdef PASSAGE_ANGLE    	
  	int bkps = m_level[ilevel].blockpair.ps();
  	int bkpe = m_level[ilevel].blockpair.pe();
  	vector<vector<int> > localcell(nodenum);
  	vector<vector<int> > to_extract_cell(nodenum);
  	vector<int> from_node;
  	vector<int> i_need_to_extract_cell; 
  	for (int i = bkps; i < bkpe; ++i)
  	{
  		int n0 = m_level[ilevel].blockpair[i].outnode.node;
  		localcell[n0].push_back(i);
  		to_extract_cell[n0].push_back(m_level[ilevel].blockpair[i].outnode.index);
  	}
  	BcastNewInfo_NodeRank(to_extract_cell, i_need_to_extract_cell, from_node, 0);
  	vector<vector<BoxIndexBC> > data_tosend(nodenum);
  	int recvnum = from_node.size();
  	for (int i = 0; i < recvnum; ++i)
  	{
  		data_tosend[from_node[i]].push_back(BoxIndexBC(
  																				m_level[ilevel].m_box[i_need_to_extract_cell[i]].lowpt,
  																				m_level[ilevel].m_geom[i_need_to_extract_cell[i]].boxcenter));
  	}
  	vector<BoxIndexBC> recvbc;
  	vector<int> recvbc_node;
  	BcastNewInfo_NodeRank(data_tosend, recvbc, recvbc_node, 10);
  	int s0 = 0;
  	for (int i = 0; i < nodenum; ++i)
  	{
  		int extract_num = to_extract_cell[i].size();
  		for (int j = 0; j < extract_num; ++j)
  		{
  			int bkp0 = localcell[i][j];
  			BoxIndexBC & gloc = recvbc[s0];
  			int in0 = m_level[ilevel].blockpair[bkp0].innode;
  			Point & in0pt = m_level[ilevel].m_box[in0].lowpt;
  			if (abs(in0pt[1] - gloc.pt[1]) == highpt.xy[1]*level_grid_ratio[ilevel][1])
  			{
  				//angle of the second point minus the first point
  				m_level[ilevel].blockpair[bkp0].theta = ComptPointAngle_Rotate_X(gloc.bc, m_level[ilevel].m_geom[in0].boxcenter);
  				if (recvbc_node[s0] != node)
  				{
  					// printf("N%dL%d local cell (%d,%d,%d) remote cell (%d,%d,%d) from node %d rotangle %f\n",
  					// 	node, ilevel, in0pt[0], in0pt[1], in0pt[2],
  					// 	gloc.pt[0], gloc.pt[1], gloc.pt[2], 
  					// 	recvbc_node[s0],
  					// 	m_level[ilevel].blockpair[bkp0].theta);
  					if (abs(m_level[ilevel].blockpair[bkp0].theta) > PASSAGE_ANGLE + 0.01)
  					{
  						printf("ghost cell (%d,%d,%d)(%f,%f,%f) localcell (%d,%d,%d)(%f,%f,%f)\n",
  							gloc.pt[0], gloc.pt[1], gloc.pt[2], gloc.bc[0], gloc.bc[1], gloc.bc[2],
  							in0pt[0], in0pt[1], in0pt[2],
  							m_level[ilevel].m_geom[in0].boxcenter[0],
  							m_level[ilevel].m_geom[in0].boxcenter[1],
  							m_level[ilevel].m_geom[in0].boxcenter[2]);
  						MPI_Abort(MPI_COMM_WORLD, 2373);
  					}
  				}
  			}
  			else
  			{
  				m_level[ilevel].blockpair[bkp0].theta = 0.0;
  			}
  			++s0;
  		}
  	}
  	MPI_Barrier(share_comm);
#endif  	
  }