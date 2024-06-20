#include "AMR.H"

bool gridisnew = false;
bool adpflag = false;

	void AMR::PerformAdaptive(vector<Body> & abody)
	{
#ifdef DEBUG		
		CheckPointIndex();
#endif		
		SetRefineFlag();
		CheckLevelRefineFlag();
		GiveAFlag("Finish Check refine flag!!!", 5);
		for (int i = init_adjust_level; i < m_mesh.cur_level_num; ++i)
		{
			TagtheBox(i, abody);
		}
		GiveAFlag("Finish TagtheBox!!!", 5);
		//printf("Finish TagtheBox!!!\n");
		for (int i = init_adjust_level; i < m_mesh.cur_level_num; ++i)
		{
			m_level[i].AdjustMesh();
			GiveAFlag("Finish AdjustMesh!!!", 5);
			m_level[i].ConstructNewDomainGhost(m_mesh.FaceNormDir());	
			GiveAFlag("Finish ConstructNewDomainGhost!!!", 5);
#ifdef DEBUG
			m_level[i].CheckLevelProGhostType();
			GiveAFlag("Finish CheckLevelProGhostType!!!", 5);
#endif					
		}
		GiveAFlag("Finish attach bodies to the new level!!!", 5);
		for (int i = init_adjust_level; i < m_mesh.cur_level_num; ++i)
		{
			ManageLevelPoints(i);
			GiveAFlag("Finish ManageLevelPoints!!!", 5);
		}
		GiveAFlag("Finish AdjustMesh for each AMRLevel in AMR.H!!!", 5);
		/*--------------------------------------------------------------*/
	
		for (int i = init_adjust_level; i < m_mesh.cur_level_num; ++i)
		{
#ifndef IMPORT_MESH				
			ComptDatainNewBox(i, m_level[i].coarsegeom, m_level[i].finegeom);
#endif			
			if (flowdataexist) 
			{
				ComptDatainNewBox(i, m_level[i].coarsedata, m_level[i].finedata);
			}
		}
		GiveAFlag("Finish compute new data in AMR.H!!!", 5);
		/*--------------------------------------------------------------*/
		/*Compute new data must be perform before Assignnewloc Function*/
		/*Because after Assignnewloc, the old array will be released!!*/
		for (int i = init_adjust_level; i < m_mesh.cur_level_num; ++i)
		{
			/*This switch is to extrude the holes in the box array!!!*/
			m_mesh.m_level[i].TagRemovedBox();
//#ifndef LOCALBODY			
			if (i == m_mesh.cur_level_num-1 && m_mesh.infectbox.size() > 0)
			{
				Body::Check_Attachbox_willberemoved(abody, m_mesh);
			}
//#endif			
			m_mesh.m_level[i].RemoveConnections();
			GiveLevelFlag("Finish EndMeshHoleConnection!!!", i, 5);
			m_mesh.m_level[i].Assignnewloc(true);
			GiveAFlag("finish Assignnewloc!!!", 5);
			m_level[i].RenewGhostIndex();	
			if (i == m_mesh.cur_level_num-1 && m_mesh.infectbox.size() > 0) 
			{
				m_mesh.AdjustInfectedArray();
				Body::Renew_Attachbox(abody, m_mesh);
			}
			m_level[i].m_box.CompressArray();
			GiveAFlag("Finish Compress box array!!!", 5);			
			AllDataInitSwitch(i);
			GiveAFlag("Finish AllDataInitSwitch!!!", 5);
			PutDataBoxintoArray(i);
			GiveAFlag("Finish PutDataBoxintoArray!!!", 5);
			/*This switch is to move the ghost cells to the end of the array*/
			m_mesh.m_level[i].Assignnewloc(false);
			//m_mesh.m_level[i].RenewFacesBox();
			GiveAFlag("Finish Assignnewloc second!!!", 5);
			m_level[i].RenewGhostIndex();
			GiveAFlag("Finish RenewGhostIndex second!!!", 5);
			GiveAFlag("Finish AdjustInfectedArray second!!!", 5);
			if (i == m_mesh.cur_level_num-1 && m_mesh.infectbox.size() > 0) 
			{	
				m_mesh.AdjustInfectedArray();
				Body::Renew_Attachbox(abody, m_mesh);
			}
			GiveAFlag("Finish Renew_Attachbox!!! second", 5);
			m_level[i].ReorderBlockpair();
			GiveAFlag("Finish RenewGhostIndex!!!",5);
			SwitchNormalGhost(i);
			//m_mesh.CheckBox(i);
			/*---------------------------------------------------------*/
			//m_level[i].RenewGhostRefDmghost(m_mesh.FaceNormDir());
			m_mesh.m_level[i].RemoveOldPoints(i);
			m_mesh.ConstructBoxFace(i);
			m_mesh.ComputeDmgFaceAngle(i);
			//m_mesh.ConstructBoxFace(i);
			m_mesh.m_level[i].RemoveOldFace(i);
#ifdef IMPORT_MESH
			ComputeNewcellGeom(i);
#endif						
			/*---------------------------------------------------------*/
			GiveAFlag("Finish switch data between the ghost and the normal!!!", 5);
			// m_level[i].RenewBlockPairState();
			// GiveLevelFlag("Finish RenewBlockPairState!!!", i, 5);
			m_mesh.DataExchangeCells(i);		
		}
		if (m_mesh.infectbox.size() > 0) 
		{
			m_mesh.DataExchangeCells_Wall();
			for (int i = init_adjust_level; i < m_mesh.cur_level_num; ++i)
			{
				ComputeNewCellDistancetoDomian(i);
			}
		}
		GiveAFlag("Finish Assignnewloc for each AMRLevel in AMR.H!!!", 5);
		//CompressPointArray();							
		if (m_level.back().toderefine && m_mesh.Canbederefined())
		{
			CheckRemoveorNot(abody);
		}
		GiveAFlag("Start to compute the new level distance...", 5);
		if (morelevel)
		{
			int bodynum = abody.size();
			int l0 = m_mesh.cur_level_num - 2;
			for (int i = 0; i < bodynum; ++i)
			{
				abody[i].Attachbox_NewLevel(m_mesh, m_level[l0].m_tag, m_level[l0].allson);
			}
			GiveAFlag("Finish computing body to mesh distance...", 5);
			if (bodynum > 0)
			{
				Body::ComptDisFromBoxtoBody_Newlevel(m_mesh, m_level[l0].m_tag, m_level[l0].allson, abody);
			}
			GiveAFlag("Finish computing mesh to body distance...", 5);
		}
		/*The data array could not be resized right now!!!*/
//#ifdef DEBUG
		// PrintAllOnesideFace();
		CheckGhostLocation();
		GiveAFlag("Finish CheckGhostLocation...", 5);
		//m_mesh.CheckBoxOnlyOne();
		for (int i = init_adjust_level; i < m_mesh.cur_level_num; ++i)
		{
			CheckBoxNeib(i);
			m_mesh.CheckBox(i);
			CheckNeib27(i);
			m_mesh.ReverseIndexCheck(i);
			m_level[i].CheckResProPair_Debug();
			m_mesh.CheckBoxPointIndex(i);
			m_mesh.m_level[i].CheckGhost(i);
			m_level[i].SynGhostTag();
			m_level[i].CheckNumerofNormalghost();
		}
		CheckPointIndex();
		m_mesh.CheckFaceBoxIndex();
#ifndef IMPORT_MESH
		m_mesh.CheckBoxCenter();
		m_mesh.CheckPointxyz();
#endif		
		GiveAFlag("Q_A_Q All Processors Finish all check!!!\n",5);
		if (srank == 0)
		{
			printf("----------------------- Mesh Information------------------\n");
			//for (int i = init_adjust_level; i < m_mesh.cur_level_num; ++i)
			for (int i = m_mesh.cur_level_num-1; i < m_mesh.cur_level_num; ++i)
			{
				printf("------------------ Level %d -----------------\n", i);
				printf("N%d Point number %d Face number %d Cell number %d with ghost %d\n", 
					node,
					m_mesh.m_level[i].m_point.size(),
					m_mesh.m_level[i].m_face.size(),
					m_mesh.m_level[i].m_box.size(),
					m_mesh.m_level[i].m_box.realsize());
				// printf("Number of prolongation pair %d\n", int(m_level[i].f_pro_ghost.size()));
				// printf("Number of restriction P1 %d P2 %d\n", 
				// 	int(m_level[i].f_res_ghost[0].size()),
				// 	int(m_level[i].f_res_ghost[1].size()));
				// printf("Ghost number is %d\n", (int)m_level[i].m_box.ghost_index.size());
				// printf("Block pair number is %d\n", int(m_level[i].blockpair.size()));
			}
			printf("---------------------------------------------------------\n");		
		}
//#endif
		gridisnew = true;
		for (int i = init_adjust_level; i < m_mesh.cur_level_num; ++i)
		{
			m_level[i].MarkAllProResPair();
			m_level[i].CreatGhostProPair();
			m_level[i].MarkallPair();
			m_level[i].AssociateCelltoBlockPair(i);
		}
		GiveAFlag("Finish tag the pair in level_pro_ghost!!!", 5);
		//ConstructGhostFace();
		//FindModifyFace(m_mesh.mdyface);
		for (int i = init_adjust_level; i < m_mesh.cur_level_num; ++i)
		{
			m_level[i].cleartemparray();
			m_mesh.m_level[i].m_point.cleartemparray();
			m_mesh.m_level[i].m_face.cleartemparray();
			m_mesh.m_level[i].m_box.cleartemparray();	
		}
	}

	void AMR::ComptAdjointData(const int & flevel, const int & fbox, BoxCellGeom & newcoarsedata)
	{
		int dflevel = flevel-1;
		newcoarsedata.v = dh[dflevel][0]*dh[dflevel][1]*dh[dflevel][2];
		newcoarsedata.boxcenter.setvalue(0.0, 0.0, 0.0);
		for (Point_iterator p(0,2); p.end(); ++p)
		{
			int dfneib = m_level[flevel].m_box[fbox].neib[p.i+1][p.j+1][p.k+1];
			newcoarsedata.boxcenter += m_mesh.m_level[flevel].m_geom[dfneib].boxcenter/8.0;
		}
		newcoarsedata.boxcenter.get_theta();
		newcoarsedata.keisa[0] = Pointxyz(1.0/dh[dflevel][0], 0.0, 0.0);
		newcoarsedata.keisa[1] = Pointxyz(0.0, 1.0/dh[dflevel][1], 0.0);
		newcoarsedata.keisa[2] = Pointxyz(0.0, 0.0, 1.0/dh[dflevel][2]);
	}
	void AMR::ComptAdjointData(const int & rflevel, const int & rfbox0, Boxson<BoxCellGeom> & newfinedata)
	{
		int flevel = rflevel+1;
		Pointxyz & mombc = m_mesh.m_level[rflevel].m_geom[rfbox0].boxcenter;
		for (Point_iterator p(0,2); p.end(); ++p)
		{
			newfinedata.son[p.i][p.j][p.k].boxcenter[0] = mombc[0] + dh[flevel][0]*(double(p.i)-0.5);
			newfinedata.son[p.i][p.j][p.k].boxcenter[1] = mombc[1] + dh[flevel][1]*(double(p.j)-0.5); 
			newfinedata.son[p.i][p.j][p.k].boxcenter[2] = mombc[2] + dh[flevel][2]*(double(p.k)-0.5);
			newfinedata.son[p.i][p.j][p.k].boxcenter.get_theta();
			newfinedata.son[p.i][p.j][p.k].v = dh[flevel][0]*dh[flevel][1]*dh[flevel][2];
			newfinedata.son[p.i][p.j][p.k].keisa[0] = Pointxyz(1.0/dh[flevel][0], 0.0, 0.0);
			newfinedata.son[p.i][p.j][p.k].keisa[1] = Pointxyz(0.0, 1.0/dh[flevel][1], 0.0);
			newfinedata.son[p.i][p.j][p.k].keisa[2] = Pointxyz(0.0, 0.0, 1.0/dh[flevel][2]);
		}
	}
	void AMR::ComptAdjointData_2d(const int & flevel, const int & fbox, BoxCellGeom & newcoarsedata)
	{
		int dflevel = flevel-1;
		newcoarsedata.v = dh[dflevel][0]*dh[dflevel][1]*dh[dflevel][2];
		newcoarsedata.boxcenter.setvalue(0.0, 0.0, 0.0);
		for (Point_iterator_2d p(0,2); p.end(); ++p)
		{
			int dfneib = m_level[flevel].m_box[fbox].neib[p.i+1][p.j+1][p.k+1];
			newcoarsedata.boxcenter += m_mesh.m_level[flevel].m_geom[dfneib].boxcenter/4.0;
		}
		newcoarsedata.boxcenter.get_theta();
		newcoarsedata.keisa[0] = Pointxyz(1.0/dh[dflevel][0], 0.0, 0.0);
		newcoarsedata.keisa[1] = Pointxyz(0.0, 1.0/dh[dflevel][1], 0.0);
		newcoarsedata.keisa[2] = Pointxyz(0.0, 0.0, 1.0/dh[dflevel][2]);
	}
	void AMR::ComptAdjointData_2d(const int & rflevel, const int & rfbox0, Boxson<BoxCellGeom> & newfinedata)
	{
		int flevel = rflevel+1;
		Pointxyz & mombc = m_mesh.m_level[rflevel].m_geom[rfbox0].boxcenter;
		for (Point_iterator_2d p(0,2); p.end(); ++p)
		{
			newfinedata.son[p.i][p.j][p.k].boxcenter[0] = mombc[0] + dh[flevel][0]*(double(p.i)-0.5);
			newfinedata.son[p.i][p.j][p.k].boxcenter[1] = mombc[1] + dh[flevel][1]*(double(p.j)-0.5); 
			newfinedata.son[p.i][p.j][p.k].boxcenter[2] = mombc[2];
			newfinedata.son[p.i][p.j][p.k].boxcenter.get_theta();
			newfinedata.son[p.i][p.j][p.k].v = dh[flevel][0]*dh[flevel][1]*dh[flevel][2];
			newfinedata.son[p.i][p.j][p.k].keisa[0] = Pointxyz(1.0/dh[flevel][0], 0.0, 0.0);
			newfinedata.son[p.i][p.j][p.k].keisa[1] = Pointxyz(0.0, 1.0/dh[flevel][1], 0.0);
			newfinedata.son[p.i][p.j][p.k].keisa[2] = Pointxyz(0.0, 0.0, 1.0/dh[flevel][2]);
		}
	}

void AMR::ConstructGhostFace()
{
	vector<Face> newface;
	newface.reserve(1000);
	Point adp[3][2] = {{Point(0,1,1),Point(2,1,1)},{Point(1,0,1),Point(1,2,1)},{Point(1,1,0),Point(1,1,2)}};
	for (int i = init_adjust_level; i < m_mesh.cur_level_num-1; ++i)
	{
		newface.resize(0);
		int bs = m_level[i].g_pro.ps();
		int be = m_level[i].g_pro.pe();
		for (int b0 = bs; b0 < be; ++b0)
		{
			int g0 = m_level[i].g_pro[b0];
			int ci0 = m_level[i].level_pro_ghost[g0].ci;
			Assert(m_level[i].m_box[ci0].type == Blockghost && m_level[i].m_tag[ci0].tag > 0, "Not a ghost pro pair!!!", 344);
			for (int di = 0; di < DIM; ++di)
			{
				for (int f0 = 0; f0 < 2; ++f0)
				{
					int anb = m_level[i].m_box[ci0].neib[adp[di][f0][0]][adp[di][f0][1]][adp[di][f0][2]];
					//Assert(anb > -1, "The neib of pro pair main box must be non-negative!!!", 347);
					if (anb > -1)
					{
						if ((m_level[i].m_tag[anb].tag == -2 || m_level[i].m_tag[anb].tag == -3) && m_level[i].m_box[anb].type == Normalcell)
						{
							/*------------Construct a new face for the two cells----------*/
							int right0 = 1-f0;
							newface.push_back(Face());
							Face & aface = newface.back();
							aface[f0] = anb;
							aface[right0] = ci0;
							aface.fnv = di;
							m_mesh.m_level[i].NewFaceArea(aface, i, di, ci0, f0);
						}
					}
				}
			}
		}
		//PRINTFinLEVEL("g_pro start %d end %d!!!", i, bs, be);
		m_mesh.m_level[i].IncludeNewFace(newface);
		newface.resize(0);
	}
	GiveAFlag("Finish creat ghost face part 1!!!", 5);
	for (int i = init_adjust_level; i < m_mesh.cur_level_num-1; ++i)
	{
		int bs = m_level[i].level_pro_ghost.ps();
		int be = m_level[i].level_pro_ghost.pe();
		newface.resize(0);
		if (!m_level[i+1].twodflag)
		{
			for (int b0 = bs; b0 < be; ++b0)
			{
				int ci0 = m_level[i].level_pro_ghost[b0].ci;
				if (m_level[i].m_tag[ci0].tag == -2 || m_level[i].m_tag[ci0].tag == -3)
				{
					for (int di = 0; di < DIM; ++di)
					{
						for (int f0 = 0; f0 < 2; ++f0)
						{
							int anb = m_level[i].m_box[ci0].neib[adp[di][f0][0]][adp[di][f0][1]][adp[di][f0][2]];
							if (anb > -1)
							{
								if (m_level[i].m_box[anb].type == Normalcell)
								{
									if (m_level[i].m_tag[anb].tag > 0)
									{
										int p0 = m_level[i].m_tag[anb].tag - 1;
										CtoFPair & pro_pair = m_level[i].f_pro_ghost[p0];
										CtoFPair & lpg_pair = m_level[i].level_pro_ghost[b0];
										Assert(pro_pair.ci == anb, "The pro tag should be same with the box in creating ghost faces!!!", 271);
										Assert(lpg_pair.ci == ci0, "The ghost pro tag should be same with the box in creating ghost faces!!!", 271);
										/*------------Construct a new face for the two cells----------*/
										int right0 = 1-f0;
										if (di == 0)
										{
											for (int nx = 0; nx < 2; ++nx)
											{
												for (int ny = 0; ny < 2; ++ny)
												{
													newface.push_back(Face());
													Face & aface = newface.back();
													aface[f0] = pro_pair.fi.son[right0][nx][ny];
													aface[right0] = lpg_pair.fi.son[f0][nx][ny];
													Assert(m_level[i+1].m_box.isghost(aface[f0]), "The ghost face x side 0 should be a ghost box!!!", 285);
													Assert(m_level[i+1].m_box.isghost(aface[right0]), "The ghost face x side 0 should be a ghost box!!!", 286);
													Assert(m_level[i+1].m_box[aface[f0]].neib[2-2*f0][1][1] == aface[right0],
														"Two x 1 ghost boxes are not adjacent!!!", 288);
													Assert(m_level[i+1].m_box[aface[right0]].neib[2-2*right0][1][1] == aface[f0],
														"Two x 2 ghost boxes are not adjacent!!!", 290);
													Assert(m_level[i+1].m_box[aface[f0]].faces[di][right0]== -1, 
														"The ghost box x 1 should not have a ghost face!!!", 292);
													Assert(m_level[i+1].m_box[aface[right0]].faces[di][f0]== -1, 
														"The ghost box x 2 should not have a ghost face!!!", 294);
													aface.fnv = di;
													m_mesh.m_level[i+1].NewFaceArea(aface, i+1, di, aface[right0], f0);
												}
											}
										}
										else if (di == 1)
										{
											for (int nx = 0; nx < 2; ++nx)
											{
												for (int ny = 0; ny < 2; ++ny)
												{
													newface.push_back(Face());
													Face & aface = newface.back();
													aface[f0] = pro_pair.fi.son[nx][right0][ny];
													aface[right0] = lpg_pair.fi.son[nx][f0][ny];
													Assert(m_level[i+1].m_box.isghost(aface[f0]), "The ghost face y side 0 should be a ghost box!!!", 302);
													Assert(m_level[i+1].m_box.isghost(aface[right0]), "The ghost face y side 0 should be a ghost box!!!", 303);
													Assert(m_level[i+1].m_box[aface[f0]].neib[1][2-2*f0][1] == aface[right0],
														"Two y 1 ghost boxes are not adjacent!!!", 305);
													Assert(m_level[i+1].m_box[aface[right0]].neib[1][2-2*right0][1] == aface[f0],
														"Two y 2 ghost boxes are not adjacent!!!", 307);
													Assert(m_level[i+1].m_box[aface[f0]].faces[di][right0]== -1, 
														"The ghost box y 1 should not have a ghost face!!!", 317);
													Assert(m_level[i+1].m_box[aface[right0]].faces[di][f0]== -1, 
														"The ghost box y 2 should not have a ghost face!!!", 319);
													aface.fnv = di;
													m_mesh.m_level[i+1].NewFaceArea(aface, i+1, di, aface[right0], f0);
												}
											}
										}
										else if (di == 2)
										{
											for (int nx = 0; nx < 2; ++nx)
											{
												for (int ny = 0; ny < 2; ++ny)
												{
													newface.push_back(Face());
													Face & aface = newface.back();
													aface[f0] = pro_pair.fi.son[nx][ny][right0];
													aface[right0] = lpg_pair.fi.son[nx][ny][f0];
													Assert(m_level[i+1].m_box.isghost(aface[f0]), "The ghost face z side 0 should be a ghost box!!!", 319);
													Assert(m_level[i+1].m_box.isghost(aface[right0]), "The ghost face z side 0 should be a ghost box!!!", 320);
													Assert(m_level[i+1].m_box[aface[f0]].neib[1][1][2-2*f0] == aface[right0],
														"Two z 1 ghost boxes are not adjacent!!!", 326);
													Assert(m_level[i+1].m_box[aface[right0]].neib[1][1][2-2*right0] == aface[f0],
														"Two z 2 ghost boxes are not adjacent!!!", 328);
													Assert(m_level[i+1].m_box[aface[f0]].faces[di][right0]== -1, 
														"The ghost box z 1 should not have a ghost face!!!", 342);
													Assert(m_level[i+1].m_box[aface[right0]].faces[di][f0]== -1, 
														"The ghost box z 2 should not have a ghost face!!!", 344);
													aface.fnv = di;
													m_mesh.m_level[i+1].NewFaceArea(aface, i+1, di, aface[right0], f0);
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
			for (int b0 = bs; b0 < be; ++b0)
			{
				int ci0 = m_level[i].level_pro_ghost[b0].ci;
				if (m_level[i].m_tag[ci0].tag == -2 || m_level[i].m_tag[ci0].tag == -3)
				{
					for (int di = 0; di < DIM; ++di)
					{
						for (int f0 = 0; f0 < 2; ++f0)
						{
							int anb = m_level[i].m_box[ci0].neib[adp[di][f0][0]][adp[di][f0][1]][adp[di][f0][2]];
							if (anb > -1)
							{
								if (m_level[i].m_box[anb].type == Normalcell)
								{
									if (m_level[i].m_tag[anb].tag > 0)
									{
										// PRINTFinLEVEL("Ghost res pair box (%d,%d,%d) neib (%d,%d,%d) tag is %d\n", i, 
										// 	m_level[i].m_box[ci0].ix(),
										// 	m_level[i].m_box[ci0].iy(),
										// 	m_level[i].m_box[ci0].iz(),
										// 	adp[di][f0][0],adp[di][f0][1],adp[di][f0][2],
										// 	m_level[i].m_tag[anb].tag);
										int p0 = m_level[i].m_tag[anb].tag - 1;
										CtoFPair & pro_pair = m_level[i].f_pro_ghost[p0];
										CtoFPair & lpg_pair = m_level[i].level_pro_ghost[b0];
										Assert(pro_pair.ci == anb, "The pro tag should be same with the box in creating ghost faces!!!", 271);
										Assert(lpg_pair.ci == ci0, "The ghost pro tag should be same with the box in creating ghost faces!!!", 271);
										/*------------Construct a new face for the two cells----------*/
										int right0 = 1-f0;
										if (di == 0)
										{
											for (int nx = 0; nx < 2; ++nx)
											{
												newface.push_back(Face());
												Face & aface = newface.back();
												aface[f0] = pro_pair.fi.son[right0][nx][0];
												aface[right0] = lpg_pair.fi.son[f0][nx][0];
												Assert(m_level[i+1].m_box.isghost(aface[f0]), "The ghost face x side 0 should be a ghost box!!!", 285);
												Assert(m_level[i+1].m_box.isghost(aface[right0]), "The ghost face x side 0 should be a ghost box!!!", 286);
												Assert(m_level[i+1].m_box[aface[f0]].neib[2-2*f0][1][1] == aface[right0],
													"Two x 1 ghost boxes are not adjacent!!!", 288);
												Assert(m_level[i+1].m_box[aface[right0]].neib[2-2*right0][1][1] == aface[f0],
													"Two x 2 ghost boxes are not adjacent!!!", 290);
												Assert(m_level[i+1].m_box[aface[f0]].faces[di][right0]== -1, 
													"The ghost box x 1 should not have a ghost face!!!", 292);
												Assert(m_level[i+1].m_box[aface[right0]].faces[di][f0]== -1, 
													"The ghost box x 2 should not have a ghost face!!!", 294);
												aface.fnv = di;
												m_mesh.m_level[i+1].NewFaceArea(aface, i+1, di, aface[right0], f0);
											}
										}
										else if (di == 1)
										{
											for (int nx = 0; nx < 2; ++nx)
											{
												newface.push_back(Face());
												Face & aface = newface.back();
												aface[f0] = pro_pair.fi.son[nx][right0][0];
												aface[right0] = lpg_pair.fi.son[nx][f0][0];
												Assert(m_level[i+1].m_box.isghost(aface[f0]), "The ghost face y side 0 should be a ghost box!!!", 302);
												Assert(m_level[i+1].m_box.isghost(aface[right0]), "The ghost face y side 0 should be a ghost box!!!", 303);
												Assert(m_level[i+1].m_box[aface[f0]].neib[1][2-2*f0][1] == aface[right0],
													"Two y 1 ghost boxes are not adjacent!!!", 305);
												Assert(m_level[i+1].m_box[aface[right0]].neib[1][2-2*right0][1] == aface[f0],
													"Two y 2 ghost boxes are not adjacent!!!", 307);
												Assert(m_level[i+1].m_box[aface[f0]].faces[di][right0]== -1, 
													"The ghost box y 1 should not have a ghost face!!!", 317);
												Assert(m_level[i+1].m_box[aface[right0]].faces[di][f0]== -1,
													"The ghost box y 2 should not have a ghost face!!!", 319);
												aface.fnv = di;
												m_mesh.m_level[i+1].NewFaceArea(aface, i+1, di, aface[right0], f0);
											}
										}
										else if (di == 2)
										{
											for (int nx = 0; nx < 2; ++nx)
											{
												for (int ny = 0; ny < 2; ++ny)
												{
													newface.push_back(Face());
													Face & aface = newface.back();
													aface[f0] = pro_pair.fi.son[nx][ny][0];
													aface[right0] = lpg_pair.fi.son[nx][ny][0];
													Assert(m_level[i+1].m_box.isghost(aface[f0]), "The ghost face z side 0 should be a ghost box!!!", 319);
													Assert(m_level[i+1].m_box.isghost(aface[right0]), "The ghost face z side 0 should be a ghost box!!!", 320);
													Assert(m_level[i+1].m_box[aface[f0]].neib[1][1][2-2*f0] == aface[right0],
														"Two z 1 ghost boxes are not adjacent!!!", 326);
													Assert(m_level[i+1].m_box[aface[right0]].neib[1][1][2-2*right0] == aface[f0],
														"Two z 2 ghost boxes are not adjacent!!!", 328);
													Assert(m_level[i+1].m_box[aface[f0]].faces[di][right0]== -1, 
														"The ghost box z 1 should not have a ghost face!!!", 342);
													Assert(m_level[i+1].m_box[aface[right0]].faces[di][f0]== -1, 
														"The ghost box z 2 should not have a ghost face!!!", 344);
													aface.fnv = di;
													m_mesh.m_level[i+1].NewFaceArea(aface, i+1, di, aface[right0], f0);
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
		m_mesh.m_level[i+1].IncludeNewFace(newface);
	}
	MPI_Barrier(share_comm);
	GiveAFlag("Finish ConstructGhostFace!!!", 5);
}

	void AMR::ComputeNewCellDistancetoDomian(const int & ilevel)
	{
		int flevel = ilevel+1;
		int clevel = ilevel-1;
		if (flevel < m_mesh.cur_level_num)
		{
			if (!level_twod_flag[flevel])
			{
				for (int j = 0; j < m_level[ilevel].rfnum_procs[srank]; ++j)
				{
					//Assert((j > -1 && j < m_level[ilevel].rfbox.size()), "j error!!!", 548);
					int i0 = m_level[ilevel].rfbox[j];
					//Assert((i0 > -1 && i0 < m_level[ilevel].m_tag.realsize()), "i0 error!!!", 550);
					int tag0 = m_level[ilevel].m_tag[i0].tag;
					//Assert((tag0 > -1 && tag0 < m_level[ilevel].allson.realsize()), "tag0 error!!!", 551);
					int son0 = m_level[ilevel].allson[tag0].son[0][0][0];
					//Assert((son0 > -1 && son0 < m_level[flevel].m_box.realsize()), "son0 error!!!", 554);
					if (m_level[flevel].m_box[son0].type == Normalcell)
					{
						for (Point_iterator p(0,2); p.end(); ++p)
						{
							int sonid = m_level[ilevel].allson[tag0].son[p.i][p.j][p.k];
							//Assert((sonid > -1 && m_level[ilevel].m_box.realsize()), "sonid error!!!", 560);
							if (m_level[flevel].m_box[sonid].type == Normalcell)
							{
								m_mesh.ComptCelldistancetoDomain(flevel, sonid);
							}
						}
					}
				}
			}
			else
			{
				for (int j = 0; j < m_level[ilevel].rfnum_procs[srank]; ++j)
				{
					int i0 = m_level[ilevel].rfbox[j];
					int tag0 = m_level[ilevel].m_tag[i0].tag;
					int son0 = m_level[ilevel].allson[tag0].son[0][0][0];
					if (m_level[flevel].m_box[son0].type == Normalcell)
					{			
						for (Point_iterator_2d p(0,2); p.end(); ++p)
						{
							int sonid = m_level[ilevel].allson[tag0].son[p.i][p.j][p.k];
							m_mesh.ComptCelldistancetoDomain(flevel, sonid);
						}
					}
				}
			}
		}
		GiveAFlag("Finish ComputeNewCellDistancetoDomain for fine cells!!!", 5);
		for (int j = 0; j < m_level[ilevel].dfnum_procs[srank]; ++j)
		{
			int dfid = m_level[ilevel].derfbox[j];
			int momid = m_level[ilevel].m_tag[dfid].detag;
			if (m_level[clevel].m_box[momid].type == Normalcell)
			{	
				m_mesh.ComptCelldistancetoDomain(clevel, momid);
			}
		}
		MPI_Barrier(share_comm);
		GiveAFlag("Finish ComputeNewCellDistancetoDomain!!!", 5);
	}

