#include "AMR.H"
#include "GhostPair.H"
#include "NS_Solver.H"

void AMR::FindModifyFace(vector<DataArray<Facepair> > & mdyface0)
{
	Point f0[6] = {Point(0,1,1),Point(2,1,1),Point(1,0,1),Point(1,2,1),Point(1,1,0),Point(1,1,2)};
	int boxface[6][2] = {{0,0},{0,1},{1,0},{1,1},{2,0},{2,1}};
	int dim_dir[3] = {2, 4, 6};
	if (mdyface0.size() != m_mesh.cur_level_num)
	{
		mdyface0.resize(m_mesh.cur_level_num);
		for (int i = 0; i < m_mesh.cur_level_num; ++i)
		{
			mdyface0[i] = DataArray<Facepair>();
		}
	}
	for (int ilevel = 0; ilevel < m_mesh.cur_level_num-1; ++ilevel)
	{
		mdyface0[ilevel].setnum_nocopy(0,0);
		vector<Facepair> mdyface;
		mdyface.reserve(1000);
		int ps0 = m_level[ilevel].f_pro_ghost.ps();
		int pe0 = m_level[ilevel].f_pro_ghost.pe();
		int gps0 = m_level[ilevel].g_pro.ps();
		int	gpe0 = m_level[ilevel].g_pro.pe();
		if (!m_level[ilevel+1].twodflag)
		{
			for (int i = ps0; i < pe0; ++i)
			{
				int ci0 = m_level[ilevel].f_pro_ghost[i].ci;
				for (int d0 = 0; d0 < dim_dir[DIM-1]; ++d0)
				{
					int aneib0 = m_level[ilevel].m_box[ci0].neib[f0[d0][0]][f0[d0][1]][f0[d0][2]];
					Assert(aneib0 > -1, "The neib of the main cell of pro ghost is -1 in flux modification!!!", 16);
					// if ((m_level[ilevel].m_box.isghost(aneib0) && m_level[ilevel].m_box[aneib0].type == Normalcell) ||
					// 	(m_level[ilevel].m_box[aneib0].type == Blockghost && m_level[ilevel].m_tag[aneib0].tag == -2))
					if (m_level[ilevel].m_tag[aneib0].tag == -2)
					{
						mdyface.push_back(Facepair(m_level[ilevel].m_box[ci0].faces[boxface[d0][0]][boxface[d0][1]]));
						if (m_level[ilevel].m_box[ci0].faces[boxface[d0][0]][boxface[d0][1]] < 0)
						{
							PRINTFinLEVEL("Flux modification face error!!! f_pro_ghost %d main box %d face (%d,%d) is %d!!!", ilevel,
								i, ci0, boxface[d0][0], boxface[d0][1],
								m_level[ilevel].m_box[ci0].faces[boxface[d0][0]][boxface[d0][1]]);
							MPI_Abort(MPI_COMM_WORLD, 33);
						}
						if (boxface[d0][0] == 0)
						{
							for (int nx = 0; nx < 2; ++nx)
							{
								for (int ny = 0; ny < 2; ++ny)
								{
									int fi0 = m_level[ilevel].f_pro_ghost[i].fi.son[boxface[d0][1]][nx][ny];
									Assert(fi0 > -1, "The son of the pro ghost in flux modification error!!!", 29);
									mdyface.back().fface[nx][ny] = m_level[ilevel+1].m_box[fi0].faces[boxface[d0][0]][boxface[d0][1]];
									Assert(mdyface.back().fface[nx][ny] > -1, "The face of the son in flux modification error!!!",31);
								}
							}
						}				
						else if (boxface[d0][0] == 1)
						{
							for (int nx = 0; nx < 2; ++nx)
							{
								for (int ny = 0; ny < 2; ++ny)
								{
									int fi0 = m_level[ilevel].f_pro_ghost[i].fi.son[nx][boxface[d0][1]][ny];
									Assert(fi0 > -1, "The son of the pro ghost in flux modification error!!!", 42);
									mdyface.back().fface[nx][ny] = m_level[ilevel+1].m_box[fi0].faces[boxface[d0][0]][boxface[d0][1]];
									Assert(mdyface.back().fface[nx][ny] > -1, "The face of the son in flux modification error!!!",44);
								}
							}
						}				
						else if (boxface[d0][0] == 2)
						{
							for (int nx = 0; nx < 2; ++nx)
							{
								for (int ny = 0; ny < 2; ++ny)
								{
									int fi0 = m_level[ilevel].f_pro_ghost[i].fi.son[nx][ny][boxface[d0][1]];
									Assert(fi0 > -1, "The son of the pro ghost in flux modification error!!!", 55);
									mdyface.back().fface[nx][ny] = m_level[ilevel+1].m_box[fi0].faces[boxface[d0][0]][boxface[d0][1]];
									Assert(mdyface.back().fface[nx][ny] > -1, "The face of the son in flux modification error!!!",57);
								}
							}
						}					
					}
				}
			}
			for (int i = gps0; i < gpe0; ++i)
			{
				int g0 = m_level[ilevel].g_pro[i];
				int ci0 = m_level[ilevel].level_pro_ghost[g0].ci;
				for (int d0 = 0; d0 < dim_dir[DIM-1]; ++d0)
				{
					int aneib0 = m_level[ilevel].m_box[ci0].neib[f0[d0][0]][f0[d0][1]][f0[d0][2]];
					Assert(aneib0 > -1, "The neib of the main cell of pro ghost is -1 in flux modification!!!", 16);
					// if ((m_level[ilevel].m_box.isghost(aneib0) && m_level[ilevel].m_box[aneib0].type == Normalcell) ||
					// 	(m_level[ilevel].m_box[aneib0].type == Blockghost && m_level[ilevel].m_tag[aneib0].tag == -2))
					if (m_level[ilevel].m_tag[aneib0].tag == -2 && m_level[ilevel].m_box[aneib0].type == Normalcell)
					{
						mdyface.push_back(Facepair(m_level[ilevel].m_box[ci0].faces[boxface[d0][0]][boxface[d0][1]]));
						if (m_level[ilevel].m_box[ci0].faces[boxface[d0][0]][boxface[d0][1]] < 0)
						{
							PRINTFinLEVEL("Flux modification face error!!! f_pro_ghost %d main box %d face (%d,%d) is %d!!!", ilevel,
								i, ci0, boxface[d0][0], boxface[d0][1],
								m_level[ilevel].m_box[ci0].faces[boxface[d0][0]][boxface[d0][1]]);
							MPI_Abort(MPI_COMM_WORLD, 33);
						}
						if (boxface[d0][0] == 0)
						{
							for (int nx = 0; nx < 2; ++nx)
							{
								for (int ny = 0; ny < 2; ++ny)
								{
									int fi0 = m_level[ilevel].level_pro_ghost[g0].fi.son[boxface[d0][1]][nx][ny];
									Assert(fi0 > -1, "The son of the pro ghost in flux modification error!!!", 29);
									mdyface.back().fface[nx][ny] = m_level[ilevel+1].m_box[fi0].faces[boxface[d0][0]][boxface[d0][1]];
									Assert(mdyface.back().fface[nx][ny] > -1, "The face of the son in flux modification error!!!",31);
								}
							}
						}				
						else if (boxface[d0][0] == 1)
						{
							for (int nx = 0; nx < 2; ++nx)
							{
								for (int ny = 0; ny < 2; ++ny)
								{
									int fi0 = m_level[ilevel].level_pro_ghost[g0].fi.son[nx][boxface[d0][1]][ny];
									Assert(fi0 > -1, "The son of the pro ghost in flux modification error!!!", 42);
									mdyface.back().fface[nx][ny] = m_level[ilevel+1].m_box[fi0].faces[boxface[d0][0]][boxface[d0][1]];
									Assert(mdyface.back().fface[nx][ny] > -1, "The face of the son in flux modification error!!!",44);
								}
							}
						}				
						else if (boxface[d0][0] == 2)
						{
							for (int nx = 0; nx < 2; ++nx)
							{
								for (int ny = 0; ny < 2; ++ny)
								{
									int fi0 = m_level[ilevel].level_pro_ghost[g0].fi.son[nx][ny][boxface[d0][1]];
									Assert(fi0 > -1, "The son of the pro ghost in flux modification error!!!", 55);
									mdyface.back().fface[nx][ny] = m_level[ilevel+1].m_box[fi0].faces[boxface[d0][0]][boxface[d0][1]];
									Assert(mdyface.back().fface[nx][ny] > -1, "The face of the son in flux modification error!!!",57);
								}
							}
						}					
					}
				}
			}
		}
		else
		{
			for (int i = ps0; i < pe0; ++i)
			{
				int ci0 = m_level[ilevel].f_pro_ghost[i].ci;				
				for (int d0 = 0; d0 < dim_dir[DIM-1]; ++d0)
				{
					int aneib0 = m_level[ilevel].m_box[ci0].neib[f0[d0][0]][f0[d0][1]][f0[d0][2]];
					Assert(aneib0 > -1, "The neib of the main cell of pro ghost is -1 in flux modification!!!", 16);
					//if ((m_level[ilevel].m_box.isghost(aneib0) && m_level[ilevel].m_box[aneib0].type == Normalcell) ||
					//	(m_level[ilevel].m_box[aneib0].type == Blockghost && m_level[ilevel].m_tag[aneib0].tag == -2))
					if (m_level[ilevel].m_tag[aneib0].tag == -2)
					{
						mdyface.push_back(Facepair(m_level[ilevel].m_box[ci0].faces[boxface[d0][0]][boxface[d0][1]]));
						if (m_level[ilevel].m_box[ci0].faces[boxface[d0][0]][boxface[d0][1]] < 0)
						{
							PRINTFinLEVEL("Flux modification face error!!! f_pro_ghost %d main box %d face (%d,%d) is %d!!!", ilevel,
								i, ci0, boxface[d0][0], boxface[d0][1],
								m_level[ilevel].m_box[ci0].faces[boxface[d0][0]][boxface[d0][1]]);
							MPI_Abort(MPI_COMM_WORLD, 33);
						}
						if (boxface[d0][0] == 0)
						{
							for (int nx = 0; nx < 2; ++nx)
							{
								int fi0 = m_level[ilevel].f_pro_ghost[i].fi.son[boxface[d0][1]][nx][0];
								Assert(fi0 > -1, "The son of the pro ghost in flux modification error!!!", 29);
								mdyface.back().fface[nx][0] = m_level[ilevel+1].m_box[fi0].faces[boxface[d0][0]][boxface[d0][1]];
								mdyface.back().fface[nx][1] = m_level[ilevel+1].m_box[fi0].faces[boxface[d0][0]][boxface[d0][1]];
								Assert(mdyface.back().fface[nx][0] > -1, "The face of the son in flux modification error!!! x",31);
							}
						}					
						else if (boxface[d0][0] == 1)
						{
							for (int nx = 0; nx < 2; ++nx)
							{
								int fi0 = m_level[ilevel].f_pro_ghost[i].fi.son[nx][boxface[d0][1]][0];
								Assert(fi0 > -1, "The son of the pro ghost in flux modification error!!!", 42);
								mdyface.back().fface[nx][0] = m_level[ilevel+1].m_box[fi0].faces[boxface[d0][0]][boxface[d0][1]];
								mdyface.back().fface[nx][1] = m_level[ilevel+1].m_box[fi0].faces[boxface[d0][0]][boxface[d0][1]];
								Assert(mdyface.back().fface[nx][0] > -1, "The face of the son in flux modification error!!! y",44);
							}
						}						
						else if (boxface[d0][0] == 2)
						{
							for (int nx = 0; nx < 2; ++nx)
							{
								for (int ny = 0; ny < 2; ++ny)
								{
									int fi0 = m_level[ilevel].f_pro_ghost[i].fi.son[nx][ny][0];
									Assert(fi0 > -1, "The son of the pro ghost in flux modification error!!!", 55);
									mdyface.back().fface[nx][ny] = m_level[ilevel+1].m_box[fi0].faces[boxface[d0][0]][boxface[d0][1]];
									if (mdyface.back().fface[nx][ny] < 0)
									{
										printf("Main box neib (%d,%d,%d) is %d is ghost %d type %d tag %d!!!\n",
											f0[d0][0],f0[d0][1],f0[d0][2],aneib0,m_level[ilevel].m_box.isghost(aneib0),
											m_level[ilevel].m_box[aneib0].type,
											m_level[ilevel].m_tag[aneib0].tag);
									}
									Assert(mdyface.back().fface[nx][ny] > -1, "The face of the son in flux modification error!!! z",57);
								}
							}
						}					
					}
				}
			}
			for (int i = gps0; i < gpe0; ++i)
			{
				int g0 = m_level[ilevel].g_pro[i];
				int ci0 = m_level[ilevel].level_pro_ghost[g0].ci;
				for (int d0 = 0; d0 < dim_dir[DIM-1]; ++d0)
				{
					int aneib0 = m_level[ilevel].m_box[ci0].neib[f0[d0][0]][f0[d0][1]][f0[d0][2]];
					Assert(aneib0 > -1, "The neib of the main cell of pro ghost is -1 in flux modification!!!", 16);
					//if ((m_level[ilevel].m_box.isghost(aneib0) && m_level[ilevel].m_box[aneib0].type == Normalcell) ||
					//	(m_level[ilevel].m_box[aneib0].type == Blockghost && m_level[ilevel].m_tag[aneib0].tag == -2))
					if (m_level[ilevel].m_tag[aneib0].tag == -2 && m_level[ilevel].m_box[aneib0].type == Normalcell)
					{
						mdyface.push_back(Facepair(m_level[ilevel].m_box[ci0].faces[boxface[d0][0]][boxface[d0][1]]));
						if (m_level[ilevel].m_box[ci0].faces[boxface[d0][0]][boxface[d0][1]] < 0)
						{
							PRINTFinLEVEL("Flux modification face error!!! f_pro_ghost %d main box %d face (%d,%d) is %d!!!", ilevel,
								i, ci0, boxface[d0][0], boxface[d0][1],
								m_level[ilevel].m_box[ci0].faces[boxface[d0][0]][boxface[d0][1]]);
							MPI_Abort(MPI_COMM_WORLD, 33);
						}
						if (boxface[d0][0] == 0)
						{
							for (int nx = 0; nx < 2; ++nx)
							{
								int fi0 = m_level[ilevel].level_pro_ghost[g0].fi.son[boxface[d0][1]][nx][0];
								Assert(fi0 > -1, "The son of the pro ghost in flux modification error!!!", 29);
								mdyface.back().fface[nx][0] = m_level[ilevel+1].m_box[fi0].faces[boxface[d0][0]][boxface[d0][1]];
								mdyface.back().fface[nx][1] = m_level[ilevel+1].m_box[fi0].faces[boxface[d0][0]][boxface[d0][1]];
								Assert(mdyface.back().fface[nx][0] > -1, "The face of the son in flux modification error!!! x 0",31);
							}
						}					
						else if (boxface[d0][0] == 1)
						{
							for (int nx = 0; nx < 2; ++nx)
							{
								int fi0 = m_level[ilevel].level_pro_ghost[g0].fi.son[nx][boxface[d0][1]][0];
								Assert(fi0 > -1, "The son of the pro ghost in flux modification error!!!", 42);
								mdyface.back().fface[nx][0] = m_level[ilevel+1].m_box[fi0].faces[boxface[d0][0]][boxface[d0][1]];
								mdyface.back().fface[nx][1] = m_level[ilevel+1].m_box[fi0].faces[boxface[d0][0]][boxface[d0][1]];
								Assert(mdyface.back().fface[nx][0] > -1, "The face of the son in flux modification error!!! y 0 ",44);
							}
						}						
						else if (boxface[d0][0] == 2)
						{
							for (int nx = 0; nx < 2; ++nx)
							{
								for (int ny = 0; ny < 2; ++ny)
								{
									int fi0 = m_level[ilevel].level_pro_ghost[g0].fi.son[nx][ny][0];
									Assert(fi0 > -1, "The son of the pro ghost in flux modification error!!!", 55);
									mdyface.back().fface[nx][ny] = m_level[ilevel+1].m_box[fi0].faces[boxface[d0][0]][boxface[d0][1]];
									if (mdyface.back().fface[nx][ny] < 0)
									{
										printf("Main box neib (%d,%d,%d) is %d is ghost %d type %d tag %d!!!\n",
											f0[d0][0],f0[d0][1],f0[d0][2],aneib0,m_level[ilevel].m_box.isghost(aneib0),
											m_level[ilevel].m_box[aneib0].type,
											m_level[ilevel].m_tag[aneib0].tag);
									}
									Assert(mdyface.back().fface[nx][ny] > -1, "The face of the son in flux modification error!!! z 0",57);
								}
							}
						}					
					}
				}
			}
		}
		mdyface0[ilevel].Addnew(mdyface);
		mdyface0[ilevel].DirectlyReduceNew();
	}
	GiveAFlag("Finish finding the faces of flux modification!!!", 5);
}

void NS_Solver::BalanceFaceFlux(vector<DataArray<Facepair> > & mdyface, vector<DataArray<FlowVec> > & fluxarray)
{
#ifndef TIME_REFINE	
	for (int ilevel = 0; ilevel < a_mesh.cur_level_num-1; ++ilevel)
	{
		int bs = mdyface[ilevel].ps();
		int be = mdyface[ilevel].pe();
		if (!level_twod_flag[ilevel+1])
		{
			for (int i = bs; i < be; ++i)
			{
				int ci0 = mdyface[ilevel][i].cface;
				Assert(ci0 > -1 && ci0 < a_mesh.m_level[ilevel].m_face.size(), "Coarse face index error in flux balance!!!",86);
				FlowVec c_old_flux = fluxarray[ilevel][ci0]*0.25;					
				fluxarray[ilevel][ci0] *= 0.5;
				for (int nx = 0; nx < 2; ++nx)
				{
					for (int ny = 0; ny < 2; ++ny)
					{
						int fi0 = mdyface[ilevel][i].fface[nx][ny];
						Assert(fi0 > -1 && fi0 < a_mesh.m_level[ilevel+1].m_face.size(), "Fine face index error in flux balance!!!",86);
						//fluxarray[ilevel+1][fi0] = (fluxarray[ilevel+1][fi0] + c_old_flux)*0.5;
						fluxarray[ilevel+1][fi0] *= 0.5;
						fluxarray[ilevel][ci0] += fluxarray[ilevel+1][fi0];
						fluxarray[ilevel+1][fi0] += c_old_flux*0.5;
					}
				}
			}						
		}
		else
		{
			for (int i = bs; i < be; ++i)
			{
				int ci0 = mdyface[ilevel][i].cface;
				Assert(ci0 > -1 && ci0 < a_mesh.m_level[ilevel].m_face.size(), "Coarse face index error in flux balance!!!",86);
				FlowVec c_old_flux = fluxarray[ilevel][ci0];
				fluxarray[ilevel][ci0] *= 0.5;
				if (a_mesh.m_level[ilevel].m_face[ci0].fnv != 2)
				{			
					for (int nx = 0; nx < 2; ++nx)
					{
						int fi0 = mdyface[ilevel][i].fface[nx][0];
						Assert(fi0 > -1 && fi0 < a_mesh.m_level[ilevel+1].m_face.size(), "Fine face index error in flux balance!!!",86);
						fluxarray[ilevel+1][fi0] *= 0.5;
						fluxarray[ilevel][ci0] += fluxarray[ilevel+1][fi0];
						fluxarray[ilevel+1][fi0] += c_old_flux*0.25;
					}
				}				
				else
				{
					for (int nx = 0; nx < 2; ++nx)
					{
						for (int ny = 0; ny < 2; ++ny)
						{
							int fi0 = mdyface[ilevel][i].fface[nx][ny];
							Assert(fi0 > -1 && fi0 < a_mesh.m_level[ilevel+1].m_face.size(), "Fine face index error in flux balance!!!",86);
							fluxarray[ilevel+1][fi0] *= 0.5;
							fluxarray[ilevel][ci0] += fluxarray[ilevel+1][fi0];
							fluxarray[ilevel+1][fi0] += c_old_flux*0.125;
						}
					}					
				}				
			}
		}
	}
	MPI_Barrier(share_comm);
#endif	
}

void NS_Solver::MiddleStepFlux(vector<DataArray<Facepair> > & mdyface)
{
	for (int ilevel = 0; ilevel < a_mesh.cur_level_num; ++ilevel)
	{
		int bs = mdyface[ilevel].ps();
		int be = mdyface[ilevel].pe();
		for (int i = bs; i < be; ++i)
		{
			int ci0 = mdyface[ilevel][i].cface;
#ifdef VISCOUSITY			
			mdyface[ilevel][i].cvflux = vsflux[ilevel][ci0];
#endif
			mdyface[ilevel][i].cinvflux = faceflux[ilevel][ci0];
			for (int nx = 0; nx < 2; ++nx)
			{
				for (int ny = 0; ny < 2; ++ny)
				{	
					int fi0 = mdyface[ilevel][i].fface[nx][ny];
#ifdef VISCOUSITY					
					mdyface[ilevel][i].fvflux[nx][ny] = vsflux[ilevel+1][fi0];
#endif					
					mdyface[ilevel][i].finvflux[nx][ny] = faceflux[ilevel+1][fi0];	
				}
			}
		}
	}
}