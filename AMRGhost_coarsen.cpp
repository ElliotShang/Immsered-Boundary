#include "AMRLevel.H"
/*-------------------Notation-----------------------------*/
/*---Effect of coarsening---------------------------------*/
//1. some block pairs of level i+1 will be removed;
//---> detag = 0

//2. some block pairs of level i+1 will be closed;
//---> for active bkp of level i+1, test whether it is needed by the active cells in the ghost region;

//3. some block pairs of level i will be reopened;
//---> for nonactive bkps of level i, test whether it is needed by the active cells in the ghost region;

//4. some closed bkp pro structures will be removed;
//---> for each pair in close bkp, test whether the fine cells will be cleared;
//---> conditions: 1. detag = 0; 2. expired prolongation pair in its block.

//5. some level_pro_ghost pro structures will be removed;
//---> for each pair in level_pro_ghost, test whether the fine cells in the ghost region will be cleared;
//---> conditions: 1. detag = 0 in the ghost region; 
//-----------------2. expired prolongation pair in the ghost region, obtained from the adjacent block;

void AMRLevel::MarkProPairBeforeSyn()
{
	int bs = f_pro_ghost.ps();
	int be = f_pro_ghost.pe();
	if (!finelevel->twodflag)
	{		
		for (int i = bs; i < be; ++i)
		{
			for (Point_iterator p(0,2); p.end(); ++p)
			{
				int fi0 = f_pro_ghost[i].fi.son[p.i][p.j][p.k];
				finelevel->m_tag[fi0].detag = -20;
			}
		}
	}
	else
	{
		for (int i = bs; i < be; ++i)
		{
			for (Point_iterator_2d p(0,2); p.end(); ++p)
			{
				int fi0 = f_pro_ghost[i].fi.son[p.i][p.j][p.k];
				finelevel->m_tag[fi0].detag = -20;
			}
		}
	}
	MPI_Barrier(share_comm);
}

void AMRLevel::TagDerefineGhost()
{
	MarkProPairBeforeSyn();
	finelevel->TellBlockPairNewDetag();
	for (int i0 = finelevel->m_tag.ps(); i0 < finelevel->m_tag.pe(); ++i0)
	{
		if ((finelevel->m_tag[i0].detag == 0 || finelevel->m_tag[i0].detag == 1) && finelevel->m_box[i0].type != Dmghost)
		{
			int sx00 = finelevel->m_box[i0].ix() & 1;
			int sy00 = finelevel->m_box[i0].iy() & 1;
			int sz00 = finelevel->m_box[i0].iz() & 1;					
			for (Point_iterator p(0,2); p.end(); ++p)
			{
				int aneib = finelevel->m_box[i0].neib[sx00+p.i][sy00+p.j][sz00+p.k];
				if (aneib > -1)
				{

					if (finelevel->m_tag[aneib].detag == -1 ||
							 finelevel->m_tag[aneib].detag == -2)
					{
						int detag0 = -1;
						//PRINTFinLEVEL("Nonzero detag %d neib %d",cur_level,i0,aneib);
						Assert(finelevel->m_box[aneib].type != Notype, "The box must have a type when TagDerefineGhost!!!", 99);
						if (finelevel->m_box[aneib].type != Dmghost)
						{
							detag0 = -2;
							if (finelevel->m_tag[i0].detag == 0)
							{
								int nss[3] = {1-sx00, 1-sy00, 1-sz00};
								int nse[3] = {3-sx00, 3-sy00, 3-sz00};
								for (int nx = nss[0]; nx < nse[0]; ++nx)
								{
									for (int ny = nss[1]; ny < nse[1]; ++ny)
									{
										for (int nz = nss[2]; nz < nse[2]; ++nz)
										{
											int sneib = finelevel->m_box[i0].neib[nx][ny][nz];
											Assert(sneib > -1, "The four box must be non-negative!!! 2", 173);
											finelevel->m_tag[sneib].detag = 1;
										}
									}
								}
							}
						}
						else
						{
							bool removeflag = true;
							if (finelevel->m_tag[aneib].detag == -1)
							{
								int sx,sy,sz;
								finelevel->BoxOddeven(aneib,sx,sy,sz);
								int nxs[3] = {1-sx,1-sy,1-sz};
								int nxe[3] = {3-sx,3-sy,3-sz};
								for (int ni = nxs[0]; ni < nxe[0]; ++ni)
								{
									for (int nj = nxs[1]; nj < nxe[1]; ++nj)
									{
										for (int nk = nxs[2]; nk < nxe[2]; ++nk)
										{
											int an1 = finelevel->m_box[aneib].neib[ni][nj][nk];
											// if (sx0 == 51 && sy0 == -1 && sz0 == 40)
											Assert(an1>-1, "The four cells must be non-negative when tag coarsen ghost!!! 102", 102);
											int sx0 = finelevel->m_box[an1].ix() & 1;
											int sy0 = finelevel->m_box[an1].iy() & 1;
											int sz0 = finelevel->m_box[an1].iz() & 1;
											for (Point_iterator p(0,2); p.end(); ++p)
											{
												// sx = 0 ni = 1 p.i = 0 snx = 0
												// sx = 0 ni = 2 p.i = 0 snx = 1
												// sx = 1 ni = 0 p.i = 0 snx = 0
												// sx = 1 ni = 1 p.i = 0 snx = 1												
												int snx = sx0+p.i;
												int sny = sy0+p.j;
												int snz = sz0+p.k;
												int an2 = finelevel->m_box[an1].neib[snx][sny][snz];
												if (an2 > -1)
												{
													if (finelevel->m_box[an2].type != Dmghost && (finelevel->m_tag[an2].detag == -1 || finelevel->m_tag[an2].detag == -2))
													{
														removeflag = false;
														goto TAGVALUE;
													}
												}
											}
										}
									}
								}
								TAGVALUE:;
								if (removeflag) detag0 = 0;
								else 						detag0 = -2;
							}
						}
						if (finelevel->m_tag[aneib].detag == -1)
						{
#ifdef DEBUG
							if (finelevel->m_box[aneib].type != Dmghost && detag0 == 0)
							{
								PRINTFinLEVEL("A non-domain ghost can not be detag 0 Box (%d,%d,%d)!!!",
									cur_level,
									finelevel->m_box[aneib].ix(),
									finelevel->m_box[aneib].iy(),
									finelevel->m_box[aneib].iz());
								MPI_Abort(MPI_COMM_WORLD,166);
							}
#endif							
							Assert(detag0 != -1, "detag0 error!!!", 166);
							int sx = finelevel->m_box[aneib].ix() & 1;
							int sy = finelevel->m_box[aneib].iy() & 1;
							int sz = finelevel->m_box[aneib].iz() & 1;
							int nss[3] = {1-sx, 1-sy, 1-sz};
							int nse[3] = {3-sx, 3-sy, 3-sz};
							for (int nx = nss[0]; nx < nse[0]; ++nx)
							{
								for (int ny = nss[1]; ny < nse[1]; ++ny)
								{
									for (int nz = nss[2]; nz < nse[2]; ++nz)
									{
										int sneib = finelevel->m_box[aneib].neib[nx][ny][nz];
										Assert(sneib > -1, "The four box must be non-negative 1!!!", 146);
										finelevel->m_tag[sneib].detag = detag0;
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
	Point dpt[6] = {Point(0,1,1),Point(2,1,1),Point(1,0,1),Point(1,2,1),Point(1,1,0),Point(1,1,2)};
	for (int i = finelevel->m_tag.ps(); i < finelevel->m_tag.pe(); ++i)
	{
		if (finelevel->m_box[i].type != Dmghost)
		{
			if (finelevel->m_tag[i].detag == 1 || finelevel->m_tag[i].detag == 0)
			{
				if ((finelevel->m_box[i].ix()&1)+(finelevel->m_box[i].iy()&1)+(finelevel->m_box[i].iz()&1) == 0)
				{
					for (int d0 = 0; d0 < 6; ++d0)
					{
						int an1 = finelevel->m_box[i].neib[dpt[d0][0]][dpt[d0][1]][dpt[d0][2]];
						if (an1 > -1)
						{
							int an2 = finelevel->m_box[an1].neib[dpt[d0][0]][dpt[d0][1]][dpt[d0][2]];
							if (an2 > -1)
							{
								int an3 = finelevel->m_box[an2].neib[dpt[d0][0]][dpt[d0][1]][dpt[d0][2]];
								if (an3 > -1)
								{
									int an4 = finelevel->m_box[an3].neib[dpt[d0][0]][dpt[d0][1]][dpt[d0][2]];
									if (an4 > -1)
									{
										if (finelevel->m_tag[an4].detag == -1)
										{
											int sx = finelevel->m_box[an4].ix()&1;
											int sy = finelevel->m_box[an4].iy()&1;
											int sz = finelevel->m_box[an4].iz()&1;
											int ss[3] = {1-sx, 1-sy, 1-sz};
											int se[3] = {3-sx, 3-sy, 3-sz};
											for (int nx = ss[0]; nx < se[0]; ++nx)
											{
												for (int ny = ss[1]; ny < se[1]; ++ny)
												{
													for (int nz = ss[2]; nz < se[2]; ++nz)
													{
														int anb = finelevel->m_box[an4].neib[nx][ny][nz];
														Assert(anb > -1, "The four cells must be non-negative during coarsen 238!!!", 238);
														Assert(finelevel->m_tag[anb].detag == -1 || finelevel->m_tag[anb].detag == -3, "The four cells must be non-negative during coarsen 239!!!", 239);
														finelevel->m_tag[anb].detag = -3;
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
		}
	}
	MPI_Barrier(share_comm);
#ifdef DEBUG	
	for (int i0 = finelevel->m_tag.ps(); i0 < finelevel->m_tag.pe(); ++i0)
	{
		if (finelevel->m_tag[i0].detag != -1)
		{
			int sx00 = finelevel->m_box[i0].ix() & 1;
			int sy00 = finelevel->m_box[i0].iy() & 1;
			int sz00 = finelevel->m_box[i0].iz() & 1;
			for (int nx = 1-sx00; nx < 3-sx00; ++nx)
			{
				for (int ny = 1-sy00; ny < 3-sy00; ++ny)
				{
					for (int nz = 1-sz00; nz < 3-sz00; ++nz)
					{
						int sneib = finelevel->m_box[i0].neib[nx][ny][nz];
						if (sneib < 0)
						{
							PRINTFinLEVEL("Box odd %d %d %d box is %d neib %d %d %d is %d\n", cur_level, sx00, sy00, sz00, i0, nx, ny, nz, sneib);
						}
						Assert(sneib > -1, "The four box must be non-negative 79!!!", 79);
						if (finelevel->m_tag[sneib].detag != finelevel->m_tag[i0].detag)
						{
							printf("Box (%d,%d,%d) detag is 0 but its neib (%d,%d,%d) detag is %d\n",
								finelevel->m_box[i0].ix(), finelevel->m_box[i0].iy(), finelevel->m_box[i0].iz(), 
								nx, ny, nz, finelevel->m_tag[sneib].detag);
							MPI_Abort(MPI_COMM_WORLD, 84);
						}
					}
				}
			}	
		}
	}
#endif	
}

void AMRLevel::TagDerefineGhost_2d()
{
	MarkProPairBeforeSyn();
	finelevel->TellBlockPairNewDetag();
	for (int i0 = finelevel->m_tag.ps(); i0 < finelevel->m_tag.pe(); ++i0)
	{
		if ((finelevel->m_tag[i0].detag == 0 || finelevel->m_tag[i0].detag == 1) && finelevel->m_box[i0].type != Dmghost)
		{
			int sx00 = finelevel->m_box[i0].ix() & 1;
			int sy00 = finelevel->m_box[i0].iy() & 1;				
			for (Point_iterator_2d p(0,2); p.end(); ++p)
			{
				for (int nz0 = 0; nz0 < 3; ++nz0)
				{
					int aneib = finelevel->m_box[i0].neib[sx00+p.i][sy00+p.j][nz0];
					if (aneib > -1)
					{
						if (finelevel->m_tag[aneib].detag == -1 ||
								 finelevel->m_tag[aneib].detag == -2)
						{
							int detag0 = -1;
							Assert(finelevel->m_box[aneib].type != Notype, "The box must have a type when TagDerefineGhost!!!", 99);
							if (finelevel->m_box[aneib].type != Dmghost)
							{
								detag0 = -2;
								if (finelevel->m_tag[i0].detag == 0)
								{
									int nss[2] = {1-sx00, 1-sy00};
									int nse[2] = {3-sx00, 3-sy00};
									for (int nx = nss[0]; nx < nse[0]; ++nx)
									{
										for (int ny = nss[1]; ny < nse[1]; ++ny)
										{
											int sneib = finelevel->m_box[i0].neib[nx][ny][1];
											Assert(sneib > -1, "The four box must be non-negative!!! 2", 173);
											finelevel->m_tag[sneib].detag = 1;
										}
									}
								}
							}
							else
							{
								bool removeflag = true;
								if (finelevel->m_tag[aneib].detag == -1)
								{
									/*Cycle the son-cells of a mom domain ghost*/
									/*to see if they has near cells that will not be coarsened*/
									int sx,sy;
									finelevel->BoxOddeven(aneib,sx,sy);
									int nxs[2] = {1-sx,1-sy};
									int nxe[2] = {3-sx,3-sy};
									for (int ni = nxs[0]; ni < nxe[0]; ++ni)
									{
										for (int nj = nxs[1]; nj < nxe[1]; ++nj)
										{
											int an1 = finelevel->m_box[aneib].neib[ni][nj][1];
											Assert(an1>-1, "The four cells must be non-negative when tag coarsen ghost!!! 102", 102);
											int sx0 = finelevel->m_box[an1].ix() & 1;
											int sy0 = finelevel->m_box[an1].iy() & 1;
											for (Point_iterator_2d p(0,2); p.end(); ++p)
											{
												// sx = 0 ni = 1 p.i = 0 snx = 0
												// sx = 0 ni = 2 p.i = 0 snx = 1
												// sx = 1 ni = 0 p.i = 0 snx = 0
												// sx = 1 ni = 1 p.i = 0 snx = 1
												for (int nz1 = 0; nz1 < 3; ++nz1)
												{									
													int snx = sx0+p.i;
													int sny = sy0+p.j;
													int an2 = finelevel->m_box[an1].neib[snx][sny][nz1];
													if (an2 > -1)
													{
														if (finelevel->m_box[an2].type != Dmghost && 
															 (finelevel->m_tag[an2].detag == -1 || finelevel->m_tag[an2].detag == -2))
														{
															removeflag = false;
															goto TAGVALUE;
														}
													}
												}
											}
										}
									}
									TAGVALUE:;
									if (removeflag) detag0 = 0;
									else 						detag0 = -2;
								}
							}
							if (finelevel->m_tag[aneib].detag == -1)
							{
#ifdef DEBUG
								if (finelevel->m_box[aneib].type != Dmghost && detag0 == 0)
								{
									PRINTFinLEVEL("A non-domain ghost can not be detag 0 Box (%d,%d,%d)!!!",
										cur_level,
										finelevel->m_box[aneib].ix(),
										finelevel->m_box[aneib].iy(),
										finelevel->m_box[aneib].iz());
									MPI_Abort(MPI_COMM_WORLD,166);
								}
#endif							
								Assert(detag0 != -1, "detag0 error!!!", 166);
								int sx = finelevel->m_box[aneib].ix() & 1;
								int sy = finelevel->m_box[aneib].iy() & 1;
								int nss[2] = {1-sx, 1-sy};
								int nse[2] = {3-sx, 3-sy};
								for (int nx = nss[0]; nx < nse[0]; ++nx)
								{
									for (int ny = nss[1]; ny < nse[1]; ++ny)
									{
										int sneib = finelevel->m_box[aneib].neib[nx][ny][1];
										Assert(sneib > -1, "The four box must be non-negative 1!!!", 146);
										finelevel->m_tag[sneib].detag = detag0;
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
	Point dpt[6] = {Point(0,1,1),Point(2,1,1),Point(1,0,1),Point(1,2,1),Point(1,1,0),Point(1,1,2)};
	for (int i = finelevel->m_tag.ps(); i < finelevel->m_tag.pe(); ++i)
	{
		if (finelevel->m_box[i].type != Dmghost)
		{
			if (finelevel->m_tag[i].detag == 1 || finelevel->m_tag[i].detag == 0)
			{
				if ((finelevel->m_box[i].ix()&1)+(finelevel->m_box[i].iy()&1) == 0)
				{
					for (int d0 = 0; d0 < 4; ++d0)
					{
						int an1 = finelevel->m_box[i].neib[dpt[d0][0]][dpt[d0][1]][dpt[d0][2]];
						if (an1 > -1)
						{
							int an2 = finelevel->m_box[an1].neib[dpt[d0][0]][dpt[d0][1]][dpt[d0][2]];
							if (an2 > -1)
							{
								int an3 = finelevel->m_box[an2].neib[dpt[d0][0]][dpt[d0][1]][dpt[d0][2]];
								if (an3 > -1)
								{
									int an4 = finelevel->m_box[an3].neib[dpt[d0][0]][dpt[d0][1]][dpt[d0][2]];
									if (an4 > -1)
									{
										if (finelevel->m_tag[an4].detag == -1)
										{
											int sx = finelevel->m_box[an4].ix()&1;
											int sy = finelevel->m_box[an4].iy()&1;
											int ss[2] = {1-sx, 1-sy};
											int se[2] = {3-sx, 3-sy};
											for (int nx = ss[0]; nx < se[0]; ++nx)
											{
												for (int ny = ss[1]; ny < se[1]; ++ny)
												{
													int anb = finelevel->m_box[an4].neib[nx][ny][1];
													Assert(anb > -1, "The four cells must be non-negative during coarsen 238!!!", 238);
													Assert(finelevel->m_tag[anb].detag == -1 || finelevel->m_tag[anb].detag == -3, "The four cells must be non-negative during coarsen 239!!!", 239);
													finelevel->m_tag[anb].detag = -3;
												}
											}
										}
									}
								}
							}
						}
					}
					for (int d0 = 4; d0 < 6; ++d0)
					{
						int an1 = finelevel->m_box[i].neib[dpt[d0][0]][dpt[d0][1]][dpt[d0][2]];
						if (an1 > -1)
						{
							int an2 = finelevel->m_box[an1].neib[dpt[d0][0]][dpt[d0][1]][dpt[d0][2]];
							if (an2 > -1)
							{
								if (finelevel->m_tag[an2].detag == -1)
								{
									int sx = finelevel->m_box[an2].ix()&1;
									int sy = finelevel->m_box[an2].iy()&1;
									int ss[2] = {1-sx, 1-sy};
									int se[2] = {3-sx, 3-sy};
									for (int nx = ss[0]; nx < se[0]; ++nx)
									{
										for (int ny = ss[1]; ny < se[1]; ++ny)
										{
											int anb = finelevel->m_box[an2].neib[nx][ny][1];
											Assert(anb > -1, "The four cells must be non-negative during coarsen 238!!!", 238);
											Assert(finelevel->m_tag[anb].detag == -1 || finelevel->m_tag[anb].detag == -3, "The four cells must be non-negative during coarsen 239!!!", 239);
											finelevel->m_tag[anb].detag = -3;
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
	MPI_Barrier(share_comm);
#ifdef DEBUG	
	for (int i0 = finelevel->m_tag.ps(); i0 < finelevel->m_tag.pe(); ++i0)
	{
		if (finelevel->m_tag[i0].detag != -1)
		{
			int sx00 = finelevel->m_box[i0].ix() & 1;
			int sy00 = finelevel->m_box[i0].iy() & 1;
			for (int nx = 1-sx00; nx < 3-sx00; ++nx)
			{
				for (int ny = 1-sy00; ny < 3-sy00; ++ny)
				{
					int sneib = finelevel->m_box[i0].neib[nx][ny][1];
					if (sneib < 0)
					{
						PRINTFinLEVEL("Two D levels Box odd %d %d box is %d neib %d %d %d is %d\n", 
							cur_level, sx00, sy00, i0, nx, ny, 1, sneib);
					}
					Assert(sneib > -1, "The four box must be non-negative 79!!!", 79);
					if (finelevel->m_tag[sneib].detag != finelevel->m_tag[i0].detag)
					{
						printf("Box (%d,%d,%d) detag is 0 but its neib (%d,%d,%d) detag is %d\n",
							finelevel->m_box[i0].ix(), finelevel->m_box[i0].iy(), finelevel->m_box[i0].iz(), 
							nx, ny, 1, finelevel->m_tag[sneib].detag);
						MPI_Abort(MPI_COMM_WORLD, 84);
					}
				}
			}	
		}
	}
#endif	
}

void AMRLevel::AdjustOldPair_Coarsen()
{
	finelevel->ClearExpiredDomainGhost();
	GiveLevelFlag("Finish clear fine domain ghost!!!", cur_level, 5);
	RemoveExpiredProlongation_Coarsen();
	GiveLevelFlag("Finish remove expired prolongation pair!!!", cur_level, 5);
	RemoveFineBlockPair_Coarsen();
	GiveLevelFlag("Finish remove fine blockpair!!!", cur_level, 5);
	RemoveExpiredLevelProGhost_Coarsen();
	GiveLevelFlag("Finish remove level_pro_ghost!!!", cur_level, 5);
	CheckExitingResPair();
	GiveLevelFlag("Finish CheckExitingResPair!!!", cur_level, 5);
}
void AMRLevel::RemoveExpiredProlongation_Coarsen()
{
	/*Remove the expired prolongation pair*/
	if (!finelevel->twodflag)
	{
		int fps = f_pro_ghost.ps();
		int fpe = f_pro_ghost.pe();
		for (int i = fps; i < fpe; ++i)
		{
			bool removeflag = true;
			for (Point_iterator p(0,2); p.end(); ++p)
			{
				int fi0 = f_pro_ghost[i].fi.son[p.i][p.j][p.k];
				int sxe = p.i+2; int sye = p.j+2; int sze = p.k+2;
				for (int nx = p.i; nx < sxe; ++nx)
				{
					for (int ny = p.j; ny < sye; ++ny)
					{
						for (int nz = p.k; nz < sze; ++nz)
						{
							int aneib = finelevel->m_box[fi0].neib[nx][ny][nz];
							if (aneib > -1)
							{
								if (finelevel->m_box[aneib].type != Dmghost)
								{
									if (finelevel->m_tag[aneib].detag == -1 || 
											finelevel->m_tag[aneib].detag == -2 ||
											finelevel->m_tag[aneib].detag == -3)
									{
										removeflag = false;
										goto OPTION0;
									}
								}
							}
						}
					}
				}
			}	
			OPTION0:;
			//printf("df pro %d removeflag is %d\n", i, removeflag);
			if (removeflag)
			{
				for (Point_iterator p(0,2); p.end(); ++p)
				{
					int fi0 = f_pro_ghost[i].fi.son[p.i][p.j][p.k];
					// if (finelevel->m_box[fi0].ix() == 31 &&
					// 		finelevel->m_box[fi0].iy() == 24 &&
					// 		finelevel->m_box[fi0].iz() == 11)
					// {
					// 	printf("Why this pro cell was removed!!!!!~~~~~~~~~~~~~~~~\n");
					// }
					// ----This cell will be deleted from both the blockpair array 
					/*----and the level_pro_ghost of the neib block*/
					finelevel->m_tag[fi0].detag = -21;
				}
				drf_remove_pro_ghost_fcell.push_back(i);
			}
		}
	}
	else
	{
		int fps = f_pro_ghost.ps();
		int fpe = f_pro_ghost.pe();
		for (int i = fps; i < fpe; ++i)
		{
			bool removeflag = true;
			for (Point_iterator_2d p(0,2); p.end(); ++p)
			{
				int fi0 = f_pro_ghost[i].fi.son[p.i][p.j][p.k];
				int sxe = p.i+2; int sye = p.j+2;
				for (int nx = p.i; nx < sxe; ++nx)
				{
					for (int ny = p.j; ny < sye; ++ny)
					{
						for (int nz = 0; nz < 3; ++nz)
						{
							int aneib = finelevel->m_box[fi0].neib[nx][ny][nz];
							if (aneib > -1)
							{
								if (finelevel->m_box[aneib].type != Dmghost)
								{
									if (finelevel->m_tag[aneib].detag == -1 || 
											finelevel->m_tag[aneib].detag == -2 ||
											finelevel->m_tag[aneib].detag == -3)
									{
										removeflag = false;
										goto OPTION1;
									}
								}
							}
						}
					}
				}
			}	
			OPTION1:;
			//printf("df pro %d removeflag is %d\n", i, removeflag);
			if (removeflag)
			{
				for (Point_iterator_2d p(0,2); p.end(); ++p)
				{
					int fi0 = f_pro_ghost[i].fi.son[p.i][p.j][p.k];
					// if (finelevel->m_box[fi0].ix() == 31 &&
					// 		finelevel->m_box[fi0].iy() == 24 &&
					// 		finelevel->m_box[fi0].iz() == 11)
					// {
					// 	printf("Why this pro cell was removed!!!!!~~~~~~~~~~~~~~~~\n");
					// }
					// ----This cell will be deleted from both the blockpair array 
					/*----and the level_pro_ghost of the neib block*/
					finelevel->m_tag[fi0].detag = -21;
				}
				drf_remove_pro_ghost_fcell.push_back(i);
			}
		}
	}
	MPI_Barrier(share_comm);
	finelevel->TellBlockPairNewDetag();
}

void AMRLevel::CheckExitingResPair()
{
	if (!finelevel->twodflag)
	{
		for (int ig = 0; ig < ighost; ++ig)
		{
			for (int i = f_res_ghost[ig].ps(); i < f_res_ghost[ig].pe(); ++i)
			{
				int fi00 = f_res_ghost[ig][i].fi.son[0][0][0];
				int ci0 = f_res_ghost[ig][i].ci;
				if (finelevel->m_tag[fi00].detag < -1)
				{
					if (finelevel->m_tag[fi00].detag == -2 && ig == 1)
					{
						drf_res_move2layer.push_back(i);
					}
					for (Point_iterator p(0,2); p.end(); ++p)
					{
						int fi0 = f_res_ghost[ig][i].fi.son[p.i][p.j][p.k];
						finelevel->m_tag[fi0].detag = -1;
					}
					drf_recover_pair[ig].push_back(i);
				}
				else if (finelevel->m_tag[fi00].detag > -1)
				{
					if (finelevel->m_tag[fi00].detag == 0)
					{
						drf_res_keep_c_remove_f[ig].push_back(i);
						drf_remove_respair.push_back(f_res_ghost[ig][i]);
						for (Point_iterator p(0,2); p.end(); ++p)
						{
							int fi0 = f_res_ghost[ig][i].fi.son[p.i][p.j][p.k];
							finelevel->m_tag[fi0].detag = -1;
						}
					}
					else if (finelevel->m_tag[fi00].detag == 1)
					{
						drf_res_switch_ghost_layer[ig].push_back(i);
						for (Point_iterator p(0,2); p.end(); ++p)
						{
							int fi0 = f_res_ghost[ig][i].fi.son[p.i][p.j][p.k];
							finelevel->m_tag[fi0].detag = -1;
						}
					}
				}
			}
		}
	}
	else
	{
		for (int ig = 0; ig < ighost; ++ig)
		{
			for (int i = f_res_ghost[ig].ps(); i < f_res_ghost[ig].pe(); ++i)
			{
				int fi00 = f_res_ghost[ig][i].fi.son[0][0][0];
				int ci0 = f_res_ghost[ig][i].ci;
				if (finelevel->m_tag[fi00].detag < -1)
				{
					if (finelevel->m_tag[fi00].detag == -2 && ig == 1)
					{
						drf_res_move2layer.push_back(i);
					}
					for (Point_iterator_2d p(0,2); p.end(); ++p)
					{
						int fi0 = f_res_ghost[ig][i].fi.son[p.i][p.j][p.k];
						finelevel->m_tag[fi0].detag = -1;
					}
					drf_recover_pair[ig].push_back(i);
				}
				else if (finelevel->m_tag[fi00].detag > -1)
				{
					if (finelevel->m_tag[fi00].detag == 0)
					{
						drf_res_keep_c_remove_f[ig].push_back(i);
						drf_remove_respair.push_back(f_res_ghost[ig][i]);
						for (Point_iterator_2d p(0,2); p.end(); ++p)
						{
							int fi0 = f_res_ghost[ig][i].fi.son[p.i][p.j][p.k];
							finelevel->m_tag[fi0].detag = -1;
						}
					}
					else if (finelevel->m_tag[fi00].detag == 1)
					{
						drf_res_switch_ghost_layer[ig].push_back(i);
						for (Point_iterator_2d p(0,2); p.end(); ++p)
						{
							int fi0 = f_res_ghost[ig][i].fi.son[p.i][p.j][p.k];
							finelevel->m_tag[fi0].detag = -1;
						}
					}
				}
			}
		}
	}
	MPI_Barrier(share_comm);
}

void AMRLevel::RemoveExpiredLevelProGhost_Coarsen()
{
	/*-------------Executed for the pro pairs in the ghost region--------*/
	if (!finelevel->twodflag)
	{
		for (int i = level_pro_ghost.ps(); i < level_pro_ghost.pe(); ++i)
		{
			int i0 = level_pro_ghost[i].fi.son[0][0][0];
			if (finelevel->m_tag[i0].detag == -1 || finelevel->m_tag[i0].detag == -20)
			{
				goto NEXTPAIR_0;
			}
			else if (finelevel->m_tag[i0].detag == -21)
			{
				drf_remove_level_pro_ghost.push_back(i);
			}
			else if (finelevel->m_tag[i0].detag == 0)
			{
				drf_remove_level_pro_ghost.push_back(i);
				drf_recover_level_pro_ghost.push_back(i);
				for (Point_iterator p(0,2); p.end(); ++p)
				{
					int fi0 = level_pro_ghost[i].fi.son[p.i][p.j][p.k];
					finelevel->m_tag[fi0].detag = -1;
				}
			}
			else // if (finelevel->m_tag[i0].detag != -20)
			{
				drf_recover_level_pro_ghost.push_back(i);
				for (Point_iterator p(0,2); p.end(); ++p)
				{
					int fi0 = level_pro_ghost[i].fi.son[p.i][p.j][p.k];
					finelevel->m_tag[fi0].detag = -1;
				}
			}
			NEXTPAIR_0:;
		}
	}
	else
	{
		for (int i = level_pro_ghost.ps(); i < level_pro_ghost.pe(); ++i)
		{
			int i0 = level_pro_ghost[i].fi.son[0][0][0];
			if (finelevel->m_tag[i0].detag == -1 || finelevel->m_tag[i0].detag == -20)
			{
				goto NEXTPAIR_1;
			}
			else if (finelevel->m_tag[i0].detag == -21)
			{
				drf_remove_level_pro_ghost.push_back(i);
			}
			else if (finelevel->m_tag[i0].detag == 0)
			{
				drf_remove_level_pro_ghost.push_back(i);
				drf_recover_level_pro_ghost.push_back(i);
				for (Point_iterator_2d p(0,2); p.end(); ++p)
				{
					int fi0 = level_pro_ghost[i].fi.son[p.i][p.j][p.k];
					finelevel->m_tag[fi0].detag = -1;
				}
			}
			else // if (finelevel->m_tag[i0].detag != -20)
			{
				drf_recover_level_pro_ghost.push_back(i);
				for (Point_iterator_2d p(0,2); p.end(); ++p)
				{
					int fi0 = level_pro_ghost[i].fi.son[p.i][p.j][p.k];
					finelevel->m_tag[fi0].detag = -1;
				}
			}
			NEXTPAIR_1:;
		}
	}
	MPI_Barrier(share_comm);
}

void AMRLevel::RemoveFineBlockPair_Coarsen()
{
	if (NULL != finelevel)
	{
		for (int i = finelevel->blockpair.ps(); i < finelevel->blockpair.pe(); ++i)
		{
			int i0 = finelevel->blockpair[i].innode;
			if (finelevel->m_tag[i0].detag == 0 || finelevel->m_tag[i0].detag == -21)
			{
				finelevel->blockpair.givehole(i);
			}
		}
	}
}

void AMRLevel::BuildNewPair_Derefine()
{
	/*construct the new prolongation and restriction pair*/
	if (!finelevel->twodflag)
	{
		int bs = finelevel->m_tag.ps();
		int be = finelevel->m_tag.pe();
		for (int i = bs; i < be; ++i)
		{
			if (finelevel->m_tag[i].detag == 0)
			{
				Assert(finelevel->m_box[i].type != Dmghost, "A domain ghost detag can not be 0 when BuildNewPair_Derefine!!!", 637);
				finelevel->m_box.givehole(i);
			}
		}
		MPI_Barrier(share_comm);
		bs = finelevel->m_box.ps();
		be = finelevel->m_box.pe();
		for (int i = bs; i < be; ++i)
		{
			if (finelevel->m_tag[i].detag == 1)
			{
				int sx = finelevel->m_box[i].ix() & 1;
				int sy = finelevel->m_box[i].iy() & 1;
				int sz = finelevel->m_box[i].iz() & 1;
				if (sx+sy+sz==0)
				{
					new_f_pro_ghost.push_back(CtoFPair());
					for (Point_iterator p(0,2); p.end(); ++p)
					{
						int fi0 = finelevel->m_box[i].neib[p.i+1][p.j+1][p.k+1];
						new_f_pro_ghost.back().fi.son[p.i][p.j][p.k] = fi0;
						finelevel->m_tag[fi0].detag = 0;
						finelevel->m_box.ghost_index.push_back(fi0);
					}
				}
			}
			else if (finelevel->m_tag[i].detag == -2)
			{
				int sx = finelevel->m_box[i].ix() & 1;
				int sy = finelevel->m_box[i].iy() & 1;
				int sz = finelevel->m_box[i].iz() & 1;
				if (sx+sy+sz==0)
				{
					new_f_res_ghost[0].push_back(CtoFPair());
					for (Point_iterator p(0,2); p.end(); ++p)
					{
						int fi0 = finelevel->m_box[i].neib[p.i+1][p.j+1][p.k+1];
						new_f_res_ghost[0].back().fi.son[p.i][p.j][p.k] = fi0;
						finelevel->m_tag[fi0].detag = -2;
					}
				}
			}
			else if (finelevel->m_tag[i].detag == -3)
			{
				int sx = finelevel->m_box[i].ix() & 1;
				int sy = finelevel->m_box[i].iy() & 1;
				int sz = finelevel->m_box[i].iz() & 1;
				if (sx+sy+sz==0)
				{
					new_f_res_ghost[1].push_back(CtoFPair());
					for (Point_iterator p(0,2); p.end(); ++p)
					{
						int fi0 = finelevel->m_box[i].neib[p.i+1][p.j+1][p.k+1];
						new_f_res_ghost[1].back().fi.son[p.i][p.j][p.k] = fi0;
						finelevel->m_tag[fi0].detag = -2;
					}
				}
			}
		}
		MPI_Barrier(share_comm);
		int gbs = finelevel->m_box.gps();
		int gbe = finelevel->m_box.gpe();
		for (int i = gbs; i < gbe; ++i)
		{
#ifdef DEBUG
			if (finelevel->m_box[i].type == Dmghost)
			{
				Assert(finelevel->m_tag[i].detag == -1, "A dmghost detag should be -1!!!",713);
			}
#endif		
			if (finelevel->m_box[i].type == Blockghost)
			{
				if (finelevel->m_tag[i].detag == -1)
				{
					goto NEXTCELL_0;
				}
				else if (finelevel->m_tag[i].detag == 0)
				{
					finelevel->m_tag[i].detag = -2;
					drf_remove_blockghost.push_back(i);
				}
				else if (finelevel->m_tag[i].detag != -20 && finelevel->m_tag[i].detag != -21)
				{
					//Assert(finelevel->m_tag[i].detag != -20 && finelevel->m_tag[i].detag != -21, "A single block ghost can not be marked by the pro pair!!!",716);
					int sx, sy, sz;
					finelevel->BoxOddeven(i,sx,sy,sz);
					if (sx+sy+sz==0) drf_new_level_pro_ghost.push_back(i);
					finelevel->m_tag[i].detag = -2;
				}
			}
			NEXTCELL_0:;
		}
		MPI_Barrier(share_comm);
		for (int i = 0; i < finelevel->ghost_target.size(); ++i)
		{
			int c0 = finelevel->ghost_target[i];
			if (finelevel->m_tag[c0].detag == -2)
			{
				Assert(finelevel->m_box[c0].type==Blockghost,"Build blockpair for non-block ghost during coarsening!!!",743);
				int sx, sy, sz;
				finelevel->BoxOddeven(c0,sx,sy,sz);
				if (sx+sy+sz==0) drf_new_coarse_bpghost.push_back(i);
			}
		}
		MPI_Barrier(share_comm);
		PRINTFinLEVEL("Coarsen new pro pair is %d", cur_level, (int)new_f_pro_ghost.size());
		PRINTFinLEVEL("Coarsen new res 1 pair is %d", cur_level, (int)new_f_res_ghost[0].size());
		PRINTFinLEVEL("Coarsen new res 2 pair is %d", cur_level, (int)new_f_res_ghost[1].size());
		PRINTFinLEVEL("Coarsen ghost number is %d before including", cur_level, (int)m_box.ghost_index.size());
		PRINTFinLEVEL("Coarsen drf_new_level_pro_ghost is %d", cur_level, (int)drf_new_level_pro_ghost.size());
		PRINTFinLEVEL("Coarsen drf_new_coarse_bpghost is %d", cur_level, (int)drf_new_coarse_bpghost.size());
#ifdef DEBUG	
		for (int i = finelevel->m_tag.ps(); i < finelevel->m_tag.pe(); ++i)
		{
			if (finelevel->m_tag[i].detag == -3)
			{
				int ix0 = finelevel->m_box[i].ix();
				int iy0	= finelevel->m_box[i].iy();
				int iz0	= finelevel->m_box[i].iz();
				PRINTFinLEVEL("There should not be -3 detag at this stage!!! Box(%d,%d,%d) type %d isghost %d",
					finelevel->cur_level, ix0, iy0, iz0,
					finelevel->m_box[i].type,
					finelevel->m_box.isghost(i));
				for (int nx = 1-(ix0&1); nx < 3-(ix0&1); ++nx)
				{
					for (int ny = 1-(iy0&1); ny < 3-(iy0&1); ++ny)
					{
						for (int nz = 1-(iz0&1); nz < 3-(iz0&1); ++nz)
						{
							int an00 = finelevel->m_box[i].neib[nx][ny][nz];
							PRINTFinLEVEL("Box (%d, %d, %d) neib (%d,%d,%d) detag is %d", finelevel->cur_level,
								ix0,iy0,iz0,nx,ny,nz,
								finelevel->m_tag[an00].detag);
						}
					}
				}
				MPI_Abort(MPI_COMM_WORLD,679);
			}
		}
#endif
	}
	/*-----------For 2d refined cases---------------------*/
	else
	{
		int bs = finelevel->m_tag.ps();
		int be = finelevel->m_tag.pe();
		for (int i = bs; i < be; ++i)
		{
			if (finelevel->m_tag[i].detag == 0)
			{
				Assert(finelevel->m_box[i].type != Dmghost, "A domain ghost detag can not be 0 when BuildNewPair_Derefine!!!", 637);
				finelevel->m_box.givehole(i);
			}
		}
		MPI_Barrier(share_comm);
		bs = finelevel->m_box.ps();
		be = finelevel->m_box.pe();
		for (int i = bs; i < be; ++i)
		{
			if (finelevel->m_tag[i].detag == 1)
			{
				int sx = finelevel->m_box[i].ix() & 1;
				int sy = finelevel->m_box[i].iy() & 1;
				if (sx+sy==0)
				{
					new_f_pro_ghost.push_back(CtoFPair());
					for (Point_iterator_2d p(0,2); p.end(); ++p)
					{
						int fi0 = finelevel->m_box[i].neib[p.i+1][p.j+1][p.k+1];
						new_f_pro_ghost.back().fi.son[p.i][p.j][p.k] = fi0;
						finelevel->m_tag[fi0].detag = 0;
						finelevel->m_box.ghost_index.push_back(fi0);
					}
				}
			}
			else if (finelevel->m_tag[i].detag == -2)
			{
				int sx = finelevel->m_box[i].ix() & 1;
				int sy = finelevel->m_box[i].iy() & 1;
				if (sx+sy==0)
				{
					new_f_res_ghost[0].push_back(CtoFPair());
					for (Point_iterator_2d p(0,2); p.end(); ++p)
					{
						int fi0 = finelevel->m_box[i].neib[p.i+1][p.j+1][p.k+1];
						new_f_res_ghost[0].back().fi.son[p.i][p.j][p.k] = fi0;
						finelevel->m_tag[fi0].detag = -2;
					}
				}
			}
			else if (finelevel->m_tag[i].detag == -3)
			{
				int sx = finelevel->m_box[i].ix() & 1;
				int sy = finelevel->m_box[i].iy() & 1;
				if (sx+sy==0)
				{
					new_f_res_ghost[1].push_back(CtoFPair());
					for (Point_iterator_2d p(0,2); p.end(); ++p)
					{
						int fi0 = finelevel->m_box[i].neib[p.i+1][p.j+1][p.k+1];
						new_f_res_ghost[1].back().fi.son[p.i][p.j][p.k] = fi0;
						finelevel->m_tag[fi0].detag = -2;
					}
				}
			}
		}
		MPI_Barrier(share_comm);
		int gbs = finelevel->m_box.gps();
		int gbe = finelevel->m_box.gpe();
		for (int i = gbs; i < gbe; ++i)
		{
#ifdef DEBUG
			if (finelevel->m_box[i].type == Dmghost)
			{
				Assert(finelevel->m_tag[i].detag == -1, "A dmghost detag should be -1!!!",713);
			}
#endif		
			if (finelevel->m_box[i].type == Blockghost)
			{
				if (finelevel->m_tag[i].detag == -1)
				{
					goto NEXTCELL_1;
				}
				else if (finelevel->m_tag[i].detag == 0)
				{
					finelevel->m_tag[i].detag = -2;
					drf_remove_blockghost.push_back(i);
				}
				else if (finelevel->m_tag[i].detag != -20 && finelevel->m_tag[i].detag != -21)
				{
					//Assert(finelevel->m_tag[i].detag != -20 && finelevel->m_tag[i].detag != -21, "A single block ghost can not be marked by the pro pair!!!",716);
					int sx, sy;
					finelevel->BoxOddeven(i,sx,sy);
					if (sx+sy==0) drf_new_level_pro_ghost.push_back(i);
					finelevel->m_tag[i].detag = -2;
				}
			}
			NEXTCELL_1:;
		}
		MPI_Barrier(share_comm);
		for (int i = 0; i < finelevel->ghost_target.size(); ++i)
		{
			int c0 = finelevel->ghost_target[i];
			if (finelevel->m_tag[c0].detag == -2)
			{
				Assert(finelevel->m_box[c0].type==Blockghost,"Build blockpair for non-block ghost during coarsening!!!",743);
				int sx, sy;
				finelevel->BoxOddeven(c0,sx,sy);
				if (sx+sy==0) drf_new_coarse_bpghost.push_back(i);
			}
		}
		MPI_Barrier(share_comm);
		PRINTFinLEVEL("Coarsen new pro pair is %d", cur_level, (int)new_f_pro_ghost.size());
		PRINTFinLEVEL("Coarsen new res 1 pair is %d", cur_level, (int)new_f_res_ghost[0].size());
		PRINTFinLEVEL("Coarsen new res 2 pair is %d", cur_level, (int)new_f_res_ghost[1].size());
		PRINTFinLEVEL("Coarsen ghost number is %d before including", cur_level, (int)m_box.ghost_index.size());
		PRINTFinLEVEL("Coarsen drf_new_level_pro_ghost is %d", cur_level, (int)drf_new_level_pro_ghost.size());
		PRINTFinLEVEL("Coarsen drf_new_coarse_bpghost is %d", cur_level, (int)drf_new_coarse_bpghost.size());
#ifdef DEBUG	
		for (int i = finelevel->m_tag.ps(); i < finelevel->m_tag.pe(); ++i)
		{
			if (finelevel->m_tag[i].detag == -3)
			{
				int ix0 = finelevel->m_box[i].ix();
				int iy0	= finelevel->m_box[i].iy();
				int iz0	= finelevel->m_box[i].iz();
				PRINTFinLEVEL("There should not be -3 detag at this stage!!! Box(%d,%d,%d) type %d isghost %d",
					finelevel->cur_level, ix0, iy0, iz0,
					finelevel->m_box[i].type,
					finelevel->m_box.isghost(i));
				for (int nx = 1-(ix0&1); nx < 3-(ix0&1); ++nx)
				{
					for (int ny = 1-(iy0&1); ny < 3-(iy0&1); ++ny)
					{
						int an00 = finelevel->m_box[i].neib[nx][ny][1];
						PRINTFinLEVEL("Box (%d, %d, %d) neib (%d,%d,%d) detag is %d", finelevel->cur_level,
							ix0,iy0,iz0,nx,ny,1,
							finelevel->m_tag[an00].detag);
					}
				}
				MPI_Abort(MPI_COMM_WORLD,679);
			}
		}
#endif
	}
}

void AMRLevel::ConstructBlockProPair_Coarsen()
{
	vector<vector<S_newcoarsebp> > senddata(nodenum);
	vector<S_newcoarsebp> recvdata;
	vector<int> recvnode;
	/*------------------------------------------------------------*/
	int new_level_pro_size = drf_new_level_pro_ghost.size();
	vector<CtoFPair> new_level_pro_ghost(new_level_pro_size);
	if (!finelevel->twodflag)
	{
		for (int i = 0; i < new_level_pro_size; ++i)
		{
			int i0 = drf_new_level_pro_ghost[i];
			new_level_pro_ghost[i].ci = finelevel->m_tag[i0].detag;
			Assert(finelevel->m_tag[i0].detag > -1, "A level pro ghost has a invalid coarse tag!!! ", 444);
			for (Point_iterator p(0,2); p.end(); ++p)
			{
				new_level_pro_ghost[i].fi.son[p.i][p.j][p.k] = finelevel->m_box[i0].neib[1+p.i][1+p.j][1+p.k];
			}
		}
	}
	else
	{
		for (int i = 0; i < new_level_pro_size; ++i)
		{
			int i0 = drf_new_level_pro_ghost[i];
			new_level_pro_ghost[i].ci = finelevel->m_tag[i0].detag;
			Assert(finelevel->m_tag[i0].detag > -1, "A level pro ghost has a invalid coarse tag!!! ", 444);
			for (Point_iterator_2d p(0,2); p.end(); ++p)
			{
				new_level_pro_ghost[i].fi.son[p.i][p.j][p.k] = finelevel->m_box[i0].neib[1+p.i][1+p.j][1+p.k];
			}
		}
	}
	/*--------------From the newly coarsened normal cells-----------*/
	int coarse_bp_size = drf_new_coarse_bpghost.size();
	for (int i = 0; i < coarse_bp_size; ++i)
	{
		int i0 = drf_new_coarse_bpghost[i];
		int g0 = finelevel->ghost_target[i0];
		int tag0 = finelevel->m_tag[g0].detag;
		Assert(tag0 > -1, "Error for the request of building a new block pair!!!", 645);
		int node_dest = finelevel->originbp[i0].node;
		senddata[node_dest].push_back(S_newcoarsebp(finelevel->originbp[i0].index, tag0));
		m_box[tag0].type = Blockghost;
	}
	/*------------------------------------------------------------*/
	//printf("start to BcastNewInfo_NodeRank!!! senddata size is %d\n", (int)senddata[0].size());
	MPI_Barrier(share_comm);
	BcastNewInfo_NodeRank(senddata, recvdata, recvnode, 6);
	//MPI_Barrier(share_comm);
	//printf("Finish BcastNewInfo_NodeRank!!!\n");
	int rd = recvdata.size();
	vector<NodePair> new_bkp(rd);
	for (int i = 0; i < rd; ++i)
	{
		int i0 = recvdata[i].remotebkp;
		int c0 = finelevel->blockpair[i0].innode;
		int t0 = finelevel->m_tag[c0].detag;
#ifdef DEBUG		
		if (!(t0 > -1))
		{
			PRINTFinLEVEL("Coarsened box (%d,%d,%d) the received tag is %d but the local tag is %d", 
				finelevel->cur_level,
				finelevel->m_box[c0].ix(),
				finelevel->m_box[c0].iy(),
				finelevel->m_box[c0].iz(),
				recvdata[i].localcell, t0);
			MPI_Abort(MPI_COMM_WORLD,850);
		}
#endif		
		Assert(t0 > -1, "error in build new block pair!!!", 75);
		new_bkp[i].innode = t0;
		new_bkp[i].outnode.node = recvnode[i];
		new_bkp[i].outnode.index = recvdata[i].localcell;
#ifdef PASSAGE_ANGLE
		new_bkp[i].theta = finelevel->blockpair[i0].theta;
#endif		
	}
	blockpair.Addnew(new_bkp);
	SynNewGhostTag();
	blockpair.DirectlyReduceNew();
	level_pro_ghost.Addnew(new_level_pro_ghost);
	level_pro_ghost.DirectlyReduceNew();
	GiveLevelFlag("finish build ghost pro ghost", cur_level, 5);
}
