#include "Body.H"

void Body::ComptNewLocation()
{
	bodycenter = bodycenter + (inlineosc_new - inlineosc);
	Pointxyz dangle = rotosc_new - rotosc;
	int bs = allpoint.ps();
	int be = allpoint.pe();
	for (int i = bs; i < be; ++i)
	{
		int i0 = i - bs;
		surfptoff[i0].rotate_z(dangle[2]);
		allpoint[i] = surfptoff[i0] + bodycenter;
	}
	MPI_Barrier(share_comm);
	inlineosc = inlineosc_new;
	rotosc = rotosc_new;
	// if (nrank == 0) printf("the body xy is (%f,%f,%f)!!!\n", 
	// 	inlineosc[0],
	// 	inlineosc[1],
	// 	inlineosc[2]);
	ComptPatchParams();
}

void Body::ComptPatchParams()
{
	int bs = patch.ps();
	int be = patch.pe();
	for (int i = bs; i < be; ++i)
	{
		ComptPatchCenter(patch[i]);
		ComptPatchNv(patch[i]);
	}
	MPI_Barrier(share_comm);
}

	void Body::MovingBody_Image(Mesh & amesh, vector<Body> & abody, double & dt00)
	{	
		lastlevel = amesh.MyCurNum()-1;
		vector<vector<int> > move_from_in_to_out(lastlevel+1);
		vector<vector<int> > move_from_out_to_in(lastlevel+1);		
		gridisnew = true;
		int bn = abody.size();
		for (int i = 0; i < bn; ++i)
		{
			abody[i].ComptWallForce(amesh);
		}
		for (int i = 0; i < bn/2; ++i)
		{
			abody[i*2].force = abody[i*2].force + abody[i*2+1].force;
			abody[i*2+1].force = abody[i*2].force;
			abody[i*2].moment = abody[i*2].moment + abody[i*2+1].moment;
			abody[i*2+1].moment = abody[i*2].moment;
		}
#ifdef SHOWTIME		
		double move_start_time = MPI_Wtime();
#endif		
		MPI_Barrier(share_comm);
		//vector<double> celldis;
		//amesh.CopyIBCellDistance(celldis);
		//GiveAFlag("Finish CopyIBCellDistance!!!", 5);
#ifndef FRAME_VIB	
		for (int i = 0; i < bn; ++i)
		{
			abody[i].Motion_Rule(amesh, dt00);
			abody[i].ModifyBoxtoPatch(amesh);
		}
#ifdef SHOWTIME
		double change_moment = MPI_Wtime();	
		double change_post_time = change_moment - move_start_time;
#endif		
		GiveAFlag("Finish computing the new location of the cylinder!!!", 5);
		// RenewIBCellsDistance(amesh, abody, celldis);
		// GiveAFlag("Finish renewing the ib cells distance!!!", 5);
// #ifdef TURBULENCE		
		for (int i = 0; i < amesh.MyCurNum(); ++i)
		{
			RenewLevelDistance(amesh, abody, i, move_from_in_to_out[i], move_from_out_to_in[i]);
		}
#ifdef SHOWTIME
		double renew_moment = MPI_Wtime();		
		double renew_dis_time = renew_moment - change_moment;
#endif		
		for (int i = 0; i < amesh.MyCurNum(); ++i)
		{
			amesh.SynNormal_BlockDistance(i, move_from_in_to_out[i], move_from_out_to_in[i]);
		}
#ifdef SHOWTIME
		double renew_syn_moment = MPI_Wtime();		
		double renew_syn_time = renew_syn_moment - renew_moment;
#endif		
// #else
// 		RenewLevelDistance(amesh, abody, lastlevel, move_from_in_to_out[lastlevel], move_from_out_to_in[lastlevel]);
// 		amesh.SynNormal_BlockDistance(lastlevel, move_from_in_to_out[lastlevel], move_from_out_to_in[lastlevel]);
// #endif		
		GiveAFlag("Finish RenewLevelDistance!!!", 5);		
		//GiveAFlag("Finish SynNormal_BlockDistance!!!", 5);
		RemoveIBCells(amesh, abody, move_from_in_to_out[lastlevel], move_from_out_to_in[lastlevel]);
		GiveAFlag("Finish locating the removed ib cells!!!", 5);
		InfectNewIBCells(amesh, abody, move_from_in_to_out[lastlevel], move_from_out_to_in[lastlevel]);
		GiveAFlag("Finish infecting new ib cells!!!", 5);
		amesh.m_dis.holeplan(true);
		amesh.m_dis.CompressArray();
		amesh.RenewInfectIndex();
		GiveAFlag("Finish renewing the infected index!!!", 5);
		FindImageCell(amesh, abody);
		amesh.DataExchangeCells_Wall();
		FindImageCell_GhostIBCell(amesh, abody);
#ifdef SHOWTIME
		double ibcell_moment = MPI_Wtime();		
		double renew_ibcell_time = ibcell_moment - renew_syn_moment;
		if (nrank == 0) printf(">>>>>>>>>>In body move: change location %f renew distance %f syn %f renew ib cells %f\n", change_post_time, renew_dis_time, renew_syn_time, renew_ibcell_time);
#endif		
#else
		for (int i = 0; i < bn; ++i)
		{
			abody[i].ComptWallForce(amesh);	
		}
		for (int i = 0; i < bn; ++i)
		{
			abody[i].force = abody[0].force + abody[1].force;
			abody[i].Motion_Rule(amesh);
		}
#endif				
	}

	void Body::RenewLevelDistance(Mesh & amesh, vector<Body> & abody, const int & level0,
		vector<int> & move_from_in_to_out,
		vector<int> & move_from_out_to_in)
	{
		int bps, bpe;
		double critic_dis = 7.6;
		amesh.m_level[level0].m_box.GlobalOrder(bps, bpe);
		// bps = amesh.m_level[level0].m_box.ps();
		// bpe = amesh.m_level[level0].m_box.pe();
		for (int i = bps; i < bpe; ++i)
		{
			//由于在求hg cell时需要判断插值点是否为solid,因此需要同时求Normal cell和Dmghost距离边界的距离
			//if (amesh.m_level[level0].m_box[i].pair.signdis < critic_dis)
			if (amesh.m_level[level0].m_box[i].pair.signdis < critic_dis && (amesh.m_level[level0].m_box[i].type != Blockghost))
			//if (amesh.m_level[level0].m_box[i].pair.signdis < 3.0 && (amesh.m_level[level0].m_box[i].type != Blockghost))
			//double dx = amesh.m_level[level0].m_geom[i].boxcenter[0]-abody[0].bodycenter[0];
			//double dy = amesh.m_level[level0].m_geom[i].boxcenter[1]-abody[0].bodycenter[1];
			//if (dx*dx/(0.5*0.5)+(dy*dy)/(1.0*1.0) < 1.0 && (amesh.m_level[level0].m_box[i].type != Blockghost))
			{
				BoxtoWall & mypair = amesh.m_level[level0].m_box[i].pair;
				Assert(mypair.body > -1 && mypair.body < abody.size(), "The cell close body index error 117!!!", 5);
				Assert(mypair.patch > -1 && mypair.patch < abody[mypair.body].patch.size(), "The cell close body patch index error 118!!!", 5);
				Pointxyz & cellcenter = amesh.m_level[level0].m_geom[i].boxcenter;
				//mypair.distance = (cellcenter - abody[mypair.body].patch[mypair.patch].pc).length();
				mypair.distance = 999.0;
				abody[mypair.body].FindClosePatchtoCell(mypair, cellcenter);
			}
		}
		MPI_Barrier(share_comm);
		vector<int> cycle_box;
		vector<int> cycle_body;
		vector<int> cycle_patch;
		for (int i = bps; i < bpe; ++i)
		{
			//由于在求hg cell时需要判断插值点是否为solid,因此需要同时求Normal cell和Dmghost距离边界的距离
			if (amesh.m_level[level0].m_box[i].pair.signdis < critic_dis && (amesh.m_level[level0].m_box[i].type != Blockghost))
			{
				BoxtoWall & mypair = amesh.m_level[level0].m_box[i].pair;
				Pointxyz & cellcenter = amesh.m_level[level0].m_geom[i].boxcenter;
#if DIM ==3
				for (Point_iterator p(0,3); p.end(); ++p)
#elif DIM == 2
				for (Point_iterator_2d p(0,3); p.end(); ++p)
#endif							
				{
#if DIM == 3				
					int an0 = amesh.m_level[level0].m_box[i].neib[p.i][p.j][p.k];
#elif DIM == 2
					int an0 = amesh.m_level[level0].m_box[i].neib[p.i][p.j][1];
#endif
					if (an0 > -1)
					{
						if (amesh.m_level[level0].m_box[an0].pair.body != amesh.m_level[level0].m_box[i].pair.body)
						{
							cycle_box.push_back(i);
							cycle_body.push_back(amesh.m_level[level0].m_box[an0].pair.body);
							cycle_patch.push_back(amesh.m_level[level0].m_box[an0].pair.patch);
						}
					}
				}
			}
		}
		MPI_Barrier(share_comm);
		int cycle_num = cycle_box.size();
		BoxtoWall neibtw;
		for (int i = 0; i < cycle_num; ++i)
		{
			int b0 = cycle_box[i];
			BoxtoWall & mypair = amesh.m_level[level0].m_box[b0].pair;
			if (mypair.body != cycle_body[i])
			{
				neibtw.body = cycle_body[i];
				neibtw.patch = cycle_patch[i];
				neibtw.distance = 999.0;
				abody[neibtw.body].FindClosePatchtoCell(neibtw, amesh.m_level[level0].m_geom[b0].boxcenter);
				bool changeflag = false;
				if (neibtw.signdis > 0.0 && mypair.signdis < 0.0)
				{
					changeflag = true;
				}
				else if (neibtw.distance < mypair.distance)
				{
					if (mypair.signdis*neibtw.signdis > 0.0)
					{
						changeflag = true;
					}
					else if (neibtw.signdis > 0.0)
					{
						changeflag = true;
					}
				}
				if (changeflag)
				{
					mypair.body = neibtw.body;
					mypair.patch = neibtw.patch;
					mypair.distance = neibtw.distance;
					mypair.signdis = neibtw.signdis;
					mypair.distance_to_body = neibtw.distance_to_body;
				}
			}
		}
		MPI_Barrier(share_comm);
		for (int i = bps; i < bpe; ++i)
		{
			//由于在求hg cell时需要判断插值点是否为solid,因此需要同时求Normal cell和Dmghost距离边界的距离
			if (amesh.m_level[level0].m_box[i].pair.signdis < critic_dis && (amesh.m_level[level0].m_box[i].type != Blockghost))
			{
				BoxtoWall & mypair = amesh.m_level[level0].m_box[i].pair;
 				mypair.signdis = min(mypair.signdis, mypair.distance_to_dm);
				if (mypair.signdis < 0.0)
				{
					if (!amesh.m_level[level0].m_box[i].solid && amesh.m_level[level0].m_box[i].type != Dmghost)
					{
						move_from_out_to_in.push_back(i);
					}
					amesh.m_level[level0].m_box[i].solid = true;
				}
				else
				{
					if (amesh.m_level[level0].m_box[i].solid && amesh.m_level[level0].m_box[i].type != Dmghost)
					{
						move_from_in_to_out.push_back(i);
					}
					amesh.m_level[level0].m_box[i].solid = false;
				}
			}
			// Pointxyz & boxbc = amesh.m_level[level0].m_geom[i].boxcenter;
			// Pointxyz dxyz = boxbc - abody[0].bodycenter;
// 			if (abs(dxyz.length_2d()-0.5 - amesh.m_level[level0].m_box[i].pair.signdis) > 0.2 && level0 >= max_mesh_level -3)
// 			{
// 				int bd0 = amesh.m_level[level0].m_box[i].pair.body;
// 				int pat0 = amesh.m_level[level0].m_box[i].pair.patch;
// 				printf("N%dL%dBox(%d,%d,%d)(%f,%f,%f) close body %d patch %d (%f,%f,%f) distance %20.14f signdis %f\n",
// 									node,level0,
// 									amesh.m_level[level0].m_box[i].ix(),
// 									amesh.m_level[level0].m_box[i].iy(),
// 									amesh.m_level[level0].m_box[i].iz(),
// 									boxbc[0], boxbc[1], boxbc[2],
// 									amesh.m_level[level0].m_box[i].pair.body,
// 									amesh.m_level[level0].m_box[i].pair.patch,
// 									abody[bd0].patch[pat0].pc[0],
// 									abody[bd0].patch[pat0].pc[1],
// 									abody[bd0].patch[pat0].pc[2],
// 									amesh.m_level[level0].m_box[i].pair.distance,
// 									amesh.m_level[level0].m_box[i].pair.signdis);
// 				for (int pd = 0; pd < abody[bd0].patneib[pat0].nbnum; ++pd)
// 				{
// 					int pnb = abody[bd0].patneib[pat0].nb[pd];
// #ifndef PASSAGE_ANGLE		
// 					Pointxyz dxyz = amesh.m_level[level0].m_geom[i].boxcenter - abody[bd0].patch[pnb].pc;
// 					PeriodicLength(dxyz);
// #else
// 					Pointxyz newbcxyz, dxyz;
// 					PeriodicAnnulaLength(amesh.m_level[level0].m_geom[i].boxcenter, abody[bd0].patch[pnb].pc, newbcxyz);
// 					dxyz = newbcxyz - abody[bd0].patch[pnb].pc;
// #endif
// 					double dis_pnb = dxyz.length();
// 					int lessthanthepair = (dis_pnb < amesh.m_level[level0].m_box[i].pair.distance);
// 					printf("Body%dPatch%d has %d neibs, the %d one is %d (%f,%f,%f) box dis to the patch is %20.14f lessthanpair %d!!!\n",
// 						bd0, pat0, abody[bd0].patneib[pat0].nbnum, pd,
// 						pnb,
// 						abody[bd0].patch[pnb].pc[0],
// 						abody[bd0].patch[pnb].pc[1],
// 						abody[bd0].patch[pnb].pc[2],
// 						dis_pnb, lessthanthepair);
// 				}
// 				MPI_Abort(MPI_COMM_WORLD, 140);
// 			}
		}
// 		for (int i = bps; i < bpe; ++i)
// 		{
// 			//由于在求hg cell时需要判断插值点是否为solid,因此需要同时求Normal cell和Dmghost距离边界的距离
// 			if (amesh.m_level[level0].m_box[i].pair.signdis < 3.0 && (amesh.m_level[level0].m_box[i].type != Blockghost))
// 			{
// 				BoxtoWall & mypair = amesh.m_level[level0].m_box[i].pair;
// 				Pointxyz & cellcenter = amesh.m_level[level0].m_geom[i].boxcenter;
// #if DIM ==3
// 				for (Point_iterator p(0,3); p.end(); ++p)
// #elif DIM == 2
// 				for (Point_iterator_2d p(0,3); p.end(); ++p)
// #endif							
// 				{
// #if DIM == 3				
// 					int an0 = amesh.m_level[level0].m_box[i].neib[p.i][p.j][p.k];
// #elif DIM == 2
// 					int an0 = amesh.m_level[level0].m_box[i].neib[p.i][p.j][1];
// #endif
// 					if (an0 > -1)
// 					{
// 						if (amesh.m_level[level0].m_box[an0].pair.body != amesh.m_level[level0].m_box[i].pair.body)
// 						{
// 							BoxtoWall neibtw = amesh.m_level[level0].m_box[an0].pair;
// 							Assert(neibtw.body > -1 && neibtw.body < abody.size(), "The cell close body index error 139!!!", 5);
// 							if (!(neibtw.patch > -1 && neibtw.patch < abody[neibtw.body].patch.size()))
// 							{
// 								printf("The cell close body patch index error!!! N%dL%dBox(%d,%d,%d) close body %d patch %d (%f,%f,%f) but its neib box (%d,%d,%d) close body %d patch %d\n",
// 									node,level0,
// 									amesh.m_level[level0].m_box[i].ix(),
// 									amesh.m_level[level0].m_box[i].iy(),
// 									amesh.m_level[level0].m_box[i].iz(),
// 									amesh.m_level[level0].m_box[i].pair.body,
// 									amesh.m_level[level0].m_box[i].pair.patch,
// 									abody[mypair.body].patch[mypair.patch].pc[0],
// 									abody[mypair.body].patch[mypair.patch].pc[1],
// 									abody[mypair.body].patch[mypair.patch].pc[2],
// 									amesh.m_level[level0].m_box[an0].ix(),
// 									amesh.m_level[level0].m_box[an0].iy(),
// 									amesh.m_level[level0].m_box[an0].iz(),
// 									amesh.m_level[level0].m_box[an0].pair.body,
// 									amesh.m_level[level0].m_box[an0].pair.patch);
// 								MPI_Abort(MPI_COMM_WORLD, 140);
// 							}
// 							Assert(neibtw.patch > -1 && neibtw.patch < abody[neibtw.body].patch.size(), "The cell close body patch index error 140!!!", 5);	
// 							neibtw.distance = 999.0;			
// 							abody[neibtw.body].FindClosePatchtoCell(neibtw, cellcenter);
// 							if (neibtw.distance < mypair.distance)
// 							{
// 								bool changeflag = false;
// 								if (mypair.signdis*neibtw.signdis > 0.0)
// 								{
// 									changeflag = true;
// 								}
// 								else if (neibtw.signdis > 0.0)
// 								{
// 									changeflag = true;
// 								}
// 								if (changeflag)
// 								{
// 									mypair.body = neibtw.body;
// 									mypair.patch = neibtw.patch;
// 									mypair.distance = neibtw.distance;
// 									mypair.signdis = neibtw.signdis;
// 									break;
// 								}
// 							}
// 						}
// 					}								
// 				}
// 				mypair.signdis = min(mypair.signdis, mypair.distance_to_dm);
// 				if (mypair.signdis < 0.0)
// 				{
// 					if (!amesh.m_level[level0].m_box[i].solid && amesh.m_level[level0].m_box[i].type != Dmghost)
// 					{
// 						move_from_out_to_in.push_back(i);
// 					}
// 					amesh.m_level[level0].m_box[i].solid = true;
// 				}
// 				else
// 				{
// 					if (amesh.m_level[level0].m_box[i].solid && amesh.m_level[level0].m_box[i].type != Dmghost)
// 					{
// 						move_from_in_to_out.push_back(i);
// 					}
// 					amesh.m_level[level0].m_box[i].solid = false;
// 				}
// 			}
// 		}
// 		MPI_Barrier(share_comm);
		for (int i = bps; i < bpe; ++i)
		{
			if (amesh.m_level[level0].m_box[i].ix() == 564 &&
				amesh.m_level[level0].m_box[i].iy() == -2 &&
				amesh.m_level[level0].m_box[i].iz() == 85)
			{
				printf("N%d Box (403,-3,79) after renew close body %d patch %d signdis %f\n",node,
					amesh.m_level[level0].m_box[i].pair.body,
					amesh.m_level[level0].m_box[i].pair.patch,
					amesh.m_level[level0].m_box[i].pair.signdis);
			}
		}
	}

	void Body::RemoveIBCells(Mesh & amesh, vector<Body> & abody, 
		vector<int> & move_from_in_to_out,
		vector<int> & move_from_out_to_in)
	{
		int in_to_outnum = move_from_in_to_out.size();
		for (int i = 0; i < in_to_outnum; ++i)
		{
			int b0 = move_from_in_to_out[i];
#if DIM == 3				
			for (Point_iterator p(0,3); p.end(); ++p)
			{
				int an0 = amesh.m_level[lastlevel].m_box[b0].neib[p.i][p.j][p.k];
#elif DIM == 2
      		for (Point_iterator_2d p(0,3); p.end(); ++p)
      		{
      			int an0 = amesh.m_level[lastlevel].m_box[b0].neib[p.i][p.j][1];
#endif      								
				if (an0 > -1)
				{
					if (!amesh.m_level[lastlevel].m_box[an0].solid && amesh.m_level[lastlevel].m_box[an0].type != Dmghost)
					{
#if DIM == 3				
						for (Point_iterator q(0,3); q.end(); ++q)
						{
							int an1 = amesh.m_level[lastlevel].m_box[an0].neib[q.i][q.j][q.k];
#elif DIM == 2
      					for (Point_iterator_2d q(0,3); q.end(); ++q)
      					{
      						int an1 = amesh.m_level[lastlevel].m_box[an0].neib[q.i][q.j][1];
#endif
							if (an1 > -1)
							{
								if (amesh.infectbox[an1] > -1 && amesh.m_level[lastlevel].m_box[an1].type != Dmghost)
								{
									// Assert((!amesh.m_level[lastlevel].m_box[an1].solid && amesh.m_level[lastlevel].m_box[an1].type != Dmghost),
									// 	"IB cell Error 200!!!", 200);
// #ifdef DEBUG
// 									if (amesh.m_level[lastlevel].m_box[an1].solid || 
// 										 amesh.m_level[lastlevel].m_box[an1].type == Dmghost)
// 									{
// 										printf("N%d IB Box(%d,%d,%d) signdis %f neib's neib (%d,%d,%d) signdis %f\n",
// 											node, amesh.m_level[lastlevel].m_box[b0].ix(),
// 											amesh.m_level[lastlevel].m_box[b0].iy(),
// 											amesh.m_level[lastlevel].m_box[b0].iz(),
// 											amesh.m_level[lastlevel].m_box[b0].pair.signdis,
// 											amesh.m_level[lastlevel].m_box[an1].ix(),
// 											amesh.m_level[lastlevel].m_box[an1].iy(),
// 											amesh.m_level[lastlevel].m_box[an1].iz(),
// 											amesh.m_level[lastlevel].m_box[an1].pair.signdis);
// 										MPI_Abort(MPI_COMM_WORLD, 213);
// 									}

// #endif																		
#if DIM == 3				
									for (Point_iterator p1(0,3); p1.end(); ++p1)
									{
										int an2 = amesh.m_level[lastlevel].m_box[an1].neib[p1.i][p1.j][p1.k];
#elif DIM == 2
      								for (Point_iterator_2d p1(0,3); p1.end(); ++p1)
      								{
      									int an2 = amesh.m_level[lastlevel].m_box[an1].neib[p1.i][p1.j][1];
#endif      									
										if (an2 > -1)
										{
											if (amesh.m_level[lastlevel].m_box[an2].type != Dmghost)
											{
												if (amesh.m_level[lastlevel].m_box[an2].solid)
												{
													goto NEXT_an1;
												}
#if DIM == 3				
												for (Point_iterator p2(0,3); p2.end(); ++p2)
												{
													int an3 = amesh.m_level[lastlevel].m_box[an2].neib[p2.i][p2.j][p2.k];
#elif DIM == 2
      											for (Point_iterator_2d p2(0,3); p2.end(); ++p2)
      											{
      												int an3 = amesh.m_level[lastlevel].m_box[an2].neib[p2.i][p2.j][1];
#endif
													if (an3 > -1)
													{
														if (amesh.m_level[lastlevel].m_box[an3].type != Dmghost && 
															amesh.m_level[lastlevel].m_box[an3].solid)
														{
															goto NEXT_an1;
														}
													}											
												}
											}
										}
									}
									amesh.infectbox[an1] = -3;
									NEXT_an1:;
								}
							}
						}
					}
				}
			}
		}
		MPI_Barrier(share_comm);
		int mps = amesh.m_dis.ps();
		int mpe = amesh.m_dis.pe();
		for (int i = mps; i < mpe; ++i)
		{
			int ci0 = amesh.m_dis[i].ci;
			if (amesh.infectbox[ci0] == -3)
			{
				amesh.m_dis.givehole(i);
				amesh.infectbox[ci0] = i;
			}
		}
		int out_to_innum = move_from_out_to_in.size();
		for (int i = 0; i < out_to_innum; ++i)
		{
			int i0 = move_from_out_to_in[i];
			if (amesh.infectbox[i0] < 0)
			{
				bool caseistrue = false;
				if (amesh.m_level[lastlevel].m_box[i0].type == Blockghost)
				{
					if (!amesh.InTransferRange(lastlevel, i0))
					{
						caseistrue = true;
					}
					else
					{
						for (Point_iterator pd(0,3); pd.end(); ++pd)
						{
							int anpd = amesh.m_level[lastlevel].m_box[i0].neib[pd.i][pd.j][pd.k];
							if (anpd == -1)
							{
								caseistrue = true;
								break;
							}
						}
					}
				}
				if (!caseistrue)
				{
					printf("N%d A box (%d,%d,%d) (%f,%f,%f) signdis %f close body %d moved from outside into body is not infected!!! Infected index is %d\n", 
						node,
						amesh.m_level[lastlevel].m_box[i0].ix(),
						amesh.m_level[lastlevel].m_box[i0].iy(),
						amesh.m_level[lastlevel].m_box[i0].iz(),
						amesh.m_level[lastlevel].m_geom[i0].boxcenter[0],
						amesh.m_level[lastlevel].m_geom[i0].boxcenter[1],
						amesh.m_level[lastlevel].m_geom[i0].boxcenter[2],
						amesh.m_level[lastlevel].m_box[i0].pair.signdis,
						amesh.m_level[lastlevel].m_box[i0].pair.body,
						amesh.infectbox[i0]);
					MPI_Abort(MPI_COMM_WORLD, 263);
				}			
			}
			// else
			if (amesh.infectbox[i0] > -1)
			{
				amesh.m_dis.givehole(amesh.infectbox[i0]);
			}
		}
		MPI_Barrier(share_comm);
		int hnum = amesh.m_dis.hole.size();
		MPI_Allreduce(MPI_IN_PLACE, &hnum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#ifdef DEBUG		
		if (nrank == 0)
		{
			printf("All blocks removed m_dis cells number is %d\n", hnum);
		}
#endif		
		// for (int i = 0; i < amesh.m_dis.hole.size(); ++i)
		// {
		// 	int i0 = amesh.m_dis.hole[i];
		// 	int ci0 = amesh.m_dis[i0].ci;
		// 	printf("N%dR%d remove ib cell %d is box (%d,%d,%d)(%f,%f,%f) signdis to wall is %f\n", node, srank,
		// 		i, 
		// 		amesh.m_level[lastlevel].m_box[ci0].ix(),
		// 		amesh.m_level[lastlevel].m_box[ci0].iy(),
		// 		amesh.m_level[lastlevel].m_box[ci0].iz(),
		// 		amesh.m_level[lastlevel].m_geom[ci0].boxcenter[0],
		// 		amesh.m_level[lastlevel].m_geom[ci0].boxcenter[1],
		// 		amesh.m_level[lastlevel].m_geom[ci0].boxcenter[2],
		// 		amesh.m_level[lastlevel].m_box[ci0].pair.signdis);
		// }			
	}

	void Body::InfectNewIBCells(Mesh & amesh, vector<Body> & abody,
		vector<int> & move_from_in_to_out,
		vector<int> & move_from_out_to_in)
	{
		int in_to_outnum = move_from_in_to_out.size();
		for (int i = 0; i < in_to_outnum; ++i)
		{
			bool infectflag = false;
			int i0 = move_from_in_to_out[i];
			if (amesh.infectbox[i0] > -1)
			{
				printf("A old negative-distance cell was infected!!!\n");
				MPI_Abort(MPI_COMM_WORLD, 289);
			}
#if DIM == 3			
			for (Point_iterator p(0,3); p.end(); ++p)
			{
				int an0 = amesh.m_level[lastlevel].m_box[i0].neib[p.i][p.j][p.k];
#else      
			for (Point_iterator_2d p(0,3); p.end(); ++p)
			{
				int an0 = amesh.m_level[lastlevel].m_box[i0].neib[p.i][p.j][1];
#endif								
				if (an0 > -1)
				{
					if (amesh.m_level[lastlevel].m_box[an0].type != Dmghost)
					{
						if (amesh.m_level[lastlevel].m_box[an0].pair.signdis < 0.0)
						{
							amesh.InfectABox(move_from_in_to_out[i]);
							infectflag = true;
							break;
						}
					}
				}
			}
// 			if (!infectflag)
// 			{
// #if DIM == 3			
// 				for (Point_iterator p(0,3); p.end(); ++p)
// 				{
// 					int an0 = amesh.m_level[lastlevel].m_box[i0].neib[p.i][p.j][p.k];
// #else      
// 				for (Point_iterator_2d p(0,3); p.end(); ++p)
// 				{
// 					int an0 = amesh.m_level[lastlevel].m_box[i0].neib[p.i][p.j][1];
// #endif								
// 					if (an0 > -1)
// 					{
// 						if (amesh.m_level[lastlevel].m_box[an0].type != Dmghost)
// 						{
// #if DIM == 3			
// 							for (Point_iterator q(0,3); q.end(); ++q)
// 							{
// 								int an1 = amesh.m_level[lastlevel].m_box[i0].neib[q.i][q.j][q.k];
// #else      
// 							for (Point_iterator_2d q(0,3); q.end(); ++q)
// 							{
// 								int an1 = amesh.m_level[lastlevel].m_box[i0].neib[q.i][q.j][1];
// #endif
// 								if (an1 > -1)
// 								{
// 									if (amesh.m_level[lastlevel].m_box[an1].type != Dmghost)
// 									{
// 										if (amesh.m_level[lastlevel].m_box[an1].pair.signdis < 0.0)
// 										{
// 											amesh.InfectABox(move_from_in_to_out[i]);
// 											infectflag = true;
// 											goto NEXTCELL;
// 										}
// 									}
// 								}
// 							}
// 						}
// 					}
// 				}
// 				NEXTCELL:;
// 			}
			if (!infectflag)
			{
				bool caseistrue = false;
				if (amesh.m_level[lastlevel].m_box[i0].type == Blockghost)
				{
					if (!amesh.InTransferRange(lastlevel, i0))
					{
						caseistrue = true;
					}
					else
					{
						for (Point_iterator pd(0,3); pd.end(); ++pd)
						{
							int anpd = amesh.m_level[lastlevel].m_box[i0].neib[pd.i][pd.j][pd.k];
							if (anpd == -1)
							{
								caseistrue = true;
								break;
							}
						}
					}
				}
				if (!caseistrue)
				{
					printf("A old negative-distance cell does not have a negative-distance neib. Maybe the move step is too large!!!\n");
					int b0 = amesh.m_level[lastlevel].m_box[i0].pair.body;
					int p0 = amesh.m_level[lastlevel].m_box[i0].pair.patch;
					printf("The box is (%d,%d,%d) close body %d close patch %d (%f,%f,%f) distance %f signdis %f\n",
						amesh.m_level[lastlevel].m_box[i0].ix(),
						amesh.m_level[lastlevel].m_box[i0].iy(),
						amesh.m_level[lastlevel].m_box[i0].iz(),
						b0, p0, 
						abody[b0].patch[p0].pc[0],
						abody[b0].patch[p0].pc[1],
						abody[b0].patch[p0].pc[2],
						amesh.m_level[lastlevel].m_box[i0].pair.distance,
						amesh.m_level[lastlevel].m_box[i0].pair.signdis);
					MPI_Abort(MPI_COMM_WORLD, 314);
				}
			}
		}
		int out_to_innum = move_from_out_to_in.size();
		for (int i = 0; i < out_to_innum; ++i)
		{
			int i0 = move_from_out_to_in[i];
#if DIM == 3				
			for (Point_iterator p(0,3); p.end(); ++p)
			{
				int an1 = amesh.m_level[lastlevel].m_box[i0].neib[p.i][p.j][p.k];
#elif DIM == 2	
			for (Point_iterator_2d p(0,3); p.end(); ++p)
			{
				int an1 = amesh.m_level[lastlevel].m_box[i0].neib[p.i][p.j][1];
#endif								
				if (an1 > -1)
				{
					if (amesh.m_level[lastlevel].m_box[an1].type != Dmghost && amesh.m_level[lastlevel].m_box[an1].pair.signdis >= 0.0)
					{
						// if (amesh.infectbox[an1] == -1)
						// {
							
#if DIM == 3				
							for (Point_iterator p1(0,3); p1.end(); ++p1)
							{
								int an2 = amesh.m_level[lastlevel].m_box[an1].neib[p1.i][p1.j][p1.k];
#elif DIM == 2	
							for (Point_iterator_2d p1(0,3); p1.end(); ++p1)
							{
								int an2 = amesh.m_level[lastlevel].m_box[an1].neib[p1.i][p1.j][1];
#endif	
								if (an2 > -1)
								{
									if (amesh.m_level[lastlevel].m_box[an2].type != Dmghost && amesh.m_level[lastlevel].m_box[an2].pair.signdis >= 0.0)
									{
										if (amesh.infectbox[an2] < 0)
										{
											Assert(amesh.infectbox[an2] == -1 || amesh.infectbox[an2] == -2, "The box infected index can only be -1", 360);
											amesh.InfectABox(an2);
										}
									}
								}
							}
						// }
					}
				}
			}
		}
		MPI_Barrier(share_comm);
		int oldinfectnum = amesh.m_dis.size()-1;
		amesh.CountInfectedBox();
		int ifs = amesh.infectbox.ps();
		int ife = amesh.infectbox.pe();
		for (int i = ifs; i < ife; ++i)
		{
			if (amesh.infectbox[i] > oldinfectnum)
			{
				BoxtoWall & abtw = amesh.m_level[lastlevel].m_box[i].pair;
				int i0 = amesh.infectbox[i];
				amesh.m_dis[i0].ci = i;
				amesh.m_dis[i0].tp = amesh.m_level[lastlevel].m_geom[i].boxcenter-abody[abtw.body].patch[abtw.patch].nv*abtw.signdis;
				amesh.m_dis[i0].hgc.closecell = i;
				amesh.m_dis[i0].hgc.InitIntparray();
			}
		}
		MPI_Barrier(share_comm);
	}
