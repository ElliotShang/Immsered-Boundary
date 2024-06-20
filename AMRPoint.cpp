#include "AMR.H"
	/*Write in 2020-1-9*/
	/*This function*/
	/*Results in multi processor may give multi new points for the same box conner*/
	/*Consider the loop order, this kind of possibility exists, but appears to be small*/
	/*Although the newbox is used as the input, it mush be included in the box array*/
	/*When this function is executed*/
	void AMR::MakeNewPoints(vector<int> & newbox, 
										 const int & aimlevel, 
										 vector<newpointloc> & newpts)
	{
		newpts.resize(0);
		int neibpts[4] = {1,0,1,0};
		PRINTFinLEVEL("new box size is %d", aimlevel, (int)newbox.size());
		for (int i = 0; i < newbox.size(); ++i)
		{
			int sonbox = newbox[i];
			/*Loop for the [2][2][2] corner points of the box*/
			for (Point_iterator q(0,2); q.end(); ++q)
			{
				int ptindex = m_level[aimlevel].m_box[sonbox].pts[q.i][q.j][q.k];
				// printf("Rank is %d refine box is %d point(%d,%d,%d) is %d\n", nrank, sonbox.index,q.i,q.j,q.k,ptindex);
				if (ptindex == -1)
				{
					bool ptflag = true;
					for (int ni = q.i; ni <= q.i+1; ++ni)
					{
						for (int nj = q.j; nj <= q.j+1; ++nj)
						{
#if DIM == 3							
							for (int nk = q.k; nk <= q.k+1; ++nk)
							{
#else 	
								int nk = 1;
#endif								
								int ptneib = m_mesh.BoxNeib(aimlevel, sonbox, ni, nj, nk);
								if (ptneib > -1)
								{
#ifdef DEBUG
									int ptx = m_mesh.m_level[aimlevel].m_box[ptneib].ix();
									int pty = m_mesh.m_level[aimlevel].m_box[ptneib].iy();
									int ptz = m_mesh.m_level[aimlevel].m_box[ptneib].iz();
									if ((m_level[aimlevel].m_box[sonbox].ix() - ptx != 1-ni && 
											 	m_level[aimlevel].m_box[sonbox].ix() + ptx != m_mesh.Boxnum()[0]*pow(2,aimlevel)-1) ||
											(m_level[aimlevel].m_box[sonbox].iy() - pty != 1-nj && 
											 	m_level[aimlevel].m_box[sonbox].iy() + pty != m_mesh.Boxnum()[1]*pow(2,aimlevel)-1) ||
											(m_level[aimlevel].m_box[sonbox].iz() - ptz != 1-nk && 
											 	m_level[aimlevel].m_box[sonbox].iz() + ptz != m_mesh.Boxnum()[2]*m_level[aimlevel].power_ratio-1))
									{
										PRINTFinLEVEL("Error: When creating point, box %d (%d,%d,%d) neib (%d,%d,%d) is %d (%d,%d,%d)",
											aimlevel, sonbox,
											m_level[aimlevel].m_box[sonbox].ix(),
											m_level[aimlevel].m_box[sonbox].iy(),
											m_level[aimlevel].m_box[sonbox].iz(),
											ni,nj,nk, ptneib,
											ptx,pty,ptz);
										MPI_Abort(MPI_COMM_WORLD, 583);
									}
#endif
#if DIM == 3									
									if (abs(m_mesh.LevelBox(aimlevel, ptneib).ix()-m_mesh.LevelBox(aimlevel, sonbox).ix()) > 1 ||
										abs(m_mesh.LevelBox(aimlevel, ptneib).iy()-m_mesh.LevelBox(aimlevel, sonbox).iy()) > 1 ||
										abs(m_mesh.LevelBox(aimlevel, ptneib).iz()-m_mesh.LevelBox(aimlevel, sonbox).iz()) > 1)
#else 
									if (abs(m_mesh.LevelBox(aimlevel, ptneib).ix()-m_mesh.LevelBox(aimlevel, sonbox).ix()) > 1 ||
										abs(m_mesh.LevelBox(aimlevel, ptneib).iy()-m_mesh.LevelBox(aimlevel, sonbox).iy()) > 1)
#endif										
									{}
									else
									{
										int neib_pti = m_level[aimlevel].m_box[ptneib].
											pts[neibpts[q.i+ni]][neibpts[q.j+nj]][neibpts[q.k+nk]];
										if (neib_pti > -1)
										{
											m_level[aimlevel].m_box[sonbox].pts[q.i][q.j][q.k] = neib_pti;
											ptflag = false;
											goto FOUNDAPOINT;
										}
									}
								}
#if DIM == 3								
							}
#endif							
						}
					}
					FOUNDAPOINT:;
					if (ptflag)
					{
						m_level[aimlevel].m_box[sonbox].pts[q.i][q.j][q.k] = 0;
						newpts.push_back(newpointloc(sonbox, q.i, q.j, q.k, i));
					}
				}
			}
		}
		int ptnum_procs = newpts.size();
		vector<int> newptsnum(sprocs,0);
		vector<int> ptnum_arrstart(sprocs,0);
		int totnewptnum;
		MPI_Barrier(share_comm);
		//ShowAllRankData("New point number", ptnum_procs, 5);
		MPI_Allgather(&ptnum_procs, 1, MPI_INT, &newptsnum[0], 1, MPI_INT, share_comm);
		CountTotalNum(newptsnum, totnewptnum);
		ArrayProcsStart(newptsnum, ptnum_arrstart);
		int oldtotptnum = m_mesh.m_level[aimlevel].m_point.size();
		for (int ptn = 0; ptn < ptnum_procs; ++ptn)
		{
			newpts[ptn].arrindex = ptn+oldtotptnum+ptnum_arrstart[srank];
			for (int ni = newpts[ptn].boxpointx; ni <= newpts[ptn].boxpointx+1; ++ni)
			{
				for (int nj = newpts[ptn].boxpointy; nj <= newpts[ptn].boxpointy+1; ++nj)
				{
#if DIM == 3					
					for (int nk = newpts[ptn].boxpointz; nk <= newpts[ptn].boxpointz+1; ++nk)
					{
#else 
						int nk = 1;
#endif						
						int ptneib = m_mesh.BoxNeib(aimlevel, newpts[ptn].boxindex, ni, nj, nk);
						if (ptneib > -1)
						{
#if DIM == 3							
							if (abs(m_mesh.LevelBox(aimlevel,ptneib).ix()-m_mesh.LevelBox(aimlevel, newpts[ptn].boxindex).ix()) > 1 ||
								abs(m_mesh.LevelBox(aimlevel,ptneib).iy()-m_mesh.LevelBox(aimlevel, newpts[ptn].boxindex).iy()) > 1 ||
								abs(m_mesh.LevelBox(aimlevel,ptneib).iz()-m_mesh.LevelBox(aimlevel, newpts[ptn].boxindex).iz()) > 1)
#else 
							if (abs(m_mesh.LevelBox(aimlevel,ptneib).ix()-m_mesh.LevelBox(aimlevel, newpts[ptn].boxindex).ix()) > 1 ||
								abs(m_mesh.LevelBox(aimlevel,ptneib).iy()-m_mesh.LevelBox(aimlevel, newpts[ptn].boxindex).iy()) > 1)
#endif																
							{}
							else
							{
								m_level[aimlevel].m_box[ptneib].
									pts[neibpts[newpts[ptn].boxpointx+ni]]
								   	[neibpts[newpts[ptn].boxpointy+nj]]
								   	[neibpts[newpts[ptn].boxpointz+nk]] = newpts[ptn].arrindex;
								// printf("Level %d Box %d point(%d,%d,%d) is %d\n", aimlevel, ptneib,
								// 	neibpts[newpts[ptn].boxpointx+ni],neibpts[newpts[ptn].boxpointy+nj], 
								// 	neibpts[newpts[ptn].boxpointz+nk],newpts[ptn].arrindex);
							}
						}
#if DIM == 3
					}
#endif					
				}
			}
		}
		MPI_Barrier(share_comm);
	}

	void AMR::ManageLevelPoints(const int & ilevel)
	{
		/*To compute the new point coordinate*/
		int pt_range_s[3] = {0,0,1};
		int pt_range_e[3] = {1,2,2};
		/*-------------------------------------------*/
		if (NULL != m_level[ilevel].coarselevel)
		{
			MakeNewPoints(m_level[ilevel].newcoarsebox, ilevel-1, m_level[ilevel].newcoarsept);
			GiveAFlag("Create points for coarse level", 5);
			if (!m_level[ilevel].twodflag)
			{
				for (int cpt = 0; cpt < m_level[ilevel].newcoarsept.size(); ++cpt)
				{
					int dfbox0 = m_level[ilevel].newcoarsept[cpt].vectindex;
					int i0 = m_level[ilevel].derfbox[dfbox0];
					m_level[ilevel].tocoarse.push_back(mPoint());
					m_level[ilevel].tocoarse.back().index = m_level[ilevel].newcoarsept[cpt].arrindex;
					ComptDerefinePtxyz(ilevel, i0, m_level[ilevel].tocoarse.back().xyz, m_level[ilevel].newcoarsept[cpt].boxpointx,
						m_level[ilevel].newcoarsept[cpt].boxpointy, m_level[ilevel].newcoarsept[cpt].boxpointz);
				}
			}
			else
			{
				for (int cpt = 0; cpt < m_level[ilevel].newcoarsept.size(); ++cpt)
				{
					int dfbox0 = m_level[ilevel].newcoarsept[cpt].vectindex;
					int i0 = m_level[ilevel].derfbox[dfbox0];
					m_level[ilevel].tocoarse.push_back(mPoint());
					m_level[ilevel].tocoarse.back().index = m_level[ilevel].newcoarsept[cpt].arrindex;
					ComptDerefinePtxyz_2d(ilevel, i0, m_level[ilevel].tocoarse.back().xyz, m_level[ilevel].newcoarsept[cpt].boxpointx,
						m_level[ilevel].newcoarsept[cpt].boxpointy, m_level[ilevel].newcoarsept[cpt].boxpointz);
				}
			}
			m_mesh.m_level[ilevel-1].m_point.Addnew(m_level[ilevel].tocoarse);
			m_mesh.m_level[ilevel-1].m_point.DirectlyReduceNew();
		}
		
		if (NULL != m_level[ilevel].finelevel)
		{
			MakeNewPoints(m_level[ilevel].newrefinebox, ilevel+1, m_level[ilevel].newfinept);
			GiveAFlag("Create points for fine level", 5);
			if (!m_level[ilevel].finelevel->twodflag)
			{
				for (int fpt = 0; fpt < m_level[ilevel].newfinept.size(); ++fpt)
				{
					m_level[ilevel].tofine.push_back(mPoint());
					m_level[ilevel].tofine.back().index = m_level[ilevel].newfinept[fpt].arrindex;
					int sonid = m_level[ilevel].newfinept[fpt].vectindex%8;
					int momid = (m_level[ilevel].newfinept[fpt].vectindex-sonid)/8;
					int sonz = sonid%2;
					int sony = ((sonid - sonz)/2)%2;
					int sonx = (sonid-sony*2-sonz)/4;
					int i0 = m_level[ilevel].rfbox[momid];
#ifndef IMPORT_MESH				
					ComptRefinePtxyz(ilevel, i0, m_level[ilevel].tofine.back().xyz, sonx, sony, sonz, 
						m_level[ilevel].newfinept[fpt].boxpointx,
						m_level[ilevel].newfinept[fpt].boxpointy, 
						m_level[ilevel].newfinept[fpt].boxpointz);
#endif
#ifdef IMPORT_MESH
					ComptPtxyz_Refine(ilevel, i0, m_level[ilevel].tofine.back().xyz, sonx, sony, sonz,
						m_level[ilevel].newfinept[fpt].boxpointx,
						m_level[ilevel].newfinept[fpt].boxpointy, 
						m_level[ilevel].newfinept[fpt].boxpointz);
#endif
				}								
			}
			else
			{
				for (int fpt = 0; fpt < m_level[ilevel].newfinept.size(); ++fpt)
				{
					m_level[ilevel].tofine.push_back(mPoint());
					m_level[ilevel].tofine.back().index = m_level[ilevel].newfinept[fpt].arrindex;
					int sonid = m_level[ilevel].newfinept[fpt].vectindex%4;
					int momid = (m_level[ilevel].newfinept[fpt].vectindex-sonid)/4;
					int sony = sonid%2;
					int sonx = (sonid - sony)/2;
					int sonz = 0;
					int i0 = m_level[ilevel].rfbox[momid];
#ifndef IMPORT_MESH				
					ComptRefinePtxyz_2d(ilevel, i0, m_level[ilevel].tofine.back().xyz, sonx, sony, sonz, 
						m_level[ilevel].newfinept[fpt].boxpointx,
						m_level[ilevel].newfinept[fpt].boxpointy, 
						m_level[ilevel].newfinept[fpt].boxpointz);
#endif
#ifdef IMPORT_MESH
					ComptPtxyz_Refine_2d(ilevel, i0, m_level[ilevel].tofine.back().xyz, sonx, sony, sonz,
						m_level[ilevel].newfinept[fpt].boxpointx,
						m_level[ilevel].newfinept[fpt].boxpointy, 
						m_level[ilevel].newfinept[fpt].boxpointz);
#endif
				}	
			}
			m_mesh.m_level[ilevel+1].m_point.Addnew(m_level[ilevel].tofine);
			m_mesh.m_level[ilevel+1].m_point.DirectlyReduceNew();
		}
		PRINTFinLEVEL("New point number for the coarse level: %d",ilevel,(int)m_level[ilevel].newcoarsept.size());
		PRINTFinLEVEL("New point number for the fine level: %d",ilevel,(int)m_level[ilevel].newfinept.size());
		//ExcludeOldPoints(ilevel);
		MPI_Barrier(share_comm);
	}