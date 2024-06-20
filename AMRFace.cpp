#include "AMR.H"
Point Face::unitnv[3] = {Point(1,0,0), Point(0,1,0), Point(0,0,1)};
	// void AMR::FaceAdjust()
	// {
	// 	vector<vector<Face> > newfineface(m_mesh.cur_level_num);
	// 	vector<vector<Face> > newcface(m_mesh.cur_level_num);
	// 	for (int i = init_adjust_level; i < m_mesh.cur_level_num; ++i)
	// 	{
	// 		MakeNewFace(i, newfineface[i], newcface[i]);
	// 		printf("R%dL%d new fine face is %d new coarse face is %d\n", 
	// 			nrank,i, (int)newfineface[i].size(), (int)newcface[i].size());
	// 	}
	// 	MPI_Barrier(MPI_COMM_WORLD);
	// 	GiveAFlag("Finish Make new faces!!!", 5);
	// 	for (int i = init_adjust_level; i < m_mesh.cur_level_num; ++i)
	// 	{
	// 		if (NULL != m_level[i].coarselevel)
	// 		{
	// 			m_mesh.m_level[i-1].m_face.Addnew(newcface[i]);
	// 		}
	// 		if (NULL != m_level[i].finelevel)
	// 		{
	// 			m_mesh.m_level[i+1].m_face.Addnew(newfineface[i]);
	// 		}
	// 	}
	// 	MPI_Barrier(MPI_COMM_WORLD);
	// 	GiveAFlag("Finish put the new faces into the array!!!", 5);
	// 	for (int i = init_adjust_level; i < m_mesh.cur_level_num; ++i)
	// 	{
	// 		//printf("Face Level %d\n", i);
	// 		FaceArrayHole(i);
	// 		GiveAFlag("Finish finding the face array holes!!!", 29);
	// 		m_mesh.m_level[i].m_face.holeplan(true);
	// 		GiveAFlag("Finish hole plan for the face array!!!", 31);
	// 		m_mesh.m_level[i].m_face.CompressArray();
	// 	}
	// 	MPI_Barrier(MPI_COMM_WORLD);
	// 	GiveAFlag("Finish compress the face array!!!", 5);
	// }
	// void AMR::RenewFaceIndex(const int & ilevel)
	// {
	// 	/*Inverse renew is not allowed because the m_geom array has not been expanded yet*/
	// 	for (int i = m_mesh.m_level[ilevel].m_face.ps(); i < m_mesh.m_level[ilevel].m_face.pe(); ++i)
	// 	{
	// 		for (int fn0 = 0; fn0 < 2; ++fn0)
	// 		{
	// 			int bn = m_mesh.m_level[ilevel].m_face[i][fn0];
	// 			// PRINTFinLEVEL("Face %d side %d is box %d (%d,%d,%d) new index is %d",ilevel,i,fn0,bn,
	// 			// 	m_level[ilevel].m_box[bn].ix(),m_level[ilevel].m_box[bn].iy(),m_level[ilevel].m_box[bn].iz(),
	// 			// 	m_mesh.m_level[ilevel].m_box[bn].neib[1][1][1]);
	// 			//Assert(bn > -1, "Negative index in AMR::RenewFaceIndex!!!", 45);
	// 			if (bn > -1)
	// 			{
	// 				m_mesh.m_level[ilevel].m_face[i][fn0] = 
	// 					m_mesh.m_level[ilevel].m_box[bn].neib[1][1][1];
	// 			}
	// 			//printf("Level %d face %d side %d box %d\n", ilevel,i,fn0,m_mesh.m_level[ilevel].m_box[bn].neib[1][1][1].index);
	// 			// Assert(m_mesh.m_level[ilevel].m_face[i][fn0] > -1, 
	// 			// 	"The new side box must be positive", 49);				
	// 		}
	// 	}
	// 	MPI_Barrier(MPI_COMM_WORLD);
	// }
// 	void AMR::MakeNewFace(const int & i, vector<Face> & newfineface, vector<Face> & newcface)
// 	{
// 		for (int j = 0; j < m_level[i].rfbox.size(); ++j)
// 		{
// 			int j0 = m_level[i].rfbox[j];
// 			if (m_level[i].m_tag[j0].tag_noghost > -1)
// 			{
// 				for (int f1 = 0; f1 < 2; ++f1)
// 				{
// 					for (int f2 = 0; f2 < 2; ++f2)
// 					{
// 						newfineface.push_back(Face());
// 						newfineface.back()[0] = m_level[i].allson[m_level[i].m_tag[j0].tag].son[0][f1][f2];
// 						newfineface.back()[1] = m_level[i].allson[m_level[i].m_tag[j0].tag].son[1][f1][f2];
// 						newfineface.back().fnv = 0;
// 						ComptFaceArea(i+1, newfineface.back());
// 					}
// 				}
// 				for (int f1 = 0; f1 < 2; ++f1)
// 				{
// 					for (int f2 = 0; f2 < 2; ++f2)
// 					{
// 						newfineface.push_back(Face());
// 						newfineface.back()[0] = m_level[i].allson[m_level[i].m_tag[j0].tag].son[f1][0][f2];
// 						newfineface.back()[1] = m_level[i].allson[m_level[i].m_tag[j0].tag].son[f1][1][f2];
// 						newfineface.back().fnv = 1;
// 						ComptFaceArea(i+1, newfineface.back());
// 					}
// 				}
// 				for (int f1 = 0; f1 < 2; ++f1)
// 				{
// 					for (int f2 = 0; f2 < 2; ++f2)
// 					{
// 						newfineface.push_back(Face());
// 						newfineface.back()[0] = m_level[i].allson[m_level[i].m_tag[j0].tag].son[f1][f2][0];
// 						newfineface.back()[1] = m_level[i].allson[m_level[i].m_tag[j0].tag].son[f1][f2][1];
// 						newfineface.back().fnv = 2;
// 						ComptFaceArea(i+1, newfineface.back());
// 					}
// 				}
// 			}
// 		}
// 		int twoboxnum = 0;
// 		//printf("MakeNewFace for level %d start is %d end is %d\n",i,m_mesh.m_level[i].m_face.ps(),
// 		//	m_mesh.m_level[i].m_face.pe());
// 		for (int j = m_mesh.m_level[i].m_face.ps(); j < m_mesh.m_level[i].m_face.pe(); ++j)
// 		{
// 			int bn[2] = {m_mesh.m_level[i].m_face[j][0], m_mesh.m_level[i].m_face[j][1]};
// 			/*make new fine faces*/
// 			for (int fn = 0; fn < 2; ++fn)
// 			{
// 				if (m_mesh.m_level[i].m_box.isnormal(bn[fn]))
// 				{
// 					if (m_level[i].m_tag[bn[fn]].tag_noghost > -1)
// 					{
// 						for (int f1 = 0; f1 < 4; ++f1)
// 						{
// 								newfineface.push_back(Face());
// 								newfineface.back().fnv = m_mesh.m_level[i].m_face[j].fnv;
// 						}
// 						FaceSideBox(i, newfineface, fn, bn[fn]);
// 						if (m_level[i].m_box.isghost(bn[1-fn]))
// 						{
// 							m_level[i].remove_my_face.push_back(j);
// 						}
// 						else
// 						{
// 							if (m_level[i].m_tag[bn[1-fn]].tag_noghost > -1)
// 							{
// 								m_level[i].remove_my_face.push_back(j);
// 							}
// 						}
// 						goto FINISHCYCLE;
// 					}															
// 					else if (m_level[i].m_tag[bn[fn]].detag_noghost > -1)
// 					{
// 						Point va(0,0,0);
// 						int * va0 = &va[0];
// 						va0[m_mesh.m_level[i].m_face[j].fnv] = 1-fn;
// 						if (m_level[i].m_box[bn[fn]].ix()%2 == va[0] &&
// 								m_level[i].m_box[bn[fn]].iy()%2 == va[1] &&
// 								m_level[i].m_box[bn[fn]].iz()%2 == va[2])
// 						{
// 							Point apt(1,1,1);
// 							newcface.push_back(Face());
// 							newcface.back().fnv = m_mesh.m_level[i].m_face[j].fnv;
// 							newcface.back()[fn] = m_level[i].m_tag[bn[fn]].detag;
// 							apt[newcface.back().fnv] = 2-2*fn;
// 							newcface.back()[1-fn] = m_level[i-1].m_box[newcface.back()[fn]].neib[apt[0]][apt[1]][apt[2]];
// 							//printf("new coarse face %d side detag is %d %d\n", (int)newcface.size()-1,
// 							//	m_level[i].m_tag[bn[fn]].detag, 
// 							//	m_level[i].m_tag[m_level[i].m_box[bn[fn]].neib[apt[0]][apt[1]][apt[2]].index].detag);
// 							//printf("face first box is %d neib (%d,%d,%d) is %d\n", newcface.back()[0], apt[0], apt[1], apt[2], newcface.back()[1]);
// 							ComptFaceArea(i-1, newcface.back());
// 						}
// 						if (m_level[i].m_box.isghost(bn[1-fn]))
// 						{
// 							m_level[i].remove_my_face.push_back(j);
// 						}
// 						else
// 						{
// 							if (m_level[i].m_tag[bn[1-fn]].detag_noghost > -1)
// 							{
// 								m_level[i].remove_my_face.push_back(j);
// 							}
// 						}
// 						goto FINISHCYCLE;
// 					}
// 				}
// 			}
// 			/*Important!!!!!!!*/
// 			/*the faces between the normal and ghost cells in the finer level have been replace*/
// 			/*so they must be removed*/
// 			/*make new coarse faces*/			
// 			FINISHCYCLE:;
// 		}
// // #ifdef DEBUG
// // 			else if (m_mesh.m_level[i].m_face[j].haszerobox())
// // 			{
// // 				printf("###A face without box index is detected!!!\n");
// // 				MPI_Abort(MPI_COMM_WORLD, 741);
// // 			}
// // #endif
// 		//}
// 		// Point apd(1,1,1);
// 		// int * apddir = &apd[0];
// 		// Boxloc * anb0
// 		// for (int j = 0; j < m_level[i].f_pro_ghost.size(); ++j)
// 		// {
// 		// 	int ci0 = m_level[i].f_pro_ghost[j].ci;
// 		// 	if (m_level[i].m_tag[ci0].tag_noghost > -1)
// 		// 	{
// 		// 		for (Point_iterator p(0,2); p,end(); ++p)
// 		// 		{
// 		// 			int * piptr = &p.i;
// 		// 			int fi0 = m_level[i].f_pro_ghost[j].fi[p.i][p.j][p.k];
// 		// 			Assert(m_level[i+1].m_box.isghost(fi0), "Error in One side face ghost!!!", 191);
// 		// 			for (int di = 0; di < 3; ++di)
// 		// 			{
// 		// 				apddir[di] = 2*p.i;
// 		// 				anb0 = &m_level[i+1].m_box[fi0].neib[apd[0]][apd[1]][apd[2]];
// 		// 				if (m_level[i+1].m_box.isnormal(anb0->index))
// 		// 				{
// 		// 					int finefaceindex = m_level[i+1].m_geom[anb0->index].faces[di][1-piptr[dir]];
// 		// 					m_level[i+1].m_face[finefaceindex][]
// 		// 				}
// 		// 			}
// 		// 		}
// 		// 	}
// 		// }
// 		// 	/*In this option, the face will only has one box*/
// 		// 	else
// 		// 	{
// 		// 		int bn = -1;
// 		// 		int fn = -1;
// 		// 		for (int fn0 = 0; fn0 < 2; ++fn0)
// 		// 		{
// 		// 			if (m_mesh.m_level[i].m_face[j][fn0] > -1)
// 		// 			{
// 		// 				bn = m_mesh.m_level[i].m_face[j][fn0];
// 		// 				fn = fn0;
// 		// 			}
// 		// 		}
// 		// 		Assert((bn>-1&&(fn>-1&&fn<2)), "One-side face error!!!", 234);					
// 		// 		Assert(m_mesh.m_level[i].m_face[j][1-fn] == -1, "One side face error!!!", 802);
// 		// 		if (m_level[i].m_tag[bn].tag > -1 || m_level[i].m_tag[bn].detag > -1)
// 		// 		{
// 		// 			m_mesh.m_level[i].m_face[j][fn] = -1;
// 		// 			Boxloc myneib(-1,-1);
// 		// 			if (m_mesh.m_level[i].m_face[j].fnv == 0)
// 		// 			{
// 		// 				myneib = m_mesh.m_level[i].m_box[bn].neib[2-2*fn][1][1];
// 		// 			}
// 		// 			else if (m_mesh.m_level[i].m_face[j].fnv == 1)
// 		// 			{
// 		// 				myneib = m_mesh.m_level[i].m_box[bn].neib[1][2-2*fn][1];
// 		// 			}
// 		// 			else if (m_mesh.m_level[i].m_face[j].fnv == 2)
// 		// 			{
// 		// 				myneib = m_mesh.m_level[i].m_box[bn].neib[1][1][2-2*fn];
// 		// 			}
// 		// 			/*This is due to the AMR cycle order!!!*/
// 		// 			Level i is first refined and the box will connected with the fine level!
// 		// 			/*Derefined box in level i+1 will not be connected with all box of level i*/
// 		// 			if (myneib.level > -1)
// 		// 			{
// 		// 				if (myneib.level == i)
// 		// 				{
// 		// 					m_mesh.m_level[i].m_face[j][1-fn] = myneib.index;
// 		// 				}
// 		// 				else if (myneib.level == i+1)
// 		// 				{
// 		// 					if (myneib.index < m_level[i+1].m_box.size())
// 		// 					{
// 		// 						if (m_level[i+1].m_tag[myneib.index].detag > -1)
// 		// 						{
// 		// 							m_mesh.m_level[i].m_face[j][1-fn] = m_level[i+1].m_tag[myneib.index].detag;
// 		// 						}
// 		// 					}
// 		// 				}
// 		// 			}
// 		// 			else if (myneib.level < -1)
// 		// 			{
// 		// 				for (int f1 = 0; f1 < 4; ++f1)
// 		// 				{
// 		// 					newfineface.push_back(Face());
// 		// 					newfineface.back().fnv = m_mesh.m_level[i].m_face[j].fnv;
// 		// 				}
// 		// 				FaceSideBox(i, newfineface, fn, bn);
// 		// 			}
// 		// 		}
// 		// 		// printf("P%dL%d One side face in MakeNewFace box %d(%d,%d,%d) neib (L%dB%d)\n", nrank, i, bn,
// 		// 		// 	m_mesh.m_level[i].m_box[bn].ix(), m_mesh.m_level[i].m_box[bn].iy(),
// 		// 		// 	m_mesh.m_level[i].m_box[bn].iz(), myneib.level, myneib.index);
// 		// 	}
// 		// }
// 		MPI_Barrier(MPI_COMM_WORLD);
// #ifdef DEBUG
// 		for (int i0 = 0; i0 < newfineface.size(); ++i0)
// 		{
// 			for (int fn = 0; fn < 2; ++fn)
// 			{
// 				int bn = newfineface[i0][fn];
// 				if (bn < 0)
// 				{
// 					printf("[%d]New fine face %d dir is %d side %d box index is %d\n",
// 						nrank,i0,newfineface[i0].fnv,fn,bn);
// 					MPI_Abort(MPI_COMM_WORLD,279);
// 				}
// 			}
// 		}
// 		for (int i0 = 0; i0 < newcface.size(); ++i0)
// 		{
// 			for (int fn = 0; fn < 2; ++fn)
// 			{
// 				int bn = newcface[i0][fn];
// 				// printf("new coarse face %d side %d box is %d (%d,%d,%d)\n", i0,fn,bn,
// 				// 	m_level[i-1].m_box[bn].ix(),
// 				// 	m_level[i-1].m_box[bn].iy(),
// 				// 	m_level[i-1].m_box[bn].iz());
// 				if (bn < 0)
// 				{
// 					printf("[%d]New coarse face %d dir is %d side %d box index is %d\n",
// 						nrank,i0,newcface[i0].fnv,fn,bn);
// 					MPI_Abort(MPI_COMM_WORLD,292);
// 				}
// 			}
// 		}
// #endif			
// 	}

	// void AMR::FaceArrayHole(const int & ilevel)
	// {
	// 	// int holenum = 0;
	// 	//printf("[%d] Level %d face removed %d\n", nrank, ilevel, 
	// 	//	(int)m_level[ilevel].remove_my_face.size());
	// 	for (int i = 0; i < m_level[ilevel].remove_my_face.size(); ++i)
	// 	{
	// 		m_mesh.m_level[ilevel].m_face.givehole(m_level[ilevel].remove_my_face[i]);
	// 		// ++holenum;
	// 		// int i0 = m_level[ilevel].remove_my_face[i];
	// 		// int bn[2] = {m_mesh.m_level[ilevel].m_face[i0][0],m_mesh.m_level[ilevel].m_face[i0][1]};
	// 		// Pointxyz left = m_mesh.bc(ilevel,bn[0]);
	// 		// Pointxyz right = m_mesh.bc(ilevel,bn[1]);
	// 		// printf("extrude face %d left box is (%f,%f,%f) right box is (%f,%f,%f)\n", 
	// 		// 	i0,left[0],left[1],left[2],right[0],right[1],right[2]);
	// 		// printf("Rank %d Level %d hole %d is face %d dir %d\n", nrank, ilevel, holenum, i,
	// 		// 	m_mesh.m_level[ilevel].m_face[i].fnv);
	// 	}
	// 	if (NULL != m_level[ilevel].coarselevel)
	// 	{
	// 		for (int i = 0; i < m_level[ilevel].coarselevel->remove_finelevel_face.size(); ++i)
	// 		{
	// 			m_mesh.m_level[ilevel].m_face.givehole(m_level[ilevel].coarselevel->remove_finelevel_face[i]);
	// 		}
	// 	}
	// }

	// void AMR::PrintAllOnesideFace()
	// {
	// 	for (int i = 0; i < m_mesh.cur_level_num; ++i)
	// 	{
	// 		for (int j = m_mesh.m_level[i].m_face.ps(); j < m_mesh.m_level[i].m_face.pe(); ++j)
	// 		{
	// 			int b1 = m_mesh.m_level[i].m_face[j][0];
	// 			int b2 = m_mesh.m_level[i].m_face[j][1];
	// 			int bn;
	// 			if (b1 > -1) bn = b1;
	// 			else if (b2 > -1) bn = b2;			
	// 			if ((b1 == -1 && b2 > -1) || (b1 > -1 && b2 == -1))
	// 			{
	// 				printf("P%dL%d One side box %d (%d,%d,%d)\n", nrank, i, bn,
	// 					m_mesh.m_level[i].m_box[bn].ix(), m_mesh.m_level[i].m_box[bn].iy(),
	// 					m_mesh.m_level[i].m_box[bn].iz());
	// 			}
	// 		}
	// 	}
	// }

// 	void AMR::FaceSideBox(const int & ilevel, vector<Face> & newfineface, 
// 		const int & fn, const int & bn)
// 	{
// 		int s = 0;
// 		int nfsize = newfineface.size();
// 		Boxson<int> * finebox = &m_level[ilevel].allson[m_level[ilevel].m_tag[bn].tag];
// 		int fn0 = 2-fn*2;
// 		int rf = 1-fn;
// 		if (newfineface.back().fnv == 0)
// 		{
// 			for (int f1 = 0; f1 < 2; ++f1)
// 			{
// 				for (int f2 = 0; f2 < 2; ++f2)
// 				{
// 					int nfi = nfsize-4+s;
// 					newfineface[nfi][fn] = finebox->son[rf][f1][f2];
// 					newfineface[nfi][rf] = m_level[ilevel+1].m_box[newfineface[nfi][fn]].neib[fn0][1][1];
// 					ComptFaceArea(ilevel+1, newfineface[nfi]);
// 					++s;
// 				}
// 			}
// 		}
// 		else if (newfineface.back().fnv == 1)
// 		{
// 			for (int f1 = 0; f1 < 2; ++f1)
// 			{
// 				for (int f2 = 0; f2 < 2; ++f2)
// 				{
// 					int nfi = nfsize-4+s;
// 					newfineface[nfi][fn] = finebox->son[f1][rf][f2];
// 					newfineface[nfi][rf] = m_level[ilevel+1].m_box[newfineface[nfi][fn]].neib[1][fn0][1];
// 					ComptFaceArea(ilevel+1, newfineface[nfi]);
// 					++s;
// 				}
// 			}
// 		}
// 		else if (newfineface.back().fnv == 2)
// 		{
// 			for (int f1 = 0; f1 < 2; ++f1)
// 			{
// 				for (int f2 = 0; f2 < 2; ++f2)
// 				{
// 					int nfi = nfsize-4+s;
// 					newfineface[nfi][fn] = finebox->son[f1][f2][rf];
// 					newfineface[nfi][rf] = m_level[ilevel+1].m_box[newfineface[nfi][fn]].neib[1][1][fn0];
// 					ComptFaceArea(ilevel+1, newfineface[nfi]);
// 					++s;
// 				}
// 			}
// 		}
// 	}

// void AMR::GetRemoveFaceInfoRefine()
// {
// 	Point apd(1,1,1);
// 	int * apddir = &apd[0];
// 	for (int i = init_adjust_level; i < m_mesh.cur_level_num; ++i)
// 	{
// 		for (int i0 = 0; i0 < m_level[i].f_pro_ghost.size(); ++i0)
// 		{
// 			if (m_level[i].m_tag[m_level[i].f_pro_ghost[i0].ci].tag > -1)
// 			{
// 				/*remove the faces between the normal and ghost cells in the fine level*/
// 				for (Point_iterator p(0,2); p.end(); ++p)
// 				{
// 					int * piptr = &p.i;
// 		 			int fi0 = m_level[i].f_pro_ghost[i0].fi.son[p.i][p.j][p.k];
// 		 			for (int di = 0; di < 3; ++di)
// 					{
// 						apddir[di] = 2*p.i;
// 						int anb0 = m_level[i+1].m_box[fi0].neib[apd[0]][apd[1]][apd[2]];
// 						if (m_level[i+1].m_box.isnormal(anb0))
// 						{
// 							if (m_level[i+1].m_tag[anb0].detag == -1)
// 							{
// 								m_level[i].remove_finelevel_face.push_back(m_mesh.m_level[i+1].m_box[anb0].faces[di][1-piptr[di]]);
// 							}
// 						}
// 						apddir[di] = 1;
// 					}
// 				}
// 			}
// 		}
// 	}
// 	MPI_Barrier(MPI_COMM_WORLD);
// }

// void AMR::GetRemoveFaceInfoDerefine()
// {
// 	Point apd(1,1,1);
// 	int * apddir = &apd[0];
// 	for (int ilevel = init_adjust_level; ilevel < m_mesh.cur_level_num; ++ilevel)
// 	{
// 		for (int i = 0; i < m_level[ilevel].f_res_ghost[0].size(); ++i)
// 		{
// 			int ci0 = m_level[ilevel].f_res_ghost[0][i].ci;
// 			int fi0 = m_level[ilevel].f_res_ghost[0][i].fi.son[0][0][0];
// 			if (m_level[ilevel+1].m_tag[fi0].detag > -1)
// 			{
// 				for (int di = 0; di < 3; ++di)
// 				{
// 					for (int f0 = 0; f0 < 2; ++f0)
// 					{
// 						apddir[di] = 2*f0;
// 						int aneib = m_level[ilevel].m_box[ci0].neib[apd[0]][apd[1]][apd[2]];
// 						if (m_level[ilevel].m_box.isnormal(aneib))
// 						{
// 							if (m_level[ilevel].m_tag[aneib].tag == -1)
// 							{
// 								m_level[ilevel].remove_my_face.push_back(m_mesh.m_level[ilevel].m_box[aneib].faces[di][1-f0]);
// 							}
// 						}
// 					}
// 					apddir[di] = 1;
// 				}
// 			}
// 		}
// 		printf("Level %d Removed face due to fine level derefine is %d\n", ilevel, (int)m_level[ilevel].remove_my_face.size());
// 	}
// }
