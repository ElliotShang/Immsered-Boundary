#include "Mesh.H"
#include "AMRmpi.H"
#include "AMRLevel.H"

	void Mesh::ReverseIndexCheck(const int & level_n)
	{
		bool right = true;
		int boxstart, boxend;
		ArrayOrder_s(0, m_level[level_n].m_box.realsize(), boxstart, boxend, sprocs, srank);
		for (int i = boxstart; i < boxend; ++i)
		{
			for (Point_iterator p(0,3); p.end() ; ++p)
			{
				int neib0 = m_level[level_n].m_box[i].neib[p.i][p.j][p.k];
				if (neib0 > -1)
				{
					Assert(!m_level[level_n].m_box.outrange(neib0), "Box index out of range in ReverseIndexCheck!!!", 912);
					if (m_level[level_n].m_box[neib0].neib[2-p.i][2-p.j][2-p.k] != 
						m_level[level_n].m_box[i].neib[1][1][1])
					{
						int s0 = m_level[level_n].m_box[neib0].neib[2-p.i][2-p.j][2-p.k];
						printf("Box (L%d, B%d) (%d,%d,%d) neib (%d, %d, %d) is (%d, %d) (%d,%d,%d) but the neib box's neib is (L%d, B%d) (%d,%d,%d) \n",
							level_n, i, 
							m_level[level_n].m_box[i].ix(),m_level[level_n].m_box[i].iy(),m_level[level_n].m_box[i].iz(),
							p.i, p.j, p.k, 
							level_n, neib0, 
							m_level[level_n].m_box[neib0].ix(),m_level[level_n].m_box[neib0].iy(),m_level[level_n].m_box[neib0].iz(),
							level_n, s0,
							m_level[level_n].m_box[s0].ix(),m_level[level_n].m_box[s0].iy(),m_level[level_n].m_box[s0].iz());
						right = false;
					}
				}
			}
		}
		boxstart = m_level[level_n].m_box.ps();
		boxend = m_level[level_n].m_box.pe();
		for (int i = boxstart; i < boxend; ++i)
		{
			for (Point_iterator p(0,3); p.end() ; ++p)
			{
				int neib0 = m_level[level_n].m_box[i].neib[p.i][p.j][p.k];
				if (m_level[level_n].m_box.outrange(neib0))
				{
					printf("N%dL%dB%d(%d,%d,%d) neib(%d,%d,%d) is %d\n",node, level_n, i,
						m_level[level_n].m_box[i].ix(),
						m_level[level_n].m_box[i].iy(),
						m_level[level_n].m_box[i].iz(),
						p.i,p.j,p.k,neib0);
				}
				Assert(!m_level[level_n].m_box.outrange(neib0), "Box index out of range in ReverseIndexCheck!!!", 912);
				if (m_level[level_n].m_box[neib0].neib[2-p.i][2-p.j][2-p.k] != 
					m_level[level_n].m_box[i].neib[1][1][1])
				{
					int s0 = m_level[level_n].m_box[neib0].neib[2-p.i][2-p.j][2-p.k];
					printf("Box (L%d, B%d) (%d,%d,%d) neib (%d, %d, %d) is (%d, %d) (%d,%d,%d) but the neib box's neib is (L%d, B%d) (%d,%d,%d) \n",
							level_n, i, 
							m_level[level_n].m_box[i].ix(),m_level[level_n].m_box[i].iy(),m_level[level_n].m_box[i].iz(),
							p.i, p.j, p.k, 
							level_n, neib0, 
							m_level[level_n].m_box[neib0].ix(),m_level[level_n].m_box[neib0].iy(),m_level[level_n].m_box[neib0].iz(),
							level_n, s0,
							m_level[level_n].m_box[s0].ix(),m_level[level_n].m_box[s0].iy(),m_level[level_n].m_box[s0].iz());
						right = false;
				}
			}
		}
		if (!right)
		{
			printf("###Rank %d did not pass the ReverseIndexCheck!!!\n", nrank);
			MPI_Abort(MPI_COMM_WORLD, 328);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		if (nrank == 0)
		{
			printf("$$$Level %d pass the ReverseIndexCheck!!!\n", level_n);
		}		
	}

	void Mesh::CheckBox(const int & level_n)
	{
		int nc[3];
		int myc[3];
		bool right = true;
		// int bs = m_level[level_n].m_box.ps();
		// int be = m_level[level_n].m_box.pe();
		// for (int i = bs; i < be; ++i)
		// {
		// 	for (Point_iterator p(0,3); p.end(); ++p)
		// 	{
		// 		int myneib = m_level[level_n].m_box[i].neib[p.i][p.j][p.k];
		// 		if (myneib > -1)
		// 		{
		// 			nc[0] = m_level[level_n].m_box[myneib].ix();
		// 			nc[1] = m_level[level_n].m_box[myneib].iy();
		// 			nc[2] = m_level[level_n].m_box[myneib].iz();
		// 			myc[0] = m_level[level_n].m_box[i].ix();
		// 			myc[1] = m_level[level_n].m_box[i].iy();
		// 			myc[2] = m_level[level_n].m_box[i].iz();
		// 			if (myc[0]-nc[0] != 1-p.i)
		// 			{
		// 				if (!PeriodicX() || 
		// 					(PeriodicX() && 
		// 						myc[0]+nc[0] != Boxnum().ix()*pow(2, level_n)-1))
		// 				{
		// 					PRINTFinLEVEL("Normal Box %d (%d,%d,%d) but its neib (%d, %d, %d)"
		// 						" is Box [%d] and its x is %d",
		// 						level_n, i, 
		// 						m_level[level_n].m_box[i].ix(), 
		// 						m_level[level_n].m_box[i].iy(),
		// 						m_level[level_n].m_box[i].iz(),
		// 						p.i, p.j, p.k, 
		// 						myneib, nc[0]);
		// 					right = false;
		// 					printf("###Level Check Normal Box failed!!!\n");
		// 					MPI_Abort(MPI_COMM_WORLD, 78);
		// 				}
		// 			}
		// 			if (myc[1]-nc[1] != 1-p.j)
		// 			{
		// 				if (!PeriodicY() || 
		// 					(PeriodicY() && 
		// 						myc[1]+nc[1] != Boxnum().iy()*pow(2, level_n)-1))
		// 				{
		// 					PRINTFinLEVEL("Normal Box %d (%d,%d,%d) but its neib (%d, %d, %d)"
		// 						" is Box [%d] and its y is %d",
		// 						level_n, i, 
		// 						m_level[level_n].m_box[i].ix(), 
		// 						m_level[level_n].m_box[i].iy(),
		// 						m_level[level_n].m_box[i].iz(), 
		// 						p.i, p.j, p.k, 
		// 						myneib, nc[1]);
		// 					right = false;
		// 					printf("###Level Check Normal Box failed!!!\n");
		// 					MPI_Abort(MPI_COMM_WORLD, 94);
		// 				}						
		// 			}
		// 			if (myc[2]-nc[2] != 1-p.k)
		// 			{
		// 				if (!PeriodicZ() || 
		// 					(PeriodicZ() && 
		// 						myc[2]+nc[2] !=Boxnum().iz()*level_power_ratio[level_n]-1))
		// 				{
		// 					PRINTFinLEVEL("Normal Box %d (%d,%d,%d) but its neib (%d, %d, %d)"
		// 						" is Box [%d] and its z is %d",
		// 						level_n, i, 
		// 						m_level[level_n].m_box[i].ix(), 
		// 						m_level[level_n].m_box[i].iy(),
		// 						m_level[level_n].m_box[i].iz(), 
		// 						p.i, p.j, p.k, 
		// 							myneib, nc[2]);
		// 					right = false;
		// 					printf("###Level Check Normal Box failed!!!\n");
		// 					MPI_Abort(MPI_COMM_WORLD, 109);
		// 				}
		// 			}
		// 		}
		// 		else
		// 		{
		// 			printf("Level %d Box (%d,%d,%d) neib (%d,%d,%d) is %d\n", level_n, 
		// 				m_level[level_n].m_box[i].ix(), 
		// 				m_level[level_n].m_box[i].iy(),
		// 				m_level[level_n].m_box[i].iz(), 
		// 				p.i, p.j, p.k, 
		// 				myneib);
		// 			printf("###Level Check Normal Box failed!!!\n");
		// 			MPI_Abort(MPI_COMM_WORLD, 136);
		// 		}
		// 	}
		// }
		// MPI_Barrier(MPI_COMM_WORLD);
		// //printf("Finish normal cells neib check!!!\n");
		// if (nrank == 0)
		// {
		// 	printf("$$$All ranks finished level %d check!!!\n", level_n);
		// }
		/*----------------------------------------------------------*/
		int boxstart, boxend;
		m_level[level_n].m_box.GlobalOrder(boxstart, boxend);
		right = true;
		for (int i = boxstart; i < boxend; ++i)
		{
			for (Point_iterator p(0,3); p.end(); ++p)
			{
				int myneib = m_level[level_n].m_box[i].neib[p.i][p.j][p.k];
				if (myneib > -1)
				{
					nc[0] = m_level[level_n].m_box[myneib].ix();
					nc[1] = m_level[level_n].m_box[myneib].iy();
					nc[2] = m_level[level_n].m_box[myneib].iz();
					myc[0] = m_level[level_n].m_box[i].ix();
					myc[1] = m_level[level_n].m_box[i].iy();
					myc[2] = m_level[level_n].m_box[i].iz();
					// if (myc[0] == 32 && myc[1] == 25 && myc[2] == 12)
					// {
					// 	PRINTFinLEVEL("My neib (%d,%d,%d) is (%d,%d,%d)",level_n,p.i,p.j,p.k,nc[0],nc[1],nc[2]);
					// }
					if (myc[0]-nc[0] != 1-p.i)
					{
						if (!PeriodicX() || 
							(PeriodicX() && 
								myc[0]+nc[0] != Boxnum().ix()*pow(2, level_n)-1))
						{
							PRINTFinLEVEL("Box %d (%d,%d,%d) but its neib (%d, %d, %d)"
								" is Box [%d] and its x is %d",
								level_n, i, 
								m_level[level_n].m_box[i].ix(), 
								m_level[level_n].m_box[i].iy(),
								m_level[level_n].m_box[i].iz(),
								p.i, p.j, p.k, 
								myneib, nc[0]);
							right = false;
							printf("###Level CheckBox failed!!!\n");
							MPI_Abort(MPI_COMM_WORLD, 78);
						}
					}
					if (myc[1]-nc[1] != 1-p.j)
					{
						if (!PeriodicY() || 
							(PeriodicY() && 
								myc[1]+nc[1] != Boxnum().iy()*pow(2, level_n)-1))
						{
							PRINTFinLEVEL("Box %d (%d,%d,%d) but its neib (%d, %d, %d)"
								" is Box [%d] and its y is %d",
								level_n, i, 
								m_level[level_n].m_box[i].ix(), 
								m_level[level_n].m_box[i].iy(),
								m_level[level_n].m_box[i].iz(), 
								p.i, p.j, p.k, 
								myneib, nc[1]);
							right = false;
							printf("###Level CheckBox failed!!!\n");
							MPI_Abort(MPI_COMM_WORLD, 94);
						}						
					}
					if (myc[2]-nc[2] != 1-p.k)
					{
						if (!PeriodicZ() || 
							(PeriodicZ() && 
								myc[2]+nc[2] !=Boxnum().iz()*level_power_ratio[level_n]-1))
						{
							PRINTFinLEVEL("Box %d (%d,%d,%d) but its neib (%d, %d, %d)"
								" is Box [%d] and its z is %d boxnum z %d level_power_ratio %d",
								level_n, i, 
								m_level[level_n].m_box[i].ix(), 
								m_level[level_n].m_box[i].iy(),
								m_level[level_n].m_box[i].iz(), 
								p.i, p.j, p.k, 
									myneib, nc[2],
									Boxnum().iz(), level_power_ratio[level_n]);
							right = false;
							printf("###Level CheckBox failed!!!\n");
							MPI_Abort(MPI_COMM_WORLD, 109);
						}
					}
				}
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
		if (nrank == 0)
		{
			printf("$$$All ranks finished level %d check!!!\n", level_n);
		}		
	}

	void Mesh::CheckBoxOnlyOne()
	{
		for (int i = 0; i < cur_level_num; ++i)
		{
			int bs, be;
			int tbsize = m_level[i].m_box.realsize();
			m_level[i].m_box.GlobalOrder(bs, be);
			for (int b0 = bs; b0 < be; ++b0)
			{
				for (int bn = b0+1; bn < tbsize; ++bn)
				{
					if (m_level[i].m_box[b0].ix() == m_level[i].m_box[bn].ix() &&
						m_level[i].m_box[b0].iy() == m_level[i].m_box[bn].iy() &&
						m_level[i].m_box[b0].iz() == m_level[i].m_box[bn].iz())
					{
						printf("N%dL%d Box %d and %d are both (%d,%d,%d)!!!\n",
							node, i, b0, bn, m_level[i].m_box[b0].ix(),
							m_level[i].m_box[b0].iy(),
							m_level[i].m_box[b0].iz());
						MPI_Abort(MPI_COMM_WORLD, 283);
					}
				}
			}
		}
	}

	void Mesh::CheckPointxyz()
	{
		for (int i = 0; i < cur_level_num; ++i)
		{
			bool right = true;
			int boxstart, boxend;
			m_level[i].m_box.GlobalOrder(boxstart, boxend);
			for (int bn = boxstart; bn < boxend; ++bn)
			{
				int p0i = m_level[i].m_box[bn].pts[0][0][0];
				double x0 = m_level[i].m_point[p0i][0];
				double y0 = m_level[i].m_point[p0i][1];
				double z0 = m_level[i].m_point[p0i][2];
				for (Point_iterator q(0,2); q.end(); ++q)
				{
					int pij = m_level[i].m_box[bn].pts[q.i][q.j][q.k];
					double xij = m_level[i].m_point[pij][0];
					double yij = m_level[i].m_point[pij][1];
					double zij = m_level[i].m_point[pij][2];
					if (abs(xij-x0) > dh[i][0]*double(q.i)*1.001+0.0000001)
					{
						printf("@@XE N%dP%dL%dB%d point 0 is (%f,%f,%f) but point(%d,%d,%d) is (%f,%f,%f)"
							" dh[0] is %f\n",
							node,srank,i,bn,x0,y0,z0,q.i,q.j,q.k,xij,yij,zij, dh[i][0]);
						right = false;
					}
					if (abs(yij-y0) > dh[i][1]*double(q.j)*1.001+0.0000001)
					{
						printf("@@YE N%dP%dL%dB%d point 0 is (%f,%f,%f) but point(%d,%d,%d) is (%f,%f,%f)"
							" dh[1] is %f\n",
							node,srank,i,bn,x0,y0,z0,q.i,q.j,q.k,xij,yij,zij,dh[i][1]);
						right = false;
					}
					if (abs(zij-z0) > dh[i][2]*double(q.k)*1.001+0.0000001)
					{
						printf("@@ZE N%dP%dL%dB%d point 0 is (%f,%f,%f) but point(%d,%d,%d) is (%f,%f,%f)"
							" dh[2] is %f\n",
							node,srank,i,bn,x0,y0,z0,q.i,q.j,q.k,xij,yij,zij,dh[i][2]);
						right = false;
					}
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (!right)
			{
				MPI_Abort(MPI_COMM_WORLD, 935);
			}
			else
			{
				if (nrank == 0)
				{
					printf("All processors have passed Level %d the point xyz check!!!\n", i);
				}				
			}
		}
	}

	void Mesh::CheckFaceBoxIndex()
	{
		bool right = true;
		int xface = 0; int yface = 0; int zface = 0; 
		for (int ilevel = 0; ilevel < cur_level_num; ++ilevel)
		{
			int xface = 0; int yface = 0; int zface = 0;
			for (int i = m_level[ilevel].m_face.ps(); i < m_level[ilevel].m_face.pe(); ++i)
			{
				if (m_level[ilevel].m_face[i].fnv == 0) ++xface;
				else if (m_level[ilevel].m_face[i].fnv == 1) ++yface;
				else if (m_level[ilevel].m_face[i].fnv == 2) ++zface;
				for (int fn = 0; fn < 2; ++fn)
				{
					int bn = m_level[ilevel].m_face[i][fn];
					if (m_level[ilevel].m_box.outrange(bn))
					{
						printf("[%d] Level %d Face %d Side %d Box %d is out of range [0, %d]\n",
							nrank, ilevel, i, fn, bn, m_level[ilevel].m_box.realsize());
						right = false;
						printf("[%d] Level %d face box index check failed!!!\n", nrank, ilevel);
						MPI_Abort(MPI_COMM_WORLD,324);
					}
					else if (m_level[ilevel].m_box.isghost(bn))
					{
						/*This box is a ghost of the domain*/
						int bn1 = m_level[ilevel].m_face[i][1-fn];
						if (m_level[ilevel].m_box.isghost(bn1))
						{
							printf("[%d]Level %d face %d has two ghost box [%d, %d]\n", 
								nrank, ilevel, i, bn, bn1);
							right = false;
							printf("[%d] Level %d face box index check failed!!!\n", nrank, ilevel);
							MPI_Abort(MPI_COMM_WORLD,324);
						}
					}
					else if (m_level[ilevel].m_box.isnormal(bn))
					{
						int thisface = m_level[ilevel].m_box[bn].
							faces[m_level[ilevel].m_face[i].fnv][1-fn];
						if (thisface != i)
						{
							printf("Rank %d Level %d Face %d side %d box is %d but the box face index is %d\n", 
								nrank, ilevel, i, fn, bn, thisface);
							printf("face side box center is (%f,%f,%f)\n", bc(ilevel,bn)[0], bc(ilevel,bn)[1], bc(ilevel,bn)[2]);
							right = false;
							printf("[%d] Level %d face box index check failed!!!\n", nrank, ilevel);
							MPI_Abort(MPI_COMM_WORLD,324);
						}
					}					
				}
				if (!right)
				{
					printf("[%d] Level %d face box index check failed!!!\n", nrank, ilevel);
					MPI_Abort(MPI_COMM_WORLD,324);
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (right)
			{
				if (nrank == 0)
				{
					printf("WOW~~All processors have passed face to box index check!!!\n");
				}
			}
			int mbps = m_level[ilevel].m_box.ps();
			int mbpe = m_level[ilevel].m_box.pe();
			for (int i = mbps; i < mbpe; ++i)
			{
				for (int fi = 0; fi < 3; ++fi)
				{
					for (int fn0 = 0; fn0 < 2; ++fn0)
					{
						int faceindex = m_level[ilevel].m_box[i].faces[fi][fn0];
						if (m_level[ilevel].m_face.outrange(faceindex))
						{
							printf("[%d] Level %d Box %d(%d,%d,%d) Face(%d,%d) is %d not in range!!!\n",
								nrank, ilevel, i, m_level[ilevel].m_box[i].ix(),
								m_level[ilevel].m_box[i].iy(), m_level[ilevel].m_box[i].iz(), 
								fi, fn0, faceindex);
							right = false;
							printf("[%d] Level %d face box index check R2 failed!!!\n", nrank, ilevel);
							MPI_Abort(MPI_COMM_WORLD, 406);
						}
						else
						{
							int boxindex = m_level[ilevel].m_face[faceindex][1-fn0];
							if (boxindex != i)
							{
								printf("[%d] Level %d Box %d(%d,%d,%d) Face(%d,%d) is %d but the "
									"face side box is %d!!!\n",nrank,ilevel,i,m_level[ilevel].m_box[i].ix(),
									m_level[ilevel].m_box[i].iy(), m_level[ilevel].m_box[i].iz(),
									fi, fn0, faceindex, boxindex);
								right = false;
								if (boxindex > -1)
								{
									PRINTFinLEVEL("The face size box %d is (%d,%d,%d)",ilevel, boxindex,
										m_level[ilevel].m_box[boxindex].ix(),
										m_level[ilevel].m_box[boxindex].iy(),
										m_level[ilevel].m_box[boxindex].iz());
								}
								printf("[%d] Level %d face box index check R2 failed!!!\n", nrank, ilevel);
								MPI_Abort(MPI_COMM_WORLD, 406);
							}
						}
					}
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (right)
			{
				if (nrank == 0)
				{
					printf("WOW~~All processors have passed box to face index check!!!\n");
				}
			}
			else
			{
				printf("[%d] Level %d face box index check R2 failed!!!\n", nrank, ilevel);
				MPI_Abort(MPI_COMM_WORLD, 406);
			}
		}
	}

	void Mesh::CheckBoxPointIndex(const int & i)
	{
		int boxstart, boxend;
		m_level[i].m_box.GlobalOrder(boxstart, boxend);
		for (int bn0 = boxstart; bn0 < boxend; ++bn0)
		{
			for (Point_iterator p(0,2); p.end(); ++p)
			{
				int pti = m_level[i].m_box[bn0].pts[p.i][p.j][p.k];
				if (pti < 0 || pti > m_level[i].m_point.size()-1)
				{
					printf("Node %d Rank %d Level %d Box %d Point(%d,%d,%d) is %d but the up limit is %d box bc is (%f,%f,%f)\n", 
						node,srank,i,bn0,p.i,p.j,p.k,pti,
						(int)m_level[i].m_point.size(), bc(i,bn0)[0], bc(i,bn0)[1], bc(i,bn0)[2]);
					MPI_Abort(MPI_COMM_WORLD, 993);
				}
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}

	void Mesh::CheckBoxCenter()
	{
		for (int i = 0; i < cur_level_num; ++i)
		{
			for (int bn = m_level[i].m_box.ps(); bn < m_level[i].m_box.pe(); ++bn)
			{
				Pointxyz & bc = m_level[i].m_geom[bn].boxcenter;
				if (abs(bc[0] - (double(m_level[i].m_box[bn].ix())+0.5)*dh[i][0]) > 0.1*dh[i][0] ||
						abs(bc[1] - (double(m_level[i].m_box[bn].iy())+0.5)*dh[i][1]) > 0.1*dh[i][1] ||
						abs(bc[2] - (double(m_level[i].m_box[bn].iz())+0.5)*dh[i][2]) > 0.1*dh[i][2])
				{
					printf("Box center Error!!!L%dB%d is (%d,%d,%d) but its box center is (%f,%f,%f)\n",
						i,bn,m_level[i].m_box[bn].ix(),m_level[i].m_box[bn].iy(),m_level[i].m_box[bn].iz(),
						bc[0],bc[1],bc[2]);
					MPI_Abort(MPI_COMM_WORLD, 509);
				}
			}
			for (int bn = m_level[i].m_box.gps(); bn < m_level[i].m_box.gpe(); ++bn)
			{
				Pointxyz & bc = m_level[i].m_geom[bn].boxcenter;
				if (abs(bc[0] - (double(m_level[i].m_box[bn].ix())+0.5)*dh[i][0]) > 0.1*dh[i][0] ||
						abs(bc[1] - (double(m_level[i].m_box[bn].iy())+0.5)*dh[i][1]) > 0.1*dh[i][1] ||
						abs(bc[2] - (double(m_level[i].m_box[bn].iz())+0.5)*dh[i][2]) > 0.1*dh[i][2])
				{
					printf("Ghost Box center Error!!!L%dB%d is (%d,%d,%d) but its box center is (%f,%f,%f)\n",
						i,bn,m_level[i].m_box[bn].ix(),m_level[i].m_box[bn].iy(),m_level[i].m_box[bn].iz(),
						bc[0],bc[1],bc[2]);
					MPI_Abort(MPI_COMM_WORLD, 517);
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
	}

	void Meshlevel::CheckGhost(const int & ilevel)
	{
		int gn = m_box.ghost_index.size();
		MPI_Allreduce(MPI_IN_PLACE, &gn, 1, MPI_INT, MPI_SUM, share_comm);
		if (gn != m_box.arrayghostnum())
		{
			PRINTFinLEVEL("Error of array ghost number!!! Array ghost number is %d but sum of index is %d", ilevel, m_box.arrayghostnum(), gn);
			MPI_Abort(MPI_COMM_WORLD, 428);
		}
		for (int i = 0; i < m_box.ghost_index.size(); ++i)
		{
			int i0 = m_box.ghost_index[i];
			m_box[i0].neib[1][1][1] = -1;
		}
		MPI_Barrier(share_comm);
		for (int i = m_box.gps(); i < m_box.gpe(); ++i)
		{
			if (m_box[i].neib[1][1][1] != -1)
			{
				PRINTFinLEVEL("Error of array ghost index!!! Box (%d,%d,%d) is not a ghost!!!",ilevel,
					m_box[i].ix(),
					m_box[i].iy(),
					m_box[i].iz());
				MPI_Abort(MPI_COMM_WORLD, 444);
			}
		}
		MPI_Barrier(share_comm);
		int mgisize = m_box.ghost_index.size();
		for (int i = 0; i < mgisize; ++i)
		{
			int i0 = m_box.ghost_index[i];
			m_box[i0].neib[1][1][1] = i0;
		}
		MPI_Barrier(share_comm);
	}

	void AMRLevel::CheckNumerofNormalghost()
	{
		m_tag.setnum_nocopy(m_box.realsize(), 0);
		int mps = m_tag.ps();
		int mpe = m_tag.pe();
		for (int i = mps; i < mpe; ++i)
		{
			m_tag[i].tag = -1;
		}
		MPI_Barrier(MPI_COMM_WORLD);
		if (NULL != coarselevel)
		{
			int fps = coarselevel->f_pro_ghost.ps();
			int fpe = coarselevel->f_pro_ghost.pe();
			if (!twodflag)
			{
				for (int i = fps; i < fpe; ++i)
				{
					for (Point_iterator p(0,2); p.end(); ++p)
					{
						int fi0 = coarselevel->f_pro_ghost[i].fi.son[p.i][p.j][p.k];
						m_tag[fi0].tag = 0;
					}
				}
			}
			else
			{
				for (int i = fps; i < fpe; ++i)
				{
					for (Point_iterator_2d p(0,2); p.end(); ++p)
					{
						int fi0 = coarselevel->f_pro_ghost[i].fi.son[p.i][p.j][p.k];
						m_tag[fi0].tag = 0;
					}
				}
			}
		}
		for (int ig = 0; ig < ighost; ++ig)
		{
			for (int i = f_res_ghost[ig].ps(); i < f_res_ghost[ig].pe(); ++i)
			{
				int ci0 = f_res_ghost[ig][i].ci;
				m_tag[ci0].tag = 0;
			}
		}
		for (int i = 0; i < ghost_target.size(); ++i)
		{
			int i0 = ghost_target[i];
			Assert(i0 > -1 && i0 < m_box.realsize(), "The ghost target should be in range!!!", 597);
			m_tag[i0].tag = 2;
		}
		for (int i = dmghost.ps(); i < dmghost.pe(); ++i)
		{
			int i0 = dmghost[i].cell;
			m_tag[i0].tag = 3;
		}
		MPI_Barrier(MPI_COMM_WORLD);
		for (int i = mps; i < mpe; ++i)
		{
			if (m_box[i].type == Normalcell)
			{
				if (m_box.isghost(i))
				{
					if (m_tag[i].tag != 0)
					{
						printf("N %d Level %d Box (%d,%d,%d) is a normal ghost but not marked by pro or res 0 1 type is %d!!!", 
							node, cur_level, m_box[i].ix(), m_box[i].iy(), m_box[i].iz(),m_tag[i].tag);
						MPI_Abort(MPI_COMM_WORLD, 608);
					}
				}
				else
				{
					Assert(m_box.isnormal(i), "The box is not a normal or ghost!!!", 607);
					if (m_tag[i].tag != -1)
					{
						PRINTFinLEVEL("Box (%d,%d,%d) is a normal cell type is %d!!!", cur_level, m_box[i].ix(), m_box[i].iy(), m_box[i].iz(),m_tag[i].tag);
						MPI_Abort(MPI_COMM_WORLD, 617);
					}
				}
			}
			else if (m_box[i].type == Dmghost)
			{
				if (m_tag[i].tag != 3)
				{
					PRINTFinLEVEL("Box (%d,%d,%d) is a domain ghost but type is %d!!!", cur_level, m_box[i].ix(), m_box[i].iy(), m_box[i].iz(),m_tag[i].tag);
					MPI_Abort(MPI_COMM_WORLD, 626);
				}
			}
			else if (m_box[i].type == Blockghost)
			{
				if (m_tag[i].tag != 2)
				{
					PRINTFinLEVEL("Box (%d,%d,%d) is a block ghost but type is %d!!!", cur_level, m_box[i].ix(), m_box[i].iy(), m_box[i].iz(),m_tag[i].tag);
					MPI_Abort(MPI_COMM_WORLD, 634);
				}
			}
			else 
			{
				PRINTFinLEVEL("Box (%d,%d,%d) has no type!!!", cur_level, m_box[i].ix(), m_box[i].iy(), m_box[i].iz());
				MPI_Abort(MPI_COMM_WORLD, 640);
			}

		}
	}

	void AMRLevel::CheckLevelProGhostType()
	{
		int bs = level_pro_ghost.ps();
		int be = level_pro_ghost.pe();
		if (NULL != finelevel)
		{
			if (!finelevel->twodflag)
			{
				for (int i = bs; i < be; ++i)
				{
					int ci0 = level_pro_ghost[i].ci;
					for (Point_iterator p(0,2); p.end(); ++p)
					{
						int fi0 = level_pro_ghost[i].fi.son[p.i][p.j][p.k];
						if (finelevel->m_box[fi0].type != m_box[ci0].type)
						{
							PRINTFinLEVEL("Level pro ghost type error!!! main type is %d son type is %d!!!",
								cur_level,m_box[ci0].type,finelevel->m_box[fi0].type);
							MPI_Abort(MPI_COMM_WORLD,1338);
						}
					}
				}
			}
			else
			{
				for (int i = bs; i < be; ++i)
				{
					int ci0 = level_pro_ghost[i].ci;
					for (Point_iterator_2d p(0,2); p.end(); ++p)
					{
						int fi0 = level_pro_ghost[i].fi.son[p.i][p.j][p.k];
						if (finelevel->m_box[fi0].type != m_box[ci0].type)
						{
							PRINTFinLEVEL("Level pro ghost type error!!! main type is %d son type is %d!!!",
								cur_level,m_box[ci0].type,finelevel->m_box[fi0].type);
							MPI_Abort(MPI_COMM_WORLD,1338);
						}
					}
				}
			}
		}
	}

