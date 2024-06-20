#include "Body.H"
#include "NS_Solver.H"

void Body::FindWallfaces(vector<Body> & abody, Mesh & amesh)
{
#ifdef TURBULENCE	
	int lastlevel = amesh.cur_level_num;
	Point adp[3][2] = {{Point(0,1,1),Point(2,1,1)},{Point(1,0,1),Point(1,2,1)},{Point(1,1,0),Point(1,1,2)}};
	Point start_pt[3][2] = {{Point(0,0,0),Point(1,0,0)},{Point(0,0,0),Point(0,1,0)},{Point(0,0,0),Point(0,0,1)}};
	Pointxyz dxyz[3][2] = {{Pointxyz(0.0, dh[lastlevel][1], 0.0)/double(segnum), Pointxyz(0.0, 0.0, dh[lastlevel][2]/double(segnum))},
						  {Pointxyz(dh[lastlevel][0]/double(segnum), 0.0, 0.0), Pointxyz(0.0, 0.0, dh[lastlevel][2]/double(segnum))},
						  {Pointxyz(dh[lastlevel][0]/double(segnum), 0.0, 0.0), Pointxyz(0.0, dh[lastlevel][1]/double(segnum), 0.0)}};
	int mps = amesh.m_dis.ps();
	int mpe = amesh.m_dis.pe();
	vector<wallface> localwf;
	for (int i = mps; i < mpe; ++i)
	{
		int ci0 = amesh.m_dis[i].ci;
		if (!amesh.m_level[lastlevel].m_box[ci0].solid && amesh.m_level[lastlevel].m_box.isnormal(ci0))
		{	
			for (int di = 0; di < DIM; ++di)
			{
				for (int f0 = 0; f0 < 2; ++f0)
				{
					int n0 = amesh.m_level[lastlevel].m_box[ci0].neib[adp[di][f0][0]][adp[di][f0][1]][adp[di][f0][2]];
					if (n0 < -1)
					{
						printf("Error in FindWallfaces!!! The box should have a neib!!!\n");
						MPI_Abort(MPI_COMM_WORLD, 34);
					}
					if (amesh.infectbox[n0] == -1)
					{
						localwf.push_back(wallface());
						localwf.back().ci = i;
						localwf.back().dir = di;
						localwf.back().side = f0;
						int corpt = amesh.m_level[lastlevel].m_box[ci0].pts[start_pt[di][f0][0]][start_pt[di][f0][1]][start_pt[di][f0][2]];
						Pointxyz & ptxyz = amesh.m_level[lastlevel].m_point[corpt].xyz;
						Pointxyz & boxbc = amesh.m_level[lastlevel].m_geom[ci0].boxcenter;
						BoxtoWall & boxpair = amesh.m_level[lastlevel].m_box[ci0].pair;
						Pointxyz & ptnmv = abody[boxpair.body].patch[boxpair.patch].nv;
						for (int s1 = 0; s1 < segnum; ++s1)
						{
							for (int s2 = 0; s2 < segnum; ++s2)
							{
								localwf.back().segcenter[s1][s2] = ptxyz + dxyz[di][0]*(double(s1)+0.5) + dxyz[di][1]*(double(s2)+0.5);
								localwf.back().segdistowall[s1][s2] = boxpair.signdis + (localwf.back().segcenter[s1][s2]-boxbc).dot(ptnmv);
								if (localwf.back().segdistowall[s1][s2] < 0.0)
								{
									printf("The wall face segment distance to the wall should be positive!!!\n");
									MPI_Abort(MPI_COMM_WORLD, 61);
								}
							}
						}
					}
				}
			}
		}
	}
	amesh.m_wallface.Addnew(localwf);
#endif	
}

void Body::ModifyWallfaceflux(vector<Body> & abody, Mesh & amesh)
{
#ifdef TURBULENCE
	int lastlevel = amesh.cur_level_num;
	int bps = amesh.m_wallface.ps();
	int bpe = amesh.m_wallface.pe();
	FlowVec segflux(0.0);
	for (int i = bps; i < bpe; ++i)
	{
		wallface & awf = amesh.m_wallface[i];
		awf.newflux = FlowVec(0.0);
		double patut = amesh.m_dis[awf.ci].ut;
		int b1 = amesh.m_dis[awf.ci].ci;
		BoxtoWall & btw = amesh.m_level[lastlevel].m_box[b1].pair;

		int fid = amesh.m_level[lastlevel].m_box[b1].faces[awf.dir][awf.side];
		int b2 = amesh.m_level[lastlevel].m_face[fid][awf.side];
		Pointxyz vn = amesh.m_dis[awf.ci].vn;
		for (int s1 = 0; s1 < segnum; ++s1)
		{
			for (int s2 = 0; s2 < segnum; ++s2)
			{
				double seg_yplus = patut*awf.segdistowall[s1][s2]/amesh.m_dis[awf.ci].hg_viseddy*Re;
				double seg_uplus;
				uplus_as_function_of_yplus(seg_uplus, seg_yplus);
				double seg_vt = seg_uplus*patut;
				double seg_vn = vn.length()/btw.signdis*awf.segdistowall[s1][s2];
				Pointxyz seg_vel = amesh.m_dis[awf.ci].tangdir*seg_vt + abody[btw.body].patch[btw.patch].nv*seg_vn + amesh.m_dis[awf.ci].patv;
				double seg_vel_shad = seg_vel.dot(amesh.m_level[lastlevel].m_face[fid].keisa);
				segflux[0] = seg_vel_shad*amesh.m_level[lastlevel].m_data[b1].roe;
				segflux[1] = seg_vel_shad*amesh.m_level[lastlevel].m_data[b1].roe*seg_vel[0];
				segflux[2] = seg_vel_shad*amesh.m_level[lastlevel].m_data[b1].roe*seg_vel[1];
				segflux[3] = seg_vel_shad*amesh.m_level[lastlevel].m_data[b1].roe*seg_vel[2];
				segflux[4] = seg_vel_shad*(amesh.m_level[lastlevel].m_data[b1].roe*amesh.m_level[lastlevel].m_data[b1].e+
										   amesh.m_level[lastlevel].m_data[b1].p);
				awf.newflux += segflux/double(segnum*segnum);
			}
		}
		awf.newflux *= amesh.m_level[lastlevel].m_face[fid].area;
	}
#endif	
}

void NS_Solver::ConfineWallface()
{
#ifdef TURBULENCE	
	int lastlevel = a_mesh.cur_level_num;
	int bps = a_mesh.m_wallface.ps();
	int bpe = a_mesh.m_wallface.pe();
	for (int i = bps; i < bpe; ++i)
	{
		wallface & awf = a_mesh.m_wallface[i];
		int b1 = a_mesh.m_dis[awf.ci].ci;
		int fid = a_mesh.m_level[lastlevel].m_box[b1].faces[awf.dir][awf.side];
		faceflux[lastlevel][b1] = awf.newflux;
	}
#endif	
}