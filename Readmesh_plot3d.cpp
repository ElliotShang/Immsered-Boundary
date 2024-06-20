#include "Mesh.H"
#include <iostream>
#include <fstream>
#include <cstdio>
#include "AMRmpi.H"

void Mesh::ReadPlot3dFile(vector<vector<vector<Pointxyz> > > & meshpts, Point & ptnums_g, double * dh0)
{
	ifstream pmesh(meshfile.c_str());
	if (!pmesh.good())
	{
		printf("Mesh file does not exist!!!\n");
		MPI_Abort(MPI_COMM_WORLD, 9);
	}
	if (srank == 0)
	{
		int meshblock;
		pmesh >> meshblock;
		Assert(meshblock > 0, "Mesh block number error!!!", 15);
		vector<Point> block_node_num(meshblock);
		vector<Point> block_origin(meshblock);
		Point m_low_pt(0,0,0), m_high_pt(0,0,0);
		for (int i = 0; i < meshblock; ++i)
		{
			pmesh >> block_node_num[i][0];
			pmesh >> block_node_num[i][1];
			pmesh >> block_node_num[i][2];
			pmesh >> block_origin[i][0];
			pmesh >> block_origin[i][1];
			pmesh >> block_origin[i][2];
			for (int j = 0; j < 3; ++j)
			{
				m_low_pt[j] = min(m_low_pt[j], block_origin[i][j]);
				m_high_pt[j] = max(m_high_pt[j], block_origin[i][j]+block_node_num[i][j]-1);				
			}
		}
#if DIM == 2
		Assert(m_high_pt[2]-m_low_pt[2]==1, "The mesh is not 2-D!!!", 38);
#endif				
		for (int i = 0; i < 3; ++i)
		{
			Assert(m_low_pt[i]==lowpt.xy[i], "Read mesh low point error!!!", 35);
			Assert(m_high_pt[i]==highpt.xy[i], "Read mesh high point error!!!", 35);
		}
		if (nrank == 0) printf("start to read the point coordinates in the plot3d file...\n");
		//int dt0 = 0;
		for (int i = 0; i < meshblock; ++i)
		{
			for (int dir = 0; dir < 3; ++dir)
			{
				for (int nz = 0; nz < block_node_num[i][2]; ++nz)
				{
					int mz = nz + block_origin[i][2]+ighost;
					Assert(mz >= ighost && mz < m_high_pt[2]+1+ighost, "Point z in Plot3D file index must in range!!!", 54);
					for (int ny = 0; ny < block_node_num[i][1]; ++ny)
					{
						int my = ny + block_origin[i][1]+ighost;
						Assert(my >= ighost && my < m_high_pt[1]+1+ighost, "Point y in Plot3D file index must in range!!!", 54);
						for (int nx = 0; nx < block_node_num[i][0]; ++nx)
						{
							int mx = nx + block_origin[i][0]+ighost;
							//printf("start to read Block %d dir %d (%d,%d,%d)", i, dir, mx,my,mz);
							Assert(mx >= ighost && mx < m_high_pt[0]+1+ighost, "Point x in Plot3D file index must in range!!!", 54);
							pmesh >> meshpts[mx][my][mz][dir];
							//++dt0;
							//printf("Read Block %d dir %d (%d,%d,%d) data number %d data is %20.16f\n", i, dir, mx,my,mz, dt0, meshpts[mx][my][mz][dir]);
							meshpts[mx][my][mz][dir] *= meshscale[dir];
						}
					}
				}
			}
		}
		if (nrank == 0) printf("Finish read the point coordinates in the plot3d file...\n");
		for (int i = ighost; i < ptnums_g[0]-ighost-1; ++i)
		{
			for (int j = ighost; j < ptnums_g[1]-ighost-1; ++j)
			{
				for (int k = ighost; k < ptnums_g[2]-ighost-1; ++k)
				{
					dh0[0] = min(dh0[0],abs(meshpts[i+1][j][k][0]-meshpts[i][j][k][0]));
					dh0[1] = min(dh0[1],abs(meshpts[i][j+1][k][1]-meshpts[i][j][k][1]));
					dh0[2] = min(dh0[2],abs(meshpts[i][j][k+1][2]-meshpts[i][j][k][2]));
					//printf("Left is %f right is %f\n", meshpts[i][j][k][0], meshpts[i+1][j][k][0]);
				}
			}
		}
		if (nrank == 0) printf("Finish Read Point coordinates from the file dh0 is %f,%f,%f!!!\n",0,dh0[0], dh0[1], dh0[2]);
		pmesh.close();
		Point nxyz;
		for (nxyz[0] = 0; nxyz[0] < ptnums_g[0]; ++nxyz[0])
		{
			for (nxyz[1] = 0; nxyz[1] < ptnums_g[1]; ++nxyz[1])
			{
				for (nxyz[2] = 0; nxyz[2] < ptnums_g[2]; ++nxyz[2])
				{
					bool newpoint = false;
					Point ref1(nxyz[0],nxyz[1],nxyz[2]), ref2(nxyz[0],nxyz[1],nxyz[2]);
					Pointxyz dr0(0.0,0.0,0.0);
					for (int di = 0; di < 3; ++di)
					{
						if (nxyz[di] == 0) 											  
							{ref1[di] = ighost, ref2[di] = ighost+1, dr0[di] = 2.0; newpoint = true;}

						else if (nxyz[di] == 1) 								  
							{ref1[di] = ighost, ref2[di] = ighost+1, dr0[di] = 1.0; newpoint = true;}

						else if (nxyz[di] == ptnums_g[di]-ighost) 
							{ref1[di] = nxyz[di]-1, ref2[di] = nxyz[di]-2, dr0[di] = 1.0; newpoint = true;}

						else if (nxyz[di] == ptnums_g[di]-1) 		  
							{ref1[di] = nxyz[di]-2, ref2[di] = nxyz[di]-3, dr0[di] = 2.0; newpoint = true;}
					}
					if (newpoint)
					{
						for (int di = 0; di < 3; ++di)
						{
							meshpts[nxyz[0]][nxyz[1]][nxyz[2]][di] = meshpts[ref1[0]][ref1[1]][ref1[2]][di] + 
								(meshpts[ref1[0]][ref1[1]][ref1[2]][di] - meshpts[ref2[0]][ref2[1]][ref2[2]][di])*dr0[di];
						}
					}
				}
			}
		}
		if (nrank == 0) printf("Finish construct the ghost points!!!%d\n",0,-1);
	}
}

void Mesh::ReadAlphaFlowGrid(vector<vector<vector<Pointxyz> > > & meshpts, Point & ptnums_g, double * dh0)
{
	ifstream pmesh(meshfile.c_str());
	if (!pmesh.good())
	{
		printf("Mesh file error!!!\n");
		MPI_Abort(MPI_COMM_WORLD, 9);
	}
	// if (srank == 0)
	// {
		Point block_node_num;	
		pmesh >> block_node_num[0];
		pmesh >> block_node_num[1];
		pmesh >> block_node_num[2];
#if DIM == 2
		if (block_node_num[2] != 2)
		{
			printf("The mesh has %d nodes in the Z direction!!!\n", block_node_num[2]);
			MPI_Abort(MPI_COMM_WORLD, 146);
		}
#endif				
		for (int i = 0; i < 3; ++i)
		{
			if (block_node_num[i] != highpt.xy[i]+1)
			{
				printf("Mesh direction %d has %d nodes but the imported mesh has %d nodes!!!\n", i, highpt.xy[i]+1, block_node_num[i]);
				MPI_Abort(MPI_COMM_WORLD, 154);
			}
		}
		printf("start to read the point coordinates in the plot3d file...\n");
		for (int dir = 0; dir < 3; ++dir)
		{
			for (int nz = 0; nz < block_node_num[2]; ++nz)
			{
				int mz = nz + ighost;
				for (int ny = 0; ny < block_node_num[1]; ++ny)
				{
					int my = ny + ighost;
					for (int nx = 0; nx < block_node_num[0]; ++nx)
					{
						int mx = nx + ighost;
						pmesh >> meshpts[mx][my][mz][dir];
						meshpts[mx][my][mz][dir] *= meshscale[dir];
					}
				}
			}
		}				
		printf("Finish read the point coordinates in the plot3d file...\n");
		for (int i = ighost; i < ptnums_g[0]-ighost-1; ++i)
		{
			for (int j = ighost; j < ptnums_g[1]-ighost-1; ++j)
			{
				for (int k = ighost; k < ptnums_g[2]-ighost-1; ++k)
				{
					dh0[0] = min(dh0[0],abs(meshpts[i+1][j][k][0]-meshpts[i][j][k][0]));
					dh0[1] = min(dh0[1],abs(meshpts[i][j+1][k][1]-meshpts[i][j][k][1]));
					dh0[2] = min(dh0[2],abs(meshpts[i][j][k+1][2]-meshpts[i][j][k][2]));
					//printf("Left is %f right is %f\n", meshpts[i][j][k][0], meshpts[i+1][j][k][0]);
				}
			}
		}
		PRINTFinLEVEL("Finish Read Point coordinates from the file dh0 is %20.10f,%20.10f,%20.10f!!!\n",0,dh0[0], dh0[1], dh0[2]);
		pmesh.close();
		ConstructGhostPoint(meshpts, ptnums_g);
	//}
}

void Mesh::ConstructBoxPoint_FromImportMesh()
{
	Point ptnums_g(highpt.xy[0]-lowpt.xy[0]+1+2*ighost,
								 highpt.xy[1]-lowpt.xy[1]+1+2*ighost,
								 highpt.xy[2]-lowpt.xy[2]+1+2*ighost);
	PRINTFinLEVELRANK0("Point number with ghost is %d %d %d\n", 0, ptnums_g[0], ptnums_g[1], ptnums_g[2]);
	m_level[0].m_point.setnum_nocopy(ptnums_g[0]*ptnums_g[1]*ptnums_g[2],0);
	double dh0[3] = {999.9, 999.9, 999.9};
	// if (srank == 0)
	// {
		vector<vector<vector<Pointxyz> > > meshpts(ptnums_g[0],
			     vector<vector<Pointxyz> > (ptnums_g[1],
			     	      vector<Pointxyz> (ptnums_g[2])));
		GiveAFlag("start to read the imported mesh!!!", 5);
#ifdef PASSAGE_ANGLE		
		ReadUnstructGrid(meshpts, ptnums_g, dh0);
#else	
		ReadAlphaFlowGrid(meshpts, ptnums_g, dh0);
#endif				
		GiveAFlag("Finish reading the imported mesh!!!", 5);
		for (int i = 0; i < m_level[0].m_box.realsize(); ++i)
		{
			for (Point_iterator p(0,2); p.end(); ++p)
			{
				int px = m_level[0].m_box[i].ix()+p.i+ighost;
				int py = m_level[0].m_box[i].iy()+p.j+ighost;
				int pz = m_level[0].m_box[i].iz()+p.k+ighost;
				Assert(px > -1 && px < ptnums_g[0], "Import mesh box point x out of range!!!",118);
				Assert(py > -1 && py < ptnums_g[1], "Import mesh box point y out of range!!!",119);
				Assert(pz > -1 && pz < ptnums_g[2], "Import mesh box point z out of range!!!",120);
				m_level[0].m_box[i].pts[p.i][p.j][p.k] = px*ptnums_g[1]*ptnums_g[2] + py*ptnums_g[2] + pz;
				Assert(m_level[0].m_box[i].pts[p.i][p.j][p.k] < m_level[0].m_point.size(), "Point index out of range!!!",124);
			}
		}
		for (int i = 0; i < ptnums_g[0]; ++i)
		{
			for (int j = 0; j < ptnums_g[1]; ++j)
			{
				for (int k = 0; k < ptnums_g[2]; ++k)
				{
					int i0 = i*ptnums_g[1]*ptnums_g[2]+j*ptnums_g[2]+k;
					m_level[0].m_point[i0].index = i0;
					m_level[0].m_point[i0].xyz = meshpts[i][j][k];
				}
			}
		}
	// }
	MPI_Barrier(share_comm);
	//MPI_Bcast(&dh0[0], 3, MPI_DOUBLE, 0, share_comm);
	ComptLevelgridsize(dh0);
	if (nrank == 0) printf("dh0 is (%20.10f,%20.10f,%20.10f)\n", dh0[0], dh0[1], dh0[2]);
	ConstructDomainCellFace(meshpts);
}

void Mesh::ConstructDomainCellFace(vector<vector<vector<Pointxyz> > > & meshpts)
{
	double nmvdir[2] = {1.0, -1.0};
	int section[3][2] = {{ighost, ighost+highpt.xy[0]},{ighost, ighost+highpt.xy[1]},{ighost, ighost+highpt.xy[2]}};
	for (int j = 0; j < 2; ++j)
	{
		facecenter[0][j] = vector<vector<Dmcellface> > (highpt.xy[1], vector<Dmcellface>(highpt.xy[2]));
		facecenter[1][j] = vector<vector<Dmcellface> > (highpt.xy[0], vector<Dmcellface>(highpt.xy[2]));
		facecenter[2][j] = vector<vector<Dmcellface> > (highpt.xy[0], vector<Dmcellface>(highpt.xy[1]));
		for (int d1 = 0; d1 < highpt.xy[1]; ++d1)
		{
			for (int d2 = 0; d2 < highpt.xy[2]; ++d2)
			{
				facecenter[0][j][d1][d2].fc = (meshpts[section[0][j]][d1+ighost][d2+ighost]+
									meshpts[section[0][j]][d1+ighost][d2+ighost+1]+
									meshpts[section[0][j]][d1+ighost+1][d2+ighost]+
									meshpts[section[0][j]][d1+ighost+1][d2+ighost+1])/4.0;
				Pointxyz dxyz1 = meshpts[section[0][j]][d1+ighost+1][d2+ighost] - meshpts[section[0][j]][d1+ighost][d2+ighost];
				Pointxyz dxyz2 = meshpts[section[0][j]][d1+ighost][d2+ighost+1] - meshpts[section[0][j]][d1+ighost][d2+ighost];
				facecenter[0][j][d1][d2].fcnmv = dxyz1.cross(dxyz2)*nmvdir[j];
				facecenter[0][j][d1][d2].fcnmv.normalize();
			}
		}
		for (int d1 = 0; d1 < highpt.xy[0]; ++d1)
		{
			for (int d2 = 0; d2 < highpt.xy[2]; ++d2)
			{
				facecenter[1][j][d1][d2].fc = (meshpts[d1+ighost][section[1][j]][d2+ighost]+
									meshpts[d1+ighost][section[1][j]][d2+ighost+1]+
									meshpts[d1+ighost+1][section[1][j]][d2+ighost]+
									meshpts[d1+ighost+1][section[1][j]][d2+ighost+1])/4.0;
				Pointxyz dxyz1 = meshpts[d1+ighost][section[1][j]][d2+ighost+1] - meshpts[d1+ighost][section[1][j]][d2+ighost];
				Pointxyz dxyz2 = meshpts[d1+ighost+1][section[1][j]][d2+ighost] - meshpts[d1+ighost][section[1][j]][d2+ighost];
				facecenter[1][j][d1][d2].fcnmv = dxyz1.cross(dxyz2)*nmvdir[j];
				facecenter[1][j][d1][d2].fcnmv.normalize();
			}
		}
		for (int d1 = 0; d1 < highpt.xy[0]; ++d1)
		{
			for (int d2 = 0; d2 < highpt.xy[1]; ++d2)
			{
				facecenter[2][j][d1][d2].fc = (meshpts[d1+ighost][d2+ighost][section[2][j]]+
									meshpts[d1+ighost][d2+ighost+1][section[2][j]]+
									meshpts[d1+ighost+1][d2+ighost][section[2][j]]+
									meshpts[d1+ighost+1][d2+ighost+1][section[2][j]])/4.0;
				Pointxyz dxyz1 = meshpts[d1+ighost+1][d2+ighost][section[2][j]] - meshpts[d1+ighost][d2+ighost][section[2][j]];
				Pointxyz dxyz2 = meshpts[d1+ighost][d2+ighost+1][section[2][j]] - meshpts[d1+ighost][d2+ighost][section[2][j]];
				facecenter[2][j][d1][d2].fcnmv = dxyz1.cross(dxyz2)*nmvdir[j];
				facecenter[2][j][d1][d2].fcnmv.normalize();
			}
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
}

struct Boxpoint
{
	int pts[2][2][2];

	Boxpoint()
	{
		for (int i = 0; i < 2; ++i)
		{
			for (int j = 0; j < 2; ++j)
			{
				for (int k = 0; k < 2; ++k)
				{
					pts[i][j][k] = -1;
				}
			}
		}
	}
};

struct Pointneib
{
	int neib[3][3][3];

	Pointneib()
	{
		for (int i = 0; i < 3; ++i)
		{
			for (int j = 0; j < 3; ++j)
			{
				for (int k = 0; k < 3; ++k)
				{
					neib[i][j][k] = -1;
				}
			}
		}
	}
};
void Mesh::ReadUnstructGrid(vector<vector<vector<Pointxyz> > > & meshpts, Point & ptnums_g, double * dh0)
{
	ifstream pmesh(meshfile.c_str());
	if (!pmesh.good())
	{
		printf("Mesh file error!!!\n");
		MPI_Abort(MPI_COMM_WORLD, 9);
	}
	int mptnum, melenum;
	pmesh >> mptnum;
	pmesh >> melenum;
	vector<Pointxyz> uns_meshpts(mptnum);
	vector<Boxpoint> meshcell(melenum);
	vector<Pointneib> meshpoints(mptnum);
	for (int i = 0; i < mptnum; ++i)
	{
		pmesh >> uns_meshpts[i][0];
		pmesh >> uns_meshpts[i][1];
		pmesh >> uns_meshpts[i][2];
		for (int di = 0; di < 3; ++di)
		{
			uns_meshpts[i][di] *= meshscale[di];
		}
	}
	for (int i = 0; i < melenum; ++i)
	{
#if DIM == 2
		pmesh >> meshcell[i].pts[0][0][0];
		pmesh >> meshcell[i].pts[1][0][0];
		pmesh >> meshcell[i].pts[1][1][0];
		pmesh >> meshcell[i].pts[0][1][0];
		pmesh >> meshcell[i].pts[0][0][1];
		pmesh >> meshcell[i].pts[1][0][1];
		pmesh >> meshcell[i].pts[1][1][1];
		pmesh >> meshcell[i].pts[0][1][1];
#else					
		pmesh >> meshcell[i].pts[0][1][0];
		pmesh >> meshcell[i].pts[1][1][0];
		pmesh >> meshcell[i].pts[1][1][1];
		pmesh >> meshcell[i].pts[0][1][1];
		pmesh >> meshcell[i].pts[0][0][0];
		pmesh >> meshcell[i].pts[1][0][0];
		pmesh >> meshcell[i].pts[1][0][1];
		pmesh >> meshcell[i].pts[0][0][1];
#endif		
	}
	pmesh.close();
	GiveAFlag("The mesh file has been closed!!!", 5);
	for (int i = 0; i < melenum; ++i)
	{
		for (Point_iterator p(0,2); p.end(); ++p)
		{
			meshcell[i].pts[p.i][p.j][p.k] -= 1;
		}
		for (Point_iterator p(0,2); p.end(); ++p)
		{
			int an0 = meshcell[i].pts[p.i][p.j][p.k];
			for (Point_iterator q(0,2); q.end(); ++q)
			{
				int nx = 1 - p.i+q.i;
				int ny = 1 - p.j+q.j;
				int nz = 1 - p.k+q.k;
				if (meshpoints[an0].neib[nx][ny][nz] == -1)
				{
					meshpoints[an0].neib[nx][ny][nz] = meshcell[i].pts[q.i][q.j][q.k];
					// if (nrank == 0 && an0 == 5095)
					// {
					// 	printf("cell 5095 neib %d %d %d is %d (%f,%f,%f)!!!\n", nx,ny,nz,
					// 		meshcell[i].pts[q.i][q.j][q.k],
					// 		uns_meshpts[meshcell[i].pts[q.i][q.j][q.k]][0],
					// 		uns_meshpts[meshcell[i].pts[q.i][q.j][q.k]][1],
					// 		uns_meshpts[meshcell[i].pts[q.i][q.j][q.k]][2]);
					// }
				}
				else
				{
					if (meshpoints[an0].neib[nx][ny][nz] != meshcell[i].pts[q.i][q.j][q.k])
					{
						printf("Read unstruct mesh error 422!!! Old pt is %d but the new is %d!!!\n", meshpoints[an0].neib[nx][ny][nz], meshcell[i].pts[q.i][q.j][q.k]);
						MPI_Abort(MPI_COMM_WORLD, 423);
					}
				}
			}
		}
	}
	GiveAFlag("Find the mesh points neib!!!", 5);
	vector<Point> ptsxyz(mptnum);
	Point adp[3] = {Point(0,1,1), Point(1,0,1), Point(1,1,0)};
	Point re_adp[3] = {Point(2,1,1), Point(1,2,1), Point(1,1,2)};
	for (int i = 0; i < mptnum; ++i)
	{
		for (int dir = 0; dir < 3; ++dir)
		{
			if (meshpoints[i].neib[adp[dir][0]][adp[dir][1]][adp[dir][2]] == -1)
			{
				ptsxyz[i][dir] = 0;
				int start_index = 0;
				int nextcell = i;
				while (meshpoints[nextcell].neib[re_adp[dir][0]][re_adp[dir][1]][re_adp[dir][2]] > -1)
				{
					++start_index;
					nextcell = meshpoints[nextcell].neib[re_adp[dir][0]][re_adp[dir][1]][re_adp[dir][2]];
					ptsxyz[nextcell][dir] = start_index;
				}
				if (start_index != highpt.xy[dir])
				{
					printf("Read unstruct mesh error 448!!! Direction %d start point is %d(%f,%f,%f) ends at index %d cell %d(%f,%f,%f) This direction should have %d points!!!\n",
						dir, i, uns_meshpts[i][0], uns_meshpts[i][1], uns_meshpts[i][2], start_index, nextcell,
						uns_meshpts[nextcell][0], uns_meshpts[nextcell][1], uns_meshpts[nextcell][2],
						highpt.xy[dir]);
					MPI_Abort(MPI_COMM_WORLD, 450);
				}
			}
		}
	}
	GiveAFlag("Finish constructing the mesh structure!!!", 5);
	for (int i = 0; i < mptnum; ++i)
	{
		for (int dir = 0; dir < 3; ++dir)
		{
			if (ptsxyz[i][dir] < 0 || ptsxyz[i][dir] > highpt.xy[dir])
			{
				printf("Mesh point %d(%f,%f,%f) dir %d index is %d out of range!!!\n", i,
					uns_meshpts[i][0],
					uns_meshpts[i][1],
					uns_meshpts[i][2],
					dir, ptsxyz[i][dir]);
				MPI_Abort(MPI_COMM_WORLD, 466);
			}
		}
		meshpts[ighost+ptsxyz[i][0]][ighost+ptsxyz[i][1]][ighost+ptsxyz[i][2]] = uns_meshpts[i];
	}
	GiveAFlag("Finish constructing the point coordinates in the plot3d file...", 5);
	for (int i = ighost; i < ptnums_g[0]-ighost-1; ++i)
	{
		for (int j = ighost; j < ptnums_g[1]-ighost-1; ++j)
		{
			for (int k = ighost; k < ptnums_g[2]-ighost-1; ++k)
			{
				dh0[0] = min(dh0[0],abs(meshpts[i+1][j][k][0]-meshpts[i][j][k][0]));
				dh0[1] = min(dh0[1],abs(meshpts[i][j+1][k][1]-meshpts[i][j][k][1]));
				dh0[2] = min(dh0[2],abs(meshpts[i][j][k+1][2]-meshpts[i][j][k][2]));
				// if (nrank == 0)
				// {
				// 	printf("Point (%d,%d,%d) is (%f,%f,%f)\n", i,j,k, meshpts[i][j][k][0],meshpts[i][j][k][1],meshpts[i][j][k][2]);
				// }
					//printf("Left is %f right is %f\n", meshpts[i][j][k][0], meshpts[i+1][j][k][0]);
			}
		}
	}
	PRINTFinLEVEL("Finish Read Point coordinates from the file dh0 is %20.10f,%20.10f,%20.10f!!!\n",0,dh0[0], dh0[1], dh0[2]);
	ConstructGhostPoint(meshpts, ptnums_g);
	
}

void Mesh::ConstructGhostPoint(vector<vector<vector<Pointxyz> > > & meshpts, Point & ptnums_g)
{
	Point nxyz;
	for (nxyz[0] = 0; nxyz[0] < ptnums_g[0]; ++nxyz[0])
	{
		for (nxyz[1] = 0; nxyz[1] < ptnums_g[1]; ++nxyz[1])
		{
			for (nxyz[2] = 0; nxyz[2] < ptnums_g[2]; ++nxyz[2])
			{
				bool newpoint = false;
				Point ref1(nxyz[0],nxyz[1],nxyz[2]), ref2(nxyz[0],nxyz[1],nxyz[2]);
				Pointxyz dr0(0.0,0.0,0.0);
				for (int di = 0; di < 3; ++di)
				{
					if (nxyz[di] == 0) 											  
						{ref1[di] = ighost, ref2[di] = ighost+1, dr0[di] = 2.0; newpoint = true;}

					else if (nxyz[di] == 1) 								  
						{ref1[di] = ighost, ref2[di] = ighost+1, dr0[di] = 1.0; newpoint = true;}

					else if (nxyz[di] == ptnums_g[di]-ighost) 
						{ref1[di] = nxyz[di]-1, ref2[di] = nxyz[di]-2, dr0[di] = 1.0; newpoint = true;}

					else if (nxyz[di] == ptnums_g[di]-1) 		  
						{ref1[di] = nxyz[di]-2, ref2[di] = nxyz[di]-3, dr0[di] = 2.0; newpoint = true;}
				}
				if (newpoint)
				{
					for (int di = 0; di < 3; ++di)
					{
						meshpts[nxyz[0]][nxyz[1]][nxyz[2]][di] = meshpts[ref1[0]][ref1[1]][ref1[2]][di] + 
							(meshpts[ref1[0]][ref1[1]][ref1[2]][di] - meshpts[ref2[0]][ref2[1]][ref2[2]][di])*dr0[di];
					}
				}
			}
		}
	}
	if (periodic[1])
	{
		int ny0;
		for (nxyz[1] = 0; nxyz[1] < ptnums_g[1]; ++nxyz[1])
		{
			int changeflag = -1;
			if (nxyz[1] < 2)
			{
				ny0 = nxyz[1] + highpt.xy[1];
				changeflag = 0;
			}
			else if (nxyz[1] > ptnums_g[1]-ighost-1)
			{
				ny0 = nxyz[1] - highpt.xy[1];
				changeflag = 1;
			}
			if (changeflag > -1)
			{
				for (nxyz[0] = 0; nxyz[0] < ptnums_g[0]; ++nxyz[0])
				{
					for (nxyz[2] = 0; nxyz[2] < ptnums_g[2]; ++nxyz[2])
					{
#ifndef PASSAGE_ANGLE
						meshpts[nxyz[0]][nxyz[1]][nxyz[2]][0] = meshpts[nxyz[0]][ny0][nxyz[2]][0];
						meshpts[nxyz[0]][nxyz[1]][nxyz[2]][2] = meshpts[nxyz[0]][ny0][nxyz[2]][2];
						if (changeflag == 0)
						{
							meshpts[nxyz[0]][nxyz[1]][nxyz[2]][1] = meshpts[nxyz[0]][ny0][nxyz[2]][1] - dmlength[1];
						}
						else
						{
							meshpts[nxyz[0]][nxyz[1]][nxyz[2]][1] = meshpts[nxyz[0]][ny0][nxyz[2]][1] + dmlength[1];
						}
#endif
#ifdef PASSAGE_ANGLE
						meshpts[nxyz[0]][nxyz[1]][nxyz[2]][0] = meshpts[nxyz[0]][ny0][nxyz[2]][0];
						double newy0, newz0;
						if (changeflag == 0)
						{
							newy0 = meshpts[nxyz[0]][ny0][nxyz[2]][1]*cos(PASSAGE_ANGLE) - meshpts[nxyz[0]][ny0][nxyz[2]][2]*sin(PASSAGE_ANGLE);
							newz0 = meshpts[nxyz[0]][ny0][nxyz[2]][2]*cos(PASSAGE_ANGLE) + meshpts[nxyz[0]][ny0][nxyz[2]][1]*sin(PASSAGE_ANGLE);
						}
						else
						{
							newy0 = meshpts[nxyz[0]][ny0][nxyz[2]][1]*cos(-PASSAGE_ANGLE) - meshpts[nxyz[0]][ny0][nxyz[2]][2]*sin(-PASSAGE_ANGLE);
							newz0 = meshpts[nxyz[0]][ny0][nxyz[2]][2]*cos(-PASSAGE_ANGLE) + meshpts[nxyz[0]][ny0][nxyz[2]][1]*sin(-PASSAGE_ANGLE);
						}
						meshpts[nxyz[0]][nxyz[1]][nxyz[2]][1] = newy0;
						meshpts[nxyz[0]][nxyz[1]][nxyz[2]][2] = newz0;
#endif								
					}
				}
			}
		}
	}
	PRINTFinLEVEL("Finish construct the ghost points!!!%d\n",0,-1);	
}
