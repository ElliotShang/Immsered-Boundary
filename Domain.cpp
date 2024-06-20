#include "Domain.H"
	void Domain::splitdomain_tonode(DataArray<Box> & mbox, DataArray<NodePair> & mbkp)
	{
		vector<Box> b0;
		vector<int> b0_ghosttag;
		vector<NodePair> np_ghost;
		if (srank == 0 ) printf("N%dR%d domain x [%d,%d] y [%d,%d] z [%d,%d] length [%d,%d,%d]\n", 
			node, srank, dmbound[0][0], dmbound[0][1], dmbound[1][0], dmbound[1][1],
			dmbound[2][0], dmbound[2][1], dmbound[0][2], dmbound[1][2], dmbound[2][2]);
		if (srank == 0)
		{	
			vector<vector<vector<int> > > 
				bxloc(dmbound[0][2]+2*ighost, vector<vector<int> >
							(dmbound[1][2]+2*ighost, vector<int>(dmbound[2][2]+2*ighost, -1)));
			vector<Point> box3dloc;
			for (int i = 0; i < dmbound[0][2]+2*ighost; ++i)
			{
				for (int j = 0; j < dmbound[1][2]+2*ighost; ++j)
				{
					for (int k = 0; k < dmbound[2][2]+2*ighost; ++k)
					{
						int x0 = i - ighost + dmbound[0][0]; 
						int y0 = j - ighost + dmbound[1][0];
						int z0 = k - ighost + dmbound[2][0];
						
						if (x0 >= dmbound[0][0] && x0 < dmbound[0][1] &&
								y0 >= dmbound[1][0] && y0 < dmbound[1][1] &&
								z0 >= dmbound[2][0] && z0 < dmbound[2][1])
						{
							b0_ghosttag.push_back(-1);
							b0.push_back(Box(x0,y0,z0));
							bxloc[i][j][k] = b0.size()-1;
							box3dloc.push_back(Point(i,j,k));
							b0.back().type = Normalcell;
						}
					}
				}
			}
			for (int i = 0; i < dmbound[0][2]+2*ighost; ++i)
			{
				for (int j = 0; j < dmbound[1][2]+2*ighost; ++j)
				{
					for (int k = 0; k < dmbound[2][2]+2*ighost; ++k)
					{
						int x0 = i - ighost + dmbound[0][0]; 
						int y0 = j - ighost + dmbound[1][0];
						int z0 = k - ighost + dmbound[2][0];
						
						if (x0 < dmbound[0][0] || x0 >= dmbound[0][1] ||
								y0 < dmbound[1][0] || y0 >= dmbound[1][1] ||
								z0 < dmbound[2][0] || z0 >= dmbound[2][1])
						{
							b0_ghosttag.push_back(0);
							b0.push_back(Box(x0,y0,z0));
							bxloc[i][j][k] = b0.size()-1;
							box3dloc.push_back(Point(i,j,k));
							b0.back().type = Dmghost;
							//printf("N%dR%d box %d %d %d\n",node,srank,x0, y0, z0);
						}
					}
				}
			}
			for (int i = 0; i < b0.size(); ++i)
			{
				for (Point_iterator p(0,3); p.end(); ++p)
				{
					int nx = box3dloc[i][0]+p.i-1;
					int ny = box3dloc[i][1]+p.j-1;
					int nz = box3dloc[i][2]+p.k-1;
					if (nx > -1 && nx < dmbound[0][2]+2*ighost 
					 && ny > -1 && ny < dmbound[1][2]+2*ighost 
					 && nz > -1 && nz < dmbound[2][2]+2*ighost)
					{
						if (bxloc[nx][ny][nz] < 0)
						{
							printf("N%d back box (%d,%d,%d) is %d!!!\n", node,nx,ny,nz,bxloc[nx][ny][nz]);
							MPI_Abort(MPI_COMM_WORLD,573);
						}
						b0[i].setneib(p.i,p.j,p.k,bxloc[nx][ny][nz]);
					}
					else 
					{
						if (b0[i].type != Dmghost)
						{
							printf("Error of the domain neib at box (%d,%d,%d)!!!\n",b0[i].ix(),b0[i].iy(),b0[i].iz());
							MPI_Abort(MPI_COMM_WORLD,582);
						}
						b0[i].setneib(p.i,p.j,p.k,-1);
					}
				}
			}
			// printf("Start to construct the block interfaces...\n");
			ConstructBlockInterface(np_ghost, bxloc, b0);
		}
		//PRINTFinLEVEL("B0 size %d\n", 0, (int)b0.size());
		GiveAFlag("Try to put the back mesh into the array...", 5);
		mbox.Addnew(b0, b0_ghosttag);
		GiveAFlag("Finish add new boxes for first layer!!!", 5);
		mbox.DirectlyReduceNew();
		GiveAFlag("Finish make box array of the first layer!!!", 5);
		mbkp.Addnew(np_ghost);
		mbkp.DirectlyReduceNew();
		//PRINTFinLEVELRANK0("Back ground box number %d block pair number %d\n", 0, (int)mbox.size(), (int)mbkp.size());
		if (srank ==0) printf("N%dR%d Back ground box number %d block pair number %d\n", node, srank, (int)mbox.size(), (int)mbkp.size());
		GiveAFlag("Finish make block pair of the first layer!!!", 5);
		MPI_Barrier(MPI_COMM_WORLD);
	}
	void Domain::ConstructBlockInterface(vector<NodePair> & np_ghost, vector<vector<vector<int> > > bxloc, vector<Box> & b0)
	{
		int gssx[3] = {0, ighost, dmbound[0][2]+ighost};
		int gssy[3] = {0, ighost, dmbound[1][2]+ighost};
		int gssz[3] = {0, ighost, dmbound[2][2]+ighost};
		int gsex[3] = {ighost, dmbound[0][2]+ighost, dmbound[0][2]+2*ighost};
		int gsey[3] = {ighost, dmbound[1][2]+ighost, dmbound[1][2]+2*ighost};
		int gsez[3] = {ighost, dmbound[2][2]+ighost, dmbound[2][2]+2*ighost};
		int sign_thisnode[3] = {0,0,-1};
		int sign_neibnode[3] = {1,0,0};

		vector<vector<int> > senddata_ghost(nodenum);
		for (int bx = 0; bx < 3; ++bx)
		{
			for (int by = 0; by < 3; ++by)
			{
				for (int bz = 0; bz < 3; ++bz)
				{
					if (bx != 1 || by != 1 || bz != 1)
					{
						if (neibblocks[bx][by][bz] > -1)
						{
							senddata_ghost[neibblocks[bx][by][bz]].push_back(-99);
							senddata_ghost[neibblocks[bx][by][bz]].push_back(2-bx);
							senddata_ghost[neibblocks[bx][by][bz]].push_back(2-by);
							senddata_ghost[neibblocks[bx][by][bz]].push_back(2-bz);
							for (int i = gssx[bx]; i < gsex[bx]; ++i)
							{
								for (int j = gssy[by]; j < gsey[by]; ++j)
								{
									for (int k = gssz[bz]; k < gsez[bz]; ++k)
									{
										// int i0 = i+sign0[bx]*dmbound[0][2];
										// int j0 = j+sign0[by]*dmbound[1][2];
										// int k0 = k+sign0[bz]*dmbound[2][2];
										int i0 = i+sign_thisnode[bx]*alldmbound[node][0][2]+sign_neibnode[bx]*alldmbound[neibblocks[bx][by][bz]][0][2];
										int j0 = j+sign_thisnode[by]*alldmbound[node][1][2]+sign_neibnode[by]*alldmbound[neibblocks[bx][by][bz]][1][2];
										int k0 = k+sign_thisnode[bz]*alldmbound[node][2][2]+sign_neibnode[bz]*alldmbound[neibblocks[bx][by][bz]][2][2];
										senddata_ghost[neibblocks[bx][by][bz]].push_back(i0);
										senddata_ghost[neibblocks[bx][by][bz]].push_back(j0);
										senddata_ghost[neibblocks[bx][by][bz]].push_back(k0);
										senddata_ghost[neibblocks[bx][by][bz]].push_back(bxloc[i][j][k]);
										b0[bxloc[i][j][k]].type = Blockghost;
										//printf("N%d box(%d,%d,%d) senddata_ghost (%d,%d,%d)\n",node,i,j,k,i0,j0,k0);
									}
								}
							}
						}
					}
				}
			}
		}
		//printf("Finish construct local block pair...\n");
		vector<int> recvdata;
		vector<int> ghost_node;
		BcastNewInfo_NodeRank(senddata_ghost, recvdata, ghost_node, 0);
		int recvnum0 = recvdata.size();
		int s = 0; int bx, by, bz;
		while (s < recvnum0)
		{
			//printf("N%dR%d the received data %d is %d\n", node, srank, s, recvdata[s]);
			if (recvdata[s] == -99)
			{
				++s;
				bx = recvdata[s]; ++s;
				by = recvdata[s]; ++s;
				bz = recvdata[s]; ++s;
			}
			int i0 = recvdata[s]; ++s; 
			int j0 = recvdata[s]; ++s; 
			int k0 = recvdata[s]; ++s; 
			while ((i0 < ighost || i0 >= ighost + dmbound[0][2]) && PeriodicX())
			{
				int i00 = i0 - ighost;
				PeriodicNeib(i00, dmbound[0][2]);
				i0 = i00 + ighost;
			}
			while ((j0 < ighost || j0 >= ighost + dmbound[1][2]) && PeriodicY())
			{
				int i00 = j0 -ighost;
				PeriodicNeib(i00, dmbound[1][2]);
				j0 = i00 + ighost;
			}
			while ((k0 < ighost || k0 >= ighost + dmbound[2][2]) && PeriodicZ())
			{
				int i00 = k0 -ighost;
				PeriodicNeib(i00, dmbound[2][2]);
				k0 = i00 + ighost;
			}
			Assert(i0 >= ighost && i0 < ighost + dmbound[0][2], "Domain block target x error!!!", 67);
			Assert(j0 >= ighost && j0 < ighost + dmbound[1][2], "Domain block target y error!!!", 68);
			Assert(k0 >= ighost && k0 < ighost + dmbound[2][2], "Domain block target z error!!!", 67);
			//printf("bx by bz [%d,%d,%d] i0 %d j0 %d k0 %d\n",bx,by,bz, i0, j0, k0);
			np_ghost.push_back(NodePair(bxloc[i0][j0][k0], ghost_node[s], recvdata[s]));
			++s;
		}
	}

	void Domain::Set_Z_Bound(vector<int> & nodestart, vector<int> & nodeend)
	{
		Assert(nodestart.size()==nodenum, "Node slice start is not equal to the nodenum", 287);
		Assert(nodeend.size()==nodenum, "Node slice end is not equal to the nodenum", 288);
		int s = 0;
		for (int i = 0; i < blocknum[0]; ++i)
		{
			for (int j = 0; j < blocknum[1]; ++j)
			{
				for (int k = 0; k < blocknum[2]; ++k)
				{
					alldmbound[s][2][0] = nodestart[s];
					alldmbound[s][2][1] = nodeend[s]+1;
					alldmbound[s][2][2] = alldmbound[s][2][1]-alldmbound[s][2][0];
					if (node == s)
					{
						dmbound[0] = alldmbound[s][0];
						dmbound[1] = alldmbound[s][1];
						dmbound[2] = alldmbound[s][2];
					}
					++s;
				}
			}
		}
	}
