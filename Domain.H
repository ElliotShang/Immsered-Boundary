#ifndef _Domain_H_
#define _Domain_H_

#include <iostream>
#include <vector>
#include <math.h>

#include "Point.H"
#include "Box.H"
#include "Boxloc.H"
#include "Check.H"
#include "Constants.H"
#include "Point_iterator.H"
#include "DataArray.H"
#include "Face.H"
#include "GhostPair.H"

using namespace std;

struct Dmcellface
{
	Pointxyz fc;
	Pointxyz fcnmv;
};

class Domain
{
public:
	friend class BCValues;
	friend class Mesh;

	explicit Domain(const Point & lp, const Point & hp, bool periodic0[]):
	lowpt(lp),highpt(hp)
	{
		periodicx = periodic0[0];
		periodicy = periodic0[1];
		periodicz = periodic0[2];
		int bnx = highpt.ix() - lowpt.ix();
		int bny = highpt.iy() - lowpt.iy();
		int bnz = highpt.iz() - lowpt.iz();
		boxnum = Point(bnx, bny, bnz);
		for (int i = 0; i < 3; ++i)
		{
			for (int j = 0; j < 2; ++j)
			{
				boundaryiswall[i][j] = false;
			}
		}
		DomainBound();
	}

	Point & Domainboxnum()
	{
		return boxnum;
	}

	inline void PeriodicNeib(int & i0, const int & a_boxnum)
	{
		if (i0 < 0) i0 += a_boxnum;
		else if (i0 >= a_boxnum) i0 -= a_boxnum;
	}
	/*dv and uv are the minimum and maximum values that i0 can reach, respectively*/
	inline void PeriodicNeib(int & i0, const int & dv, const int & uv)
	{
		if (i0 < dv) i0 = uv - (dv-i0) + 1 ;
		else if (i0 > uv) i0 = i0 - uv + dv - 1;
	}

	// void CreateBox(vector<Box> & b0, vector<vector<vector<int> > > & bxloc)
	// {
	// 	int bnx = boxnum.ix();
	// 	int bny = boxnum.iy();
	// 	int bnz = boxnum.iz();
	// 	for (int i = 0; i < bnx+2*ighost; ++i)
	// 	{
	// 		for (int j = 0; j < bny+2*ighost; ++j)
	// 		{
	// 			for (int k = 0; k < bnz+2*ighost; ++k)
	// 			{
	// 				bxloc[i][j][k] = -2;
	// 			}
	// 		}
	// 	}
	// 	if (periodicx)
	// 	{
	// 		int x0[2*ighost];
	// 		for (int i = 0; i < ighost; ++i)
	// 		{
	// 			x0[i] = i;
	// 			x0[2*ighost-1-i] = bnx+2*ighost-1-i;
	// 		}
	// 		for (int i = 0; i < 2*ighost; ++i)
	// 		{
	// 			for (int j = 0; j < bny+2*ighost; ++j)
	// 			{
	// 				for (int k = 0; k < bnz+2*ighost; ++k)
	// 				{
	// 					bxloc[x0[i]][j][k] = -1;
	// 				}
	// 			}			
	// 		}
	// 	}

	// 	if (periodicy)
	// 	{
	// 		int y0[2*ighost];
	// 		for (int i = 0; i < ighost; ++i)
	// 		{
	// 			y0[i] = i;
	// 			y0[2*ighost-1-i] = bny+2*ighost-1-i;
	// 		}
	// 		for (int i = 0; i < bnx+2*ighost; ++i)
	// 		{
	// 			for (int j = 0; j < 2*ighost; ++j)
	// 			{
	// 				for (int k = 0; k < bnz+2*ighost; ++k)
	// 				{
	// 					bxloc[i][y0[j]][k] = -1;
	// 				}
	// 			}			
	// 		}
	// 	}
	// 	if (periodicz)
	// 	{
	// 		int z0[2*ighost];
	// 		for (int i = 0; i < ighost; ++i)
	// 		{
	// 			z0[i] = i;
	// 			z0[2*ighost-1-i] = bnz+2*ighost-1-i;
	// 		}
	// 		for (int i = 0; i < bnx+2*ighost; ++i)
	// 		{
	// 			for (int j = 0; j < bny+2*ighost; ++j)
	// 			{
	// 				for (int k = 0; k < 2*ighost; ++k)
	// 				{
	// 					bxloc[i][j][z0[k]] = -1;
	// 				}
	// 			}			
	// 		}
	// 	}
	// 	int s = 0;
	// 	for (int i = ighost; i < bnx+ighost; ++i)
	// 	{
	// 		for (int j = ighost; j < bny+ighost; ++j)
	// 		{
	// 			for (int k = ighost; k < bnz+ighost; ++k)
	// 			{
	// 				bxloc[i][j][k] = s;
	// 				++s;
	// 			}
	// 		}
	// 	}
	// 	domain_ghost_num = 0;
	// 	for (int i = 0; i < bnx+2*ighost; ++i)
	// 	{
	// 		for (int j = 0; j < bny+2*ighost; ++j)
	// 		{
	// 			for (int k = 0; k < bnz+2*ighost; ++k)
	// 			{
	// 				if (bxloc[i][j][k] == -2)
	// 				{
	// 					bxloc[i][j][k] = s;
	// 					++s;
	// 					++domain_ghost_num;
	// 				}
	// 			}
	// 		}
	// 	}
	// 	ShowOneRankData(0, "Total box number", s, 5);
	// 	ShowOneRankData(0, "Ghost number", domain_ghost_num, 5);
	// 	Point dp(-1,-1,-1);
	// 	Point up(bnx+2*ighost, bny+2*ighost, bnz+2*ighost);
	// 	for (int i = ighost; i < bnx+ighost; ++i)
	// 	{
	// 		for (int j = ighost; j < bny+ighost; ++j)
	// 		{
	// 			for (int k = ighost; k < bnz+ighost; ++k)
	// 			{
	// 				Assert(bxloc[i][j][k] > -1, "no box in the point!!!", 5);
	// 				AddaBox(i,j,k,b0,bxloc,dp,up);
	// 			}
	// 		}
	// 	}
	// 	for (int i = 0; i < bnx+2*ighost; ++i)
	// 	{
	// 		for (int j = 0; j < bny+2*ighost; ++j)
	// 		{
	// 			for (int k = 0; k < bnz+2*ighost; ++k)
	// 			{
	// 				if (bxloc[i][j][k] > -1 &&
	// 						(i < ighost || i > bnx+ighost-1 ||
	// 						j < ighost || j > bny+ighost-1 ||
	// 						k < ighost || k > bnz+ighost-1))
	// 				{
	// 					AddaBox(i,j,k,b0,bxloc,dp,up);
	// 				}
	// 			}
	// 		}
	// 	}
	// }
 
	// void AddaBox(const int & i, const int & j, const int & k, 
	// 	vector<Box> & b0, vector<vector<vector<int> > > & bxloc,
	// 	Point & dp, Point & up)
	// {
	// 	b0.push_back(Box(i-ighost, j-ighost, k-ighost));
	// 	for (Point_iterator p(0,3); p.end(); ++p)
	// 	{
	// 		int i0 = i+p.i-1;
	// 		int j0 = j+p.j-1;
	// 		int k0 = k+p.k-1;
	// 		Point apt = Point(i0,j0,k0);
	// 		if (apt > dp && apt < up)
	// 		{
	// 			if (bxloc[i0][j0][k0] > -1)
	// 			{
	// 				b0.back().setneib(p.i, p.j, p.k, bxloc[i0][j0][k0]);
	// 			}
	// 			else
	// 			{
	// 				if (periodicx && (i==ighost || i==boxnum[0]+ighost-1))
	// 				{
	// 					PeriodicNeib(i0, ighost, boxnum[0]+ighost-1);
	// 				}
	// 				if (periodicy && (j==ighost || j==boxnum[1]+ighost-1))
	// 				{
	// 					PeriodicNeib(j0, ighost, boxnum[1]+ighost-1);
	// 				}
	// 				if (periodicz && (k==ighost || k==boxnum[2]+ighost-1))
	// 				{
	// 					PeriodicNeib(k0, ighost, boxnum[2]+ighost-1);
	// 				}
	// 				Assert(bxloc[i0][j0][k0] > -1, "Error in find a neib!!!", 220);
	// 				if (bxloc[i0][j0][k0] < 0)
	// 				{
	// 					printf("Error in find the x neib!!! Box(%d,%d,%d) is %d neib (%d,%d,%d) is %d i0 %d j0 %d k0 %d\n", 
	// 						i,j,k,bxloc[i][j][k],p.i,p.j,p.k,bxloc[i0][j0][k0],i0,j0,k0);
	// 					MPI_Abort(MPI_COMM_WORLD, 210);
	// 				}
	// 				b0.back().setneib(p.i, p.j, p.k, bxloc[i0][j0][k0]);
	// 			}
	// 		}
	// 	}
	// }

	// void AttachPoints(vector<Box> & b0, vector<mPoint> & tp)
	// {
	// 	int bnx = boxnum[0];
	// 	int bny = boxnum[1];
	// 	int bnz = boxnum[2];

	// 	int ptloc[bnx+1][bny+1][bnz+1];		
	// 	for (int i = 0; i < bnx+1; ++i)
	// 	{
	// 		for (int j = 0; j < bny+1; ++j)
	// 		{
	// 			for (int k = 0; k < bnz+1; ++k)
	// 			{
	// 				tp.push_back(mPoint());
	// 				tp.back().index = tp.size()-1;
	// 				tp.back()[0] = double(i)*dh[0][0];
	// 				tp.back()[1] = double(j)*dh[0][1];
	// 				tp.back()[2] = double(k)*dh[0][2];
	// 				ptloc[i][j][k] = tp.back().index;
	// 			}
	// 		}
	// 	}
	// 	for (int i = 0; i < b0.size(); ++i)
	// 	{
	// 		for (Point_iterator p(0,2); p.end(); ++p)
	// 		{
	// 			int ptx = b0[i].ix()+p.i;
	// 			int pty = b0[i].iy()+p.j;
	// 			int ptz = b0[i].iz()+p.k;
	// 			if (ptx > -1 && ptx < bnx+1 &&
	// 					pty > -1 && pty < bny+1 &&
	// 					ptz > -1 && ptz < bnz+1)
	// 			{
	// 				b0[i].pts[p.i][p.j][p.k] = ptloc[ptx][pty][ptz];
	// 			}
	// 		}
	// 	}
	// }

	// void splitdomain(vector<Box> & b0, vector<mPoint> & tp, vector<Face> & tface)
	// {
	// 	int s0 = 0;
	// 	int bnx = boxnum[0];
	// 	int bny = boxnum[1];
	// 	int bnz = boxnum[2];	
	// 	vector<vector<vector<int> > > 
	// 		bxloc(bnx+2*ighost, vector<vector<int> >
	// 					(bny+2*ighost, vector<int>(bnz+2*ighost, -2)));
	// 	//int bxloc[bnx+2*ighost][bny+2*ighost][bnz+2*ighost];	

	// 	CreateBox(b0, bxloc);
	// 	GiveAFlag("Finish create box!!!", 5);
	// 	AttachPoints(b0, tp);
	// 	GiveAFlag("Finish attach points!!!", 5);

	// 	int bnx0 = bnx+1+ighost;
	// 	int bny0 = bny+1+ighost; 
	// 	int bnz0 = bnz+1+ighost;
	// 	if (periodicx) bnx0 = bnx+ighost;
	// 	if (periodicy) bny0 = bny+ighost;
	// 	if (periodicz) bnz0 = bnz+ighost;
	// 	for (int i = ighost; i < bnx0; ++i)
	// 	{
	// 		for (int j = ighost; j < bny+ighost; ++j)
	// 		{
	// 			for (int k = ighost; k < bnz+ighost; ++k)
	// 			{
	// 				int sbox[2] = {i-1, i};
	// 				tface.push_back(Face());
	// 				tface.back().fnv = 0;
	// 				tface.back().area = dh[0][1]*dh[0][2];
	// 				tface.back().keisa.setvalue(1.0, 0.0, 0.0);
	// 				for (int xx = 0; xx < 2; ++xx)
	// 				{
	// 					if (periodicx) PeriodicNeib(sbox[xx], ighost, bnx+ighost-1);
	// 					tface.back()[xx] = bxloc[sbox[xx]][j][k];
	// 					Assert(bxloc[sbox[xx]][j][k] > -1, "Error in finding x face box", 312);
	// 				}
	// 			}
	// 		}
	// 	}
	// 	for (int i = ighost; i < bnx+ighost; ++i)
	// 	{
	// 		for (int j = ighost; j < bny0; ++j)
	// 		{
	// 			for (int k = ighost; k < bnz+ighost; ++k)
	// 			{
	// 				int sbox[2] = {j-1, j};
	// 				tface.push_back(Face());
	// 				tface.back().fnv = 1;
	// 				tface.back().area = dh[0][0]*dh[0][2];
	// 				tface.back().keisa.setvalue(0.0, 1.0, 0.0);
	// 				for (int xx = 0; xx < 2; ++xx)
	// 				{
	// 					if (periodicy) PeriodicNeib(sbox[xx], ighost, bny+ighost-1);
	// 					tface.back()[xx] = bxloc[i][sbox[xx]][k];
	// 					Assert(bxloc[i][sbox[xx]][k] > -1, "Error in finding y face box", 312);
	// 				}
	// 			}
	// 		}
	// 	}
	// 	for (int i = ighost; i < bnx+ighost; ++i)
	// 	{
	// 		for (int j = ighost; j < bny+ighost; ++j)
	// 		{
	// 			for (int k = ighost; k < bnz0; ++k)
	// 			{
	// 				int sbox[2] = {k-1, k};
	// 				tface.push_back(Face());
	// 				tface.back().fnv = 2;
	// 				tface.back().area = dh[0][0]*dh[0][1];
	// 				tface.back().keisa.setvalue(0.0, 0.0, 1.0);
	// 				for (int xx = 0; xx < 2; ++xx)
	// 				{
	// 					if (periodicz) PeriodicNeib(sbox[xx], ighost, bnz+ighost-1);
	// 					tface.back()[xx] = bxloc[i][j][sbox[xx]];
	// 					Assert(bxloc[i][j][sbox[xx]] > -1, "Error in finding z face box", 312);
	// 				}
	// 			}
	// 		}
	// 	}			
	// }

	void SetFaceProperty(DataArray<Box> & backbox, DataArray<BoxCellGeom> & backgeom);

	void FaceProperty(DataArray<Box> & backbox, DataArray<BoxCellGeom> & backgeom, DataArray<Domainghost> & faceghost)
	{
		SetFaceProperty(backbox, backgeom);
		bool xyzmat;
		vector<int> matchnum(facenum, 0);
		vector<Domainghost> dmg_procs;
		//printf("facenum is %d\n", facenum);
		for (int i = backbox.gps(); i < backbox.gpe(); ++i)
		{
			if (backbox[i].type == Dmghost)
			{
				int x0 = backbox[i].ix();
				int y0 = backbox[i].iy();
				int z0 = backbox[i].iz();				
				int matchfacenum = -1;
				for (int j = 0; j < facenum; ++j)
				{
					xyzmat = false;
					if (inface(x0,y0,z0,j))
					{
						xyzmat = true;
					}
					if (xyzmat)
					{
						int refcell = backbox[i].neib[1+facedir[j][0]][1+facedir[j][1]][1+facedir[j][2]];
						if (inface(backbox[refcell].ix(), backbox[refcell].iy(), backbox[refcell].iz(), j))
						{
							refcell = backbox[refcell].neib[1+facedir[j][0]][1+facedir[j][1]][1+facedir[j][2]];
						}
						if (matchfacenum == -1)
						{
							dmg_procs.push_back(Domainghost(i, j, refcell));						
							matchfacenum = j;
						}
						else
						{
							dmg_procs.back().nbface = j;
							dmg_procs.back().refcell = refcell;
							matchfacenum = j;
						}							
						if (facetype[j] != -99)
						{
							break;
						}	
					}
				}
				//Assert(matchfacenum > -1, "A ghost box has not been matched by a face!!!",429);
				if (matchfacenum > -1)
				{
					++matchnum[matchfacenum];
				}
			}		
		}
		faceghost.Addnew(dmg_procs);
		faceghost.DirectlyReduceNew();
		MPI_Allreduce(MPI_IN_PLACE, &matchnum[0], facenum, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		bool right = true;
		for (int i = 0; i < facenum; ++i)
		{
			if (nrank == 0)
			{
				printf("N%d The face %d finds %d boxes!!!\n", node, i, matchnum[i]);
				if (matchnum[i] < 1 && facetype[i] != -99)
				{
					printf("***Warning***Node %d Face %d does not match any box!!! Please Check the face property!!!\n",node,i);
					right = false;
				}
			}
		}
		if (!right)
		{
			MPI_Abort(MPI_COMM_WORLD, 353);
		}
	}

	bool inface(const int & x0, const int & y0, const int & z0, const int & facej)
	{
		if (faceup[facej][0] > x0 && facedown[facej][0] < x0 &&
				faceup[facej][1] > y0 && facedown[facej][1] < y0 &&
				faceup[facej][2] > z0 && facedown[facej][2] < z0)
		{
			return true;
		}
		else
		{
			return false;
		}
	}

	const int & TellFaceNum() const
	{
		return facenum;
	}

	const Point & FaceNormVector(const int & i) const
	{
		return facedir[i];
	}

	const bool & PeriodicX()
	{
		return periodicx;
	}

	const bool & PeriodicY()
	{
		return periodicy;
	}

	const bool & PeriodicZ()
	{
		return periodicz;
	}

	Point & DmBoxNum()
	{
		return boxnum;
	}

	const int & DomainGhostNum()
	{
		return domain_ghost_num;
	}

	inline const int & FaceType(const int & i0)
	{
		return facetype[i0];
	}

	void splitdomain_tonode(DataArray<Box> & mbox, DataArray<NodePair> & mbkp);

	bool Split_Condition(const int & i0, const int & j0, const int & k0, 
											 const int & local_start, const int & local_end)
	{
		if (i0 >= local_start-ighost && i0 < local_end+ighost)
		{
			return true;
		}
		else
		{
			return false;
		}
	}

	void SetBlockNum();

	void Set_Z_Bound(vector<int> & nodestart, vector<int> & nodeend);

	void DomainBound()
	{
		SetBlockNum();
#ifdef PASSAGE_ANGLE
		if ((blocknum[0] != 1 || blocknum[1] != 1) && nrank == 0)
		{
			printf("For PASSAGE_ANGLE defined, the x/y block number must be 1!!!\n");
			MPI_Abort(MPI_COMM_WORLD, 535);
		}
#endif				
		if (blocknum[0]*blocknum[1]*blocknum[2] != nodenum)
		{
			if (nrank == 0) printf("The input node number is %d but the block number is %d\n", nodenum, blocknum[0]*blocknum[1]*blocknum[2]);
			MPI_Abort(MPI_COMM_WORLD,608);
		}
		alldmbound.resize(nodenum);
		for (int i = 0; i < nodenum; ++i)
		{
			alldmbound[i].resize(3);
		}
		int gblockindex[blocknum[0]][blocknum[1]][blocknum[2]];
		int s = 0;
		for (int i = 0; i < blocknum[0]; ++i)
		{
			for (int j = 0; j < blocknum[1]; ++j)
			{
				for (int k = 0; k < blocknum[2]; ++k)
				{
					gblockindex[i][j][k] = s;
					ArrayOrder_s(0, boxnum[0], alldmbound[s][0][0], alldmbound[s][0][1], blocknum[0], i);
					ArrayOrder_s(0, boxnum[1], alldmbound[s][1][0], alldmbound[s][1][1], blocknum[1], j);
					ArrayOrder_s(0, boxnum[2], alldmbound[s][2][0], alldmbound[s][2][1], blocknum[2], k);
					alldmbound[s][0][2] = alldmbound[s][0][1]-alldmbound[s][0][0];
					alldmbound[s][1][2] = alldmbound[s][1][1]-alldmbound[s][1][0];
					alldmbound[s][2][2] = alldmbound[s][2][1]-alldmbound[s][2][0];
					if (node == s)
					{
						block0 = Point(i,j,k);
						dmbound[0] = alldmbound[s][0];
						dmbound[1] = alldmbound[s][1];
						dmbound[2] = alldmbound[s][2];
					}
					++s;
				}
			}
		}

		for (int i = 0; i < 3; ++i)
		{
			for (int j = 0; j < 3; ++j)
			{
				for (int k = 0; k < 3; ++k)
				{
					int nx = block0[0] + i -1;
					int ny = block0[1] + j -1;
					int nz = block0[2] + k -1;
					if (PeriodicX()) PeriodicNeib(nx, blocknum[0]);
					if (PeriodicY()) PeriodicNeib(ny, blocknum[1]);
					if (PeriodicZ()) PeriodicNeib(nz, blocknum[2]);
					if (nx > -1 && nx < blocknum[0] &&
							ny > -1 && ny < blocknum[1] &&
							nz > -1 && nz < blocknum[2])
					{
						neibblocks[i][j][k] = gblockindex[nx][ny][nz];
					}
					else
					{
						neibblocks[i][j][k] = -1;
					}
				}
			}
		}
	}
	/*executed by one processor*/
	void ConstructBlockInterface(vector<NodePair> & np_ghost, vector<vector<vector<int> > > bxloc, vector<Box> & b0);

	void CheckPeriodicFaceType()
	{
		for (int i = 0; i < 3; ++i)
		{
			if (periodic[i])
			{
				for (int j = 0; j < facenum; ++j)
				{
					if (facedir[j][i] != 0)
					{
						if (facetype[j] != -99)
						{
							printf("Error of the periodic face setting!!! Dir %d is set periodic but the face %d type is not -99!!! "
								"Face vector is (%d,%d,%d)\n",i,j,facedir[j][0],facedir[j][1],facedir[j][2]);
							MPI_Abort(MPI_COMM_WORLD,691);
						}
					}
				}
			}
		}
		for (int i = 0; i < facenum; ++i)
		{
			if (facetype[i] == -99)
			{
				for (int j = 0; j < 3; ++j)
				{
					if (facedir[i][j] != 0)
					{
						if (!periodic[j])
						{
							printf("Error of periodic face setting!!! Face %d type is -99 but dir %d is not periodic!!!\n",i,j);
							MPI_Abort(MPI_COMM_WORLD,708);
						}
					}
				}
			}
		}
	}

	void SetFaceArray()
	{
		faceup.resize(facenum);
		facedown.resize(facenum);
		facetype.resize(facenum);
		facedir.resize(facenum);
		refboxface.resize(facenum);
		faceinput.resize(facenum);
		inputfiles.resize(facenum);
		for (int i = 0; i < facenum; ++i)
		{
			faceinput[i] = false;
		}
	}

private:

	Point lowpt;
	Point highpt;
	Point boxnum;
	bool periodicx;
	bool periodicy;
	bool periodicz;
	int facenum;
	int domain_ghost_num;
	vector<Point> facedown;
	vector<Point> faceup;
	vector<Point> facedir;
	vector<int>   facetype;
	vector<bool>  faceinput;
	vector<string> inputfiles;
	vector<Boxloc> refboxface;
	Point dmbound[3];
	Point blocknum;
	int neibblocks[3][3][3];
	Point block0;
	bool boundaryiswall[3][2];

	vector<vector<Point> > alldmbound;

};

#endif