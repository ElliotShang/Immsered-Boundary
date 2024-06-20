#ifndef _TwoLevelOp_H_
#define _TwoLevelOp_H_

#include "FlowVariables.H"
#include "Boxloc.H"
#include "Point_iterator.H"
#include "Mesh.H"
#include "BCValues.H"

extern int neibson[4];

#define ComptGrad(leftv, rightv, gradx) gradx = (rightv-leftv)*0.5

#ifndef TURBULENCE

#define GradientInaCell(tlop)\
		for (int ti0 = 0; ti0 < 5; ++ti0)\
		{\
			ComptGrad((*tlop.neibflow[0][1][1])[tlop.id0[ti0]], (*tlop.neibflow[2][1][1])[tlop.id0[ti0]], tlop.con_gradx[ti0]);\
			ComptGrad((*tlop.neibflow[1][0][1])[tlop.id0[ti0]], (*tlop.neibflow[1][2][1])[tlop.id0[ti0]], tlop.con_grady[ti0]);\
			ComptGrad((*tlop.neibflow[1][1][0])[tlop.id0[ti0]], (*tlop.neibflow[1][1][2])[tlop.id0[ti0]], tlop.con_gradz[ti0]);\
		}	

#define ComptFineCellData(coarsecell, finecell, sx, sy, sz, tlop)\
		for (int tvn0 = 0; tvn0 < 5; ++tvn0)\
		{\
			finecell[tlop.id0[tnv0]] = (*tlop.neibflow[1][1][1])[tlop.id0[tnv0]]+\
				tlop.fincrem[sx]*tlop.con_gradx[tvn0]+\
				tlop.fincrem[sy]*tlop.con_grady[tvn0]+\
				tlop.fincrem[sz]*tlop.con_gradz[tvn0];\
		}\
		Get_roe_ideal_gas(finecell);\
		Get_E(finecell)

#define ComptFineCellData_2d(coarsecell, finecell, sx, sy, sz, tlop)\
		for (int tvn0 = 0; tvn0 < 5; ++tvn0)\
		{\
			finecell[tlop.id0[tnv0]] = (*tlop.neibflow[1][1][1])[tlop.id0[tnv0]]+\
				tlop.fincrem[sx]*tlop.con_gradx[tvn0]+\
				tlop.fincrem[sy]*tlop.con_grady[tvn0];\
		}\
		Get_roe_ideal_gas(finecell);\
		Get_E(finecell)		

#define AverageCell(coarsevar, tlop)\
		coarsevar.roe = 0.0;\
		coarsevar.u = 0.0;\
		coarsevar.v = 0.0;\
		coarsevar.w = 0.0;\
		coarsevar.e = 0.0;\
		for (Point_iterator ttp(0,2); ttp.end(); ++ttp)\
		{\
			coarsevar.roe += tlop.finedataptr.son[ttp.i][ttp.j][ttp.k]->roe;\
			coarsevar.u += tlop.finedataptr.son[ttp.i][ttp.j][ttp.k]->u*tlop.finedataptr.son[ttp.i][ttp.j][ttp.k]->roe;\
			coarsevar.v += tlop.finedataptr.son[ttp.i][ttp.j][ttp.k]->v*tlop.finedataptr.son[ttp.i][ttp.j][ttp.k]->roe;\
			coarsevar.w += tlop.finedataptr.son[ttp.i][ttp.j][ttp.k]->w*tlop.finedataptr.son[ttp.i][ttp.j][ttp.k]->roe;\
			coarsevar.e += tlop.finedataptr.son[ttp.i][ttp.j][ttp.k]->e*tlop.finedataptr.son[ttp.i][ttp.j][ttp.k]->roe;\
		}\
		coarsevar.u /= coarsevar.roe;\
		coarsevar.v /= coarsevar.roe;\
		coarsevar.w /= coarsevar.roe;\
		coarsevar.e /= coarsevar.roe;\
		coarsevar.roe *= 0.125;\
		Get_p_energy(coarsevar);\
		Get_T_ideal_gas(coarsevar)

#define AverageCell_2d(coarsevar, tlop)\
		coarsevar.roe = 0.0;\
		coarsevar.u = 0.0;\
		coarsevar.v = 0.0;\
		coarsevar.w = 0.0;\
		coarsevar.e = 0.0;\
		for (Point_iterator_2d ttp(0,2); ttp.end(); ++ttp)\
		{\
			coarsevar.roe += tlop.finedataptr.son[ttp.i][ttp.j][ttp.k]->roe;\
			coarsevar.u += tlop.finedataptr.son[ttp.i][ttp.j][ttp.k]->u*tlop.finedataptr.son[ttp.i][ttp.j][ttp.k]->roe;\
			coarsevar.v += tlop.finedataptr.son[ttp.i][ttp.j][ttp.k]->v*tlop.finedataptr.son[ttp.i][ttp.j][ttp.k]->roe;\
			coarsevar.w += tlop.finedataptr.son[ttp.i][ttp.j][ttp.k]->w*tlop.finedataptr.son[ttp.i][ttp.j][ttp.k]->roe;\
			coarsevar.e += tlop.finedataptr.son[ttp.i][ttp.j][ttp.k]->e*tlop.finedataptr.son[ttp.i][ttp.j][ttp.k]->roe;\
		}\
		coarsevar.u /= coarsevar.roe;\
		coarsevar.v /= coarsevar.roe;\
		coarsevar.w /= coarsevar.roe;\
		coarsevar.e /= coarsevar.roe;\
		coarsevar.roe *= 0.25;\
		Get_p_energy(coarsevar);\
		Get_T_ideal_gas(coarsevar)

#else

#define GradientInaCell(tlop)\
		for (int ti0 = 0; ti0 < 5; ++ti0)\
		{\
			ComptGrad((*tlop.neibflow[0][1][1])[tlop.id0[ti0]], (*tlop.neibflow[2][1][1])[tlop.id0[ti0]], tlop.con_gradx[ti0]);\
			ComptGrad((*tlop.neibflow[1][0][1])[tlop.id0[ti0]], (*tlop.neibflow[1][2][1])[tlop.id0[ti0]], tlop.con_grady[ti0]);\
			ComptGrad((*tlop.neibflow[1][1][0])[tlop.id0[ti0]], (*tlop.neibflow[1][1][2])[tlop.id0[ti0]], tlop.con_gradz[ti0]);\
		}\
		ComptGrad(tlop.neibflow[0][1][1]->viseddy, tlop.neibflow[2][1][1]->viseddy, tlop.vis_grad[0]);\
		ComptGrad(tlop.neibflow[1][0][1]->viseddy, tlop.neibflow[1][2][1]->viseddy, tlop.vis_grad[1]);\
		ComptGrad(tlop.neibflow[1][1][0]->viseddy, tlop.neibflow[1][1][2]->viseddy, tlop.vis_grad[2])


#define ComptFineCellData(coarsecell, finecell, sx, sy, sz, tlop)\
		for (int tvn0 = 0; tvn0 < 5; ++tvn0)\
		{\
			finecell[tlop.id0[tvn0]] = (*tlop.neibflow[1][1][1])[tlop.id0[tvn0]]+\
				tlop.fincrem[sx]*tlop.con_gradx[tvn0]+\
				tlop.fincrem[sy]*tlop.con_grady[tvn0]+\
				tlop.fincrem[sz]*tlop.con_gradz[tvn0];\
		}\
		Get_roe_ideal_gas(finecell);\
		Get_E(finecell);\
		finecell.viseddy = tlop.neibflow[1][1][1]->viseddy+\
				tlop.fincrem[sx]*tlop.vis_grad[0]+\
				tlop.fincrem[sy]*tlop.vis_grad[1]+\
				tlop.fincrem[sz]*tlop.vis_grad[2]

#define ComptFineCellData_2d(coarsecell, finecell, sx, sy, sz, tlop)\
		for (int tvn0 = 0; tvn0 < 5; ++tvn0)\
		{\
			finecell[tlop.id0[tvn0]] = (*tlop.neibflow[1][1][1])[tlop.id0[tvn0]]+\
				tlop.fincrem[sx]*tlop.con_gradx[tvn0]+\
				tlop.fincrem[sy]*tlop.con_grady[tvn0];\
		}\
		Get_roe_ideal_gas(finecell);\
		Get_E(finecell);\
		finecell.viseddy = tlop.neibflow[1][1][1]->viseddy+\
				tlop.fincrem[sx]*tlop.vis_grad[0]+\
				tlop.fincrem[sy]*tlop.vis_grad[1]

#define AverageCell(coarsevar, tlop)\
		coarsevar.roe = 0.0;\
		coarsevar.u = 0.0;\
		coarsevar.v = 0.0;\
		coarsevar.w = 0.0;\
		coarsevar.e = 0.0;\
		coarsevar.viseddy = 0.0;\
		for (Point_iterator ttp(0,2); ttp.end(); ++ttp)\
		{\
			coarsevar.roe += tlop.finedataptr.son[ttp.i][ttp.j][ttp.k]->roe;\
			coarsevar.u += tlop.finedataptr.son[ttp.i][ttp.j][ttp.k]->u*tlop.finedataptr.son[ttp.i][ttp.j][ttp.k]->roe;\
			coarsevar.v += tlop.finedataptr.son[ttp.i][ttp.j][ttp.k]->v*tlop.finedataptr.son[ttp.i][ttp.j][ttp.k]->roe;\
			coarsevar.w += tlop.finedataptr.son[ttp.i][ttp.j][ttp.k]->w*tlop.finedataptr.son[ttp.i][ttp.j][ttp.k]->roe;\
			coarsevar.e += tlop.finedataptr.son[ttp.i][ttp.j][ttp.k]->e*tlop.finedataptr.son[ttp.i][ttp.j][ttp.k]->roe;\
			coarsevar.viseddy += tlop.finedataptr.son[ttp.i][ttp.j][ttp.k]->viseddy;\
		}\
		coarsevar.u /= coarsevar.roe;\
		coarsevar.v /= coarsevar.roe;\
		coarsevar.w /= coarsevar.roe;\
		coarsevar.e /= coarsevar.roe;\
		coarsevar.roe *= 0.125;\
		coarsevar.viseddy *= 0.125;\
		Get_p_energy(coarsevar);\
		Get_T_ideal_gas(coarsevar)

#define AverageCell_2d(coarsevar, tlop)\
		coarsevar.roe = 0.0;\
		coarsevar.u = 0.0;\
		coarsevar.v = 0.0;\
		coarsevar.w = 0.0;\
		coarsevar.e = 0.0;\
		coarsevar.viseddy = 0.0;\
		for (Point_iterator_2d ttp(0,2); ttp.end(); ++ttp)\
		{\
			coarsevar.roe += tlop.finedataptr.son[ttp.i][ttp.j][ttp.k]->roe;\
			coarsevar.u += tlop.finedataptr.son[ttp.i][ttp.j][ttp.k]->u*tlop.finedataptr.son[ttp.i][ttp.j][ttp.k]->roe;\
			coarsevar.v += tlop.finedataptr.son[ttp.i][ttp.j][ttp.k]->v*tlop.finedataptr.son[ttp.i][ttp.j][ttp.k]->roe;\
			coarsevar.w += tlop.finedataptr.son[ttp.i][ttp.j][ttp.k]->w*tlop.finedataptr.son[ttp.i][ttp.j][ttp.k]->roe;\
			coarsevar.e += tlop.finedataptr.son[ttp.i][ttp.j][ttp.k]->e*tlop.finedataptr.son[ttp.i][ttp.j][ttp.k]->roe;\
			coarsevar.viseddy += tlop.finedataptr.son[ttp.i][ttp.j][ttp.k]->viseddy;\
		}\
		coarsevar.u /= coarsevar.roe;\
		coarsevar.v /= coarsevar.roe;\
		coarsevar.w /= coarsevar.roe;\
		coarsevar.e /= coarsevar.roe;\
		coarsevar.roe *= 0.25;\
		coarsevar.viseddy *= 0.25;\
		Get_p_energy(coarsevar);\
		Get_T_ideal_gas(coarsevar)

#endif				

#define PairProlongation(a_mesh, ilevel, cbox, fbox, tlop)\
		int tflevel = ilevel+1;\
		for (int i7 = 0; i7 < 7; ++i7)\
		{\
			int tanb = a_mesh.m_level[ilevel].m_box[cbox].neib[tlop.cneib[i7][0]][tlop.cneib[i7][1]][tlop.cneib[i7][2]];\
			Assert(tanb > -1 && tanb < a_mesh.m_level[ilevel].m_box.realsize(), "Error of Prolongation neib!!!", 61);\
			tlop.neibflow[tlop.cneib[i7][0]][tlop.cneib[i7][1]][tlop.cneib[i7][2]] = &a_mesh.m_level[ilevel].m_data[tanb];\
		}\
		GradientInaCell(tlop);\
		for (Point_iterator ttp(0,2); ttp.end(); ++ttp)\
		{\
			ComptFineCellData(a_mesh.m_level[ilevel].m_data[cbox], a_mesh.m_level[tflevel].m_data[fbox.son[ttp.i][ttp.j][ttp.k]],\
				ttp.i,ttp.j,ttp.k,tlop);\
		}

#define PairProlongation_2d(a_mesh, ilevel, cbox, fbox, tlop)\
		int tflevel = ilevel+1;\
		for (int i7 = 0; i7 < 7; ++i7)\
		{\
			int tanb = a_mesh.m_level[ilevel].m_box[cbox].neib[tlop.cneib[i7][0]][tlop.cneib[i7][1]][tlop.cneib[i7][2]];\
			Assert(tanb > -1 && tanb < a_mesh.m_level[ilevel].m_box.realsize(), "Error of Prolongation neib!!!", 61);\
			tlop.neibflow[tlop.cneib[i7][0]][tlop.cneib[i7][1]][tlop.cneib[i7][2]] = &a_mesh.m_level[ilevel].m_data[tanb];\
		}\
		GradientInaCell(tlop);\
		for (Point_iterator_2d ttp(0,2); ttp.end(); ++ttp)\
		{\
			ComptFineCellData_2d(a_mesh.m_level[ilevel].m_data[cbox], a_mesh.m_level[tflevel].m_data[fbox.son[ttp.i][ttp.j][ttp.k]],\
				ttp.i,ttp.j,ttp.k,tlop);\
		}		

#define PairRestriction(a_mesh, ilevel, cbox, fbox, tlop)\
		int tflevel = ilevel+1;\
		for (Point_iterator ttp(0,2); ttp.end(); ++ttp)\
		{\
			tlop.finedataptr.son[ttp.i][ttp.j][ttp.k] = &a_mesh.m_level[tflevel].m_data[fbox.son[ttp.i][ttp.j][ttp.k]];\
		}\
		AverageCell(a_mesh.m_level[ilevel].m_data[cbox], tlop)

#define PairRestriction_2d(a_mesh, ilevel, cbox, fbox, tlop)\
		int tflevel = ilevel+1;\
		for (Point_iterator_2d ttp(0,2); ttp.end(); ++ttp)\
		{\
			tlop.finedataptr.son[ttp.i][ttp.j][ttp.k] = &a_mesh.m_level[tflevel].m_data[fbox.son[ttp.i][ttp.j][ttp.k]];\
		}\
		AverageCell_2d(a_mesh.m_level[ilevel].m_data[cbox], tlop)		

#define Restriction(a_mesh, finesonlevel, finesonindex, newcoarsedata, tlop)\
		for (Point_iterator ttp(0,2); ttp.end(); ++ttp)\
		{\
			int fneib = a_mesh.m_level[finesonlevel].m_box[finesonindex].neib[1+ttp.i][1+ttp.j][1+ttp.k];\
			tlop.finedataptr.son[ttp.i][ttp.j][ttp.k] = &a_mesh.m_level[finesonlevel].m_data[fneib];\
		}\
		AverageCell(newcoarsedata, tlop)

#define Restriction_2d(a_mesh, finesonlevel, finesonindex, newcoarsedata, tlop)\
		for (Point_iterator_2d ttp(0,2); ttp.end(); ++ttp)\
		{\
			int fneib = a_mesh.m_level[finesonlevel].m_box[finesonindex].neib[1+ttp.i][1+ttp.j][1];\
			tlop.finedataptr.son[ttp.i][ttp.j][ttp.k] = &a_mesh.m_level[finesonlevel].m_data[fneib];\
		}\
		AverageCell_2d(newcoarsedata, tlop)

#define Prolongation(a_mesh, coarseboxlevel, coarseboxindex, newfinedata, tlop)\
		for (int i7 = 0; i7 < 7; ++i7)\
		{\
			int tanb = a_mesh.m_level[coarseboxlevel].m_box[coarseboxindex].\
				neib[tlop.cneib[i7][0]][tlop.cneib[i7][1]][tlop.cneib[i7][2]];\
			if (!a_mesh.m_level[coarseboxlevel].m_box.isnorg(tanb))\
			{\
				printf("the Prolongation box (%d,%d,%d) neib (%d,%d,%d) is %d box type is %d is ghost %d\n",\
					a_mesh.m_level[coarseboxlevel].m_box[coarseboxindex].ix(),\
					a_mesh.m_level[coarseboxlevel].m_box[coarseboxindex].iy(),\
					a_mesh.m_level[coarseboxlevel].m_box[coarseboxindex].iz(),\
					tlop.cneib[i7][0],tlop.cneib[i7][1],tlop.cneib[i7][2],tanb,\
					a_mesh.m_level[coarseboxlevel].m_box[coarseboxindex].type,coarseboxindex);\
				MPI_Abort(MPI_COMM_WORLD, 101);\
			}\
			Assert(tanb > -1, "Negative index for data Prolongation!!!", 93);\
			tlop.neibflow[tlop.cneib[i7][0]][tlop.cneib[i7][1]][tlop.cneib[i7][2]] = \
				&a_mesh.m_level[coarseboxlevel].m_data[tanb];\
		}\
		GradientInaCell(tlop);\
		for (Point_iterator ttp(0,2); ttp.end(); ++ttp)\
		{\
			ComptFineCellData(a_mesh.m_level[coarseboxlevel].m_data[coarseboxindex], \
								newfinedata.son[ttp.i][ttp.j][ttp.k], ttp.i,ttp.j,ttp.k, tlop);\
		}

#define Prolongation_2d(a_mesh, coarseboxlevel, coarseboxindex, newfinedata, tlop)\
		for (int i7 = 0; i7 < 7; ++i7)\
		{\
			int tanb = a_mesh.m_level[coarseboxlevel].m_box[coarseboxindex].\
				neib[tlop.cneib[i7][0]][tlop.cneib[i7][1]][tlop.cneib[i7][2]];\
			if (!a_mesh.m_level[coarseboxlevel].m_box.isnorg(tanb))\
			{\
				printf("N%d the Prolongation 2d box (%d,%d,%d) signdis %f neib (%d,%d,%d) is %d box type is %d is ghost %d\n",\
					node, a_mesh.m_level[coarseboxlevel].m_box[coarseboxindex].ix(),\
					a_mesh.m_level[coarseboxlevel].m_box[coarseboxindex].iy(),\
					a_mesh.m_level[coarseboxlevel].m_box[coarseboxindex].iz(),\
					a_mesh.m_level[coarseboxlevel].m_box[coarseboxindex].pair.signdis,\
					tlop.cneib[i7][0],tlop.cneib[i7][1],tlop.cneib[i7][2],tanb,\
					a_mesh.m_level[coarseboxlevel].m_box[coarseboxindex].type,coarseboxindex);\
				MPI_Abort(MPI_COMM_WORLD, 101);\
			}\
			Assert(tanb > -1, "Negative index for data Prolongation!!!", 93);\
			tlop.neibflow[tlop.cneib[i7][0]][tlop.cneib[i7][1]][tlop.cneib[i7][2]] = \
				&a_mesh.m_level[coarseboxlevel].m_data[tanb];\
		}\
		GradientInaCell(tlop);\
		for (Point_iterator_2d ttp(0,2); ttp.end(); ++ttp)\
		{\
			ComptFineCellData_2d(a_mesh.m_level[coarseboxlevel].m_data[coarseboxindex], \
								newfinedata.son[ttp.i][ttp.j][ttp.k], ttp.i, ttp.j, 0, tlop);\
		}			

/*Two level operation must consider the coordinate transformation*/
/*Two level operation must consider the coordinate transformation*/
/*Two level operation must consider the coordinate transformation*/
/*Two level operation must consider the coordinate transformation*/
/*Two level operation must consider the coordinate transformation*/
/*Two level operation must consider the coordinate transformation*/
/*Two level operation must consider the coordinate transformation*/
/*Two level operation must consider the coordinate transformation*/
class TwoLevelOp
{
private:
	//Face * finemcptr[3][2][4];
	//Boxson<double> finevolume;
	//CellMCoef * neibmcoef[3][3][3];

	//vector<FlowVariables> tempfv;
	//vector<CellMCoef>     tempmc;

	

	// static vector<int> neibrefx;
	// static vector<int> neibrefy;
	// static vector<int> neibrefz;
	// static vector<int> neibptx;
	// static vector<int> neibpty;
	// static vector<int> neibptz;
	// static int off_ffneibx[3];
	// static int off_ffneiby[3];
	// static int off_ffniebz[3];

public:

	FlowVariables * neibflow[3][3][3];
	FlowVec neibflowvec[3][3][3];
	Boxson<FlowVariables *> finedataptr;
	FlowVec con_gradx, con_grady, con_gradz;
	Point cneib[7];
	double fincrem[2];
#ifdef TURBULENCE
	Pointxyz vis_grad;
#endif
	int id0[5];		

	TwoLevelOp()
	{
		fincrem[0] = -0.25;
		fincrem[1] = 0.25;
		cneib[0] = Point(1,1,0);
		cneib[1] = Point(1,1,2);
		cneib[2] = Point(1,0,1);
		cneib[3] = Point(1,2,1);
		cneib[4] = Point(0,1,1);
		cneib[5] = Point(2,1,1);
		cneib[6] = Point(1,1,1);
		id0[0] = 1;
		id0[1] = 2;
		id0[2] = 3;
		id0[3] = 4;
		id0[4] = 5;
	}

	// inline void cleartemparr()
	// {
	// 	tempfv.resize(0);
	// 	tempmc.resize(0);
	// }

	// void AverageMCoef(CellMCoef & coarsemc)
	// {
	// 	coarsemc.zero();
	// 	double tv = 0.0;
	// 	for (Point_iterator pt(0,2); pt.end(); ++pt)
	// 	{
	// 		tv += finevolume.son[pt.i][pt.j][pt.k];
	// 	}
	// 	for (int i = 0; i < 3; ++i)
	// 	{
	// 		double fa[2]  = {0.0,0.0};
	// 		for (int fi = 0; fi < 2; ++fi)
	// 		{				
	// 			for (int sfi = 0; sfi < 4; ++sfi)
	// 			{
	// 				fa[fi] += finemcptr[i][fi][sfi]->area;
	// 			}
	// 			for (int sfi = 0; sfi < 3; ++sfi)
	// 			{
	// 				coarsemc[i][sfi] += finemcptr[i][fi][0]->keisa[sfi]/2.0;
	// 			}
	// 		}
	// 		double carea = (fa[0]+fa[1])/2.0;
	// 		carea /= tv;
	// 		for (int sfi = 0; sfi < 3; ++sfi)
	// 		{
	// 			coarsemc[i][sfi] *= carea;
	// 		}
	// 	}
	// }
	// void CollectFineDataPtr(Mesh & a_mesh, const Boxloc & fineneib, 
	// 	const Point_iterator & p)
	// {
	// 	Boxloc fineneib0;
	// 	for (Point_iterator q(0,2); q.end(); ++q)
	// 	{
	// 		fineneib0 = a_mesh.BoxNeib(fineneib, 1+dn[p.i]*q.i, 1+dn[p.j]*q.j, 1+dn[p.k]*q.k);
	// 		//printf("fineneib L%dB%d other L%dB%d\n", fineneib.level, fineneib.index, fineneib0.level, fineneib0.index);
	// 		Assert(fineneib0.level == fineneib.level, "Erro in collect fine data pointer!!!", 131);
	// 		finedataptr.son[q.i][q.j][q.k] = a_mesh.LevelBoxDataPtr(fineneib0.level, fineneib0.index);
	// 	}
	// }

	// void CollectFineMcPtr(Mesh & a_mesh, const Boxloc & fineneib, 
	// 	const Point_iterator & p)
	// {
	// 	Boxloc fineneib0;
	// 	for (Point_iterator q(0,2); q.end(); ++q)
	// 	{
	// 		fineneib0 = a_mesh.BoxNeib(fineneib, 1+dn[p.i]*q.i, 1+dn[p.j]*q.j, 1+dn[p.k]*q.k);
	// 		Assert(fineneib0.level == fineneib.level, "Erro in collect fine mc pointer!!!", 170);
	// 		//finedataptr.son[q.i][q.j][q.k] = a_mesh.LevelBoxDataPtr(fineneib0.level, fineneib0.index);
	// 		finemcptr[0][q[0]][q[1]*2+q[2]] = &a_mesh.
	// 			Boxsideface(fineneib0.level, fineneib0.index, 0, q.i);
	// 		finemcptr[1][q[1]][q[0]*2+q[2]] = &a_mesh.
	// 			Boxsideface(fineneib0.level, fineneib0.index, 1, q.j);
	// 		finemcptr[2][q[2]][q[0]*2+q[1]] = &a_mesh.
	// 			Boxsideface(fineneib0.level, fineneib0.index, 2, q.j);
	// 		finevolume.son[q.i][q.j][q.k] = a_mesh.boxv(fineneib0);
	// 	}
	// }

	// void ComptFineCellMCoef(CellMCoef & coarsecell, CellMCoef & finecell, 
	// 	const int & sx, const int & sy, const int & sz)
	// {
	// 	for (int i = 0; i < 3; ++i)
	// 	{
	// 		for (int j = 0; j < 3; ++j)
	// 		{
	// 			finecell[i][j] = coarsecell[i][j]*2.0;
	// 		}
	// 	}
	// }

	/*This function is used for interpolation from a coarser box.*/
	/*There must be no less than two boxes of level i+1 between*/
	/*level i+2 and level i in one direction*/	
	/*ci cj ck is the cell location in the box*/
	// void ComptCellFlowGradient(Mesh & a_mesh, const Boxloc & coarserbox)
	// {
	// 	for (int cn0 = 0; cn0 < 7; ++cn0)
	// 	{
	// 		// p = cneib[cn0];
	// 		// const Boxloc * aimneib = &a_mesh.BoxNeib(coarserbox, p.i, p.j, p.k);
	// 		// if (aimneib->issamebox(coarserbox.level))
	// 		// {
	// 		// 	neibflow[p.i][p.j][p.k] = &a_mesh.LevelBoxData(aimneib->level, aimneib->index);
	// 		// 	//neibflow[p.i][p.j][p.k]->showdata("grad same level");
	// 		// }
	// 		// else if (aimneib->isfinebox(coarserbox.level))
	// 		// {
	// 		// 	tempfv.push_back(FlowVariables());
	// 		// 	CollectFineDataPtr(a_mesh, *aimneib, p);
	// 		// 	AverageCell(tempfv.back());
	// 		// 	neibflow[p.i][p.j][p.k] = &tempfv.back();
	// 		// 	//neibflow[p.i][p.j][p.k]->showdata("grad fine ave");
	// 		// }
	// 		// else if (aimneib->level < -1)
	// 		// {
	// 		// 	BCValues::BCBoxTreat(*aimneib, tempfv, 1, &a_mesh.LevelBoxData(coarserbox));
	// 		// 	neibflow[p.i][p.j][p.k] = &tempfv.back();
	// 		// 	//neibflow[p.i][p.j][p.k]->showdata("grad bd");
	// 		// }
	// 		// else if (aimneib->iscoarsebox(coarserbox.level))
	// 		// {
	// 		// 	printf("A coarse neib is not allowed when computing the cell gradient!!!\n");
	// 		// 	MPI_Abort(MPI_COMM_WORLD, 220);
	// 		// }
	// 		neibflowvec[cneib[cn0][0]][cneib[cn0][1]][cneib[cn0][2]] = a_mesh.LevelBoxDataPtr
	// 		neibflow[p.i][p.j][p.k]->Get_Cons_Vars(neibflowvec[p.i][p.j][p.k]);
	// 	}
	// 	GradientInaCell();
	// }

	// void McFromFineLevelBox(Mesh & a_mesh, vector<CellMCoef> & curmcbox,
	// 	const Boxloc & fineneib, const Point_iterator & p, const int & ighostnum)
	// {
	// 	tempmc.resize(0);
	// 	Boxloc prebox = fineneib;
	// 	for (int i = 0; i < ighostnum; ++i)
	// 	{
	// 		CollectFineMcPtr(a_mesh, prebox, p);
	// 		//curbox.push_back(FlowVariables());
	// 		curmcbox.push_back(CellMCoef());
	// 		//AverageCell(curbox.back());
	// 		AverageMCoef(curmcbox.back());
	// 		for (int ip0 = 0; ip0 < 2; ++ip0)
	// 		{
	// 			prebox = a_mesh.BoxNeib(prebox, p.i, p.j, p.k);
	// 		}
	// 	}
	// }

	// void DataFromFineLevelBox(Mesh & a_mesh, vector<FlowVariables> & curbox,
	// 	const Boxloc & fineneib, const Point_iterator & p, const int & ighostnum)
	// {
	// 	tempfv.resize(0);
	// 	Boxloc prebox = fineneib;
	// 	for (int i = 0; i < ighostnum; ++i)
	// 	{
	// 		curbox.push_back(FlowVariables());
	// 		if (prebox.issamebox(fineneib.level))
	// 		{
	// 			CollectFineDataPtr(a_mesh, prebox, p);				
	// 			AverageCell(curbox.back());
	// 		}
	// 		else
	// 		{
	// 			Assert(prebox.iscoarsebox(fineneib.level), 
	// 				"Error in get data from fine neib!!!", 287);
	// 			curbox.back() = a_mesh.LevelBoxData(prebox);
	// 		}
	// 		for (int ip0 = 0; ip0 < 2; ++ip0)
	// 		{
	// 			prebox = a_mesh.BoxNeib(prebox, p.i, p.j, p.k);
	// 		}
	// 	}
	// }


// 	static void DataFromSameLevelBox(GhostDataBox<FlowVariables> & curbox, GhostDataBox<CellMCoef *> & curmcbox, 
// 		const DataBox<FlowVariables> * fineneib, DataBox<CellMCoef> * finemcoef, 
// 		const Point_iterator & p)
// 	{
// 		for (int i = nis[p.i]; i < nie[p.i]; ++i)
// 		{
// 			for (int j = njs[p.j]; j < nje[p.j]; ++j)
// 			{
// 				for (int k = nks[p.k]; k < nke[p.k]; ++k)
// 				{
// 					curbox.var[i][j][k] = fineneib->var[i+myi0[p.i]][j+myj0[p.j]][k+myk0[p.k]];
// 					curmcbox.var[i][j][k] = &finemcoef->var[i+myi0[p.i]][j+myj0[p.j]][k+myk0[p.k]];
// #ifdef DEBUG								
// 					curbox.befilled[i][j][k] = true;
// #endif					
// 				}
// 			}
// 		}
// 	}

	// void McFromCoarseLevelBox(Mesh & a_mesh, vector<CellMCoef> & curmcbox,
	// 	const Boxloc & curloc, const Point_iterator & p, const Boxloc & coarsebox, const int & ighostnum)
	// {
	// 	tempmc.resize(0);
	// 	int boxx = a_mesh.LevelBox(curloc).ix()%2;
	// 	int boxy = a_mesh.LevelBox(curloc).iy()%2;
	// 	int boxz = a_mesh.LevelBox(curloc).iz()%2;
	// 	int sx0 = neibson[boxx+p.i];
	// 	int sy0 = neibson[boxy+p.j];
	// 	int sz0 = neibson[boxz+p.k];
	// 	Boxloc prebox = coarsebox;
	// 	// ComptCellGradient_Mc(a_mesh, prebox);
	// 	for (int ip = 0; ip < ighostnum; ++ip)
	// 	{
	// 		int sn0 = ip%2;
	// 		int n0 = (ip - sn0)/2;
	// 		for (int nn0 = 0; nn0 < n0; ++nn0)
	// 		{
	// 			prebox = a_mesh.BoxNeib(prebox, p.i, p.j, p.k);
	// 			// ComptCellGradient_Mc(a_mesh, prebox);
	// 		}
	// 		//curbox.push_back(FlowVariables());
	// 		curmcbox.push_back(CellMCoef());
	// 		//ComptFineCellData(a_mesh.LevelBoxData(prebox), curbox.back(), sx0, sy0, sz0);
	// 		ComptFineCellMCoef(a_mesh.boxmc(prebox.level, prebox.index), curmcbox.back(), 
	// 		sx0, sy0, sz0);
	// 		if (p.i != 1)
	// 		{
	// 			sx0 = 1-sx0;
	// 		}
	// 		else if (p.j != 1)
	// 		{
	// 			sy0 = 1-sy0;
	// 		}
	// 		else if (p.k != 1)
	// 		{
	// 			sz0 = 1-sz0;
	// 		}			
	// 	}
	// }

	// void DataFromCoarseLevelBox(Mesh & a_mesh, vector<FlowVariables> & curbox,
	// 	const Boxloc & curloc, const Point_iterator & p, const Boxloc & coarsebox, const int & ighostnum)
	// {
	// 	tempfv.resize(0);
	// 	int boxx = a_mesh.LevelBox(curloc).ix()%2;
	// 	int boxy = a_mesh.LevelBox(curloc).iy()%2;
	// 	int boxz = a_mesh.LevelBox(curloc).iz()%2;
	// 	int sx0 = neibson[boxx+p.i];
	// 	int sy0 = neibson[boxy+p.j];
	// 	int sz0 = neibson[boxz+p.k];
	// 	Boxloc prebox = coarsebox;
	// 	ComptCellFlowGradient(a_mesh, prebox);
	// 	for (int ip = 0; ip < ighostnum; ++ip)
	// 	{
	// 		int sn0 = ip%2;
	// 		int n0 = (ip - sn0)/2;
	// 		for (int nn0 = 0; nn0 < n0; ++nn0)
	// 		{
	// 			prebox = a_mesh.BoxNeib(prebox, p.i, p.j, p.k);
	// 			ComptCellFlowGradient(a_mesh, prebox);
	// 		}
	// 		curbox.push_back(FlowVariables());
	// 		ComptFineCellData(a_mesh.LevelBoxData(prebox), curbox.back(), sx0, sy0, sz0);
	// 		if (p.i != 1)
	// 		{
	// 			sx0 = 1-sx0;
	// 		}
	// 		else if (p.j != 1)
	// 		{
	// 			sy0 = 1-sy0;
	// 		}
	// 		else if (p.k != 1)
	// 		{
	// 			sz0 = 1-sz0;
	// 		}						
	// 	}		
	// }
};

#endif