#include "Body.H"
#include "Cell_iterator.H"
#include "Filt.H"

double Body::c_length = 1.0;
double Body::hg_ratio = 1.0;
double Body::force_hg_ratio[2] = {1.0, 2.0};
int Body::lastlevel = -1;
int Body::bodynum = 0;

double dt00 = 0.0;

int Cell_iterator::nbse[3][2] = {{0,2},{1,2},{1,3}};


	void Body::ComptPatchNv(Surfpatch & apatch)
	{
		Pointxyz n1 = (allpoint[apatch.corpt[1]] - allpoint[apatch.corpt[0]]);
		Pointxyz n2 = (allpoint[apatch.corpt[2]] - allpoint[apatch.corpt[0]]);
		apatch.nv = n1.cross(n2);
		apatch.nv.normalize();
		//Pointxyz stoc = (apatch.pc - bodycenter);
		if (apatch.reverse_nv_dir)
		{
			apatch.nv *= -1.0;
		}
	}

void Body::InitIBCellParams(Mesh & amesh, vector<Body> & abody)
{
	Pointxyz patv;
	IBCell * aibc;
	int pps = amesh.m_dis.ps();
	int ppe = amesh.m_dis.pe();
	MPI_Barrier(MPI_COMM_WORLD);
	for (int i = pps; i < ppe; ++i)
	{
		aibc = &amesh.ibcell(i);
		if (!amesh.m_level[lastlevel].m_box[aibc->ci].solid && amesh.m_level[lastlevel].m_box.isnormal(aibc->ci))
		{
			abody[amesh.m_level[lastlevel].m_box[aibc->ci].pair.body].Surface_Point_Velocity(aibc->tp, patv);
			amesh.LevelBoxData(lastlevel, aibc->ci).u = patv[0];
			amesh.LevelBoxData(lastlevel, aibc->ci).v = patv[1];
			amesh.LevelBoxData(lastlevel, aibc->ci).w = patv[2];
		}
	}
	MPI_Barrier(share_comm);
	GiveAFlag("Start to find the data exchange location...", 5);
	// amesh.DataExchange(amesh.m_level[lastlevel].nd, lastlevel);
	amesh.DataExchange(amesh.nd_mdis, lastlevel);
}
void Body::ComptIBHglength(Mesh & amesh, const int & ilevel, const int & ibox, double & hglength, Pointxyz & patnv)
{
	int pt0 = amesh.m_level[ilevel].m_box[ibox].pts[0][0][0];
	int pt1 = amesh.m_level[ilevel].m_box[ibox].pts[0][0][1];
	int pt2 = amesh.m_level[ilevel].m_box[ibox].pts[0][1][0];
	int pt3 = amesh.m_level[ilevel].m_box[ibox].pts[1][0][0];
	Pointxyz ps1 = amesh.m_level[ilevel].m_point[pt1].xyz - amesh.m_level[ilevel].m_point[pt0].xyz;
	Pointxyz ps2 = amesh.m_level[ilevel].m_point[pt2].xyz - amesh.m_level[ilevel].m_point[pt0].xyz;
	Pointxyz ps3 = amesh.m_level[ilevel].m_point[pt3].xyz - amesh.m_level[ilevel].m_point[pt0].xyz;
	double pl1 = ps1.length();
	double pl2 = ps2.length();
	double pl3 = ps3.length();
	ps1 /= pl1;
	ps2 /= pl2;
	ps3 /= pl3;
	double ds1 = pl1/(abs(ps1.dot(patnv)));
	double ds2 = pl2/(abs(ps2.dot(patnv)));
	double ds3 = pl3/(abs(ps3.dot(patnv)));
	hglength = min(ds1, ds2);
	hglength = min(hglength, ds3)*0.8;
	// Assert(pt0 > 0 && amesh.m_level[ilevel].m_point.size(), "The IB cell close box point 0 is not in range!!!", 39);
	// Assert(pt1 > 0 && amesh.m_level[ilevel].m_point.size(), "The IB cell close box point 1 is not in range!!!", 39);
	// Pointxyz dx(abs(amesh.m_level[ilevel].m_point[pt1][0] - amesh.m_level[ilevel].m_point[pt0][0]),
	// 			abs(amesh.m_level[ilevel].m_point[pt1][1] - amesh.m_level[ilevel].m_point[pt0][1]),
	// 			abs(amesh.m_level[ilevel].m_point[pt1][2] - amesh.m_level[ilevel].m_point[pt0][2]));
	// hglength = min(dx[0]/abs(patnv[0]), dx[1]/abs(patnv[1]));
	// hglength = min(hglength, dx[2]/abs(patnv[2]));
	// hglength = -1.0;
	// for (Point_iterator p(0,2); p.end(); ++p)
	// {
	// 	int ix = p.i*2;
	// 	int iy = p.j*2;
	// 	int iz = p.k*2;
	// 	int an0 = amesh.m_level[ilevel].m_box[ibox].neib[ix][iy][iz];
	// 	Pointxyz dr = amesh.m_level[ilevel].m_geom[an0].boxcenter - amesh.m_level[ilevel].m_geom[ibox].boxcenter;
	// 	hglength = max(dr.dot(patnv), hglength);
	// }
	// hglength *= 0.8;
}
// void Body::IBConditions(Mesh & amesh, vector<Body> & abody)
// {
// #ifdef SHOWTIME
// 	double start_time = MPI_Wtime();
// #endif	
// 	if (moveflag)
// 	{
// 		MovingBody(amesh, abody);
// 	}	
// 	Pointxyz hgpt, patv;
// 	double hglength;
// 	FlowVariables hgfv;
// 	int end0 = amesh.ibarray().pe();
// 	for (int i = amesh.ibarray().ps(); i < end0; ++i)
// 	{
// 		//printf("IB cell %d distance is %f\n", i, m_dis[i].distance);
// 		if (amesh.m_dis[i].boxpatch[0].signdis >= 0.0 && amesh.m_level[lastlevel].m_box.isnormal(amesh.m_dis[i].ci))
// 		{
// 			IBCell * aibc = &amesh.m_dis[i];
// 			ComptIBHglength(amesh, lastlevel, amesh.m_dis[i].ci, hglength, abody[aibc->boxpatch[0].bi].patch[aibc->boxpatch[0].pi].nv);
// 			//hglength = hg_ratio*dh[lastlevel][3];
// 			Pointxyz & ibpt = amesh.bc(lastlevel, aibc->ci);
// 			hgpt = ibpt + abody[aibc->boxpatch[0].bi].patch[aibc->boxpatch[0].pi].nv*hglength;
// 			//printf("start HGCellVars for IBCell %d hanging node distance is %f\n", i, hglength);
// 			HGCellVars(amesh, 
// 				aibc->ci, 
// 				hgpt,
// 				hgfv,
// 				abody[aibc->boxpatch[0].bi].patch[aibc->boxpatch[0].pi].nv,
// 				Pointxyz(amesh.m_dis[i].boxpatch[0].signdis, 0.0, 0.0),
// 				hglength);
// 			// hgfv.showdata("hgdata");
// 			//printf("Finish HGCellVars for IBCell %d hanging node distanc is %f\n", i, hglength);

// 			Assert(!(hgfv.hasnan(lastlevel, i, "HGPoint")), 
// 				"Error in HGPoint parameters!!!", 36);
// 			Assert(hgfv.roeisright("HGPoint Density"), 
// 				"Error in HGPoint Density!!!", 38);

// 			abody[aibc->boxpatch[0].bi].Surface_Point_Velocity(aibc->tp, patv);

// 			LinearDist(hglength, 
// 				hgfv, 
// 				amesh.ibcell(i), 
// 				amesh.LevelBoxData(lastlevel,aibc->ci), 
// 				patv,
// 				abody[aibc->boxpatch[0].bi].patch[aibc->boxpatch[0].pi].nv);
// 			Assert(!(aibc->fv.hasnan(lastlevel, i, "IBPoint")), 
// 				"Error in IBPoint parameters!!!", 42);
// 			Assert(aibc->fv.roeisright("IBPoint Density"), 
// 				"Error in IBPoint Density!!!", 45);
// 			//printf("Finish LinearDist for IBCell %d\n", i);
// 		}
// 	}
// 	MPI_Barrier(share_comm);
// 	GiveAFlag("Getting the IB cells' values!!!", 5);
// 	for (int i = amesh.ibarray().ps(); i < amesh.ibarray().pe(); ++i)
// 	{
// 		if (amesh.ibcell(i).boxpatch[0].signdis >= 0.0 && amesh.m_level[lastlevel].m_box.isnormal(amesh.m_dis[i].ci))
// 		{
// 			amesh.m_level[lastlevel].m_data[amesh.ibcell(i).ci] = amesh.ibcell(i).fv;
// 		}
// 		else if (amesh.ibcell(i).boxpatch[0].signdis < 0.0)
// 		{
// 			amesh.m_level[lastlevel].m_data[amesh.ibcell(i).ci] = initvar;
// 		}
// 	}
// 	MPI_Barrier(MPI_COMM_WORLD);
// 	//amesh.DataExchange(amesh.m_level[lastlevel].nd, lastlevel);
// 	amesh.DataExchange(amesh.nd_mdis, lastlevel);
// #ifdef SHOWTIME	
// 	double end_time = MPI_Wtime();
// 	step_ib_time = end_time - start_time;
// 	total_ib_time += step_ib_time;	
// #endif	
// }

// void Body::ReverseDist(const double & ibtowall, const double & hgtoib, 
// 	FlowVariables & hgfv, FlowVariables & ibfv_n1, FlowVariables & ibfv_n0, Pointxyz & wallvel)
// {
// 	if (ibtowall < 0.0000001)
// 	{
// 		ibfv.u = wallvel.u;
// 		ibfv.v = wallvel.v;
// 		ibfv.w = wallvel.w;
// 	}
// 	else
// 	{

// 	}
// }

void Body::LinearDist(FlowVariables & hgfv, IBCell & ibpt, 
	FlowVariables & ibfv_n0, Pointxyz & wallvel, Pointxyz & fnv, double & ibsigndis)
{
	double rt = ibsigndis + ibpt.hgc.hgdis;
	double r0 = ibsigndis/rt;
	Pointxyz du(hgfv.u-wallvel[0], hgfv.v-wallvel[1], hgfv.w-wallvel[2]);
	ibpt.fv.u = wallvel[0]+du[0]*r0;
	ibpt.fv.v = wallvel[1]+du[1]*r0;
	ibpt.fv.w = wallvel[2]+du[2]*r0;
	double nvel_n0 = normalvel(ibfv_n0, fnv);
	double nvel_n1 = normalvel(ibpt.fv, fnv);
	//printf("Normal velocity old is %f and new is %f\n", nvel_n0, nvel_n1);
	ibpt.pre_grad = -(nvel_n1-nvel_n0)/dt*ibfv_n0.roe;
	// ibpt.fv.p = hgfv.p - ibpt.pre_grad*hgtoib;
	ibpt.fv.p = hgfv.p;
	if (ibpt.fv.p < 0.0) ibpt.fv.p = hgfv.p;
	ibpt.fv.T = hgfv.T;
	Get_roe_ideal_gas(ibpt.fv);
	Get_E(ibpt.fv);
	//ibpt.fv.showdata("IBCell data");
	ibpt.u_grad = (du - fnv*(du.dot(fnv)))/rt;
}


// void Body::ReverseFunc(const double & hg, double & ib, const double & wv)
// {

// }

int Body::HGCellVars(Mesh & amesh, const int & ibbc, Pointxyz & hgpt, FlowVariables & hgfv, Pointxyz & nmv, const Pointxyz & atp, const double & hl)
{

	FlowVariables * nfv[2][2][2];
	double distobc[2][2][2];
	double distobc_0[2][2][2];
	int closepti[3];
	double mindis = 10.0;
	int sc; 
	//printf("HGPoint is (%f, %f, %f)\n", hgpt[0], hgpt[1], hgpt[2]);
	for (Point_iterator p(0,2); p.end(); ++p)
	{
		//double adis = amesh.distoboxpt(lastlevel, ibbc, p.i, p.j, p.k, hgpt)*rr0;
		int anb9 = amesh.m_level[lastlevel].m_box[ibbc].neib[2*p.i][2*p.j][2*p.k];
		double adis = amesh.distobc(lastlevel, anb9, hgpt);
		if (adis < mindis)
		{
			closepti[0] = p.i;
			closepti[1] = p.j;
			closepti[2] = p.k;
			mindis = adis;
		}
	}
#ifdef DEBUG
	if (mindis > 9.0)
	{
		printf("Did not find the close box point for a hanging node!!!\n");
		MPI_Abort(MPI_COMM_WORLD, 203);
	}
#endif	
	//printf("close point is (%d,%d,%d) distance is %f\n", closepti[0], closepti[1], closepti[2], mindis);
	//Assert(mindis < 0.8660255*dh[lastlevel][3], "fail to locate the ref point!!!", 154);
	double dissum = 0.0;
	int infnum_n = 0;
	int infnum_p = 0;
	int notinfnum = 0;
	for (int x0 = 0; x0 < 2; ++x0)
	{
		for (int y0 = 0; y0 < 2; ++y0)
		{
#if DIM > 2			
			for (int z0 = 0; z0 < 2; ++z0)
			{
				int nk = closepti[2]+z0;
#else 
			for (int z0 = 0; z0 < 1; ++z0)
			{
				int nk = 1;
#endif
				int ni = closepti[0]+x0;
				int nj = closepti[1]+y0;		
				int aneib = amesh.BoxNeib(lastlevel, ibbc, ni,nj,nk);
				Assert(amesh.m_level[lastlevel].m_box.isnorg(aneib), "The hg near box index must in range!!!", 209);
				Assert(amesh.m_level[lastlevel].m_box.isnormal(aneib) || (amesh.m_level[lastlevel].m_box[aneib].type == Blockghost &&
					amesh.InTransferRange(lastlevel, aneib) || amesh.m_level[lastlevel].m_box[aneib].type == Dmghost), "IBCell error!!!", 56);
				if (amesh.BoxInfectIndex(aneib) > -1)
				{
					//Assert(amesh.m_dis[amesh.infectbox[aneib]].distance >= 0.0, "Error in hg cell distance!!!", 168);
					if (amesh.m_level[lastlevel].m_box[aneib].pair.signdis >= 0.0)
					{
						VALIDPOINT:;
						nfv[x0][y0][z0] = amesh.LevelBoxDataPtr(lastlevel, aneib);
#ifdef DEBUG
						if (nfv[x0][y0][z0]->hasnan(aneib, aneib, "ref neib point"))
						{
							MPI_Abort(MPI_COMM_WORLD, 176);
						}
#endif									
						Pointxyz d0 = amesh.bc(lastlevel, aneib)-hgpt;
#if DIM > 2						
						distobc[x0][y0][z0] = d0.length();
#ifdef DEBUG						
						distobc_0[x0][y0][z0] = distobc[x0][y0][z0];
#endif						
#else
						distobc[x0][y0][z0] = d0.length_2d();
#ifdef DEBUG						
						distobc_0[x0][y0][z0] = distobc[x0][y0][z0];

#endif
#endif
						// if (distobc[x0][y0][z0] < dh[lastlevel][3]*rr0*0.01)
						// {
						// 	hgfv = *nfv[x0][y0][z0];
						// 	goto FINISH;
						// }
						double ds0 = 1.0/distobc[x0][y0][z0];
						distobc[x0][y0][z0] = ds0*ds0;
						dissum += distobc[x0][y0][z0];
						++infnum_p;
					}
					else
					{
						distobc[x0][y0][z0] = -0.05;
						++infnum_n;
#ifdef DEBUG						
						Pointxyz antobc = amesh.bc(lastlevel, aneib) - hgpt;
						distobc_0[x0][y0][z0] = antobc.length();
#endif						
					}
				}
				else if (amesh.m_level[lastlevel].m_box[aneib].type != Dmghost)
				{
					++notinfnum;
					goto VALIDPOINT;
				}												
			}
		}
	}
#ifdef DEBUG
#if DIM > 2	
	double min_tobc = 10.0;
	Point anb(-1,-1,-1);
	for (Point_iterator p(0,2); p.end(); ++p)
	{
		if (distobc_0[p.i][p.j][p.k] < min_tobc)
		{
			min_tobc = distobc_0[p.i][p.j][p.k];
			anb[0] = p.i+closepti[0];
			anb[1] = p.j+closepti[1];
			anb[2] = p.k+closepti[2];
		}
	}
	Assert(anb[0] != -1, "Error in HGCellVars 297!!!", 297);
	Assert(min_tobc > -0.00001, "Error in HGCellVars 301!!!", 301);
	int as0 = amesh.m_level[lastlevel].m_box[ibbc].neib[anb[0]][anb[1]][anb[2]];
	double mindis0 = 100.0;
	int as1 = -1;
	for (Point_iterator p(0,3); p.end(); ++p)
	{
		int aneib = amesh.m_level[lastlevel].m_box[as0].neib[p.i][p.j][p.k];
		Pointxyz abnbc = amesh.m_level[lastlevel].m_geom[aneib].boxcenter - hgpt;
		double len_abnbc = abnbc.length();
		if (len_abnbc < mindis0)
		{
			mindis0 = len_abnbc;
			as1 = aneib;
		}
	}
	if (as0 != as1)
	{
		printf("The HGCellVars mini distance is %f to box (%d,%d,%d) orient (%d,%d,%d) but a lower value is %f to box (%d,%d,%d)\n", min_tobc,
			amesh.m_level[lastlevel].m_box[as0].ix(), amesh.m_level[lastlevel].m_box[as0].iy(), amesh.m_level[lastlevel].m_box[as0].iz(),
			anb[0], anb[1], anb[2],
			mindis0,
			amesh.m_level[lastlevel].m_box[as1].ix(), amesh.m_level[lastlevel].m_box[as1].iy(), amesh.m_level[lastlevel].m_box[as1].iz());
			MPI_Abort(MPI_COMM_WORLD, 307);
	}		
#endif
#endif
	hgfv.zero();
	sc = 0;
#if DIM > 2	
	for (Point_iterator p(0,2); p.end(); ++p)
	{
		if (distobc[p.i][p.j][p.k] > 0.0)
		{
			++sc;
			double aratio = distobc[p.i][p.j][p.k]/dissum;
			hgfv.u += aratio*nfv[p.i][p.j][p.k]->u;
			hgfv.v += aratio*nfv[p.i][p.j][p.k]->v;
			hgfv.w += aratio*nfv[p.i][p.j][p.k]->w;
			hgfv.roe += aratio*nfv[p.i][p.j][p.k]->roe;
			hgfv.p += aratio*nfv[p.i][p.j][p.k]->p;
			//nfv[p.i][p.j][p.k]->showdata(p.i,p.j,p.k,distobc[p.i][p.j][p.k], dissum);
		}
	}
	if (sc < 4)
	{
		Pointxyz & bcc = amesh.bc(lastlevel, ibbc);
		printf("HG cell near box (%d,%d,%d) center is (%f,%f,%f) hg point (%f,%f,%f) close box point (%d,%d,%d) infect positive %d negative %d not infect %d\n"
			"normal vector is (%f,%f,%f) origin point (%f,%f,%f) hg length is %f\n",
			amesh.m_level[lastlevel].m_box[ibbc].ix(),
			amesh.m_level[lastlevel].m_box[ibbc].iy(),
			amesh.m_level[lastlevel].m_box[ibbc].iz(),
			bcc[0], bcc[1], bcc[2], hgpt[0], hgpt[1], hgpt[2], closepti[0], closepti[1], closepti[2], 
			infnum_p, infnum_n, notinfnum,
			nmv[0], nmv[1], nmv[2], atp.xyz[0], atp.xyz[1], atp.xyz[2], hl);
		for (Point_iterator p(0,3); p.end(); ++p)
		{
			int aneib = amesh.BoxNeib(lastlevel, ibbc, p.i,p.j,p.k);
			double distoan = amesh.distobc(lastlevel, aneib, hgpt);
			printf("HG near eight boxes (%d,%d,%d) is (%d,%d,%d) distance to the box is %f\n", p.i, p.j, p.k, 
						amesh.m_level[lastlevel].m_box[aneib].ix(),
						amesh.m_level[lastlevel].m_box[aneib].iy(),
						amesh.m_level[lastlevel].m_box[aneib].iz(),
						distoan);
		}
		int pt1 = amesh.m_level[lastlevel].m_box[ibbc].pts[1][1][1];
		int pt0 = amesh.m_level[lastlevel].m_box[ibbc].pts[0][0][0];
		double dx = amesh.m_level[lastlevel].m_point[pt1][0] - amesh.m_level[lastlevel].m_point[pt0][0];
		double dy = amesh.m_level[lastlevel].m_point[pt1][1] - amesh.m_level[lastlevel].m_point[pt0][1];
		double dz = amesh.m_level[lastlevel].m_point[pt1][2] - amesh.m_level[lastlevel].m_point[pt0][2];
		printf("the close box length is (%f,%f,%f)\n", dx, dy, dz);
		MPI_Abort(MPI_COMM_WORLD, 221);
	}
#else
	for (Point_iterator_2d p(0,2); p.end(); ++p)
	{
		if (distobc[p.i][p.j][p.k] > 0.0)
		{
			++sc;
			double aratio = distobc[p.i][p.j][p.k]/dissum;
			hgfv.u += aratio*nfv[p.i][p.j][p.k]->u;
			hgfv.v += aratio*nfv[p.i][p.j][p.k]->v;
			hgfv.w += aratio*nfv[p.i][p.j][p.k]->w;
			hgfv.roe += aratio*nfv[p.i][p.j][p.k]->roe;
			hgfv.p += aratio*nfv[p.i][p.j][p.k]->p;
			//nfv[p.i][p.j][p.k]->showdata(p.i,p.j,p.k,distobc[p.i][p.j][p.k], dissum);
		}
	}
	if (sc < 3)
	{
		Pointxyz & bcc = amesh.bc(lastlevel, ibbc);
		printf("HG cell near box center is (%f,%f,%f) hg point (%f,%f,%f) close box point (%d,%d,%d) infect positive %d negative %d not infect %d\n"
			"normal vector is (%f,%f,%f) origin point (%f,%f,%f)\n", 
			bcc[0], bcc[1], bcc[2], hgpt[0], hgpt[1], hgpt[2], closepti[0], closepti[1], closepti[2], 
			infnum_p, infnum_n, notinfnum,
			nmv[0], nmv[1], nmv[2], atp.xyz[0], atp.xyz[1], atp.xyz[2]);
		MPI_Abort(MPI_COMM_WORLD, 221);
	}
#endif		
	Get_T_ideal_gas(hgfv);
	Get_E(hgfv);
	FINISH: return 0;
}

void Body::ComptWallForce(Mesh & amesh)
{
	double start_time = MPI_Wtime(); 
	lastlevel = amesh.MyCurNum()-1;
	if (moveflag)
	{
		FindPatchHangingNode(amesh);
		//GiveAFlag("Finish FindPatchHangingNode in ComptWallForce!!!", 5);
	}
	force.zero();
	moment.zero();
	pre_force.zero();
	tau_force.zero();
	pre_moment.zero();
	tau_moment.zero();
	int patstart = patch.ps();
	int patend = patch.pe();
	Pointxyz patv;
	Pointxyz vtdir;
	Pointxyz useless[2]; 
	double useless0[2];
	for (int i = patstart; i < patend; ++i)
	{
		int ci0 = patchbox[i-patstart];
		//if (amesh.m_level[lastlevel].m_box.isnormal(ci0))
		if (patch[i].node == node)
		{
			HGCell & myhg = patch[i].hgc;
			amesh.IntpHGgcellVar(myhg, lastlevel);
#ifdef PASSAGE_ANGLE			
			if (myhg.rotangle != 0.0)
			{
				Pointxyz flowvel(myhg.fv.u, myhg.fv.v, myhg.fv.w);
				flowvel.rotate_x(-myhg.rotangle);
				myhg.fv.u = flowvel[0];
				myhg.fv.v = flowvel[1];
				myhg.fv.w = flowvel[2];
			}
#endif			
			if (myhg.fv.hasnan(-1,-1,"ComptWallForce fv nan"))
			{
				printf("N%d Wall force intp patch %d (%f,%f,%f) belongs to node %d patchbox %d (%d,%d,%d)(%f,%f,%f) pdistobox %f hg close box is (%d,%d,%d)(%f,%f,%f)\n",node,i,
					patch[i].pc[0],patch[i].pc[1],patch[i].pc[2], patch[i].node,
					ci0, 
					amesh.m_level[lastlevel].m_box[ci0].ix(),
					amesh.m_level[lastlevel].m_box[ci0].iy(),
					amesh.m_level[lastlevel].m_box[ci0].iz(),
					amesh.m_level[lastlevel].m_geom[ci0].boxcenter[0],
					amesh.m_level[lastlevel].m_geom[ci0].boxcenter[1],
					amesh.m_level[lastlevel].m_geom[ci0].boxcenter[2],
					pdistobox[i-patstart],
					amesh.m_level[lastlevel].m_box[myhg.closecell].ix(),
					amesh.m_level[lastlevel].m_box[myhg.closecell].iy(),
					amesh.m_level[lastlevel].m_box[myhg.closecell].iz(), 
					amesh.m_level[lastlevel].m_geom[myhg.closecell].boxcenter[0],
					amesh.m_level[lastlevel].m_geom[myhg.closecell].boxcenter[1],
					amesh.m_level[lastlevel].m_geom[myhg.closecell].boxcenter[2]);
				for (Point_iterator p(0,2); p.end(); ++p)
				{
					int itp0 = myhg.intpcell[p.i][p.j][p.k];
					printf("N%d Wall force patch %d intp (%d,%d,%d) is box (%d,%d,%d)\n", node, itp0, p.i,p.j,p.k,
						amesh.m_level[lastlevel].m_box[itp0].ix(),
							amesh.m_level[lastlevel].m_box[itp0].iy(),
							amesh.m_level[lastlevel].m_box[itp0].iz());
				}
				MPI_Abort(MPI_COMM_WORLD, 451);
			}
			// printf("Patch %d interpolation box (%d,%d,%d,%d) interpolation var is (%f,%f,%f,%f,%f)\n", 
			// 	i, myhg.intpcell[0], myhg.intpcell[1], myhg.intpcell[2], myhg.intpcell[3],
			// 	myhg.fv.u, myhg.fv.v, myhg.fv.w, myhg.fv.T, myhg.fv.p);
			Surface_Point_Velocity(patch[i].pc, patv);
			Pointxyz du(myhg.fv.u-patv[0], 
					myhg.fv.v-patv[1], 
					myhg.fv.w-patv[2]);
			Pointxyz dxyz_tobc = patch[i].pc - bodycenter;
#ifndef TURBULENCE			
			du = (du - patch[i].nv*(du.dot(patch[i].nv)))/myhg.hgdis;			
			Pointxyz df0 = du*mu*patch[i].area;
#else
			Patchut(patch[i], patv, patch[i].nv, vtdir);
			FiltSinglePatchPressure(patch[i]);
			// IntpHgUt(amesh, lastlevel, i, patv, patch[i].ut);
			// VelocityDecom(du, patch[i].nv, useless[0], useless0[0], vtdir, useless[1], useless0[1]);
			double wall_roe = myhg.fv.p/(Rg*patch[i].wall_T);
			patch[i].Cp = 2.0*(myhg.fv.p - initvar.p);
			patch[i].Cf = wall_roe*pow(patch[i].ut, 2);
			Pointxyz df0 = vtdir*patch[i].Cf*patch[i].area;
			patch[i].Cf *= 2.0;
#endif			
			Pointxyz df1 = patch[i].nv*(-myhg.fv.p*patch[i].area);
			Pointxyz dm0 = dxyz_tobc.cross(df0);
			Pointxyz dm1 = dxyz_tobc.cross(df1);
			pre_force += df1;
			tau_force += df0;
			pre_moment += dm1;
			tau_moment += dm0;
		}
		else
		{
			patch[i].yplus = -2.0;
		}
	}
	MPI_Allreduce(MPI_IN_PLACE, &pre_force[0], 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, &tau_force[0], 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, &pre_moment[0], 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, &tau_moment[0], 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	pre_force = pre_force * p_ref * area_ref;
	pre_moment = pre_moment * p_ref * area_ref * L_ref;
#ifndef TURBULENCE	
	tau_force = tau_force * shear_ref * area_ref;
	tau_moment = tau_moment * shear_ref * area_ref * L_ref;
#else
	tau_force = tau_force * p_ref * area_ref;
	tau_moment = tau_moment * p_ref * area_ref * L_ref;
#endif		
	force = pre_force + tau_force;
	moment = pre_moment + tau_moment;
	double end_time = MPI_Wtime();
	total_force_time += (end_time - start_time);
}
// #ifdef TURBULENCE
// void Body::IntpHgUt(Mesh & amesh, const int & lastlevel, const int & patindex, Pointxyz & patv, double & patchut)
// {
// 	Pointxyz & nmv = patch[patindex].nv;
// 	double cellut[4], tangdis[4];
// 	Pointxyz vtdir[4];
// 	for (int i = 0; i < 4; ++i)
// 	{
// 		int i0 = patch[patindex].hgc.intpcell[i];
// 		Cellut(amesh.m_level[lastlevel].m_data[i0], 
// 			   amesh.m_level[lastlevel].m_box[i0].pair.signdis,
// 			   patv, nmv, cellut[i], vtdir[i]);
// 		tangdis[i] = vtdir[i].dot(patch[patindex].hgc.pt - amesh.m_level[lastlevel].m_geom[i0].boxcenter);
// 	}
// 	double dissum = 0.0;
// 	patchut = 0.0;
// 	for (int i = 0; i < 4; ++i)
// 	{
// 		tangdis[i] = 1.0/(tangdis[i]*tangdis[i]);
// 		if (isinf(tangdis[i])) 
// 		{
// 			patchut = cellut[i];
// 			goto FINISH;
// 		}
// 		dissum += tangdis[i];
// 	}
// 	for (int i = 0; i < 4; ++i)
// 	{
// 		patchut += tangdis[i]/dissum*cellut[i];
// 	}
// 	FINISH:;
// }
// #endif
void Body::ComptDisFromBoxtoBody_Init(Mesh & amesh, vector<Body> & abody)
{
	Assert(amesh.MyCurNum()==1, "The mesh must only have one level during initial attach!!!", 527);
	int bodysign[2] = {-1,1};
	int bs, be;
	amesh.m_level[0].m_box.GlobalOrder(bs, be);
	int bn = abody.size();
	for (int i = bs; i < be; ++i)
	{
		amesh.m_level[0].m_box[i].init_boxtowall();			
		for (int b0 = 0; b0 < bn; ++b0)
		{			
			int tppe = abody[b0].patch.size();
			for (int p0 = 0; p0 < tppe; ++p0)
			{
				abody[b0].PatchDistoBox(amesh.m_level[0].m_box[i].pair, amesh.m_level[0].m_geom[i].boxcenter, p0);
			}			
		}		
	}
	MPI_Barrier(share_comm);
	GiveAFlag("Finish computing the distances from the initial level to bodies!!!", 5);
	DistanceSign_boxtowall(amesh, abody, 0);
	GiveAFlag("Finish computing the distances with sign from the initial level to bodies!!!", 5);
}

void Body::ComptDisFromBoxtoBody_Newlevel(Mesh & amesh, DataArray<rftag> & mtag, DataArray<Boxson<int> > & mson, vector<Body> & abody)
{
	int level0 = amesh.MyCurNum()-1;
	Assert(level0 > 0, "The new level index must be positive when compute the distance to a new level!!!", 549);
	int momlevel = level0-1;
	int bs = mtag.ps();
	int be = mtag.pe();
	if (!level_twod_flag[level0])
	{
		for (int i = bs; i < be; ++i)
		{
			int tag0 = mtag[i].tag;
			if (tag0 > -1)
			{				
				for (Point_iterator p(0,2); p.end(); ++p)									
				{
					int s0 = mson[tag0].son[p.i][p.j][p.k];
					Assert(s0 > -1 && s0 < amesh.m_level[level0].m_box.realsize(), "Son box index not in range when computing distance!!!", 641);
					BoxtoWall & atw = amesh.m_level[level0].m_box[s0].pair;
					atw.distance = 999.0;
// #ifdef BODY_SPHERE
// 					atw.patch = 0;
// 					Pointxyz dxyz = amesh.m_level[level0].m_geom[s0].boxcenter - abody[atw.body].bc();
// 					atw.distance = dxyz.length();
// #else								
					abody[atw.body].FindClosePatchtoCell(atw, amesh.m_level[level0].m_geom[s0].boxcenter);
// #endif					
				}
			}
		}

	}
	else
	{
		for (int i = bs; i < be; ++i)
		{
			int tag0 = mtag[i].tag;
			if (tag0 > -1)
			{									
// 				for (Point_iterator_2d p(0,2); p.end(); ++p)
// 				{
// 					int s0 = mson[tag0].son[p.i][p.j][p.k];
// 					//printf("Box %d tag is %d son (%d,%d,%d) is %d\n", i, tag0, p.i,p.j,p.k,s0);
// 					Assert(s0 > -1 && s0 < amesh.m_level[level0].m_box.realsize(), "Son box index not in range when computing distance!!!", 641);
// 					BoxtoWall & atw = amesh.m_level[level0].m_box[s0].pair;
// 					Assert(atw.body > -1 && atw.body < bodynum, "The new box should give a init close body!!!", 594);
// #if DIM==2					
// 					abody[atw.body].FindClosePatchtoCell_CycleAllPatches(atw, amesh.m_level[level0].m_geom[s0].boxcenter);
// #elif DIM == 3	
// 					abody[atw.body].FindClosePatchtoCell(atw, amesh.m_level[level0].m_geom[s0].boxcenter);
// #endif										
// 					int atb0 = atw.body;
// 					for (int sd = 0; sd < my_search_num; ++sd)
// 					{
// 						//printf("Box %d tag is %d son (%d,%d,%d) is %d search body %d is %d\n", i, tag0, p.i,p.j,p.k,s0,sd,abody[atb0].search_near_body[sd]);
// 						if (abody[atb0].search_near_body[sd] > -1)
// 						{
// 							if (abody[atb0].search_near_body[sd] > bodynum - 1)
// 							{
// 								printf("Search near body %d but total bodynum is %d!!!\n", abody[atb0].search_near_body[sd], bodynum);
// 								MPI_Abort(MPI_COMM_WORLD, 619);
// 							}
// 							int bsd = abody[atb0].search_near_body[sd];
// #if DIM == 2							
// 							abody[bsd].FindClosePatchtoCell_CycleAllPatches(atw, amesh.m_level[level0].m_geom[s0].boxcenter);
// #elif DIM == 3
// 							BoxtoWall btw;
// 							btw.body = bsd;
// 							btw.patch = abody[bsd].patch.size()/2;
// 							abody[bsd].FindClosePatchtoCell(btw, amesh.m_level[level0].m_geom[s0].boxcenter);
// 							atw.switchtwo(btw);
// #endif													
// 						}
// 					}
// 				}
// 			}
// 		}
				for (Point_iterator_2d p(0,2); p.end(); ++p)
				{
					int s0 = mson[tag0].son[p.i][p.j][p.k];
					//printf("Box %d tag is %d son (%d,%d,%d) is %d\n", i, tag0, p.i,p.j,p.k,s0);
					Assert(s0 > -1 && s0 < amesh.m_level[level0].m_box.realsize(), "Son box index not in range when computing distance!!!", 641);
					BoxtoWall & atw = amesh.m_level[level0].m_box[s0].pair;
					atw.distance = 999.0;
					Assert(atw.body > -1 && atw.body < bodynum, "The new box should give a init close body!!!", 594);	
					abody[atw.body].FindClosePatchtoCell(atw, amesh.m_level[level0].m_geom[s0].boxcenter);
				}
			}			
		}
		MPI_Barrier(share_comm);
		int cycle0;
		CHECKNEIBCELL:;
		cycle0 = 0;
		for (int i = bs; i < be; ++i)
		{
			int tag0 = mtag[i].tag;
			if (tag0 > -1)
			{
				bool bodycheck = false;
				vector<BoxtoWall> nearbodypair;
				for (Point_iterator_2d p(0,2); p.end(); ++p)									
				{
					int s0 = mson[tag0].son[p.i][p.j][p.k];
#if DIM == 3			
					for (Point_iterator q(0,2); q.end(); ++q)
					{
						int an0 = amesh.m_level[level0].m_box[s0].neib[q.i+p.i][q.j+p.j][q.k+p.k];
#elif DIM == 2
					for (Point_iterator_2d q(0,2); q.end(); ++q)
					{
						int an0 = amesh.m_level[level0].m_box[s0].neib[q.i+p.i][q.j+p.j][1];
#endif						
						if (an0 > -1)
						{
							int an0_closebody = amesh.m_level[level0].m_box[an0].pair.body;
							Assert(an0_closebody > -1 && an0_closebody < abody.size(), "A cell closest body not in range!!!", 662);
							if (an0_closebody != amesh.m_level[level0].m_box[s0].pair.body)
							{
								bodycheck = true;
								nearbodypair.push_back(amesh.m_level[level0].m_box[an0].pair);								
							}
						}
					}
				}
				if (bodycheck)
				{
					int checknum = nearbodypair.size();
					for (Point_iterator_2d p(0,2); p.end(); ++p)									
					{
						int s0 = mson[tag0].son[p.i][p.j][p.k];
						for (int bk = 0; bk < checknum; ++bk)
						{
							BoxtoWall & an0_btw = nearbodypair[bk];
							an0_btw.distance = 999.0;
							abody[an0_btw.body].FindClosePatchtoCell(an0_btw, amesh.m_level[level0].m_geom[s0].boxcenter);
							if (an0_btw.distance < amesh.m_level[level0].m_box[s0].pair.distance)
							{
								amesh.m_level[level0].m_box[s0].pair = an0_btw;
								cycle0 = 1;
							}
						}
					}
				}
			}
 		}
 		MPI_Allreduce(MPI_IN_PLACE, &cycle0, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
 		if (cycle0 > 0) goto CHECKNEIBCELL;	
	}
	// Nodedata & ad0 = amesh.m_level[level0].nd;
	// amesh.DataExchange_Pairinfo(ad0, level0);
	// int recv0 = ad0.recvpair.size()
	// for (int i = 0; i < recv0; ++i)
	// {
	// 	int dest0 = ad0.recvdestcell[i];
	// 	Pairinfo & rp = ad0.recvpair[i];
	// 	if (rp.body != amesh.m_level[level0].m_box[dest0].pair.body)
	// 	{
	// 		BoxtoWall newbtw;
	// 		newbtw.body = rp.body;
	// 		newbtw.patch = rp.patch;
	// 		newbtw.distance = 999.0;
	// 		abody[newbtw.body].FindClosePatchtoCell(newbtw, amesh.m_level[level0].m_geom[dest0].boxcenter);
	// 		if (newbtw.distance < amesh.m_level[level0].m_box[dest0].pair.distance)
	// 		{
	// 			amesh.m_level[level0].m_box[dest0].pair = newbtw;
	// 		}
	// 	}
	// }
	MPI_Barrier(share_comm);
	GiveAFlag("Finish computing the distance from a new level to bodies!!!", 5);
	DistanceSign_boxtowall(amesh, abody, level0);
	GiveAFlag("Finish computing the distance with signs from a new level to bodies!!!", 5);
}

	void Body::DistanceSign_boxtowall(Mesh & amesh, vector<Body> & abody, const int & ilevel)
	{
		int bs, be;
		amesh.m_level[ilevel].m_box.GlobalOrder(bs, be);
		Pointxyz dxyz;
		for (int i = bs; i < be; ++i)
		{ 
			Box & aibc = amesh.m_level[ilevel].m_box[i];
#ifndef PASSAGE_ANGLE								
			dxyz = amesh.bc(ilevel, i) - abody[aibc.pair.body].patch[aibc.pair.patch].pc;
			PeriodicLength(dxyz);
#else	
			Pointxyz & bcxyz = amesh.m_level[ilevel].m_geom[i].boxcenter;
			Pointxyz & pcxyz = abody[aibc.pair.body].patch[aibc.pair.patch].pc;
			Pointxyz newbcxyz;
			PeriodicAnnulaLength(bcxyz, pcxyz, newbcxyz);
			dxyz = newbcxyz - pcxyz;
#endif			
			aibc.pair.signdis = abody[aibc.pair.body].patch[aibc.pair.patch].nv.dot(dxyz);
			aibc.pair.distance_to_body = aibc.pair.signdis;
			if (aibc.pair.signdis < 0.0)
			{
				amesh.m_level[ilevel].m_box[i].solid = true;
			}
			else
			{
				amesh.m_level[ilevel].m_box[i].solid = false;
			}
		}
		MPI_Barrier(share_comm);
	}

void Body::CollectPatchNeib()
{
	int ppe = patch.size();
	patneib.setnum_nocopy(ppe,0);
	if (srank == 0)
	{
		int ppt0 = allpoint.size();
		vector<vector<int> > pthasneib(ppt0);
		for (int i = 0; i < ppe; ++i)
		{
#ifdef QUA_ELEMENT			
			for (int p0 = 0; p0 < 4; ++p0)
#endif
#ifdef TRI_ELEMENT
			for (int p0 = 0; p0 < 3; ++p0)
#endif				
			{
				int pt0 = patch[i].corpt[p0];
				pthasneib[pt0].push_back(i);
			}
		}
		vector<vector<int> > pneib(ppe);
		for (int i = 0; i < ppt0; ++i)
		{
			int ptneibnum = pthasneib[i].size();
			if (ptneibnum > 0)
			{					
				for (int j = 0; j < ptneibnum; ++j)
				{
					int pat0 = pthasneib[i][j];
					for (int j0 = 0; j0 < ptneibnum; ++j0)
					{
						if (j0 != j)
						{
							pneib[pat0].push_back(pthasneib[i][j0]);
						}
					}
				}
			}
		}
		for (int i = 0; i < ppe; ++i)
		{
			Assert(pneib[i].size() > 0, "One patch must have non-zero neibs...!!!", 640);
			for (int p0 = 1; p0 < pneib[i].size(); ++p0)
			{
				int pat0 = pneib[i][p0];
				for (int ps0 = 0; ps0 < p0; ++ps0)
				{
					int pat1 = pneib[i][ps0];
					if (pat0 == pat1)
					{
						pneib[i][p0] = pneib[i].back();
						--p0;
						pneib[i].pop_back();
						break;
					}
				}
			}
		}
		for (int i = 0; i < ppe; ++i)
		{
			patneib[i].nbnum = pneib[i].size();
#if DIM == 2			
			Assert(patneib[i].nbnum == 2 || patneib[i].nbnum == 1, "A patch neib number error...", 660);
#elif DIM == 3			
			Assert(patneib[i].nbnum >= 2, "A patch neib number must have be larger than 1...", 660);
#endif				
			for (int p0 = 0; p0 < patneib[i].nbnum; ++p0)
			{
				patneib[i].nb[p0] = pneib[i][p0];
			}
#ifdef QUA_ELEMENT			
			if (patneib[i].nbnum > 15)
#endif
#ifdef TRI_ELEMENT				
			if (patneib[i].nbnum > 20)
#endif
			{
				for (int i0 = 0; i0 < patneib[i].nbnum; ++i0)
				{
					printf("Patch %d (%f,%f,%f) (%d,%d,%d,%d) neib %d is %d (%f,%f,%f)(%d,%d,%d,%d)\n", i, 
						patch[i][0], patch[i][1], patch[i][2], 
						patch[i].corpt[0], patch[i].corpt[1], patch[i].corpt[2], patch[i].corpt[3],
						i0, patneib[i].nb[i0], patch[patneib[i].nb[i0]][0], patch[patneib[i].nb[i0]][1], patch[patneib[i].nb[i0]][2],
						patch[patneib[i].nb[i0]].corpt[0], patch[patneib[i].nb[i0]].corpt[1], patch[patneib[i].nb[i0]].corpt[2], patch[patneib[i].nb[i0]].corpt[3]);
				}
				MPI_Abort(MPI_COMM_WORLD, 722);
			}
		}
	}
	MPI_Barrier(share_comm);
}

void Body::InfectWallBox_Init_Negative(Mesh & amesh, vector<Body> & abody)
{
	int dirnum[4] = {-1,-1,4,6};
	Point adp[6] = {Point(0,1,1),Point(2,1,1),Point(1,0,1),Point(1,2,1),Point(1,1,0),Point(1,1,2)};
	amesh.InitWallInfect();
	lastlevel = amesh.MyCurNum()-1;
	int bs, be;
	amesh.m_level[lastlevel].m_box.GlobalOrder(bs, be);
	// int bs = amesh.m_level[lastlevel].m_box.ps();
	// int be = amesh.m_level[lastlevel].m_box.pe();
	for (int i = bs; i < be; ++i)
	{
		if (amesh.m_level[lastlevel].m_box[i].type != Dmghost)
		{
			if (!amesh.m_level[lastlevel].m_box[i].solid)
			{
				bool infectflag = false;
#if DIM == 3				
				for (Point_iterator p(0,3); p.end(); ++p)
#elif DIM == 2
      			for (Point_iterator_2d p(0,3); p.end(); ++p)
#endif
				{
#if DIM == 3					
					int ab0 = amesh.m_level[lastlevel].m_box[i].neib[p.i][p.j][p.k];
#elif DIM == 2				
					int ab0 = amesh.m_level[lastlevel].m_box[i].neib[p.i][p.j][1];
#endif
					if (ab0 > -1)
					{
						if (amesh.m_level[lastlevel].m_box[ab0].solid && amesh.m_level[lastlevel].m_box[ab0].type != Dmghost)
						{
							amesh.InfectABox(i);
							infectflag = true;		
							break;
						}
					}
				}				
				if (infectflag)
				{				
#if DIM == 3				
					for (Point_iterator p(0,3); p.end(); ++p)
#elif DIM == 2
      				for (Point_iterator_2d p(0,3); p.end(); ++p)
#endif
      				{
#if DIM == 3					
						int ab1 = amesh.m_level[lastlevel].m_box[i].neib[p.i][p.j][p.k];
#elif DIM == 2				
						int ab1 = amesh.m_level[lastlevel].m_box[i].neib[p.i][p.j][1];
#endif
						if (ab1 > -1)
						{
							if (!amesh.m_level[lastlevel].m_box[ab1].solid && amesh.m_level[lastlevel].m_box[ab1].type != Dmghost)
							{
								amesh.InfectABox(ab1);
							}
						}
      				}
				}										
			}
		}
	}
	MPI_Barrier(share_comm);
	amesh.CountInfectedBox();
	for (int i = bs; i < be; ++i)
	{
		int f0 = amesh.infectbox[i];
		if (f0 > -1)
		{
			BoxtoWall & ibboxtopatch = amesh.m_level[lastlevel].m_box[i].pair;
			amesh.m_dis[f0].ci = i;
			amesh.m_dis[f0].tp = amesh.m_level[lastlevel].m_geom[i].boxcenter-
				abody[ibboxtopatch.body].patch[ibboxtopatch.patch].nv*ibboxtopatch.signdis;
			amesh.m_dis[f0].hgc.closecell = i;
			amesh.m_dis[f0].hgc.InitIntparray();
			amesh.m_dis[f0].fv = initvar;
			amesh.m_dis[f0].mynode = node;
		}
	}
	MPI_Barrier(share_comm);
}

void Body::FindImageCell(Mesh & amesh, vector<Body> & abody)
{
	int mps = amesh.m_dis.ps();
	int mpe = amesh.m_dis.pe();
	lastlevel = amesh.MyCurNum()-1;
	bool rightpos;
	Pointxyz dkeisa, dxyz, dkeisa_tocc, pnmv, inc_dir, pnmv_start, cbox_rot, ib_to_pat;
	double dhglength, dhgcell;
	int hgnode;
	int ibnotinnode_num = 0;
	int rev_dir;
	for (int i = mps; i < mpe; ++i)
	{
		amesh.m_dis[i].mynode = node;
#ifdef PASSAGE_ANGLE
		amesh.m_dis[i].hgc.rotangle = 0.0;
		amesh.m_dis[i].angle_ib_to_pat = 0.0;
#endif				
		//amesh.m_dis[i].hgc.InitIntparray();
		if (amesh.m_level[lastlevel].m_box.isnormal(amesh.m_dis[i].ci))
		{
			IBCell & myib = amesh.m_dis[i];
			HGCell & myhg = myib.hgc;
			BoxtoWall & ibboxtopatch = amesh.m_level[lastlevel].m_box[myib.ci].pair;	
			int bi0 = ibboxtopatch.body;
			int pi0 = ibboxtopatch.patch;
			myhg.closecell = myib.ci;
			Pointxyz & cbox_c = amesh.m_level[lastlevel].m_geom[myib.ci].boxcenter;
			cbox_rot = cbox_c;					
			if (ibboxtopatch.signdis < 0.0)
			{
				printf("IBCell should has a non-negative distance!!!\n", 1078);
				MPI_Abort(MPI_COMM_WORLD, 1078);				
			}
			else
			{
				int dir_rev_num = 0;				
				myhg.hgdis = 0.0;			
				rightpos = false;
				pnmv = abody[bi0].patch[pi0].nv;
				inc_dir = abody[bi0].patch[pi0].nv;
#ifdef PASSAGE_ANGLE
				double dag = ComptPointAngle_Rotate_X(abody[bi0].patch[pi0].pc, amesh.m_level[lastlevel].m_geom[myib.ci].boxcenter);
				if (dag > 0.5*PASSAGE_ANGLE) 
				{
					pnmv.rotate_x(PASSAGE_ANGLE); 
					inc_dir.rotate_x(PASSAGE_ANGLE); 
					amesh.m_dis[i].angle_ib_to_pat = PASSAGE_ANGLE;
					cbox_rot.rotate_x(-PASSAGE_ANGLE);
				}
				else if (dag < -0.5*PASSAGE_ANGLE) 
				{
					pnmv.rotate_x(-PASSAGE_ANGLE); 
					inc_dir.rotate_x(-PASSAGE_ANGLE); 
					amesh.m_dis[i].angle_ib_to_pat = -PASSAGE_ANGLE;
					cbox_rot.rotate_x(PASSAGE_ANGLE);
				}
#else 
				ib_to_pat = amesh.m_level[lastlevel].m_geom[myib.ci].boxcenter - abody[bi0].patch[pi0].pc;
				PeriodicLength(ib_to_pat);
				cbox_rot = abody[bi0].patch[pi0].pc + ib_to_pat;
#endif
				myhg.attach_pt = cbox_rot - abody[bi0].patch[pi0].nv*ibboxtopatch.signdis;				
				int dmtime = 0;
				pnmv_start = pnmv;
				while (!rightpos)
				{
					dkeisa[0] = pnmv.dot(amesh.m_level[lastlevel].m_geom[myhg.closecell].keisa[0]);
					dkeisa[1] = pnmv.dot(amesh.m_level[lastlevel].m_geom[myhg.closecell].keisa[1]);
					dkeisa[2] = pnmv.dot(amesh.m_level[lastlevel].m_geom[myhg.closecell].keisa[2]);
					dhgcell = max(abs(dkeisa[0]), abs(dkeisa[1]));
					dhgcell = max(dhgcell, abs(dkeisa[2]));
					dxyz = pnmv/dhgcell;
					myhg.hgdis += dxyz.length();
					myhg.pt = cbox_c + inc_dir*myhg.hgdis;	
					amesh.LocateHGCell(lastlevel, myhg, abody[bi0].patch[pi0].pc, hgnode, dkeisa_tocc, rev_dir);
					// printf("IB cell(%d,%d,%d) signdis %f hg close cell (%d,%d,%d) hgdis %f\n",
					// 	amesh.m_level[lastlevel].m_box[myib.ci].ix(),
					// 	amesh.m_level[lastlevel].m_box[myib.ci].iy(),
					// 	amesh.m_level[lastlevel].m_box[myib.ci].iz(),
					// 	amesh.m_level[lastlevel].m_box[myib.ci].pair.signdis,
					// 	amesh.m_level[lastlevel].m_box[myhg.closecell].ix(),
					// 	amesh.m_level[lastlevel].m_box[myhg.closecell].iy(),
					// 	amesh.m_level[lastlevel].m_box[myhg.closecell].iz(),
					// 	myhg.hgdis);
					if (hgnode == node)
					{
						if (dir_rev_num > 0)
						{
							printf("The direction for the IB cell has been reversed once!!!\n");
							MPI_Abort(MPI_COMM_WORLD, 1048);
						}
#ifdef DEBUG						
						printf("N%d IB cell (%d,%d,%d) hg cell has been switched to the other side close box (%d,%d,%d)\n",
							node,
							amesh.m_level[lastlevel].m_box[myib.ci].ix(),
							amesh.m_level[lastlevel].m_box[myib.ci].iy(),
							amesh.m_level[lastlevel].m_box[myib.ci].iz(),
							amesh.m_level[lastlevel].m_box[myhg.closecell].ix(),
							amesh.m_level[lastlevel].m_box[myhg.closecell].iy(),
							amesh.m_level[lastlevel].m_box[myhg.closecell].iz());
#endif						
#ifdef PASSAGE_ANGLE						
						if (rev_dir == 1) {pnmv.rotate_x(PASSAGE_ANGLE); myhg.rotangle = PASSAGE_ANGLE;}
						else if (rev_dir == -1) {pnmv.rotate_x(-PASSAGE_ANGLE); myhg.rotangle = -PASSAGE_ANGLE;}
						else
						{
							printf("When the hg cell has been switched to the other side, the rev_dir should be 1 or -1 %d!!!\n", rev_dir);
							MPI_Abort(MPI_COMM_WORLD, 1050);
						}
						++dir_rev_num;
#endif						
					}
					else if (hgnode != -1)
					{
						amesh.m_dis[i].mynode = hgnode;
						if (amesh.IsNormalcellButNotBlockpair(lastlevel, myib.ci))
						{
							printf("N%d IB cell (%d,%d,%d) will be associated with a ghost cell in other node, but it is not InTransferRange!!!\n",
								node,
								amesh.m_level[lastlevel].m_box[myib.ci].ix(),
								amesh.m_level[lastlevel].m_box[myib.ci].iy(),
								amesh.m_level[lastlevel].m_box[myib.ci].iz());
							MPI_Abort(MPI_COMM_WORLD, 1089);
						}
						++ibnotinnode_num;
#ifdef DEBUG						
						printf("N%d IB cell (%d,%d,%d) close box %d patch %d signdis %f hg node (%f,%f,%f) hgdis %f close cell (%d,%d,%d) is not in this node and will be sent to node %d!!!\n",
							node, 
							amesh.m_level[lastlevel].m_box[myib.ci].ix(),
							amesh.m_level[lastlevel].m_box[myib.ci].iy(),
							amesh.m_level[lastlevel].m_box[myib.ci].iz(),
							amesh.m_level[lastlevel].m_box[myib.ci].pair.body,
							amesh.m_level[lastlevel].m_box[myib.ci].pair.patch,
							amesh.m_level[lastlevel].m_box[myib.ci].pair.signdis,
							myhg.pt[0], myhg.pt[1], myhg.pt[2], myhg.hgdis,
							amesh.m_level[lastlevel].m_box[myhg.closecell].ix(),
							amesh.m_level[lastlevel].m_box[myhg.closecell].iy(),
							amesh.m_level[lastlevel].m_box[myhg.closecell].iz(),
							hgnode);
#endif						
						goto NEXTCYCLE_NOCHECK;
						// MPI_Abort(MPI_COMM_WORLD, 1018);
					}
					if (amesh.m_level[lastlevel].m_box[myhg.closecell].type == Dmghost)
					{
#if DIM == 2
						printf("Hg cell is a Dmghost!!!\n");
						MPI_Abort(MPI_COMM_WORLD, 1026);
#endif												
						amesh.m_dis[i].closetowall = true;
						//amesh.m_dis[i].wallbox = myhg.closecell;
						amesh.m_dis[i].wallbox = amesh.m_dis[i].ci;
						myhg.closecell = amesh.m_dis[i].ci;
						Pointxyz dmnmv; 
						amesh.GetCellDomainNormalVector(lastlevel, amesh.m_dis[i].wallbox, dmnmv);
						// printf("Node%dIB cell (%d,%d,%d) change its normal direction (%f,%f,%f) to (%f,%f,%f)\n",
						// 	node,
						// 	amesh.m_level[lastlevel].m_box[myib.ci].ix(),
						// 	amesh.m_level[lastlevel].m_box[myib.ci].iy(),
						// 	amesh.m_level[lastlevel].m_box[myib.ci].iz(),
						// 	inc_dir[0], inc_dir[1], inc_dir[2],
						// 	dmnmv[0],dmnmv[1],dmnmv[2]);
						inc_dir = inc_dir + dmnmv;
						inc_dir.normalize();
						pnmv = pnmv_start;
						rightpos = false;
						myhg.hgdis = 0.0;
						++dmtime;
						if (dmtime > 10)
						{
							printf("N%d An IB cell (%d,%d,%d) close to domain ghost fail to locate its hg cell!!!\n", node,
								amesh.m_level[lastlevel].m_box[myib.ci].ix(),
								amesh.m_level[lastlevel].m_box[myib.ci].iy(),
								amesh.m_level[lastlevel].m_box[myib.ci].iz());
							MPI_Abort(MPI_COMM_WORLD, 1114);
						}
						goto NEXTCYCLE;
					}
#ifdef WALLINTPFLAG_NOTSOLID					
					amesh.IntpcellNotSolid(rightpos, lastlevel, myhg);
#endif					
#ifdef WALLINTPFLAG_NOTINFECTED							
					amesh.IntpcellNotSolid_NotInfect(rightpos, lastlevel, myhg);
#endif 
					NEXTCYCLE:;
				}
				amesh.CheckHgIntpcells(lastlevel, myhg);
				amesh.ComptHgIntpCoef(myhg, dkeisa_tocc, lastlevel);
			}
			NEXTCYCLE_NOCHECK:;
		}
	}
	MPI_Barrier(share_comm);
#ifdef DEBUG	
	if (ibnotinnode_num > 0) printf("N%dR%d has %d ib cells not in the node!!!\n", node, srank, ibnotinnode_num);
#endif	
}

void Body::FindImageCell_GhostIBCell(Mesh & amesh, vector<Body> & abody)
{
	lastlevel = amesh.MyCurNum()-1;
	bool rightpos;
	Pointxyz dkeisa, dxyz, dkeisa_tocc, cbox_rot;
	double dhglength, dhgcell;
	int hgnode, rev_dir;
	for (int n0 = 0; n0 < nodenum; ++n0)
	{
		int revps = amesh.rev_mdis.global_cell_extrac[n0].ps();
		int revpe = amesh.rev_mdis.global_cell_extrac[n0].pe();
		for (int i = revps; i < revpe; ++i)
		{
			int aaa = amesh.rev_mdis.global_cell_extrac[n0][i].localcell;
			int ifc = amesh.infectbox[aaa];
			if (ifc > -1)
			{
				IBCell & myib = amesh.m_dis[ifc];
				HGCell & myhg = myib.hgc;
#ifdef PASSAGE_ANGLE
				myhg.rotangle = 0.0;
				myib.angle_ib_to_pat = 0.0;
#endif								
				BoxtoWall & ibboxtopatch = amesh.m_level[lastlevel].m_box[myib.ci].pair;	
				int bi0 = ibboxtopatch.body;
				int pi0 = ibboxtopatch.patch;
				Pointxyz & cbox_c = amesh.m_level[lastlevel].m_geom[myib.ci].boxcenter;
				cbox_rot = cbox_c;
				myhg.closecell = myib.ci;					
				if (ibboxtopatch.signdis < 0.0)
				{
					printf("The ghost ib cell signdis is %f!!! Can not be negative!!!\n", amesh.m_level[lastlevel].m_box[aaa].pair.signdis);
					MPI_Abort(MPI_COMM_WORLD, 1127);				
				}
				else
				{
					int dir_rev_num = 0;
					myhg.hgdis = 0.0;			
					rightpos = false;
					Pointxyz pnmv = abody[bi0].patch[pi0].nv;
					Pointxyz inc_dir = abody[bi0].patch[pi0].nv;
#ifdef PASSAGE_ANGLE
					// double dag = ComptPointAngle_Rotate_X(abody[bi0].patch[pi0].pc, amesh.m_level[lastlevel].m_geom[myib.ci].boxcenter);
					// if (dag > 0.5*PASSAGE_ANGLE) {pnmv.rotate_x(PASSAGE_ANGLE); inc_dir.rotate_x(PASSAGE_ANGLE); myib.angle_ib_to_pat = PASSAGE_ANGLE;}
					// else if (dag < -0.5*PASSAGE_ANGLE) {pnmv.rotate_x(-PASSAGE_ANGLE); inc_dir.rotate_x(-PASSAGE_ANGLE); myib.angle_ib_to_pat = -PASSAGE_ANGLE;}	
					double dag = ComptPointAngle_Rotate_X(abody[bi0].patch[pi0].pc, amesh.m_level[lastlevel].m_geom[myib.ci].boxcenter);
					if (dag > 0.5*PASSAGE_ANGLE) 
					{
						pnmv.rotate_x(PASSAGE_ANGLE); 
						inc_dir.rotate_x(PASSAGE_ANGLE); 
						amesh.m_dis[i].angle_ib_to_pat = PASSAGE_ANGLE;
						cbox_rot.rotate_x(-PASSAGE_ANGLE);
					}
					else if (dag < -0.5*PASSAGE_ANGLE) 
					{
						pnmv.rotate_x(-PASSAGE_ANGLE); 
						inc_dir.rotate_x(-PASSAGE_ANGLE); 
						amesh.m_dis[i].angle_ib_to_pat = -PASSAGE_ANGLE;
						cbox_rot.rotate_x(PASSAGE_ANGLE);
					}
#endif
					myhg.attach_pt = cbox_rot - abody[bi0].patch[pi0].nv*ibboxtopatch.signdis;					
					int dmtime = 0;
					while (!rightpos)
					{
						dkeisa[0] = pnmv.dot(amesh.m_level[lastlevel].m_geom[myhg.closecell].keisa[0]);
						dkeisa[1] = pnmv.dot(amesh.m_level[lastlevel].m_geom[myhg.closecell].keisa[1]);
						dkeisa[2] = pnmv.dot(amesh.m_level[lastlevel].m_geom[myhg.closecell].keisa[2]);
						dhgcell = max(abs(dkeisa[0]), abs(dkeisa[1]));
						dhgcell = max(dhgcell, abs(dkeisa[2]));
						dxyz = pnmv/dhgcell;
						myhg.hgdis += dxyz.length();
						myhg.pt = cbox_c + inc_dir*myhg.hgdis;	
						amesh.LocateHGCell(lastlevel, myhg, abody[bi0].patch[pi0].pc, hgnode, dkeisa_tocc, rev_dir);
						if (amesh.m_level[lastlevel].m_box[myhg.closecell].type == Dmghost)
						{
							printf("Hg cell is a Dmghost!!!\n");
							MPI_Abort(MPI_COMM_WORLD, 1026);											
						}
						if (hgnode == node)
						{
							if (dir_rev_num > 0)
							{
								printf("The direction for the IB cell has been reversed once in FindImageCell_GhostIBCell!!!\n");
								MPI_Abort(MPI_COMM_WORLD, 1048);
							}
#ifdef DEBUG							
							printf("N%d IB cell (%d,%d,%d) hg cell has been switched to the other side close box (%d,%d,%d)\n",
								node,
								amesh.m_level[lastlevel].m_box[myib.ci].ix(),
								amesh.m_level[lastlevel].m_box[myib.ci].iy(),
								amesh.m_level[lastlevel].m_box[myib.ci].iz(),
								amesh.m_level[lastlevel].m_box[myhg.closecell].ix(),
								amesh.m_level[lastlevel].m_box[myhg.closecell].iy(),
								amesh.m_level[lastlevel].m_box[myhg.closecell].iz());
#endif							
#ifdef PASSAGE_ANGLE						
							if (rev_dir == 1) {pnmv.rotate_x(PASSAGE_ANGLE); myhg.rotangle = PASSAGE_ANGLE;}
							else if (rev_dir == -1) {pnmv.rotate_x(-PASSAGE_ANGLE); myhg.rotangle = -PASSAGE_ANGLE;}
							else
							{
								printf("When the hg cell has been switched to the other side, the rev_dir should be 1 or -1 %d!!!\n", rev_dir);
								MPI_Abort(MPI_COMM_WORLD, 1050);
							}
							++dir_rev_num;
#endif							
						}
						if (hgnode == -1)
						{
#ifdef WALLINTPFLAG_NOTSOLID
							amesh.IntpcellNotSolid(rightpos, lastlevel, myhg);
#endif
#ifdef WALLINTPFLAG_NOTINFECTED							
							amesh.IntpcellNotSolid_NotInfect(rightpos, lastlevel, myhg);
#endif
						}				
					}
					if (hgnode != -1)
					{
						printf("N%d Ghost IB cell hg node (%f,%f,%f) hgdis %f close cell (%d,%d,%d) is not in this node!!! IB cell is (%d,%d,%d) close body %d patch %d signdis %f new hg node is %d\n",
								node, myhg.pt[0], myhg.pt[1], myhg.pt[2], myhg.hgdis,
								amesh.m_level[lastlevel].m_box[myhg.closecell].ix(),
								amesh.m_level[lastlevel].m_box[myhg.closecell].iy(),
								amesh.m_level[lastlevel].m_box[myhg.closecell].iz(),
								amesh.m_level[lastlevel].m_box[myib.ci].ix(),
								amesh.m_level[lastlevel].m_box[myib.ci].iy(),
								amesh.m_level[lastlevel].m_box[myib.ci].iz(),
								amesh.m_level[lastlevel].m_box[myib.ci].pair.body,
								amesh.m_level[lastlevel].m_box[myib.ci].pair.patch,
								amesh.m_level[lastlevel].m_box[myib.ci].pair.signdis,
								hgnode);
						MPI_Abort(MPI_COMM_WORLD, 1018);
					}
					amesh.CheckHgIntpcells(lastlevel, myhg);
					amesh.ComptHgIntpCoef(myhg, dkeisa_tocc, lastlevel);
#ifdef DEBUG					
					printf("N%d constructed ghost IB cell (%d,%d,%d) from node %d!!!\n", node,
						amesh.m_level[lastlevel].m_box[myib.ci].ix(),
						amesh.m_level[lastlevel].m_box[myib.ci].iy(),
						amesh.m_level[lastlevel].m_box[myib.ci].iz(),
						n0);
#endif					
				}
			}
		}
	}
	MPI_Barrier(share_comm);
}

void Body::IBConditions_Image(Mesh & amesh, vector<Body> & abody)
{
#ifdef SHOWTIME
	double start_time = MPI_Wtime();
#endif
	dt00 += dt;
	if (moveflag && ts%1 == 0)
	{
		MovingBody_Image(amesh, abody, dt00);
		dt00 = 0.0;
	}
#ifdef SHOWTIME
	double move_time = MPI_Wtime();
#endif		
	int mps = amesh.m_dis.ps();
	int mpe = amesh.m_dis.pe();
	lastlevel = amesh.MyCurNum()-1;
	Pointxyz patv;
	Pointxyz patnv;
	for (int i = mps; i < mpe; ++i)
	{			
		if ((amesh.m_level[lastlevel].m_box.isnormal(amesh.m_dis[i].ci) && amesh.m_dis[i].mynode == node) ||
			(amesh.m_level[lastlevel].m_box.isghost(amesh.m_dis[i].ci) && amesh.m_dis[i].mynode != node))				
		{	
			int ci0 = amesh.m_dis[i].ci;
			BoxtoWall & ibboxtopatch = amesh.m_level[lastlevel].m_box[ci0].pair;
			HGCell & myhg = amesh.m_dis[i].hgc;
			amesh.IntpHGgcellVar(myhg, lastlevel);
			abody[ibboxtopatch.body].Surface_Point_Velocity(myhg.attach_pt, patv);
			patnv = abody[ibboxtopatch.body].patch[ibboxtopatch.patch].nv;
#ifdef PASSAGE_ANGLE			
			if (myhg.rotangle != 0.0)
			{
				Pointxyz ofv(myhg.fv.u, myhg.fv.v, myhg.fv.w);
				ofv.rotate_x(-myhg.rotangle);
				myhg.fv.u = ofv[0];
				myhg.fv.v = ofv[1];
				myhg.fv.w = ofv[2];
			}
			if (amesh.m_dis[i].angle_ib_to_pat != 0.0)
			{
				patv.rotate_x(amesh.m_dis[i].angle_ib_to_pat);
				patnv.rotate_x(amesh.m_dis[i].angle_ib_to_pat);
			}
#endif			
			amesh.m_dis[i].patv = patv;			
			if (ibboxtopatch.signdis < 0.0)
			{
#ifndef TURBULENCE							
				amesh.m_dis[i].fv.u = 2.0*patv[0]-myhg.fv.u;
				amesh.m_dis[i].fv.v = 2.0*patv[1]-myhg.fv.v;
				amesh.m_dis[i].fv.w = 2.0*patv[2]-myhg.fv.w;
				amesh.m_dis[i].fv.p = myhg.fv.p;
				amesh.m_dis[i].fv.T = myhg.fv.T;
				Get_roe_ideal_gas(amesh.m_dis[i].fv);
				Get_E(amesh.m_dis[i].fv);
#else
				amesh.m_dis[i].fv = initvar;
				amesh.m_dis[i].fv.u = patv[0];
				amesh.m_dis[i].fv.v = patv[1];
				amesh.m_dis[i].fv.w = patv[2];
				// amesh.m_dis[i].fv.p = myhg.fv.p;
				// amesh.m_dis[i].fv.T = myhg.fv.T;
				// Get_roe_ideal_gas(amesh.m_dis[i].fv);
				// Get_E(amesh.m_dis[i].fv);
				// amesh.m_dis[i].fv.viseddy = initvar.viseddy;
#endif												
			}
			else
			{
#ifndef TURBULENCE				
				LinearDist(myhg.fv, 
					amesh.m_dis[i], 
					amesh.m_level[lastlevel].m_data[amesh.m_dis[i].ci], 
					patv,
					patnv,
					ibboxtopatch.signdis);
#else
				ImageViseddy(amesh.m_dis[i], patv, patnv, ibboxtopatch);
#endif								
			}

			if (myhg.fv.T < 0.0)
			{	
				printf("Hg cell T is negative!!!\n");			
				// printf("My interpolation coefficients are %f,%f,%f,%f\n",
				// 	myhg.intpcoef[0],
				// 	myhg.intpcoef[1],
				// 	myhg.intpcoef[2],
				// 	myhg.intpcoef[3]);
				MPI_Abort(MPI_COMM_WORLD, 1050);
			}
			if (amesh.m_dis[i].fv.hasnan(lastlevel,amesh.m_dis[i].ci,"ib solid cell"))
			{
				printf("IB Hg cell has nan!!! Patch velocity is (%f,%f,%f)!\n", patv[0], patv[1], patv[2]);
				printf("The cell is (%d,%d,%d) ghost flag is %d, hgnode is %d!!!\n", 
					amesh.m_level[lastlevel].m_box[ci0].ix(),
					amesh.m_level[lastlevel].m_box[ci0].iy(),
					amesh.m_level[lastlevel].m_box[ci0].iz(),
					amesh.m_level[lastlevel].m_box.isghost(amesh.m_dis[i].ci),
					amesh.m_dis[i].mynode);
				printf("The patch close body is %d patch is %d signdis is %f!!!\n",
					amesh.m_level[lastlevel].m_box[ci0].pair.body,
					amesh.m_level[lastlevel].m_box[ci0].pair.patch,
					amesh.m_level[lastlevel].m_box[ci0].pair.signdis);
				printf("Attach point is (%f,%f,%f)!!!\n", myhg.attach_pt[0],myhg.attach_pt[1],myhg.attach_pt[2]);
				int ci0 = amesh.m_dis[i].ci;
				for (Point_iterator p(0,2); p.end(); ++p)
				{
					int an1 = myhg.intpcell[p.i][p.j][p.k];
					printf("IB (%d,%d,%d) hg interpolation (%d,%d,%d) is %d (%d,%d,%d) signdis %f coefficient is %f!!!\n"
						"u%f v%f w%f p%f T%f roe%f..\n",
						amesh.m_level[lastlevel].m_box[ci0].ix(),
						amesh.m_level[lastlevel].m_box[ci0].iy(),
						amesh.m_level[lastlevel].m_box[ci0].iz(),
						p.i,p.j,p.k,an1,
						amesh.m_level[lastlevel].m_box[an1].ix(),
						amesh.m_level[lastlevel].m_box[an1].iy(),
						amesh.m_level[lastlevel].m_box[an1].iz(),
						amesh.m_level[lastlevel].m_box[an1].pair.signdis,
						myhg.intpcoef[p.i][p.j][p.k],
						amesh.m_level[lastlevel].m_data[an1].u,
						amesh.m_level[lastlevel].m_data[an1].v,
						amesh.m_level[lastlevel].m_data[an1].w,
						amesh.m_level[lastlevel].m_data[an1].p,
						amesh.m_level[lastlevel].m_data[an1].T,
						amesh.m_level[lastlevel].m_data[an1].roe);
				}
#ifndef LOCA_TIME_STEPPING
				MPI_Abort(MPI_COMM_WORLD, 1076);
#else 	
				amesh.m_dis[i].fv = amesh.m_level[lastlevel].m_data[ci0];
#endif							
			}
			Assert(!amesh.m_dis[i].fv.hasnan(lastlevel,amesh.m_dis[i].ci,"ib solid cell"), "IB image cell has nan!!!", 987);
		}
	}
	MPI_Barrier(share_comm);
	for (int i = mps; i < mpe; ++i)
	{
		IBCell & myib = amesh.m_dis[i];
		if ((amesh.m_level[lastlevel].m_box.isnormal(amesh.m_dis[i].ci) && amesh.m_dis[i].mynode == node) ||
			(amesh.m_level[lastlevel].m_box.isghost(amesh.m_dis[i].ci) && amesh.m_dis[i].mynode != node))
		{
			amesh.m_level[lastlevel].m_data[myib.ci] = amesh.m_dis[i].fv;			
		}				
	}
	MPI_Barrier(share_comm);
	amesh.DataExchange(amesh.rev_mdis, lastlevel);
	amesh.DataExchange(amesh.nd_mdis, lastlevel);
	GiveAFlag("Finish exchange wall data...", 5);
#ifdef SHOWTIME	
	double end_time = MPI_Wtime();
	step_ib_time = end_time - start_time;
	step_move_time = move_time - start_time;
	total_ib_time += step_ib_time;	
#endif
	}
#ifndef LOCALBODY
	void Body::FindPatchHangingNode(Mesh & amesh)
	{
		lastlevel = amesh.MyCurNum()-1;
		int pats = patch.ps();
		int pate = patch.pe();
		Pointxyz dxyz, newbcxyz, dkeisa, hgdxyz, dkeisa_tocc;
		double dhgcell;
		int hgnode, rev_dir;
		PatchamongNodes pack0;
		Pointxyz pnmv;
		for (int i = pats; i < pate; ++i)
		{
			int i0 = i - pats;
			int ci0 = patchbox[i0];
			Assert(amesh.m_level[lastlevel].m_box.isnorg(ci0),"Patch close box error!!!", 1383);
			//if (amesh.m_level[lastlevel].m_box.isnormal(ci0))
#ifdef PASSAGE_ANGLE			
			patch[i].hgc.rotangle = 0.0;
#endif			
			if (node == patch[i].node)
			{
				int dir_rev_num = 0;
				HGCell & pathg = patch[i].hgc;				
				pathg.hgdis = 0.0;
				pathg.attach_pt = patch[i].pc;
				pathg.closecell = ci0;
				bool rightpos = false;
				pnmv = patch[i].nv;
#ifdef PASSAGE_ANGLE
				double dag = ComptPointAngle_Rotate_X(patch[i].pc, amesh.m_level[lastlevel].m_geom[pathg.closecell].boxcenter);
				if (dag > 0.5*PASSAGE_ANGLE) {pnmv.rotate_x(PASSAGE_ANGLE); pathg.rotangle = PASSAGE_ANGLE;}
				else if (dag < -0.5*PASSAGE_ANGLE) {pnmv.rotate_x(-PASSAGE_ANGLE); pathg.rotangle = -PASSAGE_ANGLE;}
#endif								
				while (!rightpos)
				{
					dkeisa[0] = pnmv.dot(amesh.m_level[lastlevel].m_geom[pathg.closecell].keisa[0]);
					dkeisa[1] = pnmv.dot(amesh.m_level[lastlevel].m_geom[pathg.closecell].keisa[1]);
					dkeisa[2] = pnmv.dot(amesh.m_level[lastlevel].m_geom[pathg.closecell].keisa[2]);
					dhgcell = max(abs(dkeisa[0]), abs(dkeisa[1]));
					dhgcell = max(dhgcell, abs(dkeisa[2]));
					hgdxyz = pnmv/dhgcell;
					pathg.hgdis += hgdxyz.length();
					pathg.pt = patch[i].pc + patch[i].nv*pathg.hgdis;
					amesh.LocatePointinComptDomain_AllDomain(pathg.pt, lastlevel, pathg.closecell, dkeisa_tocc);
					amesh.ReversePeriodicSide(lastlevel, pathg.closecell, rev_dir);
					if (rev_dir != 0)
					{
						if (dir_rev_num > 0)
						{
							printf("The direction for the IB cell has been reversed once in FindImageCell_GhostIBCell!!!\n");
							MPI_Abort(MPI_COMM_WORLD, 1048);
						}
						++dir_rev_num;

						//printf("N%d Patch%d will be switched to the other side!!!\n", node, i);
#ifdef PASSAGE_ANGLE						
						if (rev_dir == 1) {pnmv.rotate_x(PASSAGE_ANGLE); pathg.rotangle += PASSAGE_ANGLE;}
						else if (rev_dir == -1) {pnmv.rotate_x(-PASSAGE_ANGLE); pathg.rotangle += -PASSAGE_ANGLE;}
#endif						
					}									
					amesh.FindHGCellIntpCell(pathg, lastlevel, dkeisa_tocc, hgnode);
					if (hgnode != -1)
					{
						pack0.newsendmeg(i, hgnode);
						printf("N%d Patch%d will be sent to another node %d because the intp cell is out of domain!!!\n", node, i, hgnode);
						goto NEXTCYCLE;
					}
#ifdef WALLINTPFLAG_NOTSOLID
					amesh.IntpcellNotSolid(rightpos, lastlevel, pathg);
#endif					
#ifdef WALLINTPFLAG_NOTINFECTED							
					amesh.IntpcellNotSolid_NotInfect(rightpos, lastlevel, pathg);
#endif					
					if (amesh.m_level[lastlevel].m_box[pathg.closecell].type == Dmghost)
					{
						printf("Patch hg cell close box is a Dmghost!!!\n");
						MPI_Abort(MPI_COMM_WORLD, 1280);
					}	
				}
				if (amesh.m_level[lastlevel].m_box.isghost(pathg.closecell))
				{
					//Point neibdir;
					hgnode = amesh.FindNeibBlock(amesh.m_level[lastlevel].m_box[pathg.closecell], lastlevel);
					//hgnode = m_dm.neibblocks[neibdir[0]][neibdir[1]][neibdir[2]];
					if (hgnode == node)
					{
						printf("N%d P%d The patch (%f,%f,%f) hgdis %f close box (%d,%d,%d)(%f,%f,%f) hg cell close ghost cell (%d,%d,%d)(%f,%f,%f) type %d ptype %d signdis %f can not be in the same node!!!"
							"May be the patch is very close to a normal ghost, i.e. the boundary of the mesh layer!!!\n", 
							hgnode, i,
							patch[i].pc[0], patch[i].pc[1], patch[i].pc[2],
							pathg.hgdis,
							amesh.m_level[lastlevel].m_box[patchbox[i0]].ix(),
							amesh.m_level[lastlevel].m_box[patchbox[i0]].iy(),
							amesh.m_level[lastlevel].m_box[patchbox[i0]].iz(),
							amesh.m_level[lastlevel].m_geom[patchbox[i0]].boxcenter[0],
							amesh.m_level[lastlevel].m_geom[patchbox[i0]].boxcenter[1],
							amesh.m_level[lastlevel].m_geom[patchbox[i0]].boxcenter[2],
							amesh.m_level[lastlevel].m_box[pathg.closecell].ix(),
							amesh.m_level[lastlevel].m_box[pathg.closecell].iy(),
							amesh.m_level[lastlevel].m_box[pathg.closecell].iz(),
							amesh.m_level[lastlevel].m_geom[pathg.closecell].boxcenter[0],
							amesh.m_level[lastlevel].m_geom[pathg.closecell].boxcenter[1],
							amesh.m_level[lastlevel].m_geom[pathg.closecell].boxcenter[2],
							amesh.m_level[lastlevel].m_box[pathg.closecell].type,
							amesh.m_level[lastlevel].m_box[pathg.closecell].ptype,
							amesh.m_level[lastlevel].m_box[pathg.closecell].pair.signdis);
						MPI_Abort(MPI_COMM_WORLD, 1433);
					}
					pack0.newsendmeg(i, hgnode);
#ifdef DEBUG					
					printf("N%d Patch%d (%f,%f,%f) will be sent to another node %d patchbox (%d,%d,%d) the hg close cell (%d,%d,%d) hgdis %f is a ghost!!!\n", node, i, 
							patch[i].pc[0], patch[i].pc[1], patch[i].pc[2],
							hgnode,
							amesh.m_level[lastlevel].m_box[patchbox[i0]].ix(),
							amesh.m_level[lastlevel].m_box[patchbox[i0]].iy(),
							amesh.m_level[lastlevel].m_box[patchbox[i0]].iz(),
							amesh.m_level[lastlevel].m_box[pathg.closecell].ix(),
							amesh.m_level[lastlevel].m_box[pathg.closecell].iy(),
							amesh.m_level[lastlevel].m_box[pathg.closecell].iz(),
							pathg.hgdis);
					printf("N%dP%d nv is (%f,%f,%f) pnmv is (%f,%f,%f)!!!\n", node, i,
						patch[i].nv[0], patch[i].nv[1], patch[i].nv[2],
						pnmv[0], pnmv[1], pnmv[2]);
#endif					
					goto NEXTCYCLE;

				}
				amesh.CheckHgIntpcells(lastlevel, pathg);
				amesh.ComptHgIntpCoef(pathg, dkeisa_tocc, lastlevel);
				NEXTCYCLE:;					
			}
		}
		MPI_Barrier(share_comm);
		//GiveAFlag("Finish finding local patch hg cell and intpcell!!!", 5);
		int cycle_num = 0;
		int patnode;
		SENDPATCHTONODE:;
		pack0.dataexchange();
		int recvpatnum = pack0.recv_pats.size();
		if (pack0.total_exc_num > 0)
		{
			vector<BoxlocChange> mdf_pats = pack0.recv_pats;
			pack0.initrecvarray();
			pack0.initsendarray();
			for (int i = 0; i < recvpatnum; ++i)
			{
				int pat0 = mdf_pats[i].oldloc;
				Assert(pat0 >= patch.ps() && pat0 < patch.pe(), "The patch is not in this node!!!", 484);
				int newnode = mdf_pats[i].newloc;
				patch[pat0].node = newnode;
#ifdef PASSAGE_ANGLE				
				patch[pat0].hgc.rotangle = 0.0;
#endif				
				if (newnode == node)
				{
					int dir_rev_num = 0;
					int i0 = pat0 - patch.ps();
					CELLSTART:;
					amesh.LocatePointinComptDomain_AllDomain(patch[pat0].pc, lastlevel, patchbox[i0], dkeisa_tocc);
					//%%%%%%%%%%%%%%%%%%%%%%%%%%%%
					HGCell & pathg = patch[pat0].hgc;				
					pathg.hgdis = 0.0;
					pathg.attach_pt = patch[pat0].pc;
					pathg.closecell = patchbox[i0];
					pnmv = patch[pat0].nv;
#ifdef PASSAGE_ANGLE
					double dag = ComptPointAngle_Rotate_X(patch[pat0].pc, amesh.m_level[lastlevel].m_geom[pathg.closecell].boxcenter);
#ifdef DEBUG					
					printf("N%d recv pat %d dag with patchbox (%d,%d,%d) is %f pnmv is (%f,%f,%f)\n", node, pat0, 
						amesh.m_level[lastlevel].m_box[patchbox[i0]].ix(),
						amesh.m_level[lastlevel].m_box[patchbox[i0]].iy(),
						amesh.m_level[lastlevel].m_box[patchbox[i0]].iz(),
						dag, pnmv[0], pnmv[1], pnmv[2]);
#endif					
					if (dag > 0.5*PASSAGE_ANGLE) {pnmv.rotate_x(PASSAGE_ANGLE); pathg.rotangle = PASSAGE_ANGLE;}
					else if (dag < -0.5*PASSAGE_ANGLE) {pnmv.rotate_x(-PASSAGE_ANGLE); pathg.rotangle = -PASSAGE_ANGLE;}
					//printf("N%dP%d pnmv after rotate (%f,%f,%f)\n", node, pat0, pnmv[0], pnmv[1], pnmv[2]);		
#endif					
					bool rightpos = false;

					while (!rightpos)
					{
						dkeisa[0] = pnmv.dot(amesh.m_level[lastlevel].m_geom[pathg.closecell].keisa[0]);
						dkeisa[1] = pnmv.dot(amesh.m_level[lastlevel].m_geom[pathg.closecell].keisa[1]);
						dkeisa[2] = pnmv.dot(amesh.m_level[lastlevel].m_geom[pathg.closecell].keisa[2]);
						dhgcell = max(abs(dkeisa[0]), abs(dkeisa[1]));
						dhgcell = max(dhgcell, abs(dkeisa[2]));
						hgdxyz = pnmv/dhgcell;
						pathg.hgdis += hgdxyz.length();
						pathg.pt = patch[pat0].pc + patch[pat0].nv*pathg.hgdis;
						amesh.LocatePointinComptDomain_AllDomain(pathg.pt, lastlevel, pathg.closecell, dkeisa_tocc);
						amesh.ReversePeriodicSide(lastlevel, pathg.closecell, rev_dir);
						if (rev_dir != 0)
						{
							if (dir_rev_num > 0)
							{
								printf("The direction for the IB cell has been reversed once in FindImageCell_GhostIBCell!!!\n");
								MPI_Abort(MPI_COMM_WORLD, 1048);
							}
							++dir_rev_num;
#ifdef PASSAGE_ANGLE						
							if (rev_dir == 1) {pnmv.rotate_x(PASSAGE_ANGLE); pathg.rotangle += PASSAGE_ANGLE;}
							else if (rev_dir == -1) {pnmv.rotate_x(-PASSAGE_ANGLE); pathg.rotangle += -PASSAGE_ANGLE;}
#endif						
						}
						amesh.FindHGCellIntpCell(pathg, lastlevel, dkeisa_tocc, hgnode);
						if (hgnode != -1)
						{
							for (Point_iterator p1(0,3); p1.end(); ++p1)
							{
								int ahc = amesh.m_level[lastlevel].m_box[pathg.closecell].neib[p1.i][p1.j][p1.k];
								if (ahc > -1)
								{
									if (amesh.m_level[lastlevel].m_box[ahc].bkpid > -1)
									{
										int abkpid = amesh.m_level[lastlevel].m_box[ahc].bkpid;
										int bkpneib = amesh.m_level[lastlevel].m_box[abkpid].neib[2-p1.i][2-p1.j][2-p1.k];
										if (bkpneib > -1)
										{
											patchbox[i0] = bkpneib;
											goto CELLSTART;										
										}
									}
								}
							}
							printf("Patch new node must have intp cell for a received one!!!\n");
							MPI_Abort(MPI_COMM_WORLD, 1512);
						}
#ifdef WALLINTPFLAG_NOTSOLID
						amesh.IntpcellNotSolid(rightpos, lastlevel, pathg);
#endif					
#ifdef WALLINTPFLAG_NOTINFECTED							
						amesh.IntpcellNotSolid_NotInfect(rightpos, lastlevel, pathg);
#endif	
						if (amesh.m_level[lastlevel].m_box[pathg.closecell].type == Dmghost)
						{
							printf("Patch hg cell close box is a Dmghost!!!\n");
							MPI_Abort(MPI_COMM_WORLD, 1280);
						}
					}
					if (!amesh.m_level[lastlevel].m_box.isnormal(pathg.closecell))
					{
						printf("N%dR%d recv Patch %d (%f,%f,%f) close box is a ghost cell!!! patchbox (%d,%d,%d) hg close cell (%d,%d,%d) hgdis %f!!!\n", 
							node, srank, pat0, patch[pat0].pc[0], patch[pat0].pc[1], patch[pat0].pc[2],
							amesh.m_level[lastlevel].m_box[patchbox[i0]].ix(),
							amesh.m_level[lastlevel].m_box[patchbox[i0]].iy(),
							amesh.m_level[lastlevel].m_box[patchbox[i0]].iz(),
							amesh.m_level[lastlevel].m_box[pathg.closecell].ix(),
							amesh.m_level[lastlevel].m_box[pathg.closecell].iy(),
							amesh.m_level[lastlevel].m_box[pathg.closecell].iz(),
							pathg.hgdis);
						printf("N%dP%d nv is (%f,%f,%f) pnmv is (%f,%f,%f)!!!\n", node, pat0,
							patch[pat0].nv[0], patch[pat0].nv[1], patch[pat0].nv[2],
							pnmv[0], pnmv[1], pnmv[2]);
						MPI_Abort(MPI_COMM_WORLD, 1635);
					}
					amesh.CheckHgIntpcells(lastlevel, pathg);
					amesh.ComptHgIntpCoef(pathg, dkeisa_tocc, lastlevel);
#ifdef DEBUG					
					printf("N%dR%d Patch %d (%f,%f,%f) construct a new patch hg from an other node patchbox (%d,%d,%d) hg close cell (%d,%d,%d) hgdis %f!!!\n", 
						node, srank, pat0, patch[pat0].pc[0], patch[pat0].pc[1], patch[pat0].pc[2],
						amesh.m_level[lastlevel].m_box[patchbox[i0]].ix(),
						amesh.m_level[lastlevel].m_box[patchbox[i0]].iy(),
						amesh.m_level[lastlevel].m_box[patchbox[i0]].iz(),
						amesh.m_level[lastlevel].m_box[pathg.closecell].ix(),
						amesh.m_level[lastlevel].m_box[pathg.closecell].iy(),
						amesh.m_level[lastlevel].m_box[pathg.closecell].iz(),
						pathg.hgdis);
					printf("N%dP%d nv is (%f,%f,%f) pnmv is (%f,%f,%f)!!!\n", node, pat0,
						patch[pat0].nv[0], patch[pat0].nv[1], patch[pat0].nv[2],
						pnmv[0], pnmv[1], pnmv[2]);
#endif					
					NEXTCYCLE2:;
				}
			}
			++cycle_num;
			if (cycle_num > nodenum)
			{
				printf("Patch location has cycled for %d times which is larger than the node number!!!\n", cycle_num);
				MPI_Abort(MPI_COMM_WORLD, 506);
			}
			goto SENDPATCHTONODE;
		}
		MPI_Barrier(share_comm);
#ifdef DEBUG		
		if (nrank == 0) printf("Body %d FindPatchHangingNode cycle number is %d!!!\n", bodyindex, cycle_num);
#endif		
	}
#else
	void Body::FindPatchHangingNode(Mesh & amesh)
	{
// 		lastlevel = amesh.MyCurNum()-1;
// 		int pats = patch.ps();
// 		int pate = patch.pe();
// 		Pointxyz dxyz, newbcxyz, dkeisa, hgdxyz, dkeisa_tocc;
// 		double dhgcell;
// 		int hgnode, rev_dir;
// 		PatchamongNodes pack0;
// 		Pointxyz pnmv;
// 		for (int i = pats; i < pate; ++i)
// 		{
// 			int i0 = i - pats;
// 			int ci0 = patchbox[i0];
// 			Assert(amesh.m_level[lastlevel].m_box.isnorg(ci0),"Patch close box error!!!", 1383);
// 			if (amesh.m_level[lastlevel].m_box.isnormal(ci0))
// 			{
// 				HGCell & pathg = patch[i].hgc;				
// 				pathg.hgdis = 0.0;
// 				pathg.attach_pt = patch[i].pc;
// 				pathg.closecell = ci0;
// 				bool rightpos = false;
// 				pnmv = patch[i].nv;
// #ifdef PASSAGE_ANGLE
// 				double dag = ComptPointAngle_Rotate_X(patch[i].pc, amesh.m_level[lastlevel].m_geom[pathg.closecell].boxcenter);
// 				if (dag > 0.5*PASSAGE_ANGLE) pnmv.rotate_x(PASSAGE_ANGLE);
// 				else if (dag < -0.5*PASSAGE_ANGLE) pnmv.rotate_x(-PASSAGE_ANGLE);			
// #endif								
// 				while (!rightpos)
// 				{
// 					dkeisa[0] = pnmv.dot(amesh.m_level[lastlevel].m_geom[pathg.closecell].keisa[0]);
// 					dkeisa[1] = pnmv.dot(amesh.m_level[lastlevel].m_geom[pathg.closecell].keisa[1]);
// 					dkeisa[2] = pnmv.dot(amesh.m_level[lastlevel].m_geom[pathg.closecell].keisa[2]);
// 					dhgcell = max(abs(dkeisa[0]), abs(dkeisa[1]));
// 					dhgcell = max(dhgcell, abs(dkeisa[2]));
// 					hgdxyz = pnmv/dhgcell;
// 					pathg.hgdis += hgdxyz.length();
// 					pathg.pt = patch[i].pc + patch[i].nv*pathg.hgdis;
// 					amesh.LocatePointinComptDomain_AllDomain(pathg.pt, lastlevel, pathg.closecell, dkeisa_tocc);
// 					amesh.ReversePeriodicSide(lastlevel, pathg.closecell, rev_dir);
// 					if (rev_dir != 0)
// 					{
// 						//printf("N%d Patch%d will be switched to the other side!!!\n", node, i);
// #ifdef PASSAGE_ANGLE						
// 						if (rev_dir == 1) pnmv.rotate_x(PASSAGE_ANGLE);
// 						else if (rev_dir == -1) pnmv.rotate_x(-PASSAGE_ANGLE);
// #endif						
// 					}									
// 					amesh.FindHGCellIntpCell(pathg, lastlevel, dkeisa_tocc, hgnode);
// 					if (hgnode != -1)
// 					{
// 						printf("N%d Patch%d intp cell is out of domain!!!\n", node, i, hgnode);
// 						MPI_Abort(MPI_COMM_WORLD, 1715)
// 					}
// #if WALLINTPFLAG == NOTSOLID
// 					amesh.IntpcellNotSolid(rightpos, lastlevel, pathg);
// #elif WALLINTPFLAG == NOTINFECTED							
// 					amesh.IntpcellNotSolid_NotInfect(rightpos, lastlevel, pathg);
// #endif					
// 					if (amesh.m_level[lastlevel].m_box[pathg.closecell].type == Dmghost)
// 					{
// 						printf("Patch hg cell close box is a Dmghost!!!\n");
// 						MPI_Abort(MPI_COMM_WORLD, 1280);
// 					}	
// 				}
// 				if (amesh.m_level[lastlevel].m_box.isghost(pathg.closecell))
// 				{
// 					//Point neibdir;
// 					hgnode = amesh.FindNeibBlock(amesh.m_level[lastlevel].m_box[pathg.closecell], lastlevel);
// 					//hgnode = m_dm.neibblocks[neibdir[0]][neibdir[1]][neibdir[2]];
// 					if (hgnode == node)
// 					{
// 						printf("N%d P%d The patch (%f,%f,%f) hgdis %f close box (%d,%d,%d)(%f,%f,%f) hg cell close ghost cell (%d,%d,%d)(%f,%f,%f) type %d ptype %d signdis %f can not be in the same node!!!"
// 							"May be the patch is very close to a normal ghost, i.e. the boundary of the mesh layer!!!\n", 
// 							hgnode, i,
// 							patch[i].pc[0], patch[i].pc[1], patch[i].pc[2],
// 							pathg.hgdis,
// 							amesh.m_level[lastlevel].m_box[patchbox[i0]].ix(),
// 							amesh.m_level[lastlevel].m_box[patchbox[i0]].iy(),
// 							amesh.m_level[lastlevel].m_box[patchbox[i0]].iz(),
// 							amesh.m_level[lastlevel].m_geom[patchbox[i0]].boxcenter[0],
// 							amesh.m_level[lastlevel].m_geom[patchbox[i0]].boxcenter[1],
// 							amesh.m_level[lastlevel].m_geom[patchbox[i0]].boxcenter[2],
// 							amesh.m_level[lastlevel].m_box[pathg.closecell].ix(),
// 							amesh.m_level[lastlevel].m_box[pathg.closecell].iy(),
// 							amesh.m_level[lastlevel].m_box[pathg.closecell].iz(),
// 							amesh.m_level[lastlevel].m_geom[pathg.closecell].boxcenter[0],
// 							amesh.m_level[lastlevel].m_geom[pathg.closecell].boxcenter[1],
// 							amesh.m_level[lastlevel].m_geom[pathg.closecell].boxcenter[2],
// 							amesh.m_level[lastlevel].m_box[pathg.closecell].type,
// 							amesh.m_level[lastlevel].m_box[pathg.closecell].ptype,
// 							amesh.m_level[lastlevel].m_box[pathg.closecell].pair.signdis);
// 						MPI_Abort(MPI_COMM_WORLD, 1433);
// 					}
// #ifdef DEBUG					
// 					printf("N%d Patch%d (%f,%f,%f) will be sent to another node %d patchbox (%d,%d,%d) the hg close cell (%d,%d,%d) hgdis %f is a ghost!!!\n", node, i, 
// 							patch[i].pc[0], patch[i].pc[1], patch[i].pc[2],
// 							hgnode,
// 							amesh.m_level[lastlevel].m_box[patchbox[i0]].ix(),
// 							amesh.m_level[lastlevel].m_box[patchbox[i0]].iy(),
// 							amesh.m_level[lastlevel].m_box[patchbox[i0]].iz(),
// 							amesh.m_level[lastlevel].m_box[pathg.closecell].ix(),
// 							amesh.m_level[lastlevel].m_box[pathg.closecell].iy(),
// 							amesh.m_level[lastlevel].m_box[pathg.closecell].iz(),
// 							pathg.hgdis);
// 					printf("N%dP%d nv is (%f,%f,%f) pnmv is (%f,%f,%f)!!!\n", node, i,
// 						patch[i].nv[0], patch[i].nv[1], patch[i].nv[2],
// 						pnmv[0], pnmv[1], pnmv[2]);
// #endif					
// 					goto NEXTCYCLE;

// 				}
// 				amesh.CheckHgIntpcells(lastlevel, pathg);
// 				amesh.ComptHgIntpCoef(pathg, dkeisa_tocc, lastlevel);
// 				NEXTCYCLE:;					
// 			}
// 		}
// 		MPI_Barrier(share_comm);
// #ifdef DEBUG		
// 		if (nrank == 0) printf("Body %d FindPatchHangingNode cycle number is %d!!!\n", bodyindex, cycle_num);
// #endif		
	}
#endif	
// 		void Body::FindPatchHangingNode(Mesh & amesh)
// 	{
// 		lastlevel = amesh.MyCurNum()-1;
// 		int pats = patch.ps();
// 		int pate = patch.pe();
// 		Pointxyz dxyz, newbcxyz, dkeisa, hgdxyz, dkeisa_tocc;
// 		double dhgcell;
// 		int hgnode, rev_dir;
// 		PatchamongNodes pack0;
// 		Pointxyz pnmv;
// 		for (int i = pats; i < pate; ++i)
// 		{
// 			int i0 = i - pats;
// 			int ci0 = patchbox[i0];
// 			//if (amesh.m_level[lastlevel].m_box.isnormal(ci0))
// 			if (node == patch[i].node)
// 			{
// 				HGCell & pathg = patch[i].hgc;				
// 				pathg.hgdis = 0.0;
// 				pathg.attach_pt = patch[i].pc;
// 				pathg.closecell = ci0;
// 				bool rightpos = false;
// 				pnmv = patch[i].nv;
// 				bool cc_is_normal = amesh.m_level[lastlevel].m_box.isnormal(ci0);
// #ifdef PASSAGE_ANGLE
// 				double dag = ComptPointAngle_Rotate_X(patch[i].pc, amesh.m_level[lastlevel].m_geom[pathg.closecell].boxcenter);
// 				if (dag > 0.5*PASSAGE_ANGLE) pnmv.rotate_x(PASSAGE_ANGLE);
// 				else if (dag < -0.5*PASSAGE_ANGLE) pnmv.rotate_x(-PASSAGE_ANGLE);			
// #endif								
// 				while (!rightpos)
// 				{
// 					dkeisa[0] = pnmv.dot(amesh.m_level[lastlevel].m_geom[pathg.closecell].keisa[0]);
// 					dkeisa[1] = pnmv.dot(amesh.m_level[lastlevel].m_geom[pathg.closecell].keisa[1]);
// 					dkeisa[2] = pnmv.dot(amesh.m_level[lastlevel].m_geom[pathg.closecell].keisa[2]);
// 					dhgcell = max(abs(dkeisa[0]), abs(dkeisa[1]));
// 					dhgcell = max(dhgcell, abs(dkeisa[2]));
// 					hgdxyz = pnmv/dhgcell;
// 					pathg.hgdis += hgdxyz.length();
// 					pathg.pt = patch[i].pc + patch[i].nv*pathg.hgdis;
// 					amesh.LocateHGCell(lastlevel, pathg, patch[i].pc, hgnode, dkeisa_tocc, rev_dir);
// 					if (hgnode != -1 && hgnode != node)
// 					{
// 						pack0.newsendmeg(i, hgnode);
// 						printf("N%d Patch%d will be sent to another node %d!!!\n", node, i, hgnode);
// 						goto NEXTCYCLE;
// 					}
// 					if (hgnode == node)
// 					{
// 						printf("N%d Patch%d will be switched to the other side!!!\n", node, i);
// #ifdef PASSAGE_ANGLE						
// 						if (rev_dir == 1) pnmv.rotate_x(PASSAGE_ANGLE);
// 						else if (rev_dir == -1) pnmv.rotate_x(-PASSAGE_ANGLE);
// 						else
// 						{
// 							printf("During FindPatchHangingNode, when the hg close cell has been switched, the rev_dir should be 1 or -1 but it is %d\n", rev_dir);
// 						}
// #endif						
// 					}
// #if WALLINTPFLAG == NOTSOLID
// 					amesh.IntpcellNotSolid(rightpos, lastlevel, pathg);
// #elif WALLINTPFLAG == NOTINFECTED							
// 					amesh.IntpcellNotSolid_NotInfect(rightpos, lastlevel, pathg);
// #endif					
// 					if (amesh.m_level[lastlevel].m_box[pathg.closecell].type == Dmghost)
// 					{
// 						printf("Patch hg cell close box is a Dmghost!!!\n");
// 						MPI_Abort(MPI_COMM_WORLD, 1280);
// 					}
// 				}
// 				amesh.CheckHgIntpcells(lastlevel, pathg);
// 				amesh.ComptHgIntpCoef(pathg, dkeisa_tocc, lastlevel);
// 				NEXTCYCLE:;					
// 			}
// 		}
// 		MPI_Barrier(share_comm);
// 		int cycle_num = 0;
// 		int patnode;
// 		SENDPATCHTONODE:;
// 		pack0.dataexchange();
// 		int recvpatnum = pack0.recv_pats.size();
// 		if (pack0.total_exc_num > 0)
// 		{
// 			vector<BoxlocChange> mdf_pats = pack0.recv_pats;
// 			pack0.initrecvarray();
// 			pack0.initsendarray();
// 			for (int i = 0; i < recvpatnum; ++i)
// 			{
// 				int pat0 = mdf_pats[i].oldloc;
// 				Assert(pat0 >= patch.ps() && pat0 < patch.pe(), "The patch is not in this node!!!", 484);
// 				int newnode = mdf_pats[i].newloc;
// 				patch[pat0].node = newnode;
// 				if (newnode == node)
// 				{
// 					int i0 = pat0 - patch.ps();
// 					printf("N%dR%d Recv patch %d (%f,%f,%f)!!!close box is (%d,%d,%d)\n",
// 						node, srank,
// 						pat0, patch[pat0].pc[0], patch[pat0].pc[1], patch[pat0].pc[2],
// 						amesh.m_level[lastlevel].m_box[patchbox[i0]].ix(),
// 						amesh.m_level[lastlevel].m_box[patchbox[i0]].iy(),
// 						amesh.m_level[lastlevel].m_box[patchbox[i0]].iz());	
// 					amesh.LocatePointinComptDomain_AllDomain(patch[pat0].pc, lastlevel, patchbox[i0], patnode, dkeisa_tocc);
// 					if (patnode != -1)
// 					{
// 						printf("Where the patch is going? Present node is %d new node is %d!!!\n", node, patnode);
// 						MPI_Abort(MPI_COMM_WORLD, 1309);
// 					}
// 					//%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// 					HGCell & pathg = patch[pat0].hgc;				
// 					pathg.hgdis = 0.0;
// 					pathg.attach_pt = patch[pat0].pc;
// 					pathg.closecell = patchbox[i0];
// 					pnmv = patch[pat0].nv;
// #ifdef PASSAGE_ANGLE
// 					double dag = ComptPointAngle_Rotate_X(patch[pat0].pc, amesh.m_level[lastlevel].m_geom[pathg.closecell].boxcenter);
// 					if (dag > 0.5*PASSAGE_ANGLE) pnmv.rotate_x(PASSAGE_ANGLE);
// 					else if (dag < -0.5*PASSAGE_ANGLE) pnmv.rotate_x(-PASSAGE_ANGLE);			
// #endif					
// 					bool rightpos = false;
// 					while (!rightpos)
// 					{
// 						dkeisa[0] = pnmv.dot(amesh.m_level[lastlevel].m_geom[pathg.closecell].keisa[0]);
// 						dkeisa[1] = pnmv.dot(amesh.m_level[lastlevel].m_geom[pathg.closecell].keisa[1]);
// 						dkeisa[2] = pnmv.dot(amesh.m_level[lastlevel].m_geom[pathg.closecell].keisa[2]);
// 						dhgcell = max(abs(dkeisa[0]), abs(dkeisa[1]));
// 						dhgcell = max(dhgcell, abs(dkeisa[2]));
// 						hgdxyz = patch[pat0].nv/dhgcell;
// 						pathg.hgdis += hgdxyz.length();
// 						pathg.pt = patch[pat0].pc + patch[pat0].nv*pathg.hgdis;
// 						amesh.LocateHGCell(lastlevel, pathg, patch[pat0].pc, hgnode, dkeisa_tocc, rev_dir);
// 						printf("N%dR%d Patch %d hg close box is (%d,%d,%d) hgnode is %d\n", node, srank, pat0,
// 							amesh.m_level[lastlevel].m_box[pathg.closecell].ix(),
// 							amesh.m_level[lastlevel].m_box[pathg.closecell].iy(),
// 							amesh.m_level[lastlevel].m_box[pathg.closecell].iz(), hgnode);
// 						if (hgnode != -1)
// 						{
// 							if (hgnode == node)
// 							{
// #ifdef PASSAGE_ANGLE								
// 								if (rev_dir == 1) pnmv.rotate_x(PASSAGE_ANGLE);
// 								else if (rev_dir == -1) pnmv.rotate_x(-PASSAGE_ANGLE);
// 								else
// 								{
// 									printf("During FindPatchHangingNode recv from other node," 
// 										   " when the hg close cell has been switched, the rev_dir should be 1 or -1 but it is %d\n", rev_dir);
// 									MPI_Abort(MPI_COMM_WORLD, 1479);
// 								}
// #endif								
// 							}
// 							// else
// 							// {
// 							// 	printf("N%d The ghost IB cell will be switched to other node %d!!!\n", node, hgnode);
// 							// 	MPI_Abort(MPI_COMM_WORLD, 1504);
// 							// 	pack0.newsendmeg(i, hgnode);
// 							// 	printf("N%d receive patch %d but will send it to N%d\n",node, pat0, hgnode);
// 							// 	goto NEXTCYCLE2;
// 							// }
// 						}
// 						if (hgnode == -1) amesh.IntpcellNotSolid(rightpos, lastlevel, pathg);
// 						//amesh.IntpcellNotSolid_NotInfect(rightpos, lastlevel, pathg);
// 						if (amesh.m_level[lastlevel].m_box[pathg.closecell].type == Dmghost)
// 						{
// 							printf("Patch hg cell close box is a Dmghost!!!\n");
// 							MPI_Abort(MPI_COMM_WORLD, 1280);
// 						}
// 					}
// 					if (hgnode != -1 && hgnode != node)
// 					{
// 						printf("The patch close cell in new node should be a normal cell!!!\n");
// 						MPI_Abort(MPI_COMM_WORLD, 1449);
// 					}
// 					else
// 					{
// 						printf("N%dR%d Patch %d construct a new patch hg from an other node!!!\n", node, srank, pat0);
// 					}
// 					amesh.CheckHgIntpcells(lastlevel, pathg);
// 					amesh.ComptHgIntpCoef(pathg, dkeisa_tocc, lastlevel);
// 					NEXTCYCLE2:;

// 				}
// 			}
// 			++cycle_num;
// 			if (cycle_num > nodenum)
// 			{
// 				printf("Patch location has cycled for %d times which is larger than the node number!!!\n", cycle_num);
// 				MPI_Abort(MPI_COMM_WORLD, 506);
// 			}
// 			goto SENDPATCHTONODE;
// 		}
// 		MPI_Barrier(share_comm);
// 		if (nrank == 0) printf("Body %d FindPatchHangingNode cycle number is %d!!!\n", bodyindex, cycle_num);
// 	}

	void Body::FindHangingNode_Onepatch(Mesh & amesh, const int & i)
	{
		Surfpatch & onepat = patch[i];
		Pointxyz dkeisa, hgdxyz, dkeisa_tocc;
		HGCell & pathg = onepat.hgc;		
		pathg.hgdis = 0.0;
		pathg.attach_pt = patch[i].pc;
		pathg.closecell = patchbox[i-patch.ps()];
		Pointxyz pnmv = patch[i].nv;
		bool rightpos = false;
		int rev_dir, hgnode;
		int dir_rev_num = 0;
#ifdef PASSAGE_ANGLE
		pathg.rotangle = 0.0;
		double dag = ComptPointAngle_Rotate_X(onepat.pc, amesh.m_level[lastlevel].m_geom[pathg.closecell].boxcenter);
		if (dag > 0.5*PASSAGE_ANGLE) {pnmv.rotate_x(PASSAGE_ANGLE); pathg.rotangle = PASSAGE_ANGLE;}
		else if (dag < -0.5*PASSAGE_ANGLE) {pnmv.rotate_x(-PASSAGE_ANGLE); pathg.rotangle = -PASSAGE_ANGLE;}
#endif								
		while (!rightpos)
		{
			dkeisa[0] = pnmv.dot(amesh.m_level[lastlevel].m_geom[pathg.closecell].keisa[0]);
			dkeisa[1] = pnmv.dot(amesh.m_level[lastlevel].m_geom[pathg.closecell].keisa[1]);
			dkeisa[2] = pnmv.dot(amesh.m_level[lastlevel].m_geom[pathg.closecell].keisa[2]);
			double dhgcell = max(abs(dkeisa[0]), abs(dkeisa[1]));
			dhgcell = max(dhgcell, abs(dkeisa[2]));
			hgdxyz = pnmv/dhgcell;
			pathg.hgdis += hgdxyz.length();
			pathg.pt = onepat.pc + onepat.nv*pathg.hgdis;
			amesh.LocatePointinComptDomain_AllDomain(pathg.pt, lastlevel, pathg.closecell, dkeisa_tocc);
			amesh.ReversePeriodicSide(lastlevel, pathg.closecell, rev_dir);
			if (rev_dir != 0)
			{
				dir_rev_num += rev_dir;
				if (abs(dir_rev_num) > 1)
				{
					printf("The direction for the hg of a patch has been reversed once in FindHangingNode_Onepatch!!!\n");
					MPI_Abort(MPI_COMM_WORLD, 1048);
				}
#ifdef PASSAGE_ANGLE						
				if (rev_dir == 1) {pnmv.rotate_x(PASSAGE_ANGLE); pathg.rotangle += PASSAGE_ANGLE;}
				else if (rev_dir == -1) {pnmv.rotate_x(-PASSAGE_ANGLE); pathg.rotangle += -PASSAGE_ANGLE;}
#endif						
			}									
			amesh.FindHGCellIntpCell(pathg, lastlevel, dkeisa_tocc, hgnode);
			if (hgnode != -1)
			{
				//pack0.newsendmeg(i, hgnode);
				printf("N%d Patch%d will be sent to another node %d because the intp cell is out of domain!!!\n", node, i, hgnode);
				//goto NEXTCYCLE;
				MPI_Abort(MPI_COMM_WORLD, 2017);
			}
#ifdef WALLINTPFLAG_NOTSOLID
			amesh.IntpcellNotSolid(rightpos, lastlevel, pathg);
#endif			
#ifdef WALLINTPFLAG_NOTINFECTED							
			amesh.IntpcellNotSolid_NotInfect(rightpos, lastlevel, pathg);
#endif					
			if (amesh.m_level[lastlevel].m_box[pathg.closecell].type == Dmghost)
			{
				printf("Patch hg cell close box is a Dmghost!!!\n");
				MPI_Abort(MPI_COMM_WORLD, 1280);
			}	
		}
		amesh.CheckHgIntpcells(lastlevel, pathg);
		amesh.ComptHgIntpCoef(pathg, dkeisa_tocc, lastlevel);
	}

	void Body::Rotate_Axis_X(double & rot_angle, Pointxyz rot_bc)
	{
		int bps = allpoint.ps();
		int bpe = allpoint.pe();
		double y0, z0;
		Pointxyz dxyz;
		for (int i = bps; i < bpe; ++i)
		{
			dxyz = allpoint[i] - rot_bc;
			y0 = dxyz[1]*cos(rot_angle) - dxyz[2]*sin(rot_angle);
			z0 = dxyz[2]*cos(rot_angle) + dxyz[1]*sin(rot_angle);
			allpoint[i][1] = y0 + rot_bc[1];
			allpoint[i][2] = z0 + rot_bc[2];
		}
		bps = patch.ps();
		bpe = patch.pe();
		for (int i = bps; i < bpe; ++i)
		{
			patch[i].pc.pt_theta += rot_angle;
		}
		MPI_Barrier(share_comm);
	}

	void Body::Rotate_Axis_Z(double & rot_angle, Pointxyz rot_bc)
	{
		int bps = allpoint.ps();
		int bpe = allpoint.pe();
		double x0, y0;
		Pointxyz dxyz;
		for (int i = bps; i < bpe; ++i)
		{
			dxyz = allpoint[i] - rot_bc;
			x0 = dxyz[0]*cos(rot_angle) - dxyz[1]*sin(rot_angle);
			y0 = dxyz[1]*cos(rot_angle) + dxyz[0]*sin(rot_angle);
			allpoint[i][0] = rot_bc[0] + x0;
			allpoint[i][1] = rot_bc[1] + y0;
		}
		MPI_Barrier(share_comm);
	}

	void Body::Attachbox_Init(Mesh & amesh)
	{
		/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
		/*After amr, the box attached to the wall also needs to be modified*/
		/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
		patchbox.resize(patch.pe()-patch.ps()); patchbox.assign(patchbox.size(), -1);
		pdistobox.resize(patch.pe()-patch.ps()); pdistobox.assign(pdistobox.size(), 1000.0);
		lastlevel = amesh.MyCurNum()-1;
		Assert(lastlevel==0, "Attachbox_Init Error!!! The number of mesh level is larger than 1...", 343);
		//int lastlayerboxsize = amesh.LastLayerBox().size();
		int lastlayerboxsize = amesh.LastLayerBox().realsize();
		int lastlayerboxsize0 = amesh.LastLayerBox().size();
		int pps = patch.ps();
		int ppe = patch.pe();
		for (int i = pps; i < ppe; ++i)
		{
			int i0 = i-patch.ps();
			bool disflag = false;
			for (int j = 0; j < lastlayerboxsize0; ++j)
			{
				double dis = PatchDistoPoint(amesh.m_level[lastlevel].m_geom[j].boxcenter, i);
				if (dis < pdistobox[i0])
				{
					patchbox[i0] = j;
					pdistobox[i0] = dis;
					disflag = true;
				}
			}
			CycleNearCell(amesh, patchbox[i0], pdistobox[i0], i);
			// for (int j = lastlayerboxsize0; j < lastlayerboxsize; ++j)
			// {
			// 	double dis = PatchDistoPoint(amesh.m_level[lastlevel].m_geom[j].boxcenter, i);
			// 	if (dis < pdistobox[i0]-0.00000001 && amesh.m_level[lastlevel].m_box[j].type != Dmghost)
			// 	{
			// 		int is0 = patchbox[i0];
			// 		patchbox[i0] = j;
			// 		pdistobox[i0] = dis;
			// 		disflag = true;
			// 		if (i == 3634) printf("N%d patch 3634 old close box is (%d,%d,%d)(%f,%f,%f) new close box is (%d,%d,%d)(%f,%f,%f)\n", node,
			// 			amesh.m_level[lastlevel].m_box[is0].ix(),
			// 			amesh.m_level[lastlevel].m_box[is0].iy(),
			// 			amesh.m_level[lastlevel].m_box[is0].iz(),
			// 			amesh.m_level[lastlevel].m_geom[is0].boxcenter[0],
			// 			amesh.m_level[lastlevel].m_geom[is0].boxcenter[1],
			// 			amesh.m_level[lastlevel].m_geom[is0].boxcenter[2],
			// 			amesh.m_level[lastlevel].m_box[j].ix(),
			// 			amesh.m_level[lastlevel].m_box[j].iy(),
			// 			amesh.m_level[lastlevel].m_box[j].iz(),
			// 			amesh.m_level[lastlevel].m_geom[j].boxcenter[0],
			// 			amesh.m_level[lastlevel].m_geom[j].boxcenter[1],
			// 			amesh.m_level[lastlevel].m_geom[j].boxcenter[2]);
			// 	}
			// }
			if (!disflag)
			{
				printf("Attachbox_Init Error!!! Patch distance error...\n");
				MPI_Abort(MPI_COMM_WORLD, 365);
			}
			// if (amesh.m_level[lastlevel].m_box[patchbox[i0]].type == Dmghost)
			// {
			// 	printf("The patch is (%f,%f,%f) close box is (%d,%d,%d)(%f,%f,%f) distance %f the box is a domain ghost!!!\n",
			// 		patch[i].pc[0], patch[i].pc[1], patch[i].pc[2],
			// 		amesh.m_level[lastlevel].m_box[patchbox[i0]].ix(),
			// 		amesh.m_level[lastlevel].m_box[patchbox[i0]].iy(),
			// 		amesh.m_level[lastlevel].m_box[patchbox[i0]].iz(),
			// 		amesh.m_level[lastlevel].m_geom[patchbox[i0]].boxcenter[0],
			// 		amesh.m_level[lastlevel].m_geom[patchbox[i0]].boxcenter[1],
			// 		amesh.m_level[lastlevel].m_geom[patchbox[i0]].boxcenter[2], pdistobox[i0]);
			// 	MPI_Abort(MPI_COMM_WORLD, 375);
			// }
		}
		MPI_Barrier(share_comm);
	}

	void Body::Attachbox_NewLevel(Mesh & amesh, DataArray<rftag> & mtag, DataArray<Boxson<int> > & mson)
	{
		int pps = patch.ps();
		int ppe = patch.pe();
		lastlevel = amesh.MyCurNum()-1;
		int lastlayerboxsize = amesh.LastLayerBox().realsize();
		Assert(lastlevel > 0, "Attach to a non-positive level...", 395);
		int momlevel = lastlevel - 1;
		if (!level_twod_flag[lastlevel])
		{
			for (int i = pps; i < ppe; ++i)
			{
				int i0 = i-pps;
				bool disflag = false;
				// double olddis = pdistobox[i0];
				pdistobox[i0] = 1000.0;
				Assert(patchbox[i0] > -1, "The box must be attached to a level before it is to be attached to a new level!!!", 398);
				int tag0 = mtag[patchbox[i0]].tag;
				if (tag0 > -1)
				{			
					for (Point_iterator q(0,2); q.end(); ++q)
					{
						int s0 = mson[tag0].son[q.i][q.j][q.k];
						double dis = PatchDistoPoint(amesh.m_level[lastlevel].m_geom[s0].boxcenter, i);
						if (dis < pdistobox[i0])
						{
							pdistobox[i0] = dis;
							patchbox[i0] = s0;
							disflag = true;
						}
					}
				}
				else
				{
					for (int j = amesh.m_level[lastlevel].m_box.size(); j < lastlayerboxsize; ++j)
					{
						if (amesh.m_level[lastlevel].m_box[j].type == Blockghost)
						{
							double dis = PatchDistoPoint(amesh.m_level[lastlevel].m_geom[j].boxcenter, i);
							if (dis < pdistobox[i0]-0.00000001)
							{
								patchbox[i0] = j;
								pdistobox[i0] = dis;
								disflag = true;
							}
						}
					}
				}
				if (!disflag)
				{
					printf("The body is not attached to a new 3D level. The patch is (%f,%f,%f) Old close box is (%f,%f,%f) type %d The new level is %d...\n",
						patch[i].pc[0], patch[i].pc[1], patch[i].pc[2],
						amesh.m_level[momlevel].m_geom[patchbox[i0]].boxcenter[0],
						amesh.m_level[momlevel].m_geom[patchbox[i0]].boxcenter[1],
						amesh.m_level[momlevel].m_geom[patchbox[i0]].boxcenter[2],
						amesh.m_level[momlevel].m_box[patchbox[i0]].type,
						lastlevel);
					MPI_Abort(MPI_COMM_WORLD, 425);
				}
			}			
		}
		else
		{
			for (int i = pps; i < ppe; ++i)
			{
				int i0 = i-pps;
				bool disflag = false;
				pdistobox[i0] = 1000.0;
				Assert(patchbox[i0] > -1, "The box must be attached to a level before it is to be attached to a new level!!!", 433);
				int tag0 = mtag[patchbox[i0]].tag;
				if (tag0 > -1)
				{
					for (Point_iterator_2d q(0,2); q.end(); ++q)
					{
						int s0 = mson[tag0].son[q.i][q.j][q.k];
						double dis = PatchDistoPoint(amesh.m_level[lastlevel].m_geom[s0].boxcenter, i);
						if (dis < pdistobox[i0])
						{
							pdistobox[i0] = dis;
							patchbox[i0] = s0;
							disflag = true;
						}
					}
				}
				else
				{
					for (int j = amesh.m_level[lastlevel].m_box.size(); j < lastlayerboxsize; ++j)
					{
						if (amesh.m_level[lastlevel].m_box[j].type == Blockghost)
						{
							double dis = PatchDistoPoint(amesh.m_level[lastlevel].m_geom[j].boxcenter, i);
							if (dis < pdistobox[i0]-0.00000001)
							{
								patchbox[i0] = j;
								pdistobox[i0] = dis;
								disflag = true;
							}
						}
					}
				}
				if (!disflag)
				{
					printf("The body is not attached to a new 2D level. Tag is %d the new level is %d...\n", tag0, lastlevel);
						// The patch is (%f,%f,%f) Old close box is (%f,%f,%f)\n",
						// patch[i].pc[0], patch[i].pc[1], patch[i].pc[2],
						// amesh.m_level[momlevel].m_geom[oldpat].boxcenter[0],
						// amesh.m_level[momlevel].m_geom[oldpat].boxcenter[1],
						// amesh.m_level[momlevel].m_geom[oldpat].boxcenter[2]);
					MPI_Abort(MPI_COMM_WORLD, 425);
				}
				// printf("The patch is (%f,%f,%f) close box in level %d is (%f,%f,%f)\n",
				// 	patch[i].pc[0], patch[i].pc[1], patch[i].pc[2], lastlevel,
				// 	amesh.m_level[lastlevel].m_geom[pdistobox[i0]].boxcenter[0],
				// 	amesh.m_level[lastlevel].m_geom[pdistobox[i0]].boxcenter[1],
				// 	amesh.m_level[lastlevel].m_geom[pdistobox[i0]].boxcenter[2]);
			}
		}
		MPI_Barrier(share_comm);
	}

	void Body::Check_Attachbox_willberemoved(vector<Body> & abody, Mesh & amesh)
	{
		Point adp[14] = {Point(0,1,1),Point(2,1,1),Point(1,0,1),Point(1,2,1),Point(1,1,0),Point(1,1,2),
						Point(0,0,0),Point(2,0,0),Point(0,2,0),Point(2,2,0),
						Point(0,0,2),Point(2,0,2),Point(0,2,2),Point(2,2,2)};
		lastlevel = amesh.MyCurNum()-1;
		int lastbs = amesh.m_level[lastlevel].m_box.realsize();
		for (int i = 0; i < bodynum; ++i)
		{
			int pps = abody[i].patch.ps();
			int ppe = abody[i].patch.pe();
			for (int p = pps; p < ppe; ++p)
			{
				int p0 = p - pps;
				int a0 = abody[i].patchbox[p0];
				Assert(a0 > -1, "The attached box to be renewed must be non-negative!!!", 1703);
				if (abody[i].pdistobox[p0] < patch_out_distance)
				{
					if (amesh.m_level[lastlevel].m_box[a0].neib[1][1][1] < 0)
					{
						// if (node == abody[i].patch[p].node)
						// {
						// 	// printf("N%d Patch (%f,%f,%f) belongs to node %d The new attached box is (%d,%d,%d)(%f,%f,%f) and it will be removed box type is %d!!!\n",
						// 	// 	node, 
						// 	// 	abody[i].patch[p].pc[0],
						// 	// 	abody[i].patch[p].pc[1],
						// 	// 	abody[i].patch[p].pc[2],
						// 	// 	abody[i].patch[p].node,
						// 	// 	amesh.m_level[lastlevel].m_box[a0].ix(),
						// 	// 	amesh.m_level[lastlevel].m_box[a0].iy(),
						// 	// 	amesh.m_level[lastlevel].m_box[a0].iz(),
						// 	// 	amesh.m_level[lastlevel].m_geom[a0].boxcenter[0],
						// 	// 	amesh.m_level[lastlevel].m_geom[a0].boxcenter[1],
						// 	// 	amesh.m_level[lastlevel].m_geom[a0].boxcenter[2],
						// 	// 	amesh.m_level[lastlevel].m_box[a0].type);
						// 	// MPI_Abort(MPI_COMM_WORLD, 2026);
						// }
						// else
						{	
							int start_box = a0;
							int inc_dn = 1;
							while(amesh.m_level[lastlevel].m_box[start_box].neib[1][1][1] == -1 ||
								  amesh.m_level[lastlevel].m_box[start_box].type == Dmghost)
							{
								start_box += inc_dn;
								if (start_box >= lastbs)
								{
									start_box = a0-1;
									inc_dn = -1;
								}
								else if (start_box < 0)
								{
									printf("N%d The patch (%f,%f,%f) belongs to N%d close box (%d,%d,%d)(%f,%f,%f) will be removed and we do not find a new close box for the patch!!!\n",
										node,
										abody[i].patch[p].pc[0],
										abody[i].patch[p].pc[1],
										abody[i].patch[p].pc[2],
										abody[i].patch[p].node,
										amesh.m_level[lastlevel].m_box[a0].ix(),
										amesh.m_level[lastlevel].m_box[a0].iy(),
										amesh.m_level[lastlevel].m_box[a0].iz(),
										amesh.m_level[lastlevel].m_geom[a0].boxcenter[0],
										amesh.m_level[lastlevel].m_geom[a0].boxcenter[1],
										amesh.m_level[lastlevel].m_geom[a0].boxcenter[2]);
									MPI_Abort(MPI_COMM_WORLD, 2052);
								}
							}
							abody[i].patchbox[p0] = start_box;
						// abody[i].patchbox[p0] = -1;
						// for (Point_iterator di(0,3); di.end(); ++di)
						// {
						// 	int start_box = a0;
						// 	while (amesh.m_level[lastlevel].m_box[start_box].neib[1][1][1] == -1)
						// 	{
						// 		start_box = amesh.m_level[lastlevel].m_box[start_box].neib[di.i][di.j][di.k];
						// 		if (start_box == -1)
						// 		{
						// 			break;
						// 		}
						// 	}
						// 	if (start_box != -1)
						// 	{
						// 		abody[i].patchbox[p0] = start_box;
						// 		break;
						// 	}
						// }
						// if (abody[i].patchbox[p0] == -1)
						// {
						// 	printf("N%d The patch (%f,%f,%f) belongs to N%d close box (%d,%d,%d)(%f,%f,%f) will be removed and we do not find a new close box for the patch!!!\n",
						// 			node,
						// 			abody[i].patch[p].pc[0],
						// 			abody[i].patch[p].pc[1],
						// 			abody[i].patch[p].pc[2],
						// 			abody[i].patch[p].node,
						// 			amesh.m_level[lastlevel].m_box[a0].ix(),
						// 			amesh.m_level[lastlevel].m_box[a0].iy(),
						// 			amesh.m_level[lastlevel].m_box[a0].iz(),
						// 			amesh.m_level[lastlevel].m_geom[a0].boxcenter[0],
						// 			amesh.m_level[lastlevel].m_geom[a0].boxcenter[1],
						// 			amesh.m_level[lastlevel].m_geom[a0].boxcenter[2]);
						// 	MPI_Abort(MPI_COMM_WORLD, 2052);
						// }
						}
					}					
				}
			}
		}
		MPI_Barrier(share_comm);
	}

	void Body::Renew_Attachbox(vector<Body> & abody, Mesh & amesh)
	{		
		lastlevel = amesh.MyCurNum()-1;
		for (int i = 0; i < bodynum; ++i)
		{
			int pps = abody[i].patch.ps();
			int ppe = abody[i].patch.pe();
			for (int p = pps; p < ppe; ++p)
			{
				int p0 = p - pps;
				int a0 = abody[i].patchbox[p0];
				Assert(a0 > -1, "The attached box to be renewed must be non-negative!!!", 1703);
				if (abody[i].pdistobox[p0] < patch_out_distance)
				{
					if (amesh.m_level[lastlevel].m_box[a0].neib[1][1][1] < 0)
					{
						printf("N%d The new attached box is (%d,%d,%d)!!!\n",
							node, amesh.m_level[lastlevel].m_box[a0].ix(),
							amesh.m_level[lastlevel].m_box[a0].iy(),
							amesh.m_level[lastlevel].m_box[a0].iz());
						MPI_Abort(MPI_COMM_WORLD, 2397);
					}
					abody[i].patchbox[p0] = amesh.m_level[lastlevel].m_box[a0].neib[1][1][1];
					Assert(abody[i].patchbox[p0] > -1, "The renewed attached box becomes negative!!!", 1705);
				// if (node == abody[i].patch[p].node)
				// {
				// 	a0 = abody[i].patch[p].hgc.closecell;
				// 	Assert(a0 > -1, "The hg close cell to be renewed must be non-negative!!!", 1703);
				// 	abody[i].patch[p].hgc.closecell = amesh.m_level[lastlevel].m_box[a0].neib[1][1][1];
				// 	Assert(abody[i].patch[p].hgc.closecell > -1, "The hg close cell must be non-negative!!!", 1718);
				// 	for (Point_iterator pit(0,2); pit.end(); ++pit)
				// 	{
				// 		int an0 = abody[i].patch[p].hgc.intpcell[pit.i][pit.j][pit.k];
				// 		Assert(an0 > -1, "The hg intp cell to be renewed must be non-negative!!!", 1703);
				// 		abody[i].patch[p].hgc.intpcell[pit.i][pit.j][pit.k] = amesh.m_level[lastlevel].m_box[an0].neib[1][1][1];
				// 	}
				// }
				}
			}
		}
		MPI_Barrier(share_comm);
	}

	void Body::InitPatchsNode(Mesh & amesh)
  	{
  		
  		int pps = patch.ps();
  		int ppe = patch.pe();
#ifndef LOCALBODY  		
  		int patnum = ppe-pps;
  		vector<int> patnodeflag(patnum, 0);
  		vector<int> patnode(patnum,0);
  		for (int i = pps; i < ppe; ++i)
  		{
  			int i0 = i -pps;
  			pdistobox[i0] = 1000.0;
			CycleNearCell(amesh, patchbox[i0], pdistobox[i0], i);
  			if (amesh.m_level[lastlevel].m_box.isnormal(patchbox[i0]))
  			{
  				patnodeflag[i0] = 1;
  				patnode[i0] = node+1;
  			}
  		}
  		MPI_Allreduce(MPI_IN_PLACE, &patnodeflag[0], patnum, MPI_INT, MPI_SUM, nodecomm);
  		MPI_Allreduce(MPI_IN_PLACE, &patnode[0], patnum, MPI_INT, MPI_SUM, nodecomm);
  		for (int i = 0; i < patnum; ++i)
  		{
  			if (patnodeflag[i] != 1)
  			{
  				printf("N%d Patch %d (%f,%f,%f) was caught by %d node!!! close box is (%d,%d,%d) signdis is %f\n", 
  					node,
  					pps+i, 
  					patch[pps+i].pc[0],patch[pps+i].pc[1],patch[pps+i].pc[2],patnodeflag[i],
  					amesh.m_level[lastlevel].m_box[patchbox[i]].ix(),
  					amesh.m_level[lastlevel].m_box[patchbox[i]].iy(),
  					amesh.m_level[lastlevel].m_box[patchbox[i]].iz(),
  					amesh.m_level[lastlevel].m_box[patchbox[i]].pair.signdis);
  				if (amesh.m_level[lastlevel].m_box.isnormal(patchbox[i]))
  				{
  					printf("...Patch %d (%f,%f,%f) was caught by node %d box (%d,%d,%d)(%f,%f,%f)!!!\n",
  						pps+i,
  						patch[pps+i].pc[0],patch[pps+i].pc[1],patch[pps+i].pc[2],
  						node,
  						amesh.m_level[lastlevel].m_box[patchbox[i]].ix(),
  						amesh.m_level[lastlevel].m_box[patchbox[i]].iy(),
  						amesh.m_level[lastlevel].m_box[patchbox[i]].iz(),
  						amesh.m_level[lastlevel].m_geom[patchbox[i]].boxcenter[0],
  						amesh.m_level[lastlevel].m_geom[patchbox[i]].boxcenter[1],
  						amesh.m_level[lastlevel].m_geom[patchbox[i]].boxcenter[2]);
  				}
  				MPI_Abort(MPI_COMM_WORLD, 907);
  			}
  			patch[i+pps].node = patnode[i]-1;
  		}
  		MPI_Barrier(share_comm);
#else 
  		for (int i = pps; i < ppe; ++i)
  		{
  			int i0 = i -pps;
  			pdistobox[i0] = 1000.0;
			CycleNearCell(amesh, patchbox[i0], pdistobox[i0], i);
			patch[i].node = node;
		}
#endif	
  	}

  	void Body::ModifyBoxtoPatch(Mesh & amesh)
	{
#ifndef LOCALBODY		
		int cycle_num = 0;
		int pps = patch.ps();
		int ppe = patch.pe();
		/*Find the close box to the patch*/
		// double crit_range = dh[lastlevel][3]*0.4999;
		// double crit_range0 = dh[lastlevel][3]*0.85;
		//printf("Body %d patch start %d patch end %d\n", bodyindex, patch.ps(), patch.pe());
		PatchamongNodes pack0;
		int patnode;
		Pointxyz dkeisa_tocc;
		int rev_dir;
		for (int i = pps; i < ppe; ++i)
		{
			int i0 = i - pps;
			if (patch[i].node < 0 || patch[i].node >= nodenum)
			{
				printf("Error!!!Patch node is %d!!!\n", patch[i].node);
				MPI_Abort(MPI_COMM_WORLD, 461);
			}
			else if (patch[i].node == node)
			{
				int oldpat = patchbox[i0];
				amesh.LocatePointinComptDomain_AllDomain(patch[i].pc, lastlevel, patchbox[i0], dkeisa_tocc);
				amesh.ReversePeriodicSide(lastlevel, patchbox[i0], rev_dir);
			}
			else if (pdistobox[i0] < patch_out_distance)
			{	
				pdistobox[i0] = 1000.0;
				CycleNearCell(amesh, patchbox[i0], pdistobox[i0], i);
				amesh.ReversePeriodicSide(lastlevel, patchbox[i0], rev_dir);
			}
		}
		MPI_Barrier(share_comm);
#ifdef DEBUG		
		if (nrank == 0) printf("Body %d ModifyBoxtoPatch cycle number is %d!!!\n", bodyindex, cycle_num);
#endif
#endif		
	}

	void Body::ModifyBoxtoPatch_Somepatches(Mesh & amesh, vector<Boxloc> & bdloc, const int & locnum, vector<int> & bdloc_tag)
	{
		int pps = patch.ps();
		int ppe = patch.pe();
		
		Pointxyz dkeisa_tocc;
		int rev_dir;
		for (int i = 0; i < locnum; ++i)
		{
			if (bodyindex == bdloc[i][0] && (bdloc[i][1] >= pps && bdloc[i][1] < ppe))
			{
				if (pdistobox[bdloc[i][1] - pps] < patch_out_distance)
				{
					int ip0 = bdloc[i][1];
					int i0 = ip0 - pps;
					amesh.LocatePointinComptDomain_AllDomain(patch[ip0].pc, lastlevel, patchbox[i0], dkeisa_tocc);
					amesh.ReversePeriodicSide(lastlevel, patchbox[i0], rev_dir);
					if (amesh.m_level[lastlevel].m_box.isnormal(patchbox[i0]))
					{
						bdloc_tag[i] = 1;
					}
				}
			}
		}	
	}

	void Body::CompressPointArray()
	{
		DataArray<int> pttag;
		pttag.setnum_nocopy(allpoint.size(), 0);
		for (int i = allpoint.ps(); i < allpoint.pe(); ++i)
		{
			pttag[i] = i;
		}
		MPI_Barrier(share_comm);
		for (int i = allpoint.ps(); i < allpoint.pe(); ++i)
		{
			if (pttag[i] == i)
			{
				for (int j = i+1; j < allpoint.realsize(); ++j)
				{
					Pointxyz dxyz = Pointxyz(allpoint[i][0]-allpoint[j][0], allpoint[i][1]-allpoint[j][1], allpoint[i][2]-allpoint[j][2]);
					double dif_dxyz = dxyz.length();
					if (dif_dxyz < 1e-9)
					{
							pttag[j] = i;
#ifdef DEBUG							
							printf("Body %d Point %d (%f,%f,%f) is very close to point %d (%f,%f,%f), distance is %f!!!\n",
								bodyindex,
								i, allpoint[i][0], allpoint[i][1], allpoint[i][2],
								j, allpoint[j][0], allpoint[j][1], allpoint[j][2],
								dif_dxyz);
#endif							
							// MPI_Abort(MPI_COMM_WORLD, 166);
					}
				}
			}
		}
		MPI_Barrier(share_comm);
		for (int i = allpoint.ps(); i < allpoint.pe(); ++i)
		{
			if (pttag[i] != i)
			{
				while (pttag[pttag[i]] != pttag[i])
				{
					pttag[i] = pttag[pttag[i]];
				}
			}
		}
		MPI_Barrier(share_comm);
		for (int i = allpoint.ps(); i < allpoint.pe(); ++i)
		{
			if (pttag[i] != i)
			{
				allpoint.givehole(i);
			}
		}
		//printf("SR%d hold number %d!!!\n", srank, (int)allpoint.hole.size());
		MPI_Barrier(share_comm);
		for (int i = patch.pe(); i < patch.pe(); ++i)
		{
#ifdef QUA_ELEMENT				
			for (int j = 0; j < 4; ++j)
#endif
#ifdef TRI_ELEMENT
			for (int j = 0; j < 3; ++j)
#endif
			{
				int corpt0 = patch[i].corpt[j];
				patch[i].corpt[j] = pttag[corpt0];
			}
		}
		allpoint.holeplan(true);
		vector<BoxlocChange> & ptchange = allpoint.newloctag;
		int changenum = ptchange.size();
		for (int i = 0; i < changenum; ++i)
		{
			pttag[ptchange[i].oldloc] = ptchange[i].newloc;
		}
		MPI_Barrier(share_comm);
		for (int i = patch.ps(); i < patch.pe(); ++i)
		{
#ifdef QUA_ELEMENT				
			for (int j = 0; j < 4; ++j)
#endif
#ifdef TRI_ELEMENT
			for (int j = 0; j < 3; ++j)
#endif
			{
				int corpt0 = patch[i].corpt[j];
				patch[i].corpt[j] = pttag[corpt0];
			}
		}
		MPI_Barrier(share_comm);
		allpoint.CompressArray();
		if (nrank == 0) printf("Body %d point number after compress is %d\n", bodyindex, allpoint.size());
	}

	void Body::CheckBodyPatchCenter()
	{
		for (int i = patch.ps(); i < patch.pe(); ++i)
		{
			double ptlen = patch[i].pc.length();
			if (ptlen < 1e-6)
			{
				printf("Body %d Patch %d center is zero!!! Please check it!!!\n", bodyindex, i);
				MPI_Abort(MPI_COMM_WORLD, 1955);
			}
		}
	}

	void Body::SplitBodytoNode(Mesh & amesh, Body & localbody)
	{
		DataArray<int> patchtag;
		DataArray<int> pttag;
		patchtag.setnum_nocopy(patch.size(), 0);
		pttag.setnum_nocopy(allpoint.size(), 0);
		for (int i = patchtag.ps(); i < patchtag.pe(); ++i)
		{
			patchtag[i] = 0;
		}
		for (int i = pttag.ps(); i < pttag.pe(); ++i)
		{
			pttag[i] = 0;
		}
		MPI_Barrier(share_comm);
		int bps = patch.ps();
		int bpe = patch.pe();
		double criticdis = dh[0][0] + dh[0][1] + dh[0][2];
		for (int i = bps; i < bpe; ++i)
		{
			int i0 = i - bps;
			if (pdistobox[i0] < criticdis)
			{
				patchtag[i] = 1;
#ifdef QUA_ELEMENT				
				for (int p0 = 0; p0 < 4; ++p0)
#endif
#ifdef TRI_ELEMENT
				for (int p0 = 0; p0 < 3; ++p0)
#endif									
				{
					int pt0 = patch[i].corpt[p0];
					pttag[pt0] = 1;
				}
			}
		}
		MPI_Barrier(share_comm);
		int patnum = 0;
		int ptnum = 0;
		if (srank == 0)
		{
			for (int i = 0; i < patch.size(); ++i)
			{
				if (patchtag[i] != 0)
				{
					++patnum;
					patchtag[i] = patnum;
				}
			}
			for (int i = 0; i < allpoint.size(); ++i)
			{
				if (pttag[i] != 0)
				{
					++ptnum;
					pttag[i] = ptnum;
				}
			}
		}
		MPI_Barrier(share_comm);
		MPI_Bcast(&patnum, 1, MPI_INT, 0, share_comm);
		MPI_Bcast(&ptnum, 1, MPI_INT, 0, share_comm);
		// if (srank == 0)
		// {
		// 	printf("N%d localbody %d patnum %d ptnum %d\n", node, localbody.bodyindex, patnum, ptnum);
		// }
		MPI_Barrier(MPI_COMM_WORLD);
		localbody.patch.setnum_nocopy(patnum, 0);
		localbody.allpoint.setnum_nocopy(ptnum, 0);
		for (int i = patch.ps(); i < patch.pe(); ++i)
		{
			if (patchtag[i] != 0)
			{
				localbody.patch[patchtag[i]-1] = patch[i];
			}
		}
		for (int i = allpoint.ps(); i < allpoint.pe(); ++i)
		{
			if (pttag[i] != 0)
			{
				localbody.allpoint[pttag[i]-1] = allpoint[i];
			}
		}
		MPI_Barrier(share_comm);
		//GiveAFlag("Finish get localbody patch and points!!!", 2507);
		localbody.bodycenter = bodycenter;
		localbody.bcoffset = bcoffset;
		localbody.scale = scale;
		localbody.InitVelZero();
		localbody.surfptoff.resize(localbody.allpoint.pe()-localbody.allpoint.ps());
		for (int i = 0; i < localbody.surfptoff.size(); ++i)
		{
			int pt0 = i + localbody.allpoint.ps();
			localbody.surfptoff[i][0] = localbody.allpoint[pt0][0]-localbody.bodycenter[0];
			localbody.surfptoff[i][1] = localbody.allpoint[pt0][1]-localbody.bodycenter[1];
			localbody.surfptoff[i][2] = localbody.allpoint[pt0][2]-localbody.bodycenter[2];
		}
		for (int i = localbody.patch.ps(); i < localbody.patch.pe(); ++i)
		{
#ifdef QUA_ELEMENT
			for (int p0 = 0; p0 < 4; ++p0)
#endif
#ifdef TRI_ELEMENT
			for (int p0 = 0; p0 < 3; ++p0)
#endif			
			{
				int t0 = localbody.patch[i].corpt[p0];
				localbody.patch[i].corpt[p0] = pttag[t0]-1;
			}
		}
		MPI_Barrier(share_comm);
		//GiveAFlag("Start to CollectPatchNeib...", 5);
		localbody.CollectPatchNeib();
		//GiveAFlag("Finish CollectPatchNeib!!!", 5);
	}

	void Body::SortPatchCenterLocation(const int & left, const int & right, const int & dir)
	{
		if (nrank == 0)
		{
			if (left >= right)
			{
				return;
			}
			int ileft = left;
			int iright = right;
			//printf("Left is %d, right is %d!!!\n", left, right);
			while (ileft < iright)
			{
				double key = patch[ileft].pc[dir];
				while (ileft < iright && patch[iright].pc[dir] > key)
				{
					--iright;
				}
				patch.switcholdnew(ileft, iright);
				while (ileft < iright && patch[ileft].pc[dir] < key)
				{
					++ileft;
				}
				patch.switcholdnew(ileft, iright);
			}
			SortPatchCenterLocation(left, ileft-1, dir);
			SortPatchCenterLocation(ileft+1, right, dir);
		}
	}

