#include "Turbulence.H"

const double cb1 = 0.1355;
const double segma = 2.0/3.0;
const double cb2 = 0.622;
const double kappa = 0.41;
const double cw2 = 0.3;
const double cw3 = 2.0;
const double cv1 = 7.1;
const double ct3 = 1.2;
const double ct4 = 0.5;
const double cw1 = cb1/(kappa*kappa)+(1.0+cb2)/segma;
const double cv1three = cv1*cv1*cv1;
const double far_vis_ratio = 5.0;
const double max_viseddy = 300.0;
const double min_viseddy = 0.0;
const double min_viseddy_wall = 0.0;

const double f_cut = 0.15;

#ifdef TURBULENCE

double Viseddy2mut(double & cellvis, double & cellroe, double & cellmu)
{
	double bigx = cellvis/(cellmu/cellroe);
	if (bigx < 0.0) bigx = 0.0;
	double bigxthree = bigx*bigx*bigx;
	double fv1 = bigxthree/(bigxthree + cv1three);
	return cellroe*fv1*cellvis;
}

void yplus_as_function_of_re(double & hg_yplus, double & hg_re)
{
	double sq_hgre = sqrt(hg_re);
	if (hg_re <= 9.0) hg_yplus = 1.002*sq_hgre - 0.001223;
	else if (hg_re <= 37.0) hg_yplus = 1.272*sq_hgre - 0.9881*sqrt(sq_hgre) + 0.9054;
	else if (hg_re <= 98.0) hg_yplus = 0.2715*sq_hgre + 1.144*sqrt(sq_hgre) + 0.04718*hg_re;
	else if (hg_re <= 250.0) hg_yplus = 0.04921*sq_hgre + 1.576*sqrt(sq_hgre) + 0.05582*hg_re;
	else if (hg_re <= 790.0) hg_yplus = 0.5382*sq_hgre + 0.2478*sqrt(sq_hgre) + 0.04593*hg_re;
	else if (hg_re <= 2800.0) hg_yplus = 1.105*sq_hgre - 1.782*sqrt(sq_hgre) + 0.03933*hg_re;
	else hg_yplus = 2.643*sq_hgre - 10.52*sqrt(sq_hgre) + 0.03207*hg_re;
}

void uplus_as_function_of_yplus(double & uplus, double & yplus)
{
	double logyplus = log(yplus);
	double sq_yplus = sqrt(yplus);
	double qua_yplus = sqrt(sq_yplus);
	if (yplus < 3.0) uplus = yplus;
	else if (yplus < 8.5) uplus = -6.645*logyplus + 25.59*sq_yplus - 21.6*qua_yplus - 1.864*yplus;
	else if (yplus < 50.0) uplus = 12.04*logyplus - 4.079*sq_yplus - 4.504*qua_yplus + 0.1715*yplus;
	else if (yplus < 300.0) uplus = 3.031*logyplus - 0.6252*sq_yplus + 2.639*qua_yplus + 0.007553*yplus;
	else if (yplus < 1000.0) uplus = -3.535*logyplus - 1.104*sq_yplus + 13.89*qua_yplus + 0.004011*yplus;
	else uplus = 16.66*logyplus + 1.408*sq_yplus - 23.82*qua_yplus - 0.002909*yplus;
}

void ImageViseddy(IBCell & a_ibcell, Pointxyz & patv, Pointxyz & nmv, BoxtoWall & ibboxtopatch)
{
	HGCell & myhg = a_ibcell.hgc;
	double rt0 = myhg.hgdis + ibboxtopatch.signdis;
	double ib_rdis = ibboxtopatch.signdis/rt0;

	Pointxyz hg_vel_vect = Pointxyz(myhg.fv.u, myhg.fv.v, myhg.fv.w);
	double hg_veln = hg_vel_vect.dot(nmv);
	Pointxyz hg_veln_vect = Pointxyz(hg_veln*nmv[0], hg_veln*nmv[1], hg_veln*nmv[2]);
	Pointxyz hg_velt_vect = hg_vel_vect - hg_veln_vect;
	double hg_velt = hg_velt_vect.length();

	double wall_veln = patv.dot(nmv);
	Pointxyz wall_veln_vect = nmv*wall_veln;
	Pointxyz wall_velt_vect = patv - wall_veln_vect;
	double wall_velt = wall_velt_vect.length();

	Pointxyz df_vect = hg_velt_vect-wall_velt_vect;
	double df_velt = df_vect.length();

	double hg_mu = sqrt(pow(myhg.fv.T, 3))*((1.0+S_over_T_ref)/(myhg.fv.T+S_over_T_ref)); 
	double hg_viseddy = hg_mu/myhg.fv.roe;
	a_ibcell.hg_viseddy = hg_viseddy;
	double hg_re = df_velt*rt0/hg_viseddy*Re;
	double hg_yplus;
	yplus_as_function_of_re(hg_yplus, hg_re);
	a_ibcell.hg_yplus = hg_yplus;
	a_ibcell.hg_vt = df_velt;
	/*Normalized friction velocity*/
	a_ibcell.ut = hg_yplus*hg_viseddy/rt0/Re;
	/*Use normalized friction velocity to obtain the yplus*/
	a_ibcell.yplus = ib_rdis*hg_yplus;
	uplus_as_function_of_yplus(a_ibcell.uplus, a_ibcell.yplus);
	double ib_velt = a_ibcell.uplus*a_ibcell.ut;
	double ib_r0;
	if (abs(df_velt) < 0.00000001)
	{
	 	ib_r0 = 0.0;
	 	a_ibcell.tangdir = Pointxyz(0.0,0.0,0.0);
	}
	else
	{
		ib_r0 = ib_velt/df_velt;
		a_ibcell.tangdir = df_vect/df_velt;
	}
	Pointxyz old_ib_vel = Pointxyz(a_ibcell.fv.u, a_ibcell.fv.v, a_ibcell.fv.w);
	Pointxyz ib_nvect = wall_veln_vect + (hg_veln_vect-wall_veln_vect)*ib_rdis;
	a_ibcell.vn = hg_veln_vect-wall_veln_vect;
	Pointxyz ib_tvect = wall_velt_vect + df_vect*ib_r0;
	a_ibcell.fv.u = ib_nvect[0]+ib_tvect[0];
	a_ibcell.fv.v = ib_nvect[1]+ib_tvect[1];
	a_ibcell.fv.w = ib_nvect[2]+ib_tvect[2];
	a_ibcell.fv.viseddy = kappa*a_ibcell.yplus;
	a_ibcell.pre_grad = -((a_ibcell.fv.u-old_ib_vel[0])*nmv[0]+
						  (a_ibcell.fv.v-old_ib_vel[1])*nmv[1]+
						  (a_ibcell.fv.w-old_ib_vel[2])*nmv[2])/dt*a_ibcell.fv.roe;
	//a_ibcell.fv.p = myhg.fv.p;
	a_ibcell.fv.p = myhg.fv.p; // - a_ibcell.pre_grad*myhg.hgdis;
	//if (a_ibcell.fv.p < myhg.fv.p*0.5 || a_ibcell.fv.p > myhg.fv.p * 2.0) a_ibcell.fv.p = myhg.fv.p;
	// if (a_ibcell.fv.p < myhg.fv.p*0.5)
	// {
	// 	a_ibcell.fv.p = myhg.fv.p*0.5;
	// }
	// else if (a_ibcell.fv.p > myhg.fv.p*2.0)
	// {
	// 	a_ibcell.fv.p = myhg.fv.p*2.0;
	// }
	//a_ibcell.fv.T = myhg.fv.T + 0.5*Pr_cubicroot/Cp*(pow(hg_velt, 2) - (ib_tvect[0]*ib_tvect[0]+ib_tvect[1]*ib_tvect[1]+ib_tvect[2]*ib_tvect[2]));
	a_ibcell.fv.T = myhg.fv.T + 0.5*Pr_cubicroot/Cp*(pow(df_velt, 2) - pow(ib_velt, 2));
	//if (a_ibcell.fv.T < myhg.fv.T*0.5 || a_ibcell.fv.T > myhg.fv.T*2.0) a_ibcell.fv.T = myhg.fv.T;
	Get_roe_ideal_gas(a_ibcell.fv);
	Get_E(a_ibcell.fv);
	// printf("IB cell signdis %f hgdis %f hg velt %f hg veln %f hg re %f hgyplus %f ut %f ib yplus %f ib velt %f\n",
	// 	a_ibcell.boxpatch[0].signdis,myhg.hgdis,hg_velt, hg_veln, hg_re, hg_yplus, ut, ib_yplus, ib_velt);
#ifdef DEBUG	
	if (ibboxtopatch.signdis < 0.0)
	{
		printf("The ib cell for viseddy signdis negative!!! %f\n", ibboxtopatch.signdis);
		MPI_Abort(MPI_COMM_WORLD, 49);
	}
#endif	
}

void ImageViseddy_Domain(FlowVariables & refv, FlowVariables & bcvar, Pointxyz & patv, Pointxyz & nmv, 
	double & refdis, double & targetdis)
{
	double ib_rdis = targetdis/refdis;

	double hg_veln = refv.u*nmv[0]+refv.v*nmv[1]+refv.w*nmv[2];// refv.dot(nmv);
	Pointxyz hg_veln_vect(hg_veln*nmv[0], hg_veln*nmv[1], hg_veln*nmv[2]);
	Pointxyz hg_velt_vect(refv.u-hg_veln_vect[0],
						  refv.v-hg_veln_vect[1],
						  refv.w-hg_veln_vect[2]);
	double hg_velt = hg_velt_vect.length();

	double wall_veln = patv.dot(nmv);
	Pointxyz wall_veln_vect = nmv*wall_veln;
	Pointxyz wall_velt_vect = patv - wall_veln_vect;
	double wall_velt = wall_velt_vect.length();

	Pointxyz df_vect = hg_velt_vect-wall_velt_vect;
	double df_velt = df_vect.length();

	double hg_mu = sqrt(pow(refv.T, 3))*((1.0+S_over_T_ref)/(refv.T+S_over_T_ref)); 
	double hg_viseddy = hg_mu/refv.roe;
	double hg_re = df_velt*refdis/hg_viseddy*Re;
	double hg_yplus;
	yplus_as_function_of_re(hg_yplus, hg_re);
	double myut = hg_yplus*hg_viseddy/refdis/Re;
	double myyplus = ib_rdis*hg_yplus;
	double myuplus;
	uplus_as_function_of_yplus(myuplus, myyplus);
	double ib_velt = myuplus*myut;
	double ib_r0;
	if (abs(df_velt) < 0.00000001)
	{
	 	ib_r0 = 0.0;
	}
	else
	{
		ib_r0 = ib_velt/df_velt;
	}
	Pointxyz ib_nvect = wall_veln_vect + (hg_veln_vect-wall_veln_vect)*ib_rdis;
	Pointxyz ib_tvect = wall_velt_vect + df_vect*ib_r0;
	bcvar.u = ib_nvect[0]+ib_tvect[0];
	bcvar.v = ib_nvect[1]+ib_tvect[1];
	bcvar.w = ib_nvect[2]+ib_tvect[2];
	bcvar.viseddy = kappa*myyplus;
	bcvar.p = refv.p; // - a_ibcell.pre_grad*myhg.hgdis;
	bcvar.T = refv.T;
	Get_roe_ideal_gas(bcvar);
	Get_E(bcvar);
	if (bcvar.T < 0.0)
	{
		bcvar.showdata("domain wall function!");
		MPI_Abort(MPI_COMM_WORLD, 182);
	}
		
}

void Patchut(Surfpatch & mypatch, Pointxyz & patv, Pointxyz & nmv, Pointxyz & vtdir)
{
	HGCell & pathg = mypatch.hgc;
	Pointxyz hg_vel_vect = Pointxyz(pathg.fv.u, pathg.fv.v, pathg.fv.w);
	Pointxyz hg_vel_nvect = nmv*hg_vel_vect.dot(nmv); 
	Pointxyz hg_vel_tvect = hg_vel_vect - hg_vel_nvect; 

	Pointxyz patv_nvect = nmv*patv.dot(nmv);
	Pointxyz patv_tvect = patv - patv_nvect; 

	Pointxyz df_vect = hg_vel_tvect - patv_tvect;
	double df_velt = df_vect.length(); 

	if (df_velt < 1e-8) vtdir = Pointxyz(0.0,0.0,0.0);
	else vtdir = df_vect/df_velt;

	double hg_mu = sqrt(pow(pathg.fv.T, 3))*((1.0+S_over_T_ref)/(pathg.fv.T+S_over_T_ref));
	double hg_re = pathg.fv.roe*df_velt*pathg.hgdis/hg_mu*Re;

	yplus_as_function_of_re(mypatch.yplus, hg_re);
	/*Normalized friction velocity*/
	mypatch.ut = mypatch.yplus*(hg_mu/pathg.fv.roe)/pathg.hgdis/Re;
	/*---*/
	//double hg_velt = hg_vel_tvect.length();
	//double patv_velt = patv_tvect.length();
	//mypatch.wall_T = pathg.fv.T+0.5*Pr_cubicroot/Cp*(pow(hg_velt,2)-pow(patv_velt,2));
	/*---*/
	mypatch.wall_T = pathg.fv.T+0.5*Pr_cubicroot/Cp*pow(df_velt,2);
	mypatch.hg_vt = df_velt;
}

void Cellut(FlowVariables & cellfv, double & celldis_to_wall, Pointxyz & patv, Pointxyz & nmv, double & patut, Pointxyz & vtdir)
{
	Pointxyz hg_vel_vect = Pointxyz(cellfv.u - patv[0], cellfv.v-patv[1], cellfv.w-patv[2]);
	// double hg_veln = hg_vel_vect.dot(nmv);
	// Pointxyz hg_veln_vect = Pointxyz(hg_veln*nmv[0], hg_veln*nmv[1], hg_veln*nmv[2]);
	// Pointxyz hg_velt_vect = hg_vel_vect - hg_veln_vect;
	// double hg_velt = hg_velt_vect.length();
	double hg_veln, hg_velt;
	Pointxyz hg_velt_vect, hg_veln_vect;
	VelocityDecom(hg_vel_vect, nmv, hg_velt_vect, hg_velt, vtdir, hg_veln_vect, hg_veln);
	double hg_mu = sqrt(pow(cellfv.T, 3))*((1.0+S_over_T_ref)/(cellfv.T+S_over_T_ref));
	double hg_re = cellfv.roe*abs(hg_velt)*celldis_to_wall/hg_mu*Re;
	double hg_yplus;
	yplus_as_function_of_re(hg_yplus, hg_re);
	/*Normalized friction velocity*/
	patut = hg_yplus*(hg_mu/cellfv.roe)/celldis_to_wall/Re;
}

void VelocityDecom(Pointxyz & inputvel, Pointxyz & nmv,
				   Pointxyz & tangvel, double & tangspeed, Pointxyz & tangdir,
				   Pointxyz & nmvel, double & nmspeed)
{
	nmspeed = inputvel.dot(nmv);
	nmvel = Pointxyz(nmspeed*nmv[0], nmspeed*nmv[1], nmspeed*nmv[2]);
	tangvel = inputvel - nmvel;
	tangspeed = tangvel.length();
	if (abs(tangspeed) < 0.00000001)
	{
		tangdir = Pointxyz(0.0, 0.0, 0.0);
	}
	else
	{
		tangdir = tangvel/tangspeed;
	}
}

#endif