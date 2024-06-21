#include <vector>
#include <string>

#include "BCValues.H"
#include "Box.H"
#include "Domain.H"
#include "AMRLevel.H"
#include "FlowVariables.H"
#include "AMRSpaceTime.H"
#include "Pltoutput.H"
#include "AMR.H"
#include "Body.H"
#include "InputFunc.H"
#include "Runge_Kutta.H"

using namespace std;

double f_osc = .20;
double cycle_num = 100.0;

double dt = 0.0001;
double t_sum = 1000;
const int ts_sum = 10;

const double CFL = 5.0;
const double dt_ratio = 1.00;
int	dt_var_step = 500;
const int dts_plt = 500;
const int dts_time_var = 20;
const	int dts_space_var = 100000000;

const int max_mesh_level = 8;
const int init_level_num = 8;

const Point lowpt = Point(0,0,0);
// const Point	highpt = Point(200,150,1);
const Point	highpt = Point(250,188,1);

string meshfile = "xyz-1.fmt";
Pointxyz meshscale(10.0,10.0,10.0);

Pointxyz frame_rotate_vel = Pointxyz(0.0, 0.0, 0.0);
Pointxyz frame_rotate_origin = Pointxyz(0.0, 0.0, 0.0);
Pointxyz frame_trans_vel = Pointxyz(0.0, 0.0, 0.0);

const bool moveflag = true;

double roe_ref = 1.0;
double L_ref = 1.0;
double mu_ref = 1.0/400.0;
double u_ref = 1.0;
double T_ref = u_ref*u_ref/(0.1*0.1*Rg_ref*gama);

bool additional_output = false;
/*------------------------------Read file information------------------------*/
/*----If INPUT_SOLN is defined in the Makefile, this part is necessary-------*/
/*---------------------------------------------------------------------------*/
const int init_jump = 3;
const int filenodenum = 1;
const int filelevelnum = 7;
const int filebodynum = 2;
const char soln_file_name[30] = "amr_soln-720000-1.plt";
/*---------------------------------------------------------------------------*/
double taglocation = 0.0;
vector<int> adp_no(max_mesh_level, 0);
int grid_move_dir = 0;

double rot_amp = (0.0/180.0)*pai;

void AMR::ADPRule(vector<Body> & abody)
{
	adpflag = false;
	init_adjust_level = max_mesh_level-1;
	for (int i = 0; i < max_mesh_level; ++i)
	{
		m_level[i].torefine = false;
		m_level[i].toderefine = false;
	}
	double tag_diff = abody[0].inlineosc[1] - taglocation;
	//if (abs(tag_diff) > 2.0*dh[max_mesh_level-1][0])
    if (ts%int(2) == 0)
	{
		adpflag = true;
		taglocation = abody[0].inlineosc[1];
		if (tag_diff > 0.0) 
		{
			grid_move_dir = 1;
			for (int i = 0; i < max_mesh_level; ++i)
			{
				++adp_no[i];
			}
		}
		else
		{
			grid_move_dir = -1;
			for (int i = 0; i < max_mesh_level; ++i)
			{
				--adp_no[i];
			}
		}
		for (int i = 1; i < max_mesh_level; ++i)
		{
			int level_times = pow(2, max_mesh_level-i-1);
			if (adp_no[i]%level_times == 0)
			{
				m_level[i-1].torefine = true;
				m_level[i].toderefine = true;
				adp_no[i] = 0;
				init_adjust_level = min(init_adjust_level, i-1);
			}
		}
	}
	MPI_Barrier(share_comm);
}

void Body::Motion_Rule(Mesh & amesh, double & dt00)
{
    if (ts == 0)
    {
        rotosc[2] = rot_amp;
    }
    double oldloc = inlineosc[1];
    double oldrotloc = rotosc[2];
    double mvdt = dt*t_ref;
    free_tra_vel[1] *= u_ref;
    free_rot_vel[2] *= u_ref;
    inlineosc[1] *= L_ref;
    rotosc[2] *=L_ref;
    // 使用当前位置进行计算，下一步的位置和速度
    forthRK2DOF_Anal(inlineosc[1],free_tra_vel[1],
                     rotosc[2],free_rot_vel[2],
                     mass,inertia,damping,rotdamping,stiff,rotstiff,
                     force[1],moment[2],-0.1105,mvdt);
    // -0.1105
    // -0.070432
    free_tra_vel[1] = free_tra_vel[1]/u_ref;
    free_rot_vel[2] = free_rot_vel[2]/u_ref;
    inlineosc[1] /= L_ref;
    rotosc[2] /=L_ref;
    // 赋值给更新值
    inlineosc_new[1] = inlineosc[1];
    rotosc_new[2] = rotosc[2];
    // 当前步还原
    inlineosc[1] = oldloc;
    rotosc[2] = oldrotloc;
	ComptNewLocation();
	MPI_Barrier(share_comm);
}

void SetTwodFlag()
{
	level_twod_flag.resize(max_mesh_level);
	level_twod_flag[0] = false;
#if DIM > 2
	level_twod_flag[1] = false;
	level_twod_flag[2] = false;
	for (int i = 3; i < max_mesh_level; ++i)
	{
		level_twod_flag[i] = true;
	}
#else
	/*DO NOT CHANGE*/
	for (int i = 1; i < max_mesh_level; ++i)
	{
		level_twod_flag[i] = true;
	}
#endif	
}

void InitSpaceParams()
{

    double dh0[3] = {0.25, 0.25, 0.1};
    // double dh0[3] = {0.32, 0.32, 0.1};
	periodic[0] = false;
	periodic[1] = false;
	periodic[2] = true;
	ComptLevelgridsize(dh0);
}
//翼型中心位置
double foil_cord[3] = {31.25, 23.5, 1.0};
//double foil_cord[3] = {32.0, 24.0, 1.0};
double pitchx = 0.35;
double body_gap = 3.0;
double cy_cen_x = foil_cord[0]-pitchx-body_gap-0.5;
double cy_left_x = cy_cen_x - 0.5;
double foil_right_x = foil_cord[0] + 0.65;
double foil_left_x = foil_cord[0] - pitchx;
void Body::SetBodyParams(vector<Body> & abody, Mesh & amesh, vector<Body> & localbody)
{
	bodynum = 2;
	double ob[3][3] = {{foil_cord[0], foil_cord[1], foil_cord[2]},
					   {foil_cord[0], foil_cord[1], foil_cord[2]},
					   {cy_cen_x, foil_cord[1], foil_cord[2]}};
	double scale0[3] = {1.0, 1.0, 1.0};
	double cy_scale[3] = {0.01, 0.01, 0.01};
	for (int i = 0; i < bodynum; ++i)
	{
		abody.push_back(Body(i));
		abody.back().FindNearBody();
		abody.back().SetBodyCenter(ob[i][0], ob[i][1], 0.05);
		if (i < 2)
		{
            //abody.back().SetBodyOffset(ob[i][0]+(0.5-pitchx), ob[i][1], 0.0);
            abody.back().SetBodyOffset(ob[i][0]+(0.0-pitchx), ob[i][1], 0.0);
			abody.back().SetBodyScale(scale0[0], scale0[1], scale0[2]);
		}
		else
		{
			abody.back().SetBodyOffset(ob[i][0], ob[i][1], 0.0);
			abody.back().SetBodyScale(cy_scale[0], cy_scale[1], cy_scale[2]);
		}
		abody.back().InitVelZero();
		GiveAFlag("Start ReadShapeFile!!!", 5);
		if (i == 0)
		{
			//abody[0].ReadShapeFile("sphere_surf.dat");
			abody[0].ReadShapeFile("naca12_2d_upsurface_z0.1.ugrid");
		}
		else if (i == 1)
		{
			abody[1].ReadShapeFile("naca12_2d_downsurface_z0.1.ugrid");
		}
		else if (i == 2)
		{
			abody[2].ReadShapeFile("cylinder-2d-onelayer-max0p001-high10.ugrid");
		}
		GiveAFlag("Finish ReadShapeFile!!!", 5);
		if (i < 2)
		{
			abody[i].Rotate_Axis_Z(rot_amp, abody[i].bodycenter);
			abody[i].ComptSurfptOff();
			abody[i].ComptPatchParams();
		}	
		//abody[i].CheckBodyPatchCenter();
		GiveAFlag("Finish CheckBodyPatchCenter 1", 5);
		abody[i].Attachbox_Init(amesh);
		//abody[i].CheckBodyPatchCenter();
		GiveAFlag("Finish CheckBodyPatchCenter 2", 5);
        //设置质量 弹簧刚性 阻尼等系数，旋转运动还需要设置转动惯量、转动刚性以及转动阻尼。
        abody[i].mass = 0.164419;  // 质量比=20 无量纲质量
        double fndamp = 0.12285;      // 振动频率
        abody[i].stiff = pow(2.0*pai*fndamp, 2.0)*abody[i].mass;
        abody[i].damping = 0.0;
        abody[i].rotdamping = 0.0;
        abody[i].inertia = 0.010440; // 与质量比相关联
        double rotndamp = 0.12285;    // 扭转固有频率
        abody[i].rotstiff = pow(2.0*pai*rotndamp, 2.0)*abody[i].inertia;
	}
		for (int i = 0; i < bodynum; ++i)
		{
			abody[i].CheckBodyPatchCenter();
		}
	MPI_Barrier(MPI_COMM_WORLD);
	GiveAFlag("Finish read body and attach!!!", 5);
	if (bodynum > 0)
	{
		Body::ComptDisFromBoxtoBody_Init(amesh, abody);
	}
	GiveAFlag("Finish ComptDisFromBoxtoBody_Init!!!", 5);

}
void Domain::SetBlockNum()
{
	blocknum = Point(1,1,1);
}
void Domain::SetFaceProperty(DataArray<Box> & backbox, DataArray<BoxCellGeom> & backgeom)
{
	/*-------------------*/
	//refboxface has two integers
	//the first is direction x:0 y:1 z:2
	//the second is the face of the box in the above direction
	facenum = 6;

	faceup.resize(facenum);
	facedown.resize(facenum);
	facetype.resize(facenum);
	facedir.resize(facenum);
	refboxface.resize(facenum);
	/*------------Inlet---------------------------*/
	// faceup[0] = Point(0, 99999, 99999);
	// facedown[0] = Point(-10, -99999, -99999);
	// facetype[0] = -10;
	// facedir[0] = Point(1,0,0);
	faceup[0] = Point(0, 99999, 99999);
	facedown[0] = Point(-10, -99999, -99999);
	facetype[0] = -13;
	facedir[0] = Point(1,0,0);
	/*------------Outlet---------------------------*/
	// faceup[1] = Point(boxnum.ix()+10, 99999, 99999);
	// facedown[1] = Point(boxnum.ix()-1, -99999, -99999);
	// facetype[1] = -50;
	// facedir[1] = Point(-1,0,0);
	faceup[1] = Point(boxnum.ix()+10, 99999, 99999);
	facedown[1] = Point(boxnum.ix()-1, -99999, -99999);
	facetype[1] = -12;
	facedir[1] = Point(-1,0,0);
	/*------------Down---------------------------*/
	// faceup[2] = Point(34, 0, 99999);
	// facedown[2] = Point(-99999, -10, -99999);
	// facetype[2] = -51;
	// facedir[2] = Point(0,1,0);
	// refboxface[2] = Boxloc(1,0);

	// faceup[3] = Point(103, 0, 99999);
	// facedown[3] = Point(33, -10, -99999);
	// facetype[3] = -51;
	// facedir[3] = Point(0,1,0);
	// refboxface[3] = Boxloc(1,0);

	// faceup[4] = Point(99999, 0, 99999);
	// facedown[4] = Point(102, -10, -99999);
	// facetype[4] = -51;
	// facedir[4] = Point(0,1,0);
	// refboxface[4] = Boxloc(1,0);
	faceup[2] = Point(99999, 0, 99999);
	facedown[2] = Point(-99999, -10, -99999);
	facetype[2] = -50;
	facedir[2] = Point(0,1,0);
	/*------------Up---------------------------*/
	// faceup[5] = Point(99999, boxnum.iy()+10, 99999);
	// facedown[5] = Point(-99999, boxnum.iy()-1, -99999);
	// facetype[5] = -51;
	// facedir[5] = Point(0,-1,0);
	// refboxface[5] = Boxloc(1,1);
	faceup[3] = Point(99999, boxnum.iy()+10, 99999);
	facedown[3] = Point(-99999, boxnum.iy()-1, -99999);
	facetype[3] = -50;
	facedir[3] = Point(0,-1,0);
	/*------------Left---------------------------*/
	// faceup[6] = Point(99999, 99999, 0);
	// facedown[6] = Point(-99999, -99999, -10);
	// facetype[6] = -99;
	// facedir[6] = Point(0,0,1);
	faceup[4] = Point(99999, 99999, 0);
	facedown[4] = Point(-99999, -99999, -10);
	facetype[4] = -99;
	facedir[4] = Point(0,0,1);
	/*------------Right---------------------------*/
	// faceup[7] = Point(99999, 99999, boxnum.iz()+10);
	// facedown[7] = Point(-99999, -99999, boxnum.iz()-1);
	// facetype[7] = -99;
	// facedir[7] = Point(0,0,-1);
	faceup[5] = Point(99999, 99999, boxnum.iz()+10);
	facedown[5] = Point(-99999, -99999, boxnum.iz()-1);
	facetype[5] = -99;
	facedir[5] = Point(0,0,-1);
	/*----------------------------------------*/
	CheckPeriodicFaceType();
}

void BCValues::SetBCType_125_Inputinletbdlayer(FlowVariables & bcvar, FlowVariables & refvar, const int & faceindex, double & z0, int & zindex, Mesh & amesh)
{

}

void BCValues::SetInitValues()
{
	inlet_flow_angle = 0.0;
	initvar.u = 1.0;
	initvar.v = 0.0;
	initvar.w = 0.0;
	double speed = pow(pow(initvar.u, 2)+pow(initvar.v, 2)+pow(initvar.w, 2), 0.5);
	if (isinf(Ma)) MPI_Abort(MPI_COMM_WORLD, 240);
	initvar.T = pow(speed/Ma, 2)/(gama*Rg);
	initvar.roe = 1.0;
	Get_p_ideal_gas(initvar);
	Get_E(initvar);
#ifdef TURBULENCE
	initvar.viseddy = far_vis_ratio*viseddy;
#endif	
}

void BCValues::UserInit(Mesh & amesh, const int & ilevel, const int & bn)
{
	if (amesh.bc(ilevel, bn)[0] <= 0.5)
	{
		amesh.LevelBoxData(ilevel, bn) = facevar[0];
	}
	else
	{
		amesh.LevelBoxData(ilevel, bn) = facevar[1];
	}
	//Get_Cons_Vars(amesh.LevelBoxData(ilevel,bn), amesh.LevelBoxData(ilevel,bn).his[0]);
	//Get_Cons_Vars(amesh.LevelBoxData(ilevel,bn), amesh.LevelBoxData(ilevel,bn).his[1]);
	//Get_Cons_Vars(amesh.LevelBoxData(ilevel,bn), amesh.LevelBoxData(ilevel,bn).his[2]);
}

// void BCValues::SetBDFaceVars()
// {
// 	facevar[0].p = initvar.p*pow((1.0+0.2*pow(Ma,2)), gama/(gama-1));
// 	facevar[0].T = initvar.T*(1.0+0.2*pow(Ma,2));

// 	facevar[1].p = initvar.p-3.0;
// }

void BCValues::SetBDFaceVars()
{
	// facevar[0].p = 1.0;
	// facevar[0].roe = 1.0;
	// facevar[0].u = 0.0;
	// facevar[0].v = 0.0;
	// facevar[0].w = 0.0;
	// Get_T_ideal_gas(facevar[0]);
	// Get_E(facevar[0]);

	// facevar[0] = initvar;

	// facevar[1].p = 0.1;
	// facevar[1].roe = 0.125;
	// facevar[1].u = 0.0;
	// facevar[1].v = 0.0;
	// facevar[1].w = 0.0;
	// Get_T_ideal_gas(facevar[1]);
	// Get_E(facevar[1]);
	// facevar[0].roe = initvar.roe;
	facevar[0].T = kineticenergy(initvar)*initvar.roe/Cp + initvar.T;
	facevar[0].p = initvar.p*pow(1.0+0.2*pow(Ma,2), 1.4/0.4);
	//facevar[0] = initvar;
	facevar[1].p = initvar.p;
}

void BCValues::SetBackPressure(Mesh & amesh)
{
	backpressure = -0.04;
}
void Pltoutput::Design_User_File()
{
	userfilenum  = 2;
	userfilename.resize(userfilenum);
	userfilename[0] = new char[32];
	userfilename[1] = new char[32];
	sprintf(userfilename[0], "force_foil.dat");
	sprintf(userfilename[1], "force_cy.dat");
	Create_User_File();
}

void Pltoutput::Print_Var_as_Func_of_Time(vector<Body> & abody, NS_Solver & asolver)
{
	double F_ref = 0.5*u_ref*u_ref*roe_ref*(L_ref*L_ref*0.1);
	double M_ref = F_ref*L_ref;
	double massflow;
	m_mesh.ComptMassFlow(massflow);
	bool printflag = true;
	int bodynum = abody.size();
	if (printflag)
	{
		if (nrank == 0)
		{
			// printf("file name is %s\n", userfilename[0]);
			userfileptr[0] = fopen(userfilename[0], "a");
            fprintf(userfileptr[0], "%16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f\n",
                    t, abody[0].force[0]/F_ref, abody[0].force[1]/F_ref,
                    abody[0].moment[2]/M_ref,
                    massflow, abody[0].inlineosc[1],abody[0].rotosc[2]);
			fclose(userfileptr[0]);
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
}

void Mesh::ComptSectionAvePre(const int & level, double & p0, double & T0, int z0)
{
	int bps = m_level[level].m_box.ps();
	int bpe = m_level[level].m_box.pe();
	p0 = 0.0;
	int p0_num = 0;
	for (int i = bps; i < bpe; ++i)
	{
		if (m_level[level].m_geom[i].boxcenter[i] > 26.49)
		{
			++p0_num;
			p0 += m_level[level].m_data[i].p;
		}
	}
	MPI_Allreduce(MPI_IN_PLACE, &p0_num, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, &p0, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	p0 /= double(p0_num);
}

void Pltoutput::Print_Var_as_Func_of_Space(vector<Body> & abody)
{
	bool printflag = false;
	if (printflag)
	{
		for (int i = 0; i < nodenum; ++i)
		{
			if (srank == 0 && node == i)
			{
				userfileptr[0] = fopen(userfilename[0], "a");
				string dataformat = ffm + ffm + ffm + ffm + ffm + "\n";	
				for (int ln = 0; ln < m_mesh.cur_level_num; ++ln)
				{
					for (int bn = 0; bn < m_mesh.m_level[ln].m_box.size(); ++bn)
					{
						if (m_mesh.m_level[ln].m_box[bn].iy() == 0 && m_mesh.m_level[ln].m_box[bn].iz() == 0)
						{
							fprintf(userfileptr[0], 
									dataformat.c_str(), 
									m_mesh.bc(ln,bn)[0]*L_ref,
									m_mesh.m_level[ln].m_data[bn].p*p_ref,
									m_mesh.m_level[ln].m_data[bn].roe*roe_ref,
									m_mesh.m_level[ln].m_data[bn].u*u_ref,
									m_mesh.m_level[ln].m_data[bn].e*e_ref);
						}
					}				
				}
				fclose(userfileptr[0]);
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
	}
}


void AMR::SetRefineFlag()
{
	if (ts == 0)
	{
		m_level.back().torefine = true;
		init_adjust_level = 0;
	}
	// else
	// {
	// 	init_adjust_level = 0;
	// 	for (int i = 0; i < m_mesh.MyCurNum(); ++i)
	// 	{
	// 		m_level[i].torefine = true;
	// 		m_level[i].toderefine = true;
	// 	}
	// }
}

void AMR::TagtheBox(const int & ilevel, vector<Body> & abody)
	{
        //对于翼型问题一共用了8层网格，前四层0-3网格层设置为矩形，0就是背景笛卡尔网格
        //后四层4-7设置为环形
        //rs re分别给出了环形网格区域最内侧和最外层距离壁面的距离
        //例如对于第4层网格，最内侧和最外侧距离壁面的距离分别为-5.0*dh[3][0]和0.5
        //dh[3][0]就是第三层网格x方向的网格尺寸，由于这里采用的xy尺寸一致的网格，所以dh[3][0]=dh[3][1],用哪一个都行
        double re[7] = {0.0,0.0,0.0, 0.75, 0.4, 0.125,0.04};
        double rs[7] = {0.0,0.0,0.0, -5.0*dh[3][0], -5.0*dh[4][0], -5.0*dh[5][0],-3.0*dh[6][0]};
        //矩形网格层左侧距离固壁左侧的距离
        double up_length[4] = {4.5, 3.7, 3.0, 1.8};
        //矩形网格层上下距离固壁上下的距离
        double up_xrange[4] = {4.0,3.0,2.2,0.0};
        //矩形网格层右侧距离固壁右侧的距离
        double down_length[4] = {15.0, 10.0, 8.0, 12.0};

		double d0 = Body::clength()/2.0;
		/*------------------------------------------------------*/
		int tagnum = 0;
		int detagnum = 0;
		GiveAFlag("Start to initialize the tag!!!", 4);
		m_level[ilevel].Inittag();
		ShowAllRankData("tag size", m_level[ilevel].m_tag.size(), 4);
		//printf("Level %d refine flag is %d\n", ilevel, m_level[ilevel].torefine);
		if (m_level[ilevel].torefine)
		{
			int bs = m_level[ilevel].m_box.ps();
			int be = m_level[ilevel].m_box.pe();
			if (ilevel < 3)
			{
				for (int i = bs; i < be; ++i)
				{
					if (m_mesh.bc(ilevel,i)[0] - foil_left_x > -up_xrange[ilevel] &&
							m_mesh.bc(ilevel,i)[0] - foil_right_x < down_length[ilevel] &&
							abs(m_mesh.bc(ilevel,i)[1] - abody[0].bc()[1]) < up_length[ilevel])
					{
						m_level[ilevel].m_tag[i].givetag();
						++tagnum;
						m_level[ilevel].taggedbox.push_back(i);
					}
				}
			}
			else if (ilevel >= 3)
			{
				for (int i = bs; i < be; ++i)
				{
					double & celldis = m_level[ilevel].m_box[i].pair.signdis;
					if (celldis > rs[ilevel] && celldis < re[ilevel])
					{
						m_level[ilevel].m_tag[i].tag = 0;
						++tagnum;
						m_level[ilevel].taggedbox.push_back(i);
					}
					/*-----------------------------------------------------------------*/
				}
			}
		}
		if (m_level[ilevel].toderefine)
		{
			int bs = m_level[ilevel].m_box.ps();
			int be = m_level[ilevel].m_box.pe();
			if (ilevel < 4)
			{
				for (int i = bs; i < be; ++i)
				{
					if (!(m_mesh.bc(ilevel,i)[0] - foil_left_x > -up_length[ilevel-1] &&
							m_mesh.bc(ilevel,i)[0] - foil_right_x < down_length[ilevel-1] &&
							abs(m_mesh.bc(ilevel,i)[1] - abody[0].bc()[1]) < up_length[ilevel-1]))
					{
						m_level[ilevel].m_tag[i].givedetag();
						++detagnum;
					}
				}
			}
			else if (ilevel >= 4)
			{
				for (int i = bs; i < be; ++i)
				{
					double & celldis = m_level[ilevel].m_box[i].pair.signdis;
					if (!(celldis > rs[ilevel-1] && celldis < re[ilevel-1]))
					{
						m_level[ilevel].m_tag[i].givedetag();
						++detagnum;
					}
					/*-----------------------------------------------------------------*/
				}
			}
		}
		GiveAFlag("Finish initial tag!!!", 5);
		// if (m_level[ilevel].torefine && ts > 0)
		// {
		// 	m_level[ilevel].MoveBox_Refine(1, grid_move_dir,tagnum);
		// }
		// GiveAFlag("Finish Tag MoveBox_Refine!!!", 5);
		// if (m_level[ilevel].toderefine && ts > 0)
		// {
		// 	m_level[ilevel].MoveBox_Derefine(1, grid_move_dir, detagnum);
		// }
		MPI_Barrier(share_comm);
		MPI_Allreduce(MPI_IN_PLACE, &tagnum, 1, MPI_INT, MPI_SUM, share_comm);
		MPI_Allreduce(MPI_IN_PLACE, &detagnum, 1, MPI_INT, MPI_SUM, share_comm);
		if (srank == 0)
		{
			PRINTFinLEVEL("tag number is %d detagnum %d", ilevel, tagnum, detagnum);
		}		
		CheckDetagLocation(ilevel);
		GiveAFlag("Finish CheckTagLocation!!!", 5);
		if (morelevel && ilevel == m_mesh.MyCurNum()-2)
		{
			int bodynum = abody.size();
			int l0 = m_mesh.MyCurNum()-2;
			for (int i = 0; i < bodynum; ++i)
			{
				abody[i].CheckAttachBoxTag(m_level[l0].m_tag, l0, m_mesh);
			}
		}
		MPI_Barrier(share_comm);
	}
