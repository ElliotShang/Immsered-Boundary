
#include <vector>
#include <mpi.h>

#include "AMRmpi.H"
#include "Point.H"
#include "AMRSpaceTime.H"
#include "FlowVariables.H"
#include "Constants.H"
#include "Mesh.H"
#include "AMRDatatype.H"

using namespace std;

int nrank, nprocs;

MPI_Comm share_comm;
int srank, sprocs;

MPI_Comm share_comm_mother;
int shm_rank;
int shm_procs;

MPI_Comm nodecomm;
int node;
int nodenum;

double total_vs_time = 0.0;
double total_inv_time = 0.0;
double total_ib_time = 0.0;
double total_force_time = 0.0;

double step_inv_time;
double step_vs_time;
double step_ib_time;
double step_total_time;
double step_tur_time;
double step_interface_time;
double step_exchange_time;
double step_move_time;

double t = 0.0;
int ts = 0;
vector<vector<double> > dh;
bool periodic[3];
double CFL_dt;
double dmlength[3];
double critic_ylength;

MPI_Datatype MPI_FV;
MPI_Datatype MPI_PAIRINFO;

vector<int> marching_step;
vector<int> average_step;
vector<int> prolongation_step;

vector<int> marching_left_step;
vector<int> average_left_step;
vector<int> prolongation_left_step;
vector<double> step_num;
vector<int> level_refine_ratio;
vector<Point> level_grid_ratio;

vector<MPI_Datatype> MY_MPI_DATATYPE;

int type_num = 11;

void MPI_CreateNode()
{
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &nrank);
  MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &share_comm_mother);
  MPI_Comm_rank(share_comm_mother, &shm_rank);
  MPI_Comm_size(share_comm_mother, &shm_procs);
  nodenum = nprocs/SHN;
  node = shm_rank/SHN;
  MPI_Comm_split(share_comm_mother, node, shm_rank, &share_comm);
  MPI_Comm_rank(share_comm, &srank);
  MPI_Comm_size(share_comm, &sprocs);
  MPI_Comm_split(MPI_COMM_WORLD, srank, nrank, &nodecomm);
  MPI_Comm_rank(nodecomm, &node);
  MPI_Barrier(MPI_COMM_WORLD);
  int node0 = node;
  MPI_Bcast(&node, 1, MPI_INT, 0, share_comm);
  if (node0 != node)
  {
  	printf("Rank index error!!!\n");
  	MPI_Abort(MPI_COMM_WORLD, 61);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  //printf("Node %d node number is %d My Share rank is %d Share comm size is %d\n", node, nodenum, srank, sprocs);
  DEFINE_MPI_FV();
  DEFINE_MPI_PAIRINFO();
  MPI_CreatDataType();
#ifdef TURBULENCE
#ifdef TWO_SOLID_GHOSTS
  if (nrank == 0)
  {
  	printf("For the turbulence simulation, TWO_SOLID_GHOSTS Flags should not be defined!!!\n");
  	MPI_Abort(MPI_COMM_WORLD, 91);
  }
#endif
#endif  
	 
}

void MPI_CreatDataType()
{
	MY_MPI_DATATYPE.resize(type_num);
	MPI_Type_contiguous(sizeof(int)/sizeof(int), MPI_INT, &MY_MPI_DATATYPE[0]);
	MPI_Type_contiguous(2*sizeof(int)/sizeof(int), MPI_INT, &MY_MPI_DATATYPE[1]);
	MPI_Type_contiguous(sizeof(S_newtag)/sizeof(int), MPI_INT, &MY_MPI_DATATYPE[2]);
	MPI_Type_contiguous(sizeof(S_tag)/sizeof(int), MPI_INT, &MY_MPI_DATATYPE[3]);
	MPI_Type_contiguous(sizeof(S_RemoteBKPInfo)/sizeof(int), MPI_INT, &MY_MPI_DATATYPE[4]);
	MPI_Type_contiguous(sizeof(S_newfinebp0)/sizeof(int), MPI_INT, &MY_MPI_DATATYPE[5]);
	MPI_Type_contiguous(sizeof(S_newcoarsebp)/sizeof(int), MPI_INT, &MY_MPI_DATATYPE[6]);
	MPI_Type_contiguous(sizeof(ExtractMesg_xyz)/sizeof(int), MPI_INT, &MY_MPI_DATATYPE[7]);
	MPI_Type_contiguous(sizeof(Block_dis)/sizeof(int), MPI_INT, &MY_MPI_DATATYPE[8]);
	MPI_Type_contiguous(sizeof(Pointxyz)/sizeof(int), MPI_INT, &MY_MPI_DATATYPE[9]);
	MPI_Type_contiguous(sizeof(BoxIndexBC)/sizeof(int), MPI_INT, &MY_MPI_DATATYPE[10]);
	for (int i = 0; i < type_num; ++i)
	{
		MPI_Type_commit(&MY_MPI_DATATYPE[i]);
	}
}

void MPI_FreeDataType()
{
	for (int i = 0; i < type_num; ++i)
	{
		MPI_Type_free(&MY_MPI_DATATYPE[i]);
	}
}
void ComptLevelgridsize(double * dh0)
{
	dh.resize(max_mesh_level);
	GiveAFlag("Start set level 2d flag!!!", 5);
	SetTwodFlag();
	SetPower_Ratio();
	SetLevelOperationIndex();
	GiveAFlag("Finish set level 2d flag!!!", 5);
	for (int i = 0; i < max_mesh_level; ++i)
	{
		dh[i].resize(4);
		dh[i][3] = 0.0;
		dh[i][0] = dh0[0]/double(pow(2, i));
		dh[i][1] = dh0[1]/double(pow(2, i));
		dh[i][2] = dh0[2]/level_power_ratio[i];
#if DIM == 3
		if (!level_twod_flag[i])
		{
			dh[i][3] = (dh[i][1]+dh[i][2]+dh[i][0])/3.0;
		}
		else
		{
			dh[i][3] = (dh[i][1]+dh[i][0])/2.0;
		}
#elif DIM == 2
		dh[i][3] = (dh[i][1]+dh[i][0])/2.0;
#elif DIM == 1
		dh[i][3] = dh[i][0];
#endif
	}
}

void SetPower_Ratio()
{
	level_power_ratio.resize(max_mesh_level);
	level_refine_ratio.resize(max_mesh_level);
	level_grid_ratio.resize(max_mesh_level);
	level_power_ratio[0] = 1;
	level_refine_ratio[0] = 1;
	level_grid_ratio[0] = Point(1,1,1);
	for (int i = 1; i < max_mesh_level; ++i)
	{
		level_refine_ratio[i] = pow(2, i);
		if (level_twod_flag[i]) level_power_ratio[i] = level_power_ratio[i-1];
		else level_power_ratio[i] = level_power_ratio[i-1]*2;

		level_grid_ratio[i][0] = pow(2,i);
		level_grid_ratio[i][1] = pow(2,i);
		level_grid_ratio[i][2] = level_power_ratio[i];
	}
}

void SetLevelOperationIndex()
{
	marching_step.resize(max_mesh_level);
	average_step.resize(max_mesh_level);
	prolongation_step.resize(max_mesh_level);
	step_num.resize(max_mesh_level);
	marching_left_step.resize(max_mesh_level);
	average_left_step.resize(max_mesh_level);
	prolongation_left_step.resize(max_mesh_level);
	for (int i = 0; i < max_mesh_level; ++i)
	{
		marching_step[i] = pow(2, max_mesh_level-1-i);
		average_step[i] = pow(2, max_mesh_level-i);
		prolongation_step[i] = pow(2, max_mesh_level-i);
		marching_left_step[i] = 0;
		average_left_step[i] = pow(2, max_mesh_level-i-1);
		prolongation_left_step[i] = pow(2, max_mesh_level-i-1);
#ifdef TEMPORAL_REFINE		
		step_num[i] = double(marching_step[i]);
#else
		step_num[i] = 1.0;
#endif				
	}
#ifdef TEMPORAL_REFINE
	dt_var_step = marching_step[0]*max(int(double(dt_var_step)/double(marching_step[0])), 1);
#endif
}
int CountNumber(vector<int> & vt, int t0)
{
	int tn = 0;
	for (int i = 0; i < vt.size(); ++i)
	{
		if (vt[i] == t0)
		{
			tn += 1;
		}
	}
	return tn;
}

void PrintTime(const double & st, const double & et, const string & tag)
{
	if (nrank == 0)
	{
		printf("%s Time: %f seconds\n", tag.c_str(), et-st);
	}
}

void PrintTime(const double & st, const string & tag)
{
	if (nrank == 0)
	{
		printf("%s Time: %f seconds\n", tag.c_str(), st);
	}
}

void ShowTimestep()
{
	if (nrank == 0)
	{
		printf("The present time step is %d\n", ts);
	}
}

void CountTotalNum(const vector<int> & becounted, int & totalnum)
{
	totalnum = 0;
	for (int i = 0; i < becounted.size(); ++i)
	{
		if (becounted[i] < 0)
		{
			printf("The number be counted must be positive!!!\n");
			MPI_Abort(MPI_COMM_WORLD, 77);
		}
		totalnum += becounted[i];
	}
}
void CountTotalNum(const int becounted[], int & totalnum)
{
	totalnum = 0;
	for (int i = 0; i < nprocs; ++i)
	{
		if (becounted[i] < 0)
		{
			printf("The number be counted must be positive!!!\n");
			MPI_Abort(MPI_COMM_WORLD, 77);
		}
		totalnum += becounted[i];
	}
}
void ArrayProcsStart(const vector<int> & array, vector<int> & arraystart, const int & a0)
{
	arraystart[0] = a0;
	for (int i = 1; i < array.size(); ++i)
	{
		arraystart[i] = arraystart[i-1] + array[i-1];
	}
}

void ArrayProcsStart(const int array[], int arraystart[], const int & a0)
{
	arraystart[0] = a0;
	for (int i = 1; i < nprocs; ++i)
	{
		arraystart[i] = arraystart[i-1] + array[i-1];
	}
}

void PeriodicLength(Pointxyz & dxyz)
{
	if (abs(dxyz[1]) > critic_ylength)
	{
		if (periodic[1])
		{
			if (dxyz[1] > critic_ylength)
			{
				dxyz[1] -= dmlength[1];
			}
			else if (dxyz[1] < -critic_ylength)
			{
				dxyz[1] += dmlength[1];
			}
		}
	}
}

double ComptPointAngle_Rotate_X(Pointxyz & p1, Pointxyz & p2)
{
  	double angle1 = atan2(p1[2], p1[1]);
		double angle2 = atan2(p2[2], p2[1]);
		return (angle2 - angle1);
}

#ifdef PASSAGE_ANGLE
void PeriodicAnnulaLength(Pointxyz & p1, Pointxyz & p2, Pointxyz & p1_new)
{
	double angle1 = atan2(p1[2], p1[1]);
	double angle2 = atan2(p2[2], p2[1]);
	double dangle = angle2 - angle1;
	p1_new = p1;
	double ang_dif = abs(dangle);
	if (ang_dif > 0.5*PASSAGE_ANGLE)
	{
		double rev_ang = PASSAGE_ANGLE;
		if (ang_dif > 1.5*PASSAGE_ANGLE)
		{
			rev_ang *= 2.0;
		}
		if (angle2 > angle1)
		{
			p1_new[1] = p1[1]*cos(rev_ang) - p1[2]*sin(rev_ang);
			p1_new[2] = p1[2]*cos(rev_ang) + p1[1]*sin(rev_ang);
		}
		else
		{
			p1_new[1] = p1[1]*cos(-rev_ang) - p1[2]*sin(-rev_ang);
			p1_new[2] = p1[2]*cos(-rev_ang) + p1[1]*sin(-rev_ang);
		}
	}
}
#endif
void ArrayOrder_s(const int & start0, const int & end0,
								int & rstart, int & rend,
								const int & partnum, const int & rank0)
{
	int num_in_procs = 0;
	int left_num = 0;
	int totnum = end0 - start0;
	num_in_procs = floor((double)totnum/(double)partnum);
	left_num = totnum - num_in_procs*partnum;
	//printf("num_in_procs %d left_num %d totnum %d\n", num_in_procs, left_num, totnum);
	if (rank0 < left_num)
	{
		rstart = (num_in_procs+1)*rank0+start0;
		rend = (num_in_procs+1)*(rank0+1)+start0;
	}
	else
	{
		rstart = rank0*num_in_procs+left_num+start0;
		rend = (rank0+1)*num_in_procs+left_num+start0;
	}
}

int str2unstr(const Point& a, const Point& t)
{
	int unstr = (a.ix()*t.iy()+a.iy())*t.iz()+a.iz();
	return unstr;
};

