#ifndef _Nodedata_H_
#define _Nodedata_H_

#include <vector>

#include "DataArray.H"
#include "AMRmpi.H"
#include "Point.H"

using namespace std;

struct ExtractMesg
{
  int remotecell;
  int localcell;
  double rotangle;
  ExtractMesg()
  {}

  ExtractMesg(const int & r0, const int & l0):
  remotecell(r0),localcell(l0)
  {}

};

struct Direct_Mesg
{
  int fromid;
  int toid;
  double rotangle;

  Direct_Mesg()
  {}

  Direct_Mesg(const int & i1, const int & i2):
  fromid(i1),toid(i2)
  {}
};

struct ExtractMesg_xyz
{
  int remotecell;
  int localcell;
  Point xyz;
#ifdef PASSAGE_ANGLE
  double rotangle;
#endif    
#ifdef DEBUG  
  double mydis;
  Pointxyz bc;
#endif
  ExtractMesg_xyz()
  {}

#ifdef DEBUG
  ExtractMesg_xyz(const int & r0, const int & l0, Point & cellbc, double & d0, Pointxyz & ap0):
  remotecell(r0),localcell(l0),xyz(cellbc),mydis(d0), bc(ap0)
  {}
#else
  ExtractMesg_xyz(const int & r0, const int & l0, Point & cellbc):
  remotecell(r0),localcell(l0),xyz(cellbc)
  {}
#endif    

};

struct Pairinfo
{
  int body;
  int patch;

  Pairinfo()
  {}

  Pairinfo(int & b0, int & p0):
  body(b0),patch(p0)
  {}
};

#define DEFINE_MPI_PAIRINFO() \
        MPI_Type_contiguous(sizeof(Pairinfo)/sizeof(int), MPI_INT, &MPI_PAIRINFO);\
        MPI_Type_commit(&MPI_PAIRINFO)

class Nodedata
{
	friend class Meshlevel;
	friend class Mesh;
  friend class Body;
private:

	//vector<FlowVariables>  flowmesg;

  vector<vector<FlowVariables> > flowmesg_s;

  vector<vector<Pairinfo> > boxpair;

  vector<vector<ExtractMesg> > cell_extrac;

  DataArray<Direct_Mesg> direct_extrac;

  vector<DataArray<ExtractMesg> > global_cell_extrac;

  vector<vector<int> >            msgdest;

  vector<FlowVariables> recvfv;

  vector<Pairinfo> recvpair;

  vector<int> recvdestcell;

  vector<int> numrecv;

  //vector<int> getstart;

  vector<int> num2ranks;

  vector<vector<ExtractMesg_xyz> > ibothernodes;
  vector<ExtractMesg_xyz> recvibcells;
  vector<int> recvib_src_node;

public:

	Nodedata()
	{}

  void Creat_direct_extract_array(vector<Direct_Mesg> & local_dmg)
  {
    direct_extrac.setnum_nocopy(0,0);
    direct_extrac.Addnew(local_dmg);
    direct_extrac.DirectlyReduceNew();
  }

  /*For RMDA mode*/
	// void CreateWin()
	// {
 //  	numrecv.resize(nodenum);
 //  	num2ranks.resize(nodenum);
 //    getstart.resize(nodenum);
 //  	vector<int> gnum2ranks(nodenum*nodenum);
 //    int tot_send_num = 0;

 //    for (int i = 0; i < nodenum; ++i)
 //    {
 //      num2ranks[i] = cell_extrac[i].size();
 //      tot_send_num += num2ranks[i];
 //      //printf("Total extract number 1 %d 2 %d\n", tot_extract_num, m_level[nlevel].nd.d_cell_extrac[i].size());
 //    }

 //    flowmesg.resize(tot_send_num);
 //    msgdest.resize(tot_send_num);
 //    if (tot_send_num == 0)
 //    {
 //      flowmesg.reserve(1);
 //      msgdest.reserve(1);
 //    }
 //    int s0 = 0;
 //    for (int i = 0; i < nodenum; ++i)
 //    {
 //      for (int r0 = 0; r0 < num2ranks[i]; ++r0)
 //      {
 //        msgdest[s0+r0] = cell_extrac[i][r0].remotecell;
 //      }
 //      s0 += num2ranks[i];
 //    }
 //    //MPI_Barrier(MPI_COMM_WORLD);
 //    GiveAFlag("Start to create the data wins...",5);
 //  	MPI_Allgather(&num2ranks[0], nodenum, MPI_INT, &gnum2ranks[0], nodenum, MPI_INT, nodecomm);
 //  	for (int i = 0; i < nodenum; ++i)
 //  	{
 //  	  numrecv[i] = gnum2ranks[nodenum*i+node];
 //      getstart[i] = 0;
 //      for (int i0 = 0; i0 < node; ++i0)
 //      {
 //        getstart[i] += gnum2ranks[nodenum*i+i0];
 //      }
 //  	  // MPI_Scatter(&num2ranks[0], 1, MPI_INT, &numrecv[i], 1, MPI_INT,
 //  	  //   i, nodecomm);
 //  	}
 //  	int tot_recv_num; 
 //  	CountTotalNum(numrecv, tot_recv_num);
 //  	// MPI_Barrier(MPI_COMM_WORLD);
 //  	//printf("tot_recv_num is %d\n", tot_recv_num);
 //  	int have_recv_num = 0;
 //  	recvfv.resize(tot_recv_num);
 //  	recvdestcell.resize(tot_recv_num);
 //  	/*--------------------------------------------------------*/	
  	
 //  	MPI_Win destwin;   
 //  	MPI_Win_create(&msgdest[0], tot_send_num*sizeof(int), sizeof(int), 
 //  	  MPI_INFO_NULL, nodecomm, &destwin);
 //  	MPI_Win_lock_all(0, destwin);
 //  	for (int i = 0; i < nodenum; ++i)
 //  	{
 //  	  if (numrecv[i] > 0)
 //  	  {
 //  	    MPI_Get(&recvdestcell[have_recv_num], numrecv[i], MPI_INT, i, getstart[i], numrecv[i], MPI_INT, destwin);
 //  	    have_recv_num += numrecv[i];
 //  	  }
 //  	}
 //  	MPI_Win_unlock_all(destwin);
 //    MPI_Win_free(&destwin);
	// }
  /*For RECV and SEND mode*/

  void init_recvarray()
  {
    numrecv.resize(nodenum);
    num2ranks.resize(nodenum);
    for (int i = 0; i < nodenum; ++i)
    {
      numrecv[i] = 0;
      num2ranks[i] = 0;
    }
  }
	void CreateWin()
	{
  	numrecv.resize(nodenum);
  	num2ranks.resize(nodenum);
    flowmesg_s.resize(nodenum);
    boxpair.resize(nodenum);
    msgdest.resize(nodenum);
    if (global_cell_extrac.size() != nodenum)
    {
      global_cell_extrac.resize(nodenum);
      for (int i = 0; i < nodenum; ++i)
      {
        global_cell_extrac[i] = DataArray<ExtractMesg>();
      }
    }
    for (int i = 0; i < nodenum; ++i)
    {
      global_cell_extrac[i].setnum_nocopy(0,0);
      global_cell_extrac[i].Addnew(cell_extrac[i]);
      global_cell_extrac[i].DirectlyReduceNew();
      //printf("global_cell_extrac ghost size is %d!\n", (int)global_cell_extrac[i].ghost_index.size());
    }
    //GiveAFlag("Finish allocate win data array!!!", 5);
    for (int i = 0; i < nodenum; ++i)
    {
      if (i != node) num2ranks[i] = global_cell_extrac[i].pe()-global_cell_extrac[i].ps();
      else num2ranks[i] = 0;
      flowmesg_s[i].resize(num2ranks[i]);
      boxpair[i].resize(num2ranks[i]);
      msgdest[i].resize(num2ranks[i]);
    }
    //GiveAFlag("Finish num2ranks!!!", 5);
    for (int i = 0; i < nodenum; ++i)
    {
      if (i != node)
      {
        int rps = global_cell_extrac[i].ps();
        int rpe = global_cell_extrac[i].pe();
        for (int r0 = rps; r0 < rpe; ++r0)
        {
          int i0 = r0 - rps;
          msgdest[i][i0] = global_cell_extrac[i][r0].remotecell;
        }
      }
    }
    //GiveAFlag("Finish msgdest!!!", 5);
    for (int i = 0; i < nodenum; ++i)
    {
      MPI_Scatter(&num2ranks[0], 1, MPI_INT, &numrecv[i], 1, MPI_INT,
        i, nodecomm);
    }
    //GiveAFlag("Finish scatter the data number to each node!!!", 5);
    int tot_recv_num;
    CountTotalNum(numrecv, tot_recv_num);
  	int have_recv_num = 0;
  	recvfv.resize(tot_recv_num);
    recvpair.resize(tot_recv_num);
  	recvdestcell.resize(tot_recv_num);
    MPI_Request recv_req[nodenum];
    MPI_Status recv_stats[nodenum];
    for (int i = 0; i < nodenum; ++i)
    {
      if (numrecv[i] > 0)
      {
//        printf("<<<<<My rank is %d recvnum from rank %d is %d\n", nrank, i, numrecv[i]);
        MPI_Irecv(&recvdestcell[have_recv_num], numrecv[i], MPI_INT, i, i, nodecomm, &recv_req[i]);
        have_recv_num += numrecv[i];
//        printf("My rank is %d I received data from rank %d data number is %d\n", nrank, i, numrecv[i]);
      }
    }
    for (int i = 0; i < nodenum; ++i)
    {
      if (num2ranks[i] > 0)
      {
        MPI_Send(&msgdest[i][0], num2ranks[i], MPI_INT, i, node, nodecomm);
      }
    }
    for (int i = 0; i < nodenum; ++i)
    {
      if (numrecv[i] > 0)
      {
        MPI_Wait(&recv_req[i], &recv_stats[i]);
      }
    }
	}
//   void CreateWin()
//   {
//     numrecv.resize(nprocs,0);
//     num2ranks.resize(nprocs,0);
//     flowmesg_s.resize(nprocs);
//     msgdest.resize(nprocs);
//     for (int i = 0; i < nodenum; ++i)
//     {
//       int aimprocs = i*sprocs+srank;
//       if (i != node) num2ranks[aimprocs] = cell_extrac[i].size();
//       else num2ranks[aimprocs] = 0;
//       flowmesg_s[aimprocs].resize(num2ranks[aimprocs]);
//       msgdest[aimprocs].resize(num2ranks[aimprocs]);
//     }
//     for (int i = 0; i < nodenum; ++i)
//     {
//       int aimprocs = i*sprocs+srank;
//       for (int r0 = 0; r0 < num2ranks[aimprocs]; ++r0)
//       {
//         msgdest[aimprocs][r0] = cell_extrac[i][r0].remotecell;
//       }
//     }
//     for (int i = 0; i < nprocs; ++i)
//     {
//       MPI_Scatter(&num2ranks[0], 1, MPI_INT, &numrecv[i], 1, MPI_INT,
//         i, MPI_COMM_WORLD);
//     }
//     int tot_recv_num;
//     CountTotalNum(numrecv, tot_recv_num);
//     int have_recv_num = 0;
//     recvfv.resize(tot_recv_num);
//     recvdestcell.resize(tot_recv_num);
//     MPI_Request recv_req[nodenum];
//     MPI_Status recv_stats[nodenum];
//     for (int i = 0; i < nodenum; ++i)
//     {
//       int aimprocs = i*sprocs+srank;
//       if (numrecv[aimprocs] > 0)
//       {
// //        printf("<<<<<My rank is %d recvnum from rank %d is %d\n", nrank, i, numrecv[i]);
//         MPI_Irecv(&recvdestcell[have_recv_num], numrecv[aimprocs], MPI_INT, aimprocs, i, MPI_COMM_WORLD, &recv_req[i]);
//         have_recv_num += numrecv[aimprocs];
// //        printf("My rank is %d I received data from rank %d data number is %d\n", nrank, i, numrecv[i]);
//       }
//     }
//     for (int i = 0; i < nodenum; ++i)
//     {
//       int aimprocs = i*sprocs+srank;
//       if (num2ranks[aimprocs] > 0)
//       {
//         MPI_Send(&msgdest[aimprocs][0], num2ranks[aimprocs], MPI_INT, aimprocs, node, MPI_COMM_WORLD);
//       }
//     }
//     for (int i = 0; i < nodenum; ++i)
//     {
//       int aimprocs = i*sprocs+srank;
//       if (numrecv[aimprocs] > 0)
//       {
//         MPI_Wait(&recv_req[i], &recv_stats[i]);
//       }
//     }
//   } 
};

#endif