#ifndef _AMRmpi_H_
#define _AMRmpi_H_

#include <vector>
#include <mpi.h>

#include "Point.H"
//#include "FlowVariable.H"

using namespace std;

extern int nrank;
extern int nprocs;

extern MPI_Comm share_comm;
extern int srank;
extern int sprocs;

extern MPI_Comm share_comm_mother;
extern int shm_rank;
extern int shm_procs;

extern MPI_Comm nodecomm;
extern int node;
extern int nodenum;

extern MPI_Datatype MPI_FV;
extern MPI_Datatype MPI_PAIRINFO;
extern vector<MPI_Datatype> MY_MPI_DATATYPE;

int CountNumber(vector<int> & vt, int t0);

void CountTotalNum(const vector<int> & becounted, int & totalnum);

void CountTotalNum(const int becounted[], int & totalnum);

void ArrayProcsStart(const vector<int> & array, vector<int> & arraystart, const int & a0 = 0);

void ArrayProcsStart(const int array[], int arraystart[], const int & a0 = 0);

int str2unstr(const Point& a, const Point& t);

void ShowTimestep();

void PrintTime(const double & st, const double & et, const string & tag);

void PrintTime(const double & st, const string & tag);
//
void MPI_CreateNode();

void MPI_CreatDataType();

void MPI_FreeDataType();

void ArrayOrder_s(const int & start0, const int & end0,
                int & rstart, int & rend,
                const int & partnum, const int & rank0);

// #define DEFINE_MPI_BOXLOC() \
// 				MPI_Datatype MPI_BOXLOC;\
// 				MPI_Type_contiguous(2, MPI_INT, &MPI_BOXLOC);\
// 				MPI_Type_commit(&MPI_BOXLOC);

// #define DEFINE_MPI_BOX()\
// 				MPI_Datatype MPI_BOX;\
// 				MPI_Type_contiguous(3, MPI_INT, &MPI_BOX);\
// 				MPI_Type_commit(&MPI_BOX);

#define DEFINE_MPI_TYPET() \
				MPI_Datatype MPI_TYPET;\
				MPI_Type_contiguous(sizeof(T)/sizeof(int), MPI_INT, &MPI_TYPET);\
				MPI_Type_commit(&MPI_TYPET)

#define FREE_MPI_TYPET() \
        MPI_Type_free(&MPI_TYPET);

// #define DEFINE_MPI_GHOST()\
//         MPI_Datatype MPI_GHOST;\
//         MPI_Type_contiguous(sizeof(GhostData)/sizeof(double), MPI_DOUBLE, &MPI_GHOST);\
//         MPI_Type_commit(&MPI_GHOST);           

// template <class T>
// void BcastNewInfo_NodeRank(vector<vector<T> >& datasend,
//      		vector<T>& datarecv,
//         vector<int> & recvrank)
// {
//   DEFINE_MPI_TYPET();
//   vector<int> numrecv(nodenum, 0);
//   vector<int> num2ranks(nodenum, 0);
//   vector<int> gnum2ranks(nodenum*nodenum);
//   for (int i = 0; i < nodenum; ++i)
//   {
//     int dss = datasend[i].size();
//   	num2ranks[i] = dss;
    
//     if (dss == 0)
//     {
//       datasend[i].reserve(1);
//     }
//   } 
//   MPI_Allgather(&num2ranks[0], nodenum, MPI_INT, &gnum2ranks[0], nodenum, MPI_INT, nodecomm);
//   for (int i = 0; i < nodenum; ++i)
//   {
//     numrecv[i] = gnum2ranks[nodenum*i+node];
//     // MPI_Scatter(&num2ranks[0], 1, MPI_INT, &numrecv[i], 1, MPI_INT,
//     //   i, nodecomm);
//   }
//   for (int i = 0; i < nodenum; ++i)
//   {
    
//     if (num2ranks[i] == 0)
//     {
//       ++num2ranks[i];
//     }
//   }
//   int tot_recv_num;
//   CountTotalNum(numrecv, tot_recv_num);
//   // MPI_Barrier(MPI_COMM_WORLD);
//   //printf("tot_recv_num is %d\n", tot_recv_num);
//   int have_recv_num = 0;
//   datarecv.resize(tot_recv_num);
//   recvrank.resize(tot_recv_num);
//   /*--------------------------------------------------------*/
//   MPI_Win nodewin[nodenum];
//   for (int i = 0; i < nodenum; ++i)
//   {   
//     MPI_Win_create(&datasend[i][0], num2ranks[i]*sizeof(T), sizeof(T), MPI_INFO_NULL, nodecomm, &nodewin[i]);
//     MPI_Win_lock_all(0, nodewin[i]);
//   }
//   for (int i = 0; i < nodenum; ++i)
//   {
//     if (numrecv[i] > 0)
//     {
//       MPI_Get(&datarecv[have_recv_num], numrecv[i], MPI_TYPET, i, 0, numrecv[i], MPI_TYPET, nodewin[node]);
//       for (int i0 = 0; i0 < numrecv[i]; ++i0)
//       {
//         recvrank[i0+have_recv_num] = i;
//       }
//       have_recv_num += numrecv[i];
//     }
//   }
//   for (int i = 0; i < nodenum; ++i)
//   {
//     MPI_Win_unlock_all(nodewin[i]);
//     MPI_Win_free(&nodewin[i]);
//   }
    // FREE_MPI_TYPET();
// };

// template <class T>
// void BcastNewInfo_Node(vector<vector<T> >& datasend,
//         vector<T>& datarecv)
// {
//   DEFINE_MPI_TYPET();
//   vector<int> numrecv(nodenum, 0);
//   vector<int> num2ranks(nodenum, 0);
//   vector<int> gnum2ranks(nodenum*nodenum);
//   for (int i = 0; i < nodenum; ++i)
//   {
//     int dss = datasend[i].size();
//     num2ranks[i] = dss;
//     if (dss == 0)
//     {
//       datasend[i].reserve(1);
//     }
//   } 
//   MPI_Allgather(&num2ranks[0], nodenum, MPI_INT, &gnum2ranks[0], nodenum, MPI_INT, nodecomm);
//   for (int i = 0; i < nodenum; ++i)
//   {
//     numrecv[i] = gnum2ranks[nodenum*i+node];
//     // MPI_Scatter(&num2ranks[0], 1, MPI_INT, &numrecv[i], 1, MPI_INT,
//     //   i, nodecomm);
//   }
//   for (int i = 0; i < nodenum; ++i)
//   {
    
//     if (num2ranks[i] == 0)
//     {
//       ++num2ranks[i];
//     }
//   }
//   int tot_recv_num;
//   CountTotalNum(numrecv, tot_recv_num);
//   // MPI_Barrier(MPI_COMM_WORLD);
//   //printf("tot_recv_num is %d\n", tot_recv_num);
//   int have_recv_num = 0;
//   datarecv.resize(tot_recv_num);
//   /*--------------------------------------------------------*/
//   MPI_Win nodewin[nodenum];
//   for (int i = 0; i < nodenum; ++i)
//   {   
//     MPI_Win_create(&datasend[i][0], num2ranks[i]*sizeof(T), sizeof(T), MPI_INFO_NULL, nodecomm, &nodewin[i]);
//     MPI_Win_lock_all(0, nodewin[i]);
//   }
//   for (int i = 0; i < nodenum; ++i)
//   {
//     if (numrecv[i] > 0)
//     {
//       MPI_Get(&datarecv[have_recv_num], numrecv[i], MPI_TYPET, i, 0, numrecv[i], MPI_TYPET, nodewin[node]);
//       have_recv_num += numrecv[i];
//     }
//   }
//   for (int i = 0; i < nodenum; ++i)
//   {
//     MPI_Win_unlock_all(nodewin[i]);
//     MPI_Win_free(&nodewin[i]);
//   }
        // FREE_MPI_TYPET();
// };

template <class T>
void BcastNewInfo_Node(vector<vector<T> > & datasend, vector<T> & datarecv, const int & tid)
{
  vector<int> numrecv(nodenum, 0);
  vector<int> num2ranks(nodenum, 0);
  for (int i = 0; i < nodenum; ++i)
  {
    num2ranks[i] = datasend[i].size();
  }
  for (int i = 0; i < nodenum; ++i)
  {
    MPI_Scatter(&num2ranks[0], 1, MPI_INT, &numrecv[i], 1, MPI_INT,
      i, nodecomm);
  }    
  int tot_recv_num;
  CountTotalNum(numrecv, tot_recv_num);
//    ShowAllRankData("tot_recv_num", tot_recv_num);
  int have_recv_num = 0;
  datarecv.resize(tot_recv_num);
  MPI_Request recv_req[nodenum];
  MPI_Status recv_stats[nodenum];
  for (int i = 0; i < nodenum; ++i)
  {
    if (numrecv[i] > 0)
    {
//        printf("<<<<<My rank is %d recvnum from rank %d is %d\n", nrank, i, numrecv[i]);
      MPI_Irecv(&datarecv[have_recv_num], numrecv[i], MY_MPI_DATATYPE[tid], i, i, nodecomm, &recv_req[i]);
      have_recv_num += numrecv[i];
//        printf("My rank is %d I received data from rank %d data number is %d\n", nrank, i, numrecv[i]);
    }
  }
  for (int i = 0; i < nodenum; ++i)
  {
    if (num2ranks[i] > 0)
    {
      MPI_Send(&datasend[i][0], num2ranks[i], MY_MPI_DATATYPE[tid], i, node, nodecomm);
    }
  }
  for (int i = 0; i < nodenum; ++i)
  {
    if (numrecv[i] > 0)
    {
      MPI_Wait(&recv_req[i], &recv_stats[i]);
    }
  }
  
};

template <class T>
void BcastNewInfo_NodeRank(vector<vector<T> >& datasend, vector<T>& datarecv, vector<int> & recvrank, const int & tid)
{
  vector<int> numrecv(nodenum, 0);
  vector<int> num2ranks(nodenum, 0);
  for (int i = 0; i < nodenum; ++i)
  {
    num2ranks[i] = datasend[i].size();
  }
  for (int i = 0; i < nodenum; ++i)
  {
    MPI_Scatter(&num2ranks[0], 1, MPI_INT, &numrecv[i], 1, MPI_INT,
      i, nodecomm);
  }    
  int tot_recv_num;
  CountTotalNum(numrecv, tot_recv_num);
//    ShowAllRankData("tot_recv_num", tot_recv_num);
  int have_recv_num = 0;
  datarecv.resize(tot_recv_num);
  recvrank.resize(tot_recv_num);
  MPI_Request recv_req[nodenum];
  MPI_Status recv_stats[nodenum];
  for (int i = 0; i < nodenum; ++i)
  {
    if (numrecv[i] > 0)
    {
//        printf("<<<<<My rank is %d recvnum from rank %d is %d\n", nrank, i, numrecv[i]);
      MPI_Irecv(&datarecv[have_recv_num], numrecv[i], MY_MPI_DATATYPE[tid], i, i, nodecomm, &recv_req[i]);
      for (int rcv0 = 0; rcv0 < numrecv[i]; ++rcv0)
      {
        recvrank[rcv0+have_recv_num] = i;
      }
      have_recv_num += numrecv[i];
//        printf("My rank is %d I recived data from rank %d data number is %d\n", nrank, i, numrecv[i]);
    }
  }
  for (int i = 0; i < nodenum; ++i)
  {
    if (num2ranks[i] > 0)
    {
      MPI_Send(&datasend[i][0], num2ranks[i], MY_MPI_DATATYPE[tid], i, node, nodecomm);
    }
  }
  for (int i = 0; i < nodenum; ++i)
  {
    if (numrecv[i] > 0)
    {
      MPI_Wait(&recv_req[i], &recv_stats[i]);
    }
  }
};

#endif