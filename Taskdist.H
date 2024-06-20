// //
// //  Taskdist.H
// //  AMR
// //
// //  Created by wangzhuo on 2019/11/23.
// //  Copyright © 2019 wangzhuo. All rights reserved.
// //
// #ifndef Taskdist_h
// #define Taskdist_h
// #include "AMRLevel.H"
// #include "AMRmpi.H"

// using namespace std;

// class Taskdist
// {

// public:
//   Taskdist()
//   {
//   }

//   template<class T> 
//   void LoadBalance(vector<T> & mytask)
//   {
//     ParaInit();
//     rank_task_amount = mytask.size();
//     ShowAllRankData("rank_task_amount", rank_task_amount, 4);
//     int tot_out_task_amount;
//     MPI_Allreduce(&rank_task_amount, &tot_out_task_amount, 1, MPI_INT, MPI_SUM, share_comm);
        
//     int rank_task_ave = ceil(double(tot_out_task_amount)/double(sprocs));
        
//     vector<int> gap_amount(nprocs,0);
//     int task_needed = 0;

//     if (rank_task_amount > rank_task_ave)
//     {
//       out_task_num = rank_task_amount - rank_task_ave;
//     }
//     else if (rank_task_amount < rank_task_ave)
//     {
//       task_needed = rank_task_ave - rank_task_amount;
//     }
//     vector<int> ranklabel(out_task_num, srank); 
//     if (out_task_num == 0) {ranklabel.reserve(1);}

//     vector<int> out_task_num_procs(nprocs, 0);
//     int tot_out_task_num;
//     MPI_Reduce(&out_task_num, &tot_out_task_num, 1, MPI_INT, MPI_SUM, 0, share_comm);
//     MPI_Bcast(&tot_out_task_num, 1, MPI_INT, 0, share_comm);
//     MPI_Allgather(&out_task_num, 1, MPI_INT, &out_task_num_procs[0], 1, MPI_INT, share_comm);
//     GiveAFlag("Finish gather out_task_num!!!", 5);
//     //printf("gap_amount size is %d\n", (int)gap_amount.size());
//     MPI_Allgather(&task_needed, 1, MPI_INT, &gap_amount[0], 1, MPI_INT, share_comm);
//     GiveAFlag("Finish gather gap_amount!!!", 5);
//     ArrayProcsStart(out_task_num_procs, outtask_numprost);
//     //printf("all_task_rank size %d\n", (int)all_task_rank.size());
//     all_task_rank.resize(tot_out_task_num);
//     all_task_destid.resize(tot_out_task_amount);
//     all_task_rank.assign(tot_out_task_num, -1);
//     all_task_destid.assign(tot_out_task_num, -1);
//     MPI_Allgatherv(&ranklabel[0], out_task_num, MPI_INT,
//       &all_task_rank[0], &out_task_num_procs[0], &outtask_numprost[0],
//       MPI_INT, share_comm);
//     int aimtask = 0;
//     for (int i = 0; i < sprocs && aimtask < tot_out_task_num; ++i)
//     {
//       while (gap_amount[i] > 0 && aimtask < tot_out_task_num)
//       {
//         all_task_destid[aimtask] = i;
//         if (srank == i)
//         {
//           ++in_task_num;
//         }
//         gap_amount[i] -= 1;
//         ++aimtask;
//       }
//     }
//     int real_out_task_num = 0;
//     for (int i = 0; i < out_task_num; ++i)
//     {
//       if (all_task_destid[outtask_numprost[srank]+i] > -1)
//       {
//         real_out_task_num += 1;
//       }
//     }
//     out_task_num = real_out_task_num;
//     printf("[%d]my task number %d ave is %d in %d out %d\n", srank, rank_task_amount, rank_task_ave, in_task_num, out_task_num);
//     Transfer_Task(mytask); 
//   }

//   template<class T> 
//   void Transfer_Task(vector<T> & mytask)
//   {
//     //enlarge array!!!
//     if (in_task_num > 0)
//     {
//       int oldmboxsize = mytask.size();
//       mytask.resize(oldmboxsize+in_task_num);
//     }
//     vector<vector<T> > tobesent(sprocs);
//     for (int i = 0; i < out_task_num; ++i)
//     {
//       int taskindex = rank_task_amount-1-i;
//       int destrank = all_task_destid[outtask_numprost[srank]+i];
//       tobesent[destrank].push_back(mytask[taskindex]);
//     }
//     vector<T> recvfbox;
//     BcastNewInfo_Node(tobesent, recvfbox);
//     mytask.resize(rank_task_amount+in_task_num);
//     for (int i = 0; i < in_task_num; ++i)
//     {
//       mytask[i+rank_task_amount] = recvfbox[i];
//     }
//     if (out_task_num > 0)
//     {
//       mytask.resize(rank_task_amount - out_task_num);
//     }
//   }

//   void ParaInit()
//   {
//     outtask_numprost.resize(nprocs);
//     rank_task_amount = 0;
//     in_task_num = 0;
//     out_task_num = 0;
//   }

// private:

//   vector<int> all_task_rank;
//   vector<int> all_task_destid;
//   vector<int> outtask_numprost;
//   int rank_task_amount;
//   int in_task_num;
//   int out_task_num;
// };

// #endif
// /* Taskdist_h */
