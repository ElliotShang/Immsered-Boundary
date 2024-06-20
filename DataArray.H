#ifndef _DataArray_H_
#define _DataArray_H_

#include <iostream>
#include "AMRmpi.H"
#include "Box.H"
#include "Check.H"

using namespace std;

const double array_exp_ratio = 1.2;

template <class T>
class DataArray
{
private:

	MPI_Win win_fv;
	bool win_alloc;
	T * local_ptr;
	T * global_ptr;
	MPI_Aint bufsize;
	int disp_unit;

	int total_num;
	int num_without_ghost;
	int cap;
	int procs_s;
	int procs_e;
	// vector<int> new_normal_index;

	int totholenum;

	int ghostnum;
	int g_procs_s;
	int g_procs_e;

	int new_cell_num;
	int new_start_index;
	// int n_procs_s;
	// int n_procs_e;

	void ArrayOrder()
	{
		int num_in_procs = 0;
		int left_num = 0;
		num_in_procs = floor((double)num_without_ghost/(double)sprocs);
		left_num = num_without_ghost - num_in_procs*sprocs;
		if (srank < left_num)
		{
			procs_s = (num_in_procs+1)*srank;
			procs_e = (num_in_procs+1)*(srank+1);
		}
		else
		{
			procs_s = srank*num_in_procs+left_num;
			procs_e = (srank+1)*num_in_procs+left_num;
		}
	}

	void GhostOrder()
	{
		int num_in_procs = 0;
		int left_num = 0;
		num_in_procs = floor((double)ghostnum/(double)sprocs);
		left_num = ghostnum - num_in_procs*sprocs;
		if (srank < left_num)
		{
			g_procs_s = num_without_ghost+(num_in_procs+1)*srank;
			g_procs_e = num_without_ghost+(num_in_procs+1)*(srank+1);
		}
		else
		{
			g_procs_s = num_without_ghost+srank*num_in_procs+left_num;
			g_procs_e = num_without_ghost+(srank+1)*num_in_procs+left_num;
		}
	}

	void NumberEachProcs(int & num_rank, const int & total_num_all_procs)
	{
		num_rank = floor((double)total_num_all_procs/(double)sprocs);
		int left_num = total_num_all_procs - num_rank*sprocs;
		if (srank < left_num)
		{
			num_rank += 1;
		}
	}

	

public:

	DataArray()
	{
		total_num = 0;
		cap = 0;
		procs_e = 0;
		procs_s = 0;

		num_without_ghost = 0;
		ghostnum = 0;
		g_procs_s = 0;
		g_procs_e = 0;

		new_start_index = 0;
		new_cell_num = 0;

		MPI_Info a_info;
		MPI_Info_create(&a_info);
		MPI_Info_set(a_info, "alloc_shared_noncontig", "true");

		MPI_Win_allocate_shared(0, sizeof(T), a_info, share_comm, &local_ptr, &win_fv);
		MPI_Win_shared_query(win_fv, 0, &bufsize, &disp_unit, &global_ptr);
		MPI_Info_free(&a_info);
		win_alloc = true;
		// n_procs_s = 0;
		// n_procs_e = 0;
		//MPI_Barrier(MPI_COMM_WORLD);
	}

	DataArray(const vector<T> & a_data, const int & ghostnum0)
	{
		total_num = a_data.size();
		ghostnum = ghostnum0;
		num_without_ghost = total_num - ghostnum;
		new_cell_num = 0;
		new_start_index = total_num;
		cap = 0;
		//cap = a_data.size();

		procs_e = 0;
		procs_s = 0;

		g_procs_s = 0;
		g_procs_e = 0;

		MPI_Info a_info;
		MPI_Info_create(&a_info);
		MPI_Info_set(a_info, "alloc_shared_noncontig", "true");
		MPI_Win_allocate_shared(0, sizeof(T), a_info, share_comm, &local_ptr, &win_fv);
		MPI_Win_shared_query(win_fv, 0, &bufsize, &disp_unit, &global_ptr);
		MPI_Info_free(&a_info);
		Adjustarray(total_num);
		ArrayOrder();
		GhostOrder();
		//MPI_Win_fence(0, win_fv);
		for (int i = procs_s; i < procs_e; ++i)
		{
			global_ptr[i] = a_data[i];
		}
		for (int i = g_procs_s; i < g_procs_e; ++i)
		{
			global_ptr[i] = a_data[i];
			ghost_index.push_back(i);
		}
		//MPI_Win_fence(0, win_fv);
		if (nrank == 0)
		{
			printf("A DataArray was allocated with input!!!\n");
		}
	}

	// ~DataArray()
	// {
	// 	if (total_num > 0)
	// 	{
	// 		MPI_Win_free(&win_fv);
	// 		MPI_Barrier(MPI_COMM_WORLD);
	// 		if (nrank == 0) printf("The array total size %d has been freed!!!\n", total_num);
	// 	}
	// }

	vector<int> ghost_index;
	vector<BoxlocChange> newloctag;
	vector<int> hole;

	void NewOrder(const int & startbox, const int & endbox)
	{
		int num_in_procs = 0;
		int left_num = 0;
		num_in_procs = floor((double)new_cell_num/(double)sprocs);
		left_num = new_cell_num - num_in_procs*sprocs;
		if (srank < left_num)
		{
			startbox = new_start_index+(num_in_procs+1)*srank;
			endbox = new_start_index+(num_in_procs+1)*(srank+1);
		}
		else
		{
			startbox = new_start_index+srank*num_in_procs+left_num;
			endbox = new_start_index+(srank+1)*num_in_procs+left_num;
		}
	}

	inline const int & ps() const
	{
		return procs_s;
	}
	inline const int & pe() const
	{
		return procs_e;
	}

	inline const int & gps() const
	{
		return g_procs_s;
	}

	inline const int & gpe() const
	{
		return g_procs_e;
	}

	inline MPI_Win & arraywin()
	{
		return win_fv;
	}

	inline int & newstart()
	{
		return new_start_index;
	}

	inline int & newend()
	{
		return total_num;
	}

	inline int & holeindex(const int & i0)
	{
		return hole[i0];
	}

	inline int & arrayghostnum()
	{
		return ghostnum;
	}

	inline void givehole(const int & i0)
	{
		hole.push_back(i0);
	}

	inline vector<int> & holevector()
	{
		return hole;
	}

	void GlobalOrder(int & global_start, int & global_end)
	{
		int num_in_procs = 0;
		int left_num = 0;
		num_in_procs = floor((double)total_num/(double)sprocs);
		left_num = total_num - num_in_procs*sprocs;
		if (srank < left_num)
		{
			global_start = (num_in_procs+1)*srank;
			global_end = (num_in_procs+1)*(srank+1);
		}
		else
		{
			global_start = srank*num_in_procs+left_num;
			global_end = (srank+1)*num_in_procs+left_num;
		}
	}

	inline bool isghost(const int & i0)
	{
		if (i0 > num_without_ghost-1 && i0 < num_without_ghost+ghostnum) return true;
		else return false;		
	}

	inline bool isnorg(const int & i0)
	{
		if (i0 > -1 && i0 < num_without_ghost + ghostnum) return true;
		else return false;
	}

	inline bool isnormal(const int & i0)
	{
// #ifdef DEBUG
// 		if (i0 < num_without_ghost && i0 > -1)
// 		{
// 			return true;
// 		}
// 		else
// 		{
// 			return false;
// 		}
// #endif
// #ifndef DEBUG
		if (i0 < num_without_ghost && i0 > -1) return true;
		else return false;
// #endif
	}

	inline bool isnew(const int & i0)
	{
		if (i0 >= new_start_index && i0 < total_num) return true;
		else return false;
	}

	inline bool outrange(const int & i0)
	{
		if (i0 > total_num-1 || i0 < 0) return true;	
		else return false;
	}

	inline bool hasnewcell()
	{
		if (new_cell_num > 0 && total_num > new_start_index)
		{
			return true;
		}
		else
		{
			printf("New cell number is %d new_start_index is %d and total array size is %d\n", 
				new_cell_num, new_start_index, total_num);
			return false;
		}
	}


	inline void assign(const T & t0)
	{
		for (int i = procs_s; i < procs_e; ++i)
		{
			global_ptr[i] = t0;
		}
		for (int i = g_procs_s; i < g_procs_e; ++i)
		{
			global_ptr[i] = t0;
		}
	}
	inline void assign(const T & t0, const int & i0)
	{
		global_ptr[i0] = t0;
	}
	inline int Holenum()
	{
		return hole.size();
	}

	inline const int & size() const
	{
		return num_without_ghost;
	}

	inline const int & realsize() const
	{
		return total_num;
	}

	inline const int & capacity() const
	{
		return cap;
	}

	inline void replace(const int & i0)
	{
		global_ptr[newloctag[i0].newloc] = global_ptr[newloctag[i0].oldloc];
	}

	inline void switchtwo(const int & loctag0)
	{
		T a0 = global_ptr[newloctag[loctag0].newloc];
		global_ptr[newloctag[loctag0].newloc] = global_ptr[newloctag[loctag0].oldloc];
		global_ptr[newloctag[loctag0].oldloc] = a0;
	}

	inline void switcholdnew(const int & oldloc, const int & newloc)
	{
		T a0 = global_ptr[newloc];
		global_ptr[newloc] = global_ptr[oldloc];
		global_ptr[oldloc] = a0;
	}

	inline vector<BoxlocChange> & locchange()
	{
		return newloctag;
	}

	inline vector<int> & arrayghost()
	{
		return ghost_index;
	}

	// void Adjustarray(const int & new_cap)
	// {
	// 	MPI_Win new_win_fv;
	// 	MPI_Info a_info;
	// 	MPI_Info_create(&a_info);
	// 	MPI_Info_set(a_info, "alloc_shared_noncontig", "true");
	// 	T * new_lptr; T * new_gptr;
	// 	int my_buf_num = 0;
	// 	NumberEachProcs(my_buf_num, new_cap);
	// 	MPI_Win_allocate_shared(my_buf_num*sizeof(T), sizeof(T), a_info, share_comm, &new_lptr, &new_win_fv);
	// 	MPI_Win_shared_query(new_win_fv, 0, &bufsize, &disp_unit, &new_gptr);
	// 	MPI_Info_free(&a_info);
	// 	//MPI_Win_fence(0, new_win_fv);
	// 	/*copy the data part*/
	// 	copydatatonewarray(new_gptr, new_cap);
	// 	if (cap > 0)
	// 	{
	// 		MPI_Win_free(&win_fv);
	// 	}
	// 	//MPI_Win_fence(0, new_win_fv);
	// 	win_fv = new_win_fv;
	// 	local_ptr = new_lptr;
	// 	global_ptr = new_gptr;
	// 	cap = new_cap;
	// 	MPI_Barrier(share_comm);			
	// }

	void Adjustarray(const int & new_cap)
	{
		//if (ts > 0 && nrank == 0) {printf("The array should not be adjusted!!!\n"); MPI_Abort(MPI_COMM_WORLD, 417);}
		int exp_cap = int(array_exp_ratio*double(new_cap));
		int bs, be;
		ArrayOrder_s(0, total_num, bs, be, sprocs, srank);
		vector<T> tempdata(be-bs);
		copydatatonewarray(tempdata, bs, be, exp_cap);
		MPI_Barrier(share_comm);
		MPI_Win_free(&win_fv);
		int my_buf_num = 0;
		NumberEachProcs(my_buf_num, exp_cap);
		MPI_Info a_info;
		MPI_Info_create(&a_info);
		MPI_Info_set(a_info, "alloc_shared_noncontig", "true");
		MPI_Win_allocate_shared(my_buf_num*sizeof(T), sizeof(T), a_info, share_comm, &local_ptr, &win_fv);
		MPI_Win_shared_query(win_fv, 0, &bufsize, &disp_unit, &global_ptr);
		MPI_Info_free(&a_info);
		cap = exp_cap;
		if (be <= cap)
		{
			for (int j = bs; j < be; ++j)
			{
				int i0 = j - bs;
				global_ptr[j] = tempdata[i0];
			}
		}
		else
		{
			for (int j = bs; j < cap; ++j)
			{
				int i0 = j - bs;
				global_ptr[j] = tempdata[i0];
			}
		}
		MPI_Barrier(share_comm);			
	}

	void setnum_nocopy(const int & n0, const int & ghnum)
	{
		if (n0 > cap)
		{
			//if (ts > 0 && nrank == 0) {printf("The array should not be adjusted!!!\n"); MPI_Abort(MPI_COMM_WORLD, 457);}
			int exp_cap = int(array_exp_ratio*double(n0));
			int my_buf_num = 0;
			NumberEachProcs(my_buf_num, exp_cap);
			MPI_Win_free(&win_fv);
			MPI_Info a_info;
			MPI_Info_create(&a_info);
			MPI_Info_set(a_info, "alloc_shared_noncontig", "true");
			MPI_Win_allocate_shared(my_buf_num*sizeof(T), sizeof(T), a_info, share_comm, &local_ptr, &win_fv);
			MPI_Win_shared_query(win_fv, 0, &bufsize, &disp_unit, &global_ptr);
			MPI_Info_free(&a_info);
			cap = exp_cap;
		}
		total_num = n0;
		ghostnum = ghnum;
		num_without_ghost = n0 - ghostnum;
		new_cell_num = 0;
		new_start_index = n0;
		ArrayOrder();
		GhostOrder();
		// NewOrder();
	}

	void setnum_copy(const int & n0, const int & ghnum)
	{
		//int my_buf_num = 0;
		if (n0 > cap)
		{
			Adjustarray(n0);	
		}
		total_num = n0;
		ghostnum = ghnum;
		num_without_ghost = n0 - ghostnum;
		new_cell_num = 0;
		new_start_index = n0;
		ArrayOrder();
		GhostOrder();
		// NewOrder();
	}


	void IdentifyNewGhost(const vector<int> & newarray_tag, const vector<int> & procs_start)
	{
		for (int i = 0; i < newarray_tag.size(); ++i)
		{
			if (newarray_tag[i] > -1)
			{
				ghost_index.push_back(i+procs_start[srank]+total_num);
				//printf("box %d ghost tag is %d\n", i+procs_start[srank]+total_num, newarray_tag[i]);
			}
		}
	}

	void copydatatonewarray(T * newptr, const int & new_capability)
	{
		int bs, be;
		ArrayOrder_s(0, total_num, bs, be, sprocs, srank);
		if (be <= new_capability)
		{
			for (int j = bs; j < be; ++j)
			{
				newptr[j] = global_ptr[j];
			}
		}
		else
		{
			for (int j = bs; j < new_capability; ++j)
			{
				newptr[j] = global_ptr[j];
			}
		}	
	}

	void copydatatonewarray(vector<T> & newptr, const int & bs, const int & be, const int & new_capability)
	{
		if (be <= new_capability)
		{
			for (int j = bs; j < be; ++j)
			{
				int i0 = j - bs;
				newptr[i0] = global_ptr[j];
			}
		}
		else
		{
			for (int j = bs; j < new_capability; ++j)
			{
				int i0 = j - bs;
				newptr[i0] = global_ptr[j];
			}
		}	
	}

	// void setnum_notrelease(const int & n0)
	// {
	// 	num = n0;
	// 	ArrayOrder();
	// }

	// void clearloctag()
	// {
	// 	newloctag.resize(0);
	// }

	void CompressArray()
	{
		for (int i = 0; i < newloctag.size(); ++i)
		{
			// if (newloctag[i].newloc == 39231)
			// {
			// 	printf("My rank is %d the %d locchange is from %d to %d\n", srank, i, newloctag[i].oldloc, newloctag[i].newloc);
			// }			
			replace(i);
		}
		total_num -= totholenum;
		int procs_ghostnum = ghost_index.size();
		MPI_Reduce(&procs_ghostnum, &ghostnum, 1, MPI_INT, MPI_SUM, 0, share_comm);
		MPI_Bcast(&ghostnum, 1, MPI_INT, 0, share_comm);
		num_without_ghost = total_num - ghostnum;
		new_cell_num = 0;
		new_start_index = total_num;
		//printf("total number is %d ghost number is %d\n", total_num, ghostnum);
		//MPI_Barrier(share_comm);
		ArrayOrder();
		GhostOrder();	
		// NewOrder();
		MPI_Barrier(share_comm);
		hole.resize(0);
	}

	void Addnew(const vector<T> & new_ele, const vector<int> & ghosttag = vector<int>(0))
	{
		vector<int> new_size(sprocs, 0);
		int mynewsize = new_ele.size();
		int totnewsize = 0;
		MPI_Allgather(&mynewsize, 1, MPI_INT, &new_size[0], 1, MPI_INT, share_comm);	
		CountTotalNum(new_size, totnewsize);
		new_cell_num += totnewsize;
		int newnum;
		if (totnewsize > 0)
		{	
			vector<int> newelestart(sprocs);
			ArrayProcsStart(new_size, newelestart);
			if (ghosttag.size() > 0)
			{
				//printf("new ele size %d ghosttag size %d\n", (int)new_ele.size(), (int)ghosttag.size());
				IdentifyNewGhost(ghosttag, newelestart);
			}
			newnum = total_num + totnewsize;
			if (newnum > cap)
			{
				Adjustarray(newnum);
			}
			MPI_Win_fence(0, win_fv);
			for (int i = 0; i < new_size[srank]; ++i)
			{
				global_ptr[total_num+i+newelestart[srank]] = new_ele[i];
			}
			MPI_Win_fence(0, win_fv);
			total_num = newnum;
			new_start_index = total_num - new_cell_num;
			//PRINTFinLEVEL("add new total number is %d new_start_index is %d", 0, total_num, new_start_index);
		}
		// NewOrder();
	}
	void DirectlyReduceNew()
	{
		int ghostnum0 = ghost_index.size();
		MPI_Allreduce(&ghostnum0, &ghostnum, 1, MPI_INT, MPI_SUM, share_comm);
		num_without_ghost = total_num - ghostnum;
		new_cell_num = 0;
		new_start_index = total_num;
		ArrayOrder();
		GhostOrder();
		// NewOrder();
	}

	void RenewOrder()
	{

	}

	// void addnewblank(const int & i0)
	// {

	// 	vector<int> new_size(sprocs, 0);
	// 	int mynewsize = i0;
	// 	int totnewsize = 0;
	// 	MPI_Allgather(&mynewsize, 1, MPI_INT, &new_size[0], 1, MPI_INT, share_comm);	
	// 	CountTotalNum(new_size, totnewsize);
	// 	int newnum;
	// 	if (totnewsize > 0)
	// 	{	
	// 		vector<int> newelestart(sprocs);
	// 		ArrayProcsStart(new_size, newelestart);
	// 		newnum = num + totnewsize;
	// 		if (newnum > cap)
	// 		{
	// 			Adjustarray(newnum);
	// 		}
	// 		num = newnum;
	// 		ArrayOrder();
	// 	}
	// }

	void cleartemparray()
	{
		//vector<BoxlocChange>(0).swap(newloctag);
		newloctag.resize(0);
		hole.resize(0);
		//new_ghost_tag.resize(0);
		//ghost_index.resize(0);
	}

	void holeplan(const bool & holeflag)
	{
		int global_start, global_end;
		GlobalOrder(global_start, global_end);
		newloctag.resize(0);
		vector<int> holenum(sprocs, 0);
		vector<int> * holearray;
		if (holeflag)
		{
			holearray = &hole;
		}
		else
		{
			holearray = &ghost_index;
		}
		int myholenum = holearray->size();
		//printf("[%d] hole number is %d\n", srank, myholenum);
		//ShowAllRankData("mesh hole number", myholenum, 5);		
		MPI_Allgather(&myholenum, 1, MPI_INT, &holenum[0], 1, MPI_INT, share_comm);
		//MPI_Allreduce(MPI_IN_PLACE, &holenum[0], sprocs, MPI_INT, MPI_SUM, share_comm);
		// MPI_Barrier(share_comm);
		// printf("finish gather hole number\n");
		CountTotalNum(holenum, totholenum);
		//printf("hole num is %d total is %d data array size is %d\n", myholenum, totholenum, total_num);
		if (totholenum > 0)
		{
			vector<int> holestart(sprocs, 0);
			ArrayProcsStart(holenum, holestart);
			vector<int> holeloc(totholenum);
			if (holenum[srank] == 0)
			{
				holearray->push_back(0);
			}
			//printf("holearray capacity is %d\n", (int)holearray->capacity());
			//cout << srank << " " << &holearray->at(0) << endl;
			//MPI_Barrier(share_comm);
			//printf("start gather hole location!!!\n");
			MPI_Allgatherv(&holearray->at(0), holenum[srank], MPI_INT,
				&holeloc[0], &holenum[0], &holestart[0], MPI_INT, share_comm);
			//MPI_Barrier(share_comm);
			//printf("finish gather hole location!!!\n");
			if (holenum[srank]==0)
			{
				holearray->pop_back();
			}
			// FillTailHole(num, n_totholenum, n_totholenum, holeloc);
			int fillstart = total_num - totholenum;
			int arrid = total_num; 
			int holeid = totholenum;
			vector<int> tailhole(totholenum, -1);
			//printf("data array size is %d\n", total_num);
			//ShowAllRankData("hole fillstart", fillstart, 5);
			//ShowAllRankData("total hole num", totholenum, 5);
			for (int i = 0; i < totholenum; ++i)
			{
				if (holeloc[i] >= fillstart)
				{
					//printf("holeloc %d is %d fillstart is %d\n", i, holeloc[i], fillstart);
					tailhole[holeloc[i]-fillstart] = i;
				}
			}
			//printf("start while cycle!!!\n");
			while(arrid > fillstart && holeid > 0)
			{
				--arrid;
				Assert((arrid-fillstart > -1), "Error in hole plan!!!", 690);
				if (tailhole[arrid - fillstart] == -1)
				{
					
					while (holeid > 0)
					{
						--holeid;
						if (holeloc[holeid] < fillstart)
						{
							if (holeloc[holeid] >= global_start && holeloc[holeid] < global_end)
							{
								//printf("A new location change oldloc is %d newloc is %d, fillstart is %d total_num %d hole number %d\n",
								// 	arrid,holeloc[holeid],fillstart, total_num, totholenum);
								newloctag.push_back(BoxlocChange(arrid, holeloc[holeid]));
							}
							break;
						}
					}
				}
			}
			//printf("newloctag size is %d\n", (int)newloctag.size());			
		}		
	}
	/*Move the last elements to fill the holes*/
	// void FillTailHole(const int & arrlastcell, 
	// 									const int & leftholenum, 
	// 									const int & totalholenum, 
	// 									const int * holeloc, 
	// 									const int & array_procs_start,
	// 									const int & array_procs_end)
	// {
	// 	int fillstart = arrlastcell - leftholenum;
	// 	int arrid = arrnum; 
	// 	int holeid = leftholenum;
	// 	vector<int> tailhole(leftholenum, -1);
	// 	ShowAllRankData("hole fillstart", fillstart, 4);
	// 	ShowAllRankData("total hole num", totalholenum, 4);
	// 	for (int i = totholenum - leftholenum; i < totholenum; ++i)
	// 	{
	// 		if (holeloc[i] >= fillstart)
	// 		{
	// 			tailhole[holeloc[i]-fillstart] = i;
	// 		}
	// 	}
	// 	while(arrid > 0 && holeid > 0)
	// 	{
	// 		--arrid;
	// 		if (tailhole[arrid - fillstart] == -1)
	// 		{				
	// 			while (holeid > 0)
	// 			{
	// 				--holeid;
	// 				if (holeloc[holeid] < fillstart)
	// 				{
	// 					if (holeloc[holeid] >= array_procs_start && holeloc[holeid] < array_procs_end)
	// 					{
	// 						newloctag.push_back(BoxlocChange(arrid, holeloc[holeid]));
	// 					}
	// 					break;
	// 				}
	// 			}
	// 		}
	// 	}
	// }

	void SwitchTwoData()
	{
		for (int i = 0; i < newloctag.size(); ++i)
		{
			switchtwo(i);
		}
	}

	inline T & operator[](const int & i)
	{
		return global_ptr[i];
	}

	void cleararray()
	{
		MPI_Win_free(&win_fv);
	}
	
};
#endif