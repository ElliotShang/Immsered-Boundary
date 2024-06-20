#ifndef _Check_H_
#define _Check_H_

#include <iostream>
#include <string>
#include <mpi.h>

using namespace std;

#define VERBOSITY 5

#ifdef DEBUG
#	define Assert(cond, label, abortnum) \
				if (!(cond))   \
				{              \
					cout << "Debug check failed! Fail label is " << label << "\n"; \
					MPI_Abort(MPI_COMM_WORLD, abortnum);   \
				}
# define ShowAllRankData(arrayname, arraydata, vb) \
				if (vb == VERBOSITY)\
				{\
					printf("<<<My rank is %d and %s is %d\n", nrank, arrayname, arraydata);\
				}				
# define ShowOneRankData(rank0, arrayname, arraydata, vb) \
				if (nrank == rank0 && vb == VERBOSITY) \
				{									 \
					printf("<<<My rank is %d and %s is %d\n", nrank, arrayname, arraydata);\
				}
#	define GiveAFlag(flag, vb)\
				MPI_Barrier(MPI_COMM_WORLD);\
				if (nrank == 0 && vb == VERBOSITY)\
				{\
					printf("<<<All processors have arrived %s>>>\n", flag);\
				}
#	define GiveLevelFlag(flag, curlevel, vb)\
				MPI_Barrier(MPI_COMM_WORLD);\
				if (nrank == 0 && vb == VERBOSITY)\
				{\
					printf("<<<All processors have arrived %s in level %d>>>\n", flag, curlevel);\
				}
# define PRINTF(stat, args...)\
				string stat0 = "N[%d]R[%d] ";\
				stat0 += stat;\
				stat0 += "\n";\
				printf(stat0.c_str(),node,srank,args)
# define PRINTFinLEVEL(stat,ilevel,args...)\
				{\
					string stat0 = "N[%d]R[%d]L[%d] ";\
					stat0 += stat;\
					stat0 += "\n";\
					printf(stat0.c_str(),node,srank,ilevel,args);\
				}
# define PRINTFinLEVELRANK0(stat,ilevel,args...)\
				if (srank == 0)\
				{\
					string stat0 = "N[%d]R[%d]L[%d] ";\
					stat0 += stat;\
					stat0 += "\n";\
					printf(stat0.c_str(),node,srank,ilevel,args);\
				}		
#else
#	define Assert(cond, label, abortnum) void(0)
# define ShowAllRankData(arrayname, arraydata, vb) void(0)
# define ShowOneRankData(rank0, arrayname, arraydata, vb) void(0)
# define GiveAFlag(flag, vb) void (0)	
#	define GiveLevelFlag(flag, vb, curlevel) void (0)
# define PRINTF(stat, args...) void (0)
# define PRINTFinLEVEL(stat,ilevel,args...) void (0)
# define PRINTFinLEVELRANK0(stat,ilevel,args...) void (0)
#endif


#endif

