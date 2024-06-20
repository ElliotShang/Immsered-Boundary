#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include "FlowVariables.H"
#include "AMRSpaceTime.H"
#include "Mesh.H"
#include "InputFunc.H"

using namespace std;

void ReadLevelSolution(Mesh & amesh, Pltoutput & pltfile)
{
	Assert(init_level_num == filelevelnum, "The input solution file level number error!!!", 46);
	Assert(nodenum == filenodenum, "The input solution file node number error!!!", 47);
	vector<int> fvnum = pltfile.GetReadorder();
	if (srank == 0)
	{
		ifstream solnstream(soln_file_name);
		if (!solnstream.good())
		{
			printf("the input solution file %s error!!!\n", soln_file_name);
			MPI_Abort(MPI_COMM_WORLD, 49);
		}
		JumpIfstreamLine(solnstream, init_jump);
		for (int b0 = 0; b0 < filelevelnum*filenodenum+filebodynum; ++b0)
		{
			char a[200];
			solnstream.getline(a, sizeof(a));
			PRINTFinLEVEL("A is %s b0 is %d", 0, a, b0);
			int n0 = GetANumAfterString(a, "N-", ' ');
			int level0 = GetANumAfterString(a, "L-", ' ');
			int ele_num = GetANumAfterString(a, "ELEMENTS=", ',');
			PRINTFinLEVEL("The input node %d level %d element number %d!!!",level0,n0,level0,ele_num);
			Assert(ele_num > -1, "The input element number must be positive!!!", 78);
			JumpIfstreamLine(solnstream, 1);
			if (n0 == -1 || level0 == -1 || node != n0)
			{
				for (int v0 = 0; v0 < fvnum.size(); ++v0)
				{
					JumpIfstreamString(solnstream, ele_num);
				}
			}
			else
			{
				if (node == n0)
				{	
					DataArray<FlowVariables> & adata = amesh.LevelDataArray(level0);
#ifdef INPUT_GHOST								
					if (ele_num != adata.realsize())
#else	
					if (ele_num != adata.size())
#endif					
					{
						printf("L%d The input data number %d does not equal to the array size %d!!!", level0, ele_num, adata.size());
						MPI_Abort(MPI_COMM_WORLD, 65);
					}
					for (int v0 = 0; v0 < fvnum.size(); ++v0)
					{
						if (fvnum[v0] > -1)
						{
#ifdef INPUT_GHOST						
							for (int i = 0; i < adata.realsize(); ++i)
#else 
							for (int i = 0; i < adata.size(); ++i)
#endif													
							{
								solnstream >> adata[i][fvnum[v0]]; 
								// if (fvnum[v0]==7 && adata[i][fvnum[v0]] < 0.00001)
								// {
								// 	printf("Level %d flow[7] is %f!!!\n", level0, adata[i][fvnum[v0]]);
								// 	MPI_Abort(MPI_COMM_WORLD, 72);
								// }
							}
							char c0[10];
							solnstream.getline(c0, sizeof(c0));
						}
						else
						{
							JumpIfstreamString(solnstream, ele_num);
						}
					}
				}
			}
		}
		solnstream.close();
	}
	MPI_Barrier(MPI_COMM_WORLD);
}

void JumpIfstreamLine(ifstream & as, const int & jl)
{
	char a[300];
	for (int i = 0; i < jl; ++i)
	{
		as.getline(a,sizeof(a));
	}
}

void JumpIfstreamString(ifstream & as, const int & js)
{
	double a;
	for (int i = 0; i < js; ++i)
	{
		as >> a;
	}
	JumpIfstreamLine(as, 1);
}

int GetANumAfterString(char * a, const string & astring, const char & endchar)
{
	char num[20];
	int clength = astring.length();
	//printf("the line with number is %s\n", a);
	int elementnum = -1;
	for (int i = 0; i < strlen(a)-clength; ++i)
	{
		string c0(a+i, a+(i+clength));
		if (c0 == astring)
		{
			int commaloc = -1;
			for (int i0 = i+clength; i0 < strlen(a); ++i0)
			{
				if (a[i0] == endchar)
				{
					commaloc = i0;
					goto FINDENDFLAG;
				}
			}
			if (commaloc == -1)
			{
				return -1;
			}
			FINDENDFLAG:;
			if (commaloc < 1 || commaloc > strlen(a)-1)
			{
				printf("Error in found the end location when getting the specific string!!!\n");
				MPI_Abort(MPI_COMM_WORLD, 154);
			}
			stringstream estream;
			estream << string(a+(i+clength), a+commaloc);
			estream >> elementnum;
			break;
		}
	}
	// if (elementnum < 0)
	// {
	// 	printf("The input number is non-positive!!! Input array is %s tag is %s\n", a, astring.c_str());
	// 	MPI_Abort(MPI_COMM_WORLD, 147);
	// }
	return elementnum;
}

void InputSoln(Mesh & amesh, BCValues & meshbc, AMR & myamr, Pltoutput & pltfile)
{
	ReadLevelSolution(amesh, pltfile);
	for (int i = 0; i < amesh.MyCurNum(); ++i)
	{
		DataArray<FlowVariables> & da = amesh.LevelDataArray(i);
		int bs = da.ps();
		int be = da.pe();
		for (int b0 = bs; b0 < be; ++b0)
		{
			Get_E(da[b0]);
		}
	}
	MPI_Barrier(share_comm);
	for (int i = 0; i < amesh.MyCurNum()-1; ++i)
	{	
		myamr.InterfaceRestriction(i);		
	}
	for (int i = 0; i < amesh.MyCurNum(); ++i)
	{
		amesh.DataExchange(amesh.LevelDataWin(i), i);		
	}
	for (int i = 0; i < amesh.MyCurNum()-1; ++i)
	{	
		myamr.InterfaceProlongation(i);
	}
	meshbc.BCBoxTreat(amesh);
}