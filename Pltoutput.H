#ifndef _Pltoutput_H_
#define _Pltoutput_H_

//#include <iostream>
#include <string>
#include <cstdio>
#include <mpi.h>
#include <sstream>

#include "Point_iterator.H"
#include "Mesh.H"
#include "Meshlevel.H"
#include "AMRmpi.H"
#include "Body.H"
#include "NS_Solver.H"

extern bool additional_output;

using namespace std;

//#define INFECT_DISTANCE

class Pltoutput
{
public:

	Pltoutput(Mesh & a_mesh):
	m_mesh(a_mesh)
	{
		filename = NULL;
		fpt = NULL;
		ffm = "%16.8f";
#ifdef TURBULENCE
		varnum = 7;
#else				
		varnum = 6;
#endif
		other_varnum = 7;
		grid_varnum = 3;
		SetVarloc();
		Filereadorder();
		timeavenum = 0;
	}

	bool isnull()
	{
		if (NULL != filename)
		{
			return true;
		}
		else
		{
			return false;
		}
	}

	void Filereadorder()
	{
		readorder.resize(varnum+other_varnum, -1);
		for (int i = 0; i < 6; ++i)
		{
			readorder[i] = i;
		}
#ifdef TURBULENCE
		readorder[6] = 7;
#endif
	}

	void SetVarloc()
	{
		varloc.resize(varnum, -1);
		for (int i = 0; i < varnum; ++i)
		{
			varloc[i] = i;
		}
#ifdef TURBULENCE
		varloc[6] = 7;
#endif				
	}

	inline vector<int> & GetReadorder()
	{
		return readorder;
	}

	void Create_a_File(string a_file)
	{
		if (NULL != filename)
		{
			delete[] filename;
			filename = NULL;
		}
		filename = new char[a_file.length()+1];		
		sprintf(filename, "%s", a_file.c_str());
	}

	inline void Print_a_Flowvar(DataArray<FlowVariables> & a_data, const int & vari)
	{
		int s = 0;
#ifdef PRINTGHOST		
		for (int i = 0; i < a_data.realsize(); ++i)
#else
		for (int i = 0; i < a_data.size(); ++i)	
#endif				
		{
			fprintf(fpt, ffm.c_str(), a_data[i][vari]);
			++s;
			if (s%8 == 0) fprintf(fpt, "\n");
		}
		if (s%8 != 0) fprintf(fpt, "\n");		
	}

	void Print_Patch_Params(string a_file, vector<Body> & abody);

	void Print_Distance(DataArray<Box> & abox, const int & ilevel)
	{
		if (m_mesh.infectbox.size() > 0)
		{
			int s = 0;
#ifdef INFECT_DISTANCE			
			if (ilevel == m_mesh.cur_level_num - 1)
			{
	#ifdef PRINTGHOST			
				for (int i = 0; i < abox.realsize(); ++i)
	#else
				for (int i = 0; i < abox.size(); ++i)
	#endif
				{
					if (m_mesh.BoxInfectIndex(i) > -1)
					{
						fprintf(fpt, ffm.c_str(), m_mesh.InfectedBoxDis(i));
					}
					else
					{
						//ZERODIS:;
						fprintf(fpt, " 0");
					}
					++s;
					if (s%8 == 0) fprintf(fpt, "\n");
				}
				if (s%8 != 0) fprintf(fpt, "\n");
			}
			else
			{
	#ifdef PRINTGHOST			
				Print_Zero_Distance(abox.realsize());
	#else
				Print_Zero_Distance(abox.size());
	#endif				
			}

#else

	#ifdef PRINTGHOST			
			for (int i = 0; i < abox.realsize(); ++i)
	#else
			for (int i = 0; i < abox.size(); ++i)
	#endif							
			{
				fprintf(fpt, ffm.c_str(), abox[i].pair.signdis);							
				++s;
				if (s%8 == 0) fprintf(fpt, "\n");
			}
			if (s%8 != 0) fprintf(fpt, "\n");
#endif			
		}		
		else
		{
#ifdef PRINTGHOST			
			Print_Zero_Distance(abox.realsize());
#else
			Print_Zero_Distance(abox.size());
#endif					
		}
	}

	void Print_Domain_Distance(DataArray<Box> & abox)
	{
		int s = 0;
#ifdef PRINTGHOST			
		for (int i = 0; i < abox.realsize(); ++i)
#else
		for (int i = 0; i < abox.size(); ++i)
#endif							
		{
			fprintf(fpt, ffm.c_str(), abox[i].pair.distance_to_dm);							
			++s;
			if (s%8 == 0) fprintf(fpt, "\n");
		}
		if (s%8 != 0) fprintf(fpt, "\n");
	}

	void Print_Zero_Distance(const int & boxnum)
	{
		int s = 0;
		for (int i = 0; i < boxnum; ++i)
		{
			fprintf(fpt, " 0");
			++s;
			if (s%8 == 0) fprintf(fpt, "\n");
		}
		if (s%8 != 0) fprintf(fpt, "\n");
	}

	void Print_Body_Nmv(Body & a_body, const int & nvi)
	{
		int s = 0;
		int patnum = a_body.patch.size();
		for (int i = 0; i < patnum; ++i)
		{
			fprintf(fpt, ffm.c_str(), a_body.patch[i].nv[nvi]);
			++s;
			if (s%8 == 0) fprintf(fpt, "\n");
		}
		if (s%8 != 0) fprintf(fpt, "\n");
	}

	void Print_Body_ReverseFlag(Body & a_body)
	{
		int s = 0;
		int patnum = a_body.patch.size();
		for (int i = 0; i < patnum; ++i)
		{
			fprintf(fpt, "%d ", a_body.patch[i].reverse_nv_dir);
			++s;
			if (s%8 == 0) fprintf(fpt, "\n");
		}
		if (s%8 != 0) fprintf(fpt, "\n");
	}

	void Print_Body_Yplus(Body & a_body)
	{
		int s = 0;
		int patnum = a_body.patch.size();
		for (int i = 0; i < patnum; ++i)
		{
			fprintf(fpt, ffm.c_str(), a_body.patch[i].yplus);
			++s;
			if (s%8 == 0) fprintf(fpt, "\n");
		}
		if (s%8 != 0) fprintf(fpt, "\n");
	}

	void Print_Body_HGVar(Body & a_body, const int & vari0)
	{
		int s = 0;
		int patnum = a_body.patch.size();
		for (int i = 0; i < patnum; ++i)
		{
			fprintf(fpt, ffm.c_str(), a_body.patch[i].hgc.fv[vari0]);
			++s;
			if (s%8 == 0) fprintf(fpt, "\n");
		}
		if (s%8 != 0) fprintf(fpt, "\n");
	}
	void Print_Body_HGdis(Body & a_body)
	{
		int s = 0;
		int patnum = a_body.patch.size();
		for (int i = 0; i < patnum; ++i)
		{
			fprintf(fpt, ffm.c_str(), a_body.patch[i].hgc.hgdis);
			++s;
			if (s%8 == 0) fprintf(fpt, "\n");
		}
		if (s%8 != 0) fprintf(fpt, "\n");
	}

	void Print_Body_Bi(Body & a_body)
	{
		int s = 0;
		int patnum = a_body.patch.size();
		for (int i = 0; i < patnum; ++i)
		{
			fprintf(fpt, "%d ", a_body.bodyindex);
			++s;
			if (s%8 == 0) fprintf(fpt, "\n");
		}
		if (s%8 != 0) fprintf(fpt, "\n");
	}

	void Print_Body_Pi(Body & a_body)
	{
		int s = 0;
		int patnum = a_body.patch.size();
		for (int i = 0; i < patnum; ++i)
		{
			fprintf(fpt, "%d ", i);
			++s;
			if (s%8 == 0) fprintf(fpt, "\n");
		}
		if (s%8 != 0) fprintf(fpt, "\n");
	}

	void Print_Ma(DataArray<FlowVariables> & a_data)
	{
		int s = 0;
#ifdef PRINTGHOST
		for (int i = 0; i < a_data.realsize(); ++i)
#else 
		for (int i = 0; i < a_data.size(); ++i)
#endif			
		{
			fprintf(fpt, ffm.c_str(), Manumber(a_data[i]));
			++s;
			if (s%8 == 0) fprintf(fpt, "\n");
		}
		if (s%8 != 0) fprintf(fpt, "\n");		
	}

	void Print_Omega_Z(DataArray<FlowVariables> & a_data, DataArray<Box> & a_box, const int & cur_level0);

	void Print_SolidFlag(DataArray<Box> & a_box)
	{
		int s = 0;
#ifdef PRINTGHOST
		int boxnum = a_box.realsize();
#else 
		int boxnum = a_box.size();	
#endif
		for (int i = 0; i < boxnum; ++i)			
		{
			fprintf(fpt, " %d", a_box[i].solid);
			++s;
			if (s%8 == 0) fprintf(fpt, "\n");
		}
		if (s%8 != 0) fprintf(fpt, "\n");
	}

	void Print_Level_Soln(Meshlevel & a_level, const int & cur_level0);

	void Print_Level_Grid(Meshlevel & a_level, const int & cur_level0);

	void Print_Level_CellParams(Meshlevel & a_level, const int & cur_level0);

	void Print_Body_Grid(Body & a_body, const int & body0);

	void Print_Body_Soln(Body & a_body, const int & body0);

	void Output_Grid(string a_file, vector<Body> & abody);

	void Output_Grid_Nobody(string a_file);

	void Output_Soln(string a_file, vector<Body> & abody);

	void Output_Soln_Nobody(string a_file);

	void Output_GridParameters(string a_file);

	void Output_Soln_Onlybody(string a_file, vector<Body> & abody);

	void Output_Grid_Onlybody(string a_file, vector<Body> & abody);

	void Print_Mesh_Header()
	{
		if (nrank == 0)
		{
			fpt = fopen(filename, "w");
			if (NULL != fpt && isgrid)
			{
				fprintf(fpt, "TITLE = \"Boxgrid\"\n");
				fprintf(fpt, "FILETYPE = GRID\n");
				fprintf(fpt, "VARIABLES = \"x\",\"y\",\"z\"\n");
			}
			else if (NULL != fpt && !isgrid)
			{
				fprintf(fpt, "TITLE = \"Flow Solution\"\n");
				fprintf(fpt, "FILETYPE = SOLUTION\n");
#ifdef TURBULENCE
				fprintf(fpt, "VARIABLES = \"roe\",\"u\",\"v\",\"w\",\"p\",\"T\",\"viseddy\",\"Ma\",\"dis\",\"dmdis\",\"i_index\",\"j_index\",\"k_index\",\"Omegaz\"\n");
#else				
				fprintf(fpt, "VARIABLES = \"roe\",\"u\",\"v\",\"w\",\"p\",\"T\",\"Ma\",\"dis\",\"dmdis\",\"i_index\",\"j_index\",\"k_index\",\"Omegaz\"\n");
#endif			
			}
			else
			{
				printf("<<<The file does not exist!!! Please check it!!!\n");
				MPI_Abort(MPI_COMM_WORLD, 78);
			}
			fclose(fpt);
			fpt = NULL;
		}	
	}

	void Print_GridParams_Header()
	{
		if (nrank == 0)
		{
			fpt = fopen(filename, "w");
			fprintf(fpt, "TITLE = \"Cell Parameters\"\n");
			fprintf(fpt, "FILETYPE = SOLUTION\n");
			fprintf(fpt, "VARIABLES = \"keisax\",\"keisay\",\"keisaz\",\"etax\",\"etay\",\"etaz\",\"zeitax\",\"zeitay\",\"zeitaz\"\n");			
			fclose(fpt);
			fpt = NULL;
		}
	}

	void clearfile(const char * fname)
	{
		if (nrank == 0)
		{
			fpt = fopen(fname, "w");
			fclose(fpt);
			fpt = NULL;
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}

	void Design_User_File();

	void Create_User_File()
	{
		userfileptr.resize(userfilenum);
		for (int i = 0; i < userfilenum; ++i)
		{
			clearfile(userfilename[i]);
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}

	void Print_Var_as_Func_of_Time(vector<Body> & abody, NS_Solver & asolver);

	void Print_Var_as_Func_of_Space(vector<Body> & abody);

	void Print_IBcell_yplus(string & bfile);

	void Print_Box_Index(DataArray<Box> & abox, const int & id0);

private:
	Mesh & m_mesh;
	char * filename;
	FILE * fpt;
	string ffm;
	bool isgrid;
	int varnum;
	int grid_varnum;
	int other_varnum;
	int userfilenum;
	vector<char *> userfilename;
	vector<FILE *> userfileptr;
	vector<int> readorder;
	vector<int> varloc;

	vector<vector<double> > timeavedata1[30];
	int timeavenum;
};

#endif
