#include "Pltoutput.H"

	void Pltoutput::Output_Grid(string a_file, vector<Body> & abody)
	{
		isgrid = true;
		MPI_Barrier(MPI_COMM_WORLD);
		char ctime[32];
		sprintf(ctime, "-%d.plt", ts);
		a_file += ctime;
		Create_a_File(a_file);
		Print_Mesh_Header();
		for (int ni = 0; ni < nodenum; ++ni)
		{
			if (srank == 0 && node == ni)
			{
				fpt = fopen(filename, "a");
				for (int i = 0; i < m_mesh.cur_level_num; ++i)
				{
					Print_Level_Grid(m_mesh.m_level[i], i);
				}
				if (node == 0)
				{
					for (int i = 0; i < abody.size(); ++i)
					{
						Print_Body_Grid(abody[i], i);
					}
				}
				fclose(fpt);
				fpt = NULL;
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
	}

	void Pltoutput::Output_Grid_Onlybody(string a_file, vector<Body> & abody)
	{
		isgrid = true;
		MPI_Barrier(MPI_COMM_WORLD);
		char ctime[32];
		sprintf(ctime, "-%d.plt", ts);
		a_file += ctime;
		Create_a_File(a_file);
		Print_Mesh_Header();
		for (int ni = 0; ni < nodenum; ++ni)
		{
			if (srank == 0 && node == ni)
			{
				fpt = fopen(filename, "a");
				for (int i = 0; i < abody.size(); ++i)
				{
					Print_Body_Grid(abody[i], i);
				}
				fclose(fpt);
				fpt = NULL;
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
	}

	void Pltoutput::Output_Grid_Nobody(string a_file)
	{
		isgrid = true;
		MPI_Barrier(MPI_COMM_WORLD);
		char ctime[32];
		sprintf(ctime, "-%d.plt", ts);
		a_file += ctime;
		Create_a_File(a_file);
		Print_Mesh_Header();
		for (int ni = 0; ni < nodenum; ++ni)
		{
			if (srank == 0 && node == ni)
			{
				fpt = fopen(filename, "a");
				for (int i = 0; i < m_mesh.cur_level_num; ++i)
				{
					Print_Level_Grid(m_mesh.m_level[i], i);
				}
				fclose(fpt);
				fpt = NULL;
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
	}

	void Pltoutput::Output_Soln(string a_file, vector<Body> & abody)
	{
		isgrid = false;		
		char ctime[32];
		sprintf(ctime, "-%d.plt", ts);
		a_file += ctime;
		MPI_Barrier(MPI_COMM_WORLD);
		Create_a_File(a_file);
		Print_Mesh_Header();
		for (int ni = 0; ni < nodenum; ++ni)
		{
			if (srank == 0 && node == ni)
			{
				fpt = fopen(filename, "a");
				for (int i = 0; i < m_mesh.cur_level_num; ++i)
				{
					Print_Level_Soln(m_mesh.m_level[i], i);
				}
				//printf("Finish print the flow solution!!!\n");
				if (node == 0)
				{
					for (int i = 0; i < abody.size(); ++i)
					{
						Print_Body_Soln(abody[i], i);
						//abody[i].Printout_soln(fpt, varnum+other_varnum);
					}
				}
				//printf("print the solution for the body\n");
				fclose(fpt);
				fpt = NULL;
			}
			MPI_Barrier(MPI_COMM_WORLD);			
		}
		//MPI_Barrier(MPI_COMM_WORLD);
		Print_Patch_Params(a_file, abody);
	}

	void Pltoutput::Output_Soln_Onlybody(string a_file, vector<Body> & abody)
	{
		isgrid = false;		
		char ctime[32];
		sprintf(ctime, "-%d.plt", ts);
		a_file += ctime;
		MPI_Barrier(MPI_COMM_WORLD);
		Create_a_File(a_file);
		Print_Mesh_Header();
		for (int ni = 0; ni < nodenum; ++ni)
		{
			if (srank == 0 && node == ni)
			{
				fpt = fopen(filename, "a");
				for (int i = 0; i < abody.size(); ++i)
				{
					Print_Body_Soln(abody[i], i);
				}
				fclose(fpt);
				fpt = NULL;
			}
			MPI_Barrier(MPI_COMM_WORLD);			
		}
	}

	void Pltoutput::Output_Soln_Nobody(string a_file)
	{
		isgrid = false;		
		char ctime[32];
		sprintf(ctime, "-%d.plt", ts);
		a_file += ctime;
		MPI_Barrier(MPI_COMM_WORLD);
		Create_a_File(a_file);
		Print_Mesh_Header();
		for (int ni = 0; ni < nodenum; ++ni)
		{
			if (srank == 0 && node == ni)
			{
				fpt = fopen(filename, "a");
				for (int i = 0; i < m_mesh.cur_level_num; ++i)
				{
					Print_Level_Soln(m_mesh.m_level[i], i);
				}
				fclose(fpt);
				fpt = NULL;
			}
			MPI_Barrier(MPI_COMM_WORLD);			
		}
	}

	void Pltoutput::Print_Patch_Params(string a_file, vector<Body> & abody)
	{
		char ctime[32];
		sprintf(ctime, "-%d-patchdata.dat", ts);
		a_file += ctime;
		Create_a_File(a_file);
		MPI_Barrier(MPI_COMM_WORLD);
		for (int i = 0; i < abody.size(); ++i)
		{
			for (int ni = 0; ni < nodenum; ++ni)
			{
				if (srank == 0 && node == ni)
				{
					fpt = fopen(filename, "a");
					int patnum = abody[i].patch.size();
					for (int p0 = 0; p0 < patnum; ++p0)
					{
						if (abody[i].patch[p0].node == ni)
						{
					 		fprintf(fpt, "%d %d %12.6f %12.6f %12.6f %12.6f\n",
					 			i,p0,
					 			abody[i].patch[p0].hgc.fv.p,
					 			abody[i].patch[p0].ut,
					 			abody[i].patch[p0].yplus,
					 			abody[i].patch[p0].Cf);
						}
					}
					fclose(fpt);
					fpt = NULL;
				}
				MPI_Barrier(MPI_COMM_WORLD);
			}			
		}
	}

	void Pltoutput::Output_GridParameters(string a_file)
	{
		char ctime[32];
		sprintf(ctime, "-%d.plt", ts);
		a_file += ctime;
		MPI_Barrier(MPI_COMM_WORLD);
		Create_a_File(a_file);
		Print_GridParams_Header();
		for (int ni = 0; ni < nodenum; ++ni)
		{
			if (srank == 0 && node == ni)
			{
				fpt = fopen(filename, "a");
				for (int i = 0; i < m_mesh.cur_level_num; ++i)
				{
					Print_Level_CellParams(m_mesh.m_level[i], i);
				}
				fclose(fpt);
				fpt = NULL;
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
	}

	void Pltoutput::Print_Level_Grid(Meshlevel & a_level, const int & cur_level0)
	{
		// string outstr;
		// for (int i = 0; i < grid_varnum; ++i)
		// {
		// 	outstr += ffm;
		// }
		// outstr += "\n";
		char zonename[32];	
		sprintf(zonename, "N-%d L-%d", node, cur_level0);
		// stringstream ss;
		// ss.clear();
		// ss << "R-" << nrank      << "-";
		// ss << "L-" << cur_level0 << "-";
		// ss << "B-" << j;
		// string zonename = ss.str();
		if (a_level.m_point.size() > 0)
		{
#ifdef PRINTGHOST			
			fprintf(fpt, "zone T = \"%s\", NODES=%d, ELEMENTS=%d, DATAPACKING=POINT, ZONETYPE=FEBRICK\n", 
				zonename, a_level.m_point.size(), a_level.m_box.realsize());
#else
			fprintf(fpt, "zone T = \"%s\", NODES=%d, ELEMENTS=%d, DATAPACKING=POINT, ZONETYPE=FEBRICK\n", 
				zonename, a_level.m_point.size(), a_level.m_box.size());
#endif						
		}
		//N=%d, E=%d, F=fepoint, ET=brick								
		for (int i = 0; i < a_level.m_point.size(); ++i)
		{
			fprintf(fpt, "%16.8f%16.8f%16.8f\n", a_level.m_point[i][0], 
				a_level.m_point[i][1], a_level.m_point[i][2]);
		}
#ifdef PRINTGHOST		
		for (int i = 0; i < a_level.m_box.realsize(); ++i)
#else
		for (int i = 0; i < a_level.m_box.size(); ++i)	
#endif				
		{
			fprintf(fpt, "%d %d %d %d %d %d %d %d\n", 
				a_level.m_box[i].pts[0][0][0]+1,
				a_level.m_box[i].pts[1][0][0]+1,
				a_level.m_box[i].pts[1][1][0]+1,
				a_level.m_box[i].pts[0][1][0]+1,
				a_level.m_box[i].pts[0][0][1]+1,
				a_level.m_box[i].pts[1][0][1]+1,
				a_level.m_box[i].pts[1][1][1]+1,
				a_level.m_box[i].pts[0][1][1]+1);
		}
	}

	void Pltoutput::Print_Level_CellParams(Meshlevel & a_level, const int & cur_level0)
	{
		// string outstr;
		// for (int i = 0; i < grid_varnum; ++i)
		// {
		// 	outstr += ffm;
		// }
		// outstr += "\n";
		char zonename[32];	
		sprintf(zonename, "N-%d L-%d", node, cur_level0);
		// stringstream ss;
		// ss.clear();
		// ss << "R-" << nrank      << "-";
		// ss << "L-" << cur_level0 << "-";
		// ss << "B-" << j;
		// string zonename = ss.str();
		if (a_level.m_point.size() > 0)
		{
#ifdef PRINTGHOST			
			fprintf(fpt, "zone T = \"%s\", NODES=%d, ELEMENTS=%d, DATAPACKING=BLOCK, ZONETYPE=FEBRICK\n", 
				zonename, a_level.m_point.size(), a_level.m_box.realsize());
#else
			fprintf(fpt, "zone T = \"%s\", NODES=%d, ELEMENTS=%d, DATAPACKING=BLOCK, ZONETYPE=FEBRICK\n", 
				zonename, a_level.m_point.size(), a_level.m_box.size());
#endif
			fprintf(fpt, "VARLOCATION = ([1-%d]=CELLCENTERED)\n", 9);						
		}
		//N=%d, E=%d, F=fepoint, ET=brick								
		for (int dx = 0; dx < 3; ++dx)
		{
			for (int dy = 0; dy < 3; ++dy)
			{
				int s = 0;
#ifdef PRINTGHOST		
				for (int i = 0; i < a_level.m_box.realsize(); ++i)
#else
				for (int i = 0; i < a_level.m_box.size(); ++i)	
#endif				
				{				
					fprintf(fpt, "%20.4f", a_level.m_geom[i].keisa[dx][dy]);
					++s;
					if (s%8 == 0) fprintf(fpt, "\n");
				}
				if (s%8 != 0) fprintf(fpt, "\n");
			}
		}
	}

	void Pltoutput::Print_Body_Grid(Body & a_body, const int & body0)
	{
		// string outstr;
		// for (int i = 0; i < grid_varnum; ++i)
		// {
		// 	outstr += ffm;
		// }
		// outstr += "\n";
		char zonename[32];	
		sprintf(zonename, "Body-%d-N-%d", body0, node);
		// stringstream ss;
		// ss.clear();
		// ss << "R-" << nrank      << "-";
		// ss << "L-" << cur_level0 << "-";
		// ss << "B-" << j;
		// string zonename = ss.str();
		if (a_body.allpoint.size() > 0)
		{
#ifdef QUA_ELEMENT					
			fprintf(fpt, "zone T = \"%s\", NODES=%d, ELEMENTS=%d, DATAPACKING=POINT, ZONETYPE=fequadrilateral\n", 
				zonename, a_body.allpoint.size(), a_body.patch.size());
#endif
#ifdef TRI_ELEMENT
			fprintf(fpt, "zone T = \"%s\", NODES=%d, ELEMENTS=%d, DATAPACKING=POINT, ZONETYPE=fetriangle\n", 
				zonename, a_body.allpoint.size(), a_body.patch.size());
#endif												
		}
		//N=%d, E=%d, F=fepoint, ET=brick								
		for (int i = 0; i < a_body.allpoint.size(); ++i)
		{
			fprintf(fpt, "%16.8f%16.8f%16.8f\n", a_body.allpoint[i][0], 
				a_body.allpoint[i][1], a_body.allpoint[i][2]);
		}		
		for (int i = 0; i < a_body.patch.size(); ++i)					
		{
#ifdef QUA_ELEMENT			
			fprintf(fpt, "%d %d %d %d\n", 
				a_body.patch[i].corpt[0]+1,
				a_body.patch[i].corpt[1]+1,
				a_body.patch[i].corpt[2]+1,
				a_body.patch[i].corpt[3]+1);
#endif
#ifdef TRI_ELEMENT
			fprintf(fpt, "%d %d %d\n", 
				a_body.patch[i].corpt[0]+1,
				a_body.patch[i].corpt[1]+1,
				a_body.patch[i].corpt[2]+1);
#endif							
		}
	}

	void Pltoutput::Print_Body_Soln(Body & a_body, const int & body0)
	{
		char zonename[32];
		sprintf(zonename, "Body-%d-N-%d", body0, node);
		if (a_body.allpoint.size() > 0)
		{
#ifdef QUA_ELEMENT					
			fprintf(fpt, "zone T = \"%s\", NODES=%d, ELEMENTS=%d, DATAPACKING=BLOCK, ZONETYPE=fequadrilateral\n", 
				zonename, a_body.allpoint.size(), a_body.patch.size());
#endif
#ifdef TRI_ELEMENT
			fprintf(fpt, "zone T = \"%s\", NODES=%d, ELEMENTS=%d, DATAPACKING=BLOCK, ZONETYPE=fetriangle\n", 
				zonename, a_body.allpoint.size(), a_body.patch.size());
#endif			
			fprintf(fpt, "VARLOCATION = ([1-%d]=CELLCENTERED)\n", varnum+other_varnum);
		}
		Print_Body_HGVar(a_body, 1);
		for (int i = 0; i < 3; ++i)
		{
			Print_Body_Nmv(a_body, i);	
			//printf("Finish print the %d VARIABLES!!!\n", i);	
		}
		/*------This is surface pressure-------------*/
		Print_Body_HGVar(a_body, 4);
		Print_Body_Yplus(a_body);
		//Print_Body_ReverseFlag(a_body);
		Print_Body_HGVar(a_body, 2);
		Print_Body_HGdis(a_body);
		Print_Body_Bi(a_body);
		Print_Body_Pi(a_body);
		for (int i = 10; i < varnum+other_varnum; ++i)
		{
			Print_Zero_Distance(a_body.patch.size());
		}
		
	}

	void Pltoutput::Print_Level_Soln(Meshlevel & a_level, const int & cur_level0)
	{
		// string outstr;
		// for (int i = 0; i < varnum; ++i)
		// {
		// 	outstr += ffm;
		// }
		// outstr += "\n";
		char zonename[32];
		sprintf(zonename, "N-%d L-%d", node, cur_level0);
		if (a_level.m_point.size() > 0)
		{ 
#ifdef PRINTGHOST
			fprintf(fpt, "zone T = \"%s\", NODES=%d, ELEMENTS=%d, DATAPACKING=BLOCK, ZONETYPE=FEBRICK\n",
				zonename, a_level.m_point.size(), a_level.m_box.realsize()); 
#else
			fprintf(fpt, "zone T = \"%s\", NODES=%d, ELEMENTS=%d, DATAPACKING=BLOCK, ZONETYPE=FEBRICK\n",
				zonename, a_level.m_point.size(), a_level.m_box.size());
#endif
			fprintf(fpt, "VARLOCATION = ([1-%d]=CELLCENTERED)\n", varnum+other_varnum);
		}
		//printf("var number %d other_varnum %d\n", varnum, other_varnum);
		for (int i = 0; i < varnum; ++i)
		{
			Print_a_Flowvar(a_level.m_data, varloc[i]);	
			//printf("Finish print the %d VARIABLES!!!\n", i);	
		}
		Print_Ma(a_level.m_data);
		//printf("Finish print Ma!!!\n");
		Print_Distance(a_level.m_box, cur_level0);
		Print_Domain_Distance(a_level.m_box);
		for (int i = 0; i < 3; ++i)
		{
			Print_Box_Index(a_level.m_box, i);
		}
		Print_Omega_Z(a_level.m_data, a_level.m_box, cur_level0);
		//printf("Finish print distance!!!\n");
		//Print_SolidFlag(a_level.m_box);
	}

	void Pltoutput::Print_IBcell_yplus(string & bfile)
	{
		Create_a_File(bfile);
		int level0 = m_mesh.cur_level_num-1;
		for (int i = 0; i < nodenum; ++i)
		{
			if (node == i && srank == 0)
			{
				fpt = fopen(filename, "a+");
				int boxnum = m_mesh.m_dis.size();
				string dataformat0 = ffm + ffm + ffm + ffm + ffm + "\n";
				for (int b0 = 0; b0 < boxnum; ++b0)
				{
					int ci0 = m_mesh.m_dis[b0].ci;
					if (m_mesh.m_level[level0].m_box.isnormal(ci0) && m_mesh.m_level[level0].m_box[ci0].pair.signdis > 0.0)
					{
						fprintf(fpt, dataformat0.c_str(),
							m_mesh.m_level[level0].m_geom[ci0].boxcenter[0],
							m_mesh.m_level[level0].m_geom[ci0].boxcenter[1],
							m_mesh.m_level[level0].m_geom[ci0].boxcenter[2],
							m_mesh.m_dis[b0].yplus,
							m_mesh.m_dis[b0].uplus);
					}
				}				
				fclose(fpt);
				fpt = NULL;
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
	}

	void Pltoutput::Print_Box_Index(DataArray<Box> & abox, const int & id0)
	{
		int s = 0;
#ifdef PRINTGHOST			
		for (int i = 0; i < abox.realsize(); ++i)
#else
		for (int i = 0; i < abox.size(); ++i)
#endif							
		{
			if (id0 == 0 && false) 
			{
				fprintf(fpt, "%d ", abox[i].ptype);
				// if (abox[i].ptype == 1)
				// {
				// 	fprintf(fpt, "%d ", 100);
				// }
				// else
				// {
				// 	fprintf(fpt, "%d ", -1);
				// }
			}
			else fprintf(fpt, "%d ", abox[i].lowpt.xy[id0]);							
			++s;
			if (s%8 == 0) fprintf(fpt, "\n");
		}
		if (s%8 != 0) fprintf(fpt, "\n");
	}

	void Pltoutput::Print_Omega_Z(DataArray<FlowVariables> & a_data, DataArray<Box> & a_box, const int & cur_level0)
	{
		int s = 0;
		for (int i = 0; i < a_box.size(); ++i)
		{
			int left0 = a_box[i].neib[0][1][1];
			int right0 = a_box[i].neib[2][1][1];
			int up0 = a_box[i].neib[1][2][1];
			int down0 = a_box[i].neib[1][0][1];
			double omegaz = (a_data[right0].v - a_data[left0].v)/(2.0*dh[cur_level0][0]) - 
			                (a_data[up0].u - a_data[down0].u)/(2.0*dh[cur_level0][1]);
			fprintf(fpt, ffm.c_str(), omegaz);
			++s;
			if (s%8 == 0) fprintf(fpt, "\n");
		}
#ifdef PRINTGHOST
		for (int i = abox.size(); i < abox.realsize(); ++i)
		{
			fprintf(fpt, ffm.c_str(), 0.0);
			++s;
			if (s%8 == 0) fprintf(fpt, "\n");
		}
#endif		
		if (s%8 != 0) fprintf(fpt, "\n");
	}
	