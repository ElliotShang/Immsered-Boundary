#include "NS_Solver.H"

	void NS_Solver::ComptNewPrimVar()
	{
		double maxcfl = 0.0;
		double mincfl = CFL;
		//int local_cfl_mdf = 0;
		Point adp[6] = {Point(0,1,1), Point(2,1,1), Point(1,0,1), Point(1,2,1), Point(1,1,0), Point(1,1,2)};
		int level0;
		FlowVariables tempfv;
		if (bodynum > 0)
		{
			level0 = a_mesh.cur_level_num - 1;
		}
		else
		{
			level0 = a_mesh.cur_level_num;
		}
		for (int i = 0; i < level0; ++i)
		{
#ifdef TEMPORAL_REFINE
			if (ts%marching_step[i] == marching_left_step[i])
			{
#endif			
			//MPI_Win_lock_all(0, a_mesh.ShowLevelDataPtr(i)->arraywin());
			int end0 = a_mesh.m_level[i].m_box.pe();
			int start0 = a_mesh.m_level[i].m_box.ps();
			for (int bn = start0; bn < end0; ++bn)
			{
				int b0 = bn - start0;

				if (a_mesh.m_level[i].m_box[bn].pair.signdis > 0.0)
				{
//#ifndef LOCAL_TIME_STEPPING					
					Get_Prim_Vars(a_mesh.m_level[i].m_data[bn], cv_new[i][b0]);
// #else				
// 					Get_Prim_Vars(tempfv, cv_new[i][b0]);
// 					if (tempfv.T < 0.0 || tempfv.roe < 0.0)
// 					{
// 						local_cfl_ismax[i][bn] = true;
// 						local_cfl[i][bn] = local_cfl[i][bn] / dt_ratio;
// 					}
// 					else
// 					{
// 						a_mesh.m_level[i].m_data[bn] = tempfv;
// 					}
// 					maxcfl = max(maxcfl, local_cfl[i][bn]);
// 					mincfl = min(mincfl, local_cfl[i][bn]);
// #endif						
				}
				else 
				{
					a_mesh.m_level[i].m_data[bn] = initvar;
				}
				if (a_mesh.m_level[i].m_data[bn].T < 0.0)
				{
					printf("Level %d Box (%d,%d,%d) signdis is %f T is %f box v is %f after Get_Prim_Vars!!!\n",
						i,
						a_mesh.m_level[i].m_box[bn].ix(),
						a_mesh.m_level[i].m_box[bn].iy(),
						a_mesh.m_level[i].m_box[bn].iz(),
						a_mesh.m_level[i].m_box[bn].pair.signdis,
						a_mesh.m_level[i].m_data[bn].T,
						a_mesh.m_level[i].m_geom[bn].v);
					MPI_Abort(MPI_COMM_WORLD, 37);
				}				
#ifdef DEBUG				
				if (a_mesh.LevelBoxData(i, bn).hasnan(i,bn,"new prim check"))
				{
					printf("Level %d Box (%d,%d,%d) has nan value!!!\n", 
						i, 
						a_mesh.m_level[i].m_box[bn].ix(),
						a_mesh.m_level[i].m_box[bn].iy(),
						a_mesh.m_level[i].m_box[bn].iz());
				}
				if (a_mesh.LevelBoxData(i, bn).roe < 0.0 || a_mesh.LevelBoxData(i, bn).T < 0.0)
				{
					printf("Level %d Box (%d,%d,%d) invalid values density %f u %f v %f w %f p %f T %f!!!\n",
						i, 
						a_mesh.m_level[i].m_box[bn].ix(),
						a_mesh.m_level[i].m_box[bn].iy(),
						a_mesh.m_level[i].m_box[bn].iz(),
						a_mesh.m_level[i].m_data[bn].roe,
						a_mesh.m_level[i].m_data[bn].u,
						a_mesh.m_level[i].m_data[bn].v,
						a_mesh.m_level[i].m_data[bn].w,
						a_mesh.m_level[i].m_data[bn].p,
						a_mesh.m_level[i].m_data[bn].T);
					for (Point_iterator p(0,3); p.end(); ++p)
					{
						int ab1 = a_mesh.m_level[i].m_box[bn].neib[p.i][p.j][p.k];
						printf("Level %d Box (%d,%d,%d) invalid values density %f u %f v %f w %f p %f T %f!!!\n",
							i, 
							a_mesh.m_level[i].m_box[ab1].ix(),
							a_mesh.m_level[i].m_box[ab1].iy(),
							a_mesh.m_level[i].m_box[ab1].iz(),
							a_mesh.m_level[i].m_data[ab1].roe,
							a_mesh.m_level[i].m_data[ab1].u,
							a_mesh.m_level[i].m_data[ab1].v,
							a_mesh.m_level[i].m_data[ab1].w,
							a_mesh.m_level[i].m_data[ab1].p,
							a_mesh.m_level[i].m_data[ab1].T);
					}
					MPI_Abort(MPI_COMM_WORLD, 46);
				}
				Assert(!a_mesh.LevelBoxData(i, bn).hasnan(i,bn,"new prim check"), 
					"new prim check", 152);
#endif											
			}
#ifdef TEMPORAL_REFINE
			}
#endif			
		}
		if (bodynum > 0)
		{
			int end0 = a_mesh.m_level[level0].m_box.pe();
			int start0 = a_mesh.m_level[level0].m_box.ps();
			for (int bn = start0; bn < end0; ++bn)
			{
				if (a_mesh.m_level[level0].m_box[bn].pair.signdis > 0.0 && a_mesh.infectbox[bn] == -1)									
				{
					int b0 = bn - start0;
// #ifndef LOCAL_TIME_STEPPING					
					Get_Prim_Vars(a_mesh.m_level[level0].m_data[bn], cv_new[level0][b0]);
// #else 
// 					Get_Prim_Vars(tempfv, cv_new[level0][b0]);
// 					if (tempfv.T < 0.0 || tempfv.roe < 0.0)
// 					{
// 						local_cfl_ismax[level0][bn] = true;
// 						//local_cfl[level0][bn] = local_cfl[level0][bn] / dt_ratio;
// 					}
// 					else
// 					{
// 						a_mesh.m_level[level0].m_data[bn] = tempfv;
// 					}
// // 					maxcfl = max(maxcfl, local_cfl[level0][bn]);
// // 					mincfl = min(mincfl, local_cfl[level0][bn]);
// #endif				
					if (a_mesh.m_level[level0].m_data[bn].T < 0.0)
					{
						printf("Level %d Box (%d,%d,%d) signdis is %f T is %f u %f v %f w %f roe %f p %f v %20.10f after Get_Prim_Vars!!!\n",
							level0,
							a_mesh.m_level[level0].m_box[bn].ix(),
							a_mesh.m_level[level0].m_box[bn].iy(),
							a_mesh.m_level[level0].m_box[bn].iz(),
							a_mesh.m_level[level0].m_box[bn].pair.signdis,
							a_mesh.m_level[level0].m_data[bn].T,
							a_mesh.m_level[level0].m_data[bn].u,
							a_mesh.m_level[level0].m_data[bn].v,
							a_mesh.m_level[level0].m_data[bn].w,
							a_mesh.m_level[level0].m_data[bn].roe,
							a_mesh.m_level[level0].m_data[bn].p,
							a_mesh.m_level[level0].m_geom[bn].v);
						MPI_Abort(MPI_COMM_WORLD, 37);
					}	
								
#ifdef DEBUG				
					if (a_mesh.LevelBoxData(level0, bn).hasnan(level0,bn,"new prim check"))
					{
						printf("Level %d Box (%d,%d,%d) has nan value!!!\n", 
							level0, 
							a_mesh.m_level[level0].m_box[bn].ix(),
							a_mesh.m_level[level0].m_box[bn].iy(),
							a_mesh.m_level[level0].m_box[bn].iz());
					}					
					Assert(!a_mesh.LevelBoxData(level0, bn).hasnan(level0,bn,"new prim check"), 
						"new prim check", 152);
#endif					
				}
				else
				{
					if (a_mesh.infectbox[bn] == -1)
					{
						a_mesh.m_level[level0].m_data[bn] = initvar;
					}	
				}
			}
			//MPI_Win_unlock_all(a_mesh.ShowLevelDataPtr(level0)->arraywin());
		}
// #ifdef LOCAL_TIME_STEPPING		
// 		//MPI_Allreduce(MPI_IN_PLACE, &local_cfl_mdf, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
// 		MPI_Allreduce(MPI_IN_PLACE, &mincfl, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
// 		MPI_Allreduce(MPI_IN_PLACE, &maxcfl, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
// 		if (nrank == 0)
// 		{
// 			printf("Times step %d MaxCFL is %f MinCFL is %f!!!\n", ts, maxcfl, mincfl);
// 		}
// 		// if (local_cfl_mdf)
// 		// {
// 		// 	InitFlow();
// 		// }
// #endif		
		MPI_Barrier(MPI_COMM_WORLD);
	}

	void NS_Solver::Comptdt()
	{
		CFL_dt = 10000.0;
#ifdef LOCAL_TIME_STEPPING
		if (ts == 0)
		{
			for (int i = 0; i < a_mesh.cur_level_num; ++i)
			{
				int start0 = a_mesh.LevelBoxStart(i);
				int end0 = a_mesh.LevelBoxEnd(i);
				for (int b0 = start0; b0 < end0; ++b0)
				{
					local_cfl[i][b0] = CFL;
					local_cfl_ismax[i][b0] = false;
				}
			}
		}
		// else
		// {
		// 	for (int i = 0; i < a_mesh.cur_level_num; ++i)
		// 	{
		// 		int start0 = a_mesh.LevelBoxStart(i);
		// 		int end0 = a_mesh.LevelBoxEnd(i);
		// 		for (int b0 = start0; b0 < end0; ++b0)
		// 		{
		// 			if (!local_cfl_ismax[i][b0])
		// 			{
		// 				local_cfl[i][b0] *= dt_ratio;
		// 			}
		// 		}
		// 	}-
		// }
#endif					
		double convc[3], localsoundspeed, avearea;
		Pointxyz avenmv;
		for (int i = 0; i < a_mesh.cur_level_num; ++i)
		{
			int start0 = a_mesh.LevelBoxStart(i);
			int end0 = a_mesh.LevelBoxEnd(i);
			for (int b0 = start0; b0 < end0; ++b0)
			{
				localsoundspeed = sqrt(gama*Rg*a_mesh.m_level[i].m_data[b0].T);
				for (int di = 0; di < 3; ++di)
				{
					int f1 = a_mesh.m_level[i].m_box[b0].faces[di][0];
					int f2 = a_mesh.m_level[i].m_box[b0].faces[di][1];
					avenmv = (a_mesh.m_level[i].m_face[f1].keisa+a_mesh.m_level[i].m_face[f2].keisa)*0.5;
					avearea = 0.5*(a_mesh.m_level[i].m_face[f1].area+a_mesh.m_level[i].m_face[f2].area);
					Pointxyz localvel = Pointxyz(a_mesh.m_level[i].m_data[b0].u,
												 a_mesh.m_level[i].m_data[b0].v,
												 a_mesh.m_level[i].m_data[b0].w);
					convc[di] = (abs(localvel.dot(avenmv))+localsoundspeed)*avearea;
				}
				double CFL_dt0 = CFL*a_mesh.m_level[i].m_geom[b0].v/(convc[0]+convc[1]+convc[2]);
#ifdef LOCAL_TIME_STEPPING				
				local_dt[i][b0] = local_cfl[i][b0]/CFL*CFL_dt0;
#endif				
				if (CFL_dt0 < 0.0)
				{
					printf("CFL_dt is %f which is negative!!! convc 1 %f convc 2 %f convc 3 %f localsoundspeed %f v %f!!!\n",
						CFL_dt0, convc[0], convc[1], convc[2], localsoundspeed, a_mesh.m_level[i].m_geom[b0].v);
					MPI_Abort(MPI_COMM_WORLD, 424);
				}
				CFL_dt = min(CFL_dt, CFL_dt0);
			}
		}
		MPI_Allreduce(MPI_IN_PLACE, &CFL_dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
		dt = min(CFL_dt, dt*dt_ratio);
		if (nrank == 0)
		{
			FILE * timefile;
			if (ts == 0)
			{
				timefile = fopen("timestep.dat", "w");
			}
			else
			{
				timefile = fopen("timestep.dat", "a");
			}
			fprintf(timefile, "%20d %20.10f %20.10f %20.10f\n",
				ts, t, dt, CFL_dt);
			fclose(timefile);
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}

	void NS_Solver::ComptSource()
	{
		FlowVec sincrem[3];
		for (int level_n = 0; level_n < a_mesh.cur_level_num; ++level_n)
		{
#ifdef TEMPORAL_REFINE
			if (ts%marching_step[level_n] == marching_left_step[level_n])
			{
#endif			
			//MPI_Win_lock_all(0, faceflux[level_n].arraywin());
			//MPI_Win_lock_all(0, vsflux[level_n].arraywin());
			int start0 = a_mesh.m_level[level_n].m_box.ps();
			int end0 = a_mesh.m_level[level_n].m_box.pe();
			for (int bn = start0; bn < end0; ++bn)
			{
				// if (a_mesh.m_level[level_n].m_box[bn].pair.signdis > 0.0)
				// {
				int b0 = bn-start0;

				Pointxyz flow_vel(a_mesh.m_level[level_n].m_data[bn].u,
								  a_mesh.m_level[level_n].m_data[bn].v,
								  a_mesh.m_level[level_n].m_data[bn].w);
				Pointxyz corlis = frame_rotate_vel.cross(flow_vel - frame_trans_vel);				
				for (int fi = 0; fi < DIM; ++fi)
				{
					const int & bf1 = a_mesh.m_level[level_n].m_box[bn].faces[fi][0];
					const int & bf2 = a_mesh.m_level[level_n].m_box[bn].faces[fi][1];
					Face & face1 = a_mesh.m_level[level_n].m_face[bf1];
					Face & face2 = a_mesh.m_level[level_n].m_face[bf2];
// #ifdef VISCOUSITY
// 					sincrem[fi] += vsflux[level_n][bf1];
// 					sincrem[fi] -= vsflux[level_n][bf2];
// #endif	
					Assert((bf1>-1 && bf2<a_mesh.Meshface(level_n).size()), "face index 0 error", 229);
					Assert((bf1>-1 && bf2<a_mesh.Meshface(level_n).size()), "face index 1 error", 229);	

#if DIM == 3									
					sincrem[fi] = faceflux[level_n][bf1]-faceflux[level_n][bf2];
#elif DIM == 2
					sincrem[fi][0] = faceflux[level_n][bf1][0]-faceflux[level_n][bf2][0];					
					sincrem[fi][1] = faceflux[level_n][bf1][1]-faceflux[level_n][bf2][1];					
					sincrem[fi][2] = faceflux[level_n][bf1][2]-faceflux[level_n][bf2][2];					
					sincrem[fi][4] = faceflux[level_n][bf1][4]-faceflux[level_n][bf2][4];
#elif DIM == 1
					sincrem[fi][0] = faceflux[level_n][bf1][0]-faceflux[level_n][bf2][0];					
					sincrem[fi][1] = faceflux[level_n][bf1][1]-faceflux[level_n][bf2][1];					
					sincrem[fi][4] = faceflux[level_n][bf1][4]-faceflux[level_n][bf2][4];
#endif

#ifdef VISCOUSITY
	#if DIM == 3
					sincrem[fi] += vsflux[level_n][bf2] - vsflux[level_n][bf1];
	#elif DIM == 2
					sincrem[fi][0] += vsflux[level_n][bf2][0] - vsflux[level_n][bf1][0];
					sincrem[fi][1] += vsflux[level_n][bf2][1] - vsflux[level_n][bf1][1];
					sincrem[fi][2] += vsflux[level_n][bf2][2] - vsflux[level_n][bf1][2];
					sincrem[fi][4] += vsflux[level_n][bf2][4] - vsflux[level_n][bf1][4];
	#elif DIM == 1
					sincrem[fi][0] += vsflux[level_n][bf2][0] - vsflux[level_n][bf1][0];
					sincrem[fi][1] += vsflux[level_n][bf2][1] - vsflux[level_n][bf1][1];
					sincrem[fi][4] += vsflux[level_n][bf2][4] - vsflux[level_n][bf1][4];
	#endif				
#endif
#ifdef DEBUG					
					if (sincrem[fi].hasnan(level_n, bn, fi)) MPI_Abort(MPI_COMM_WORLD, 323);
#endif
				}
				source[0][level_n][b0] = sincrem[0]+sincrem[1]+sincrem[2];								
				source[1][level_n][b0][0] = 0.0;
				source[1][level_n][b0][1] = -corlis[0]*a_mesh.m_level[level_n].m_data[bn].roe;
				source[1][level_n][b0][2] = -corlis[1]*a_mesh.m_level[level_n].m_data[bn].roe;
				source[1][level_n][b0][3] = -corlis[2]*a_mesh.m_level[level_n].m_data[bn].roe;
				source[1][level_n][b0][4] = 0.0;
				Assert(!source[0][level_n][b0].hasnan(level_n, bn, 0), "source 0 check", 245);
				}
			// }
#ifdef TEMPORAL_REFINE
			}
#endif			
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}

	void NS_Solver::ComptViscousStress()
	{
		int neibe[2] = {3,2};
		int ddn[2] = {0,1};
		double facevel[3];
		double dux, duy, duz, dvx, dvy, dvz, dwx, dwy, dwz;
		CellMCoef fk0;
		double duk0, duk1, duk2, dvk0, dvk1, dvk2, dwk0, dwk1, dwk2;
		double dTk0, dTk1, dTk2, dTx, dTy, dTz;
		double tauxx, tauxy, tauxz, tauyy, tauyz, tauzz;
		double qx, qy, qz, bbbx, bbby, bbbz, mu0, qk0;
#ifdef TURBULENCE
		double dvisk0, dvisk1, dvisk2;
		double dvisx, dvisy, dvisz;
#endif				
		Point facedir(0,0,0);
		int * faceptr = &facedir[0];
#ifdef SHOWTIME		
		double starttime = MPI_Wtime();
#endif		
		for (int level_n = 0; level_n < a_mesh.MyCurNum(); ++level_n)
		{
			MPI_Barrier(MPI_COMM_WORLD);
#ifdef TEMPORAL_REFINE
			if (ts%marching_step[level_n] == marching_left_step[level_n])
			{
#endif			
			//MPI_Win_lock_all(0, vsflux[level_n].arraywin());
			//MPI_Win_fence(0, vsflux[level_n].arraywin());
			int end0 = a_mesh.m_level[level_n].m_face.pe();
			//printf("N%dR%dCompute viscous flux step %d for L%d Face start %d face end %d\n", node, srank, ts, level_n, a_mesh.Facestart(level_n), end0);
			for (int bn = a_mesh.Facestart(level_n); bn < end0; ++bn)
			{
				Face & theface = a_mesh.m_level[level_n].m_face[bn];
#if DIM == 2
				if (theface.fnv != 2)
				{
#endif					
				faceptr[theface.fnv] = 1;
				const int & abox = theface[0];		
				Assert(abox > -1, "Error In identifying the face side box in compute vs!!!", 311);
				const int neibx = 1+ddn[faceptr[0]];
				const int neiby = 1+ddn[faceptr[1]];
				const int neibz = 1+ddn[faceptr[2]];

				CellMCoef * vsmcptr[2] = {&a_mesh.m_level[level_n].m_geom[theface[0]].keisa, 
										  &a_mesh.m_level[level_n].m_geom[theface[1]].keisa};
				
				for (int di = 0; di < neibe[faceptr[0]]; ++di)
				{
					for (int dj = 0; dj < neibe[faceptr[1]]; ++dj)
					{
						for (int dk = 0; dk < neibe[faceptr[2]]; ++dk)
						{
							Point rp(di+ddn[faceptr[0]], dj+ddn[faceptr[1]], dk+ddn[faceptr[2]]);
							int aneib0 = a_mesh.m_level[level_n].m_box[abox].neib[rp[0]][rp[1]][rp[2]];
							Assert(aneib0 > -1, "Negative index for viscous computation!!!", 425);
							vsflowptr[rp[0]][rp[1]][rp[2]] = &a_mesh.m_level[level_n].m_data[aneib0];
						}
					}
				}
				const int lbi = 1;
				const int rbi = 2;				
				//printf("Compute viscous flux 3.2 for L%dF%d lbi%drbi%d\n", level_n, bn, lbi, rbi);
				if (theface.fnv == 0)
				{	
					pfpkesai0(duk0, vsflowptr[lbi][1][1]->u, vsflowptr[rbi][1][1]->u);
					pfpkesai0(dvk0, vsflowptr[lbi][1][1]->v, vsflowptr[rbi][1][1]->v);
					pfpkesai0(dwk0, vsflowptr[lbi][1][1]->w, vsflowptr[rbi][1][1]->w);
					pfpkesai0(dTk0, vsflowptr[lbi][1][1]->T, vsflowptr[rbi][1][1]->T);

					pfpkesai1(duk1, vsflowptr[lbi][2][1]->u, vsflowptr[rbi][2][1]->u, 
						vsflowptr[lbi][0][1]->u, vsflowptr[rbi][0][1]->u);
					pfpkesai1(dvk1, vsflowptr[lbi][2][1]->v, vsflowptr[rbi][2][1]->v, 
						vsflowptr[lbi][0][1]->v, vsflowptr[rbi][0][1]->v);
					pfpkesai1(dwk1, vsflowptr[lbi][2][1]->w, vsflowptr[rbi][2][1]->w, 
						vsflowptr[lbi][0][1]->w, vsflowptr[rbi][0][1]->w);
					pfpkesai1(dTk1, vsflowptr[lbi][2][1]->T, vsflowptr[rbi][2][1]->T, 
						vsflowptr[lbi][0][1]->T, vsflowptr[rbi][0][1]->T);

					pfpkesai1(duk2, vsflowptr[lbi][1][2]->u, vsflowptr[rbi][1][2]->u, 
						vsflowptr[lbi][1][0]->u, vsflowptr[rbi][1][0]->u);
					pfpkesai1(dvk2, vsflowptr[lbi][1][2]->v, vsflowptr[rbi][1][2]->v, 
						vsflowptr[lbi][1][0]->v, vsflowptr[rbi][1][0]->v);
					pfpkesai1(dwk2, vsflowptr[lbi][1][2]->w, vsflowptr[rbi][1][2]->w, 
						vsflowptr[lbi][1][0]->w, vsflowptr[rbi][1][0]->w);
					pfpkesai1(dTk2, vsflowptr[lbi][1][2]->T, vsflowptr[rbi][1][2]->T, 
						vsflowptr[lbi][1][0]->T, vsflowptr[rbi][1][0]->T);
#ifdef TURBULENCE
					pfpkesai0(dvisk0, vsflowptr[lbi][1][1]->viseddy, vsflowptr[rbi][1][1]->viseddy);
					pfpkesai1(dvisk1, vsflowptr[lbi][2][1]->viseddy, vsflowptr[rbi][2][1]->viseddy, 
						vsflowptr[lbi][0][1]->viseddy, vsflowptr[rbi][0][1]->viseddy);
					pfpkesai1(dvisk2, vsflowptr[lbi][1][2]->viseddy, vsflowptr[rbi][1][2]->viseddy, 
						vsflowptr[lbi][1][0]->viseddy, vsflowptr[rbi][1][0]->viseddy);
#endif								
				}
				else if (theface.fnv == 1)
				{
					pfpkesai0(duk1, vsflowptr[1][lbi][1]->u, vsflowptr[1][rbi][1]->u);
					pfpkesai0(dvk1, vsflowptr[1][lbi][1]->v, vsflowptr[1][rbi][1]->v);
					pfpkesai0(dwk1, vsflowptr[1][lbi][1]->w, vsflowptr[1][rbi][1]->w);
					pfpkesai0(dTk1, vsflowptr[1][lbi][1]->T, vsflowptr[1][rbi][1]->T);

					pfpkesai1(duk0, vsflowptr[2][lbi][1]->u, vsflowptr[2][rbi][1]->u, 
						vsflowptr[0][lbi][1]->u, vsflowptr[0][rbi][1]->u);
					pfpkesai1(dvk0, vsflowptr[2][lbi][1]->v, vsflowptr[2][rbi][1]->v, 
						vsflowptr[0][lbi][1]->v, vsflowptr[0][rbi][1]->v);
					pfpkesai1(dwk0, vsflowptr[2][lbi][1]->w, vsflowptr[2][rbi][1]->w, 
						vsflowptr[0][lbi][1]->w, vsflowptr[0][rbi][1]->w);
					pfpkesai1(dTk0, vsflowptr[2][lbi][1]->T, vsflowptr[2][rbi][1]->T, 
						vsflowptr[0][lbi][1]->T, vsflowptr[0][rbi][1]->T);

					pfpkesai1(duk2, vsflowptr[1][lbi][2]->u, vsflowptr[1][rbi][2]->u, 
						vsflowptr[1][lbi][0]->u, vsflowptr[1][rbi][0]->u);
					pfpkesai1(dvk2, vsflowptr[1][lbi][2]->v, vsflowptr[1][rbi][2]->v, 
						vsflowptr[1][lbi][0]->v, vsflowptr[1][rbi][0]->v);
					pfpkesai1(dwk2, vsflowptr[1][lbi][2]->w, vsflowptr[1][rbi][2]->w, 
						vsflowptr[1][lbi][0]->w, vsflowptr[1][rbi][0]->w);
					pfpkesai1(dTk2, vsflowptr[1][lbi][2]->T, vsflowptr[1][rbi][2]->T, 
						vsflowptr[1][lbi][0]->T, vsflowptr[1][rbi][0]->T);

#ifdef TURBULENCE
					pfpkesai0(dvisk1, vsflowptr[1][lbi][1]->viseddy, vsflowptr[1][rbi][1]->viseddy);
					pfpkesai1(dvisk0, vsflowptr[2][lbi][1]->viseddy, vsflowptr[2][rbi][1]->viseddy, 
						vsflowptr[0][lbi][1]->viseddy, vsflowptr[0][rbi][1]->viseddy);
					pfpkesai1(dvisk2, vsflowptr[1][lbi][2]->viseddy, vsflowptr[1][rbi][2]->viseddy, 
						vsflowptr[1][lbi][0]->viseddy, vsflowptr[1][rbi][0]->viseddy);
#endif					
				}
				else if (theface.fnv == 2)
				{
					pfpkesai0(duk2, vsflowptr[1][1][lbi]->u, vsflowptr[1][1][rbi]->u);
					pfpkesai0(dvk2, vsflowptr[1][1][lbi]->v, vsflowptr[1][1][rbi]->v);
					pfpkesai0(dwk2, vsflowptr[1][1][lbi]->w, vsflowptr[1][1][rbi]->w);
					pfpkesai0(dTk2, vsflowptr[1][1][lbi]->T, vsflowptr[1][1][rbi]->T);

					pfpkesai1(duk0, vsflowptr[2][1][lbi]->u, vsflowptr[2][1][rbi]->u, 
						vsflowptr[0][1][lbi]->u, vsflowptr[0][1][rbi]->u);
					pfpkesai1(dvk0, vsflowptr[2][1][lbi]->v, vsflowptr[2][1][rbi]->v, 
						vsflowptr[0][1][lbi]->v, vsflowptr[0][1][rbi]->v);
					pfpkesai1(dwk0, vsflowptr[2][1][lbi]->w, vsflowptr[2][1][rbi]->w, 
						vsflowptr[0][1][lbi]->w, vsflowptr[0][1][rbi]->w);
					pfpkesai1(dTk0, vsflowptr[2][1][lbi]->T, vsflowptr[2][1][rbi]->T, 
						vsflowptr[0][1][lbi]->T, vsflowptr[0][1][rbi]->T);

					pfpkesai1(duk1, vsflowptr[1][2][lbi]->u, vsflowptr[1][2][rbi]->u, 
						vsflowptr[1][0][lbi]->u, vsflowptr[1][0][rbi]->u);
					pfpkesai1(dvk1, vsflowptr[1][2][lbi]->v, vsflowptr[1][2][rbi]->v, 
						vsflowptr[1][0][lbi]->v, vsflowptr[1][0][rbi]->v);
					pfpkesai1(dwk1, vsflowptr[1][2][lbi]->w, vsflowptr[1][2][rbi]->w, 
						vsflowptr[1][0][lbi]->w, vsflowptr[1][0][rbi]->w);
					pfpkesai1(dTk1, vsflowptr[1][2][lbi]->T, vsflowptr[1][2][rbi]->T, 
						vsflowptr[1][0][lbi]->T, vsflowptr[1][0][rbi]->T);
#ifdef TURBULENCE
					pfpkesai0(dvisk2, vsflowptr[1][1][lbi]->viseddy, vsflowptr[1][1][rbi]->viseddy);
					pfpkesai1(dvisk0, vsflowptr[2][1][lbi]->viseddy, vsflowptr[2][1][rbi]->viseddy, 
						vsflowptr[0][1][lbi]->viseddy, vsflowptr[0][1][rbi]->viseddy);
					pfpkesai1(dvisk1, vsflowptr[1][2][lbi]->viseddy, vsflowptr[1][2][rbi]->viseddy, 
						vsflowptr[1][0][lbi]->viseddy, vsflowptr[1][0][rbi]->viseddy);
#endif					

				}
				fk0[0] = (vsmcptr[0]->var[0]+vsmcptr[1]->var[0])*0.5;
				fk0[1] = (vsmcptr[0]->var[1]+vsmcptr[1]->var[1])*0.5;
				fk0[2] = (vsmcptr[0]->var[2]+vsmcptr[1]->var[2])*0.5;
				VarFromKesaitoX();
				double ftemp = (vsflowptr[1][1][1]->T+vsflowptr[neibx][neiby][neibz]->T)*0.5;
				mu0 = sqrt(pow(ftemp, 3))*((1.0+S_over_T_ref)/(ftemp+S_over_T_ref));
				qk0 = qk*mu0; // qk0 = qk/mu*mu0 mu = 1.0;
#ifdef TURBULENCE
				double fvis = (vsflowptr[1][1][1]->viseddy+vsflowptr[neibx][neiby][neibz]->viseddy)*0.5;
				double froe = (vsflowptr[1][1][1]->roe+vsflowptr[neibx][neiby][neibz]->roe)*0.5;
				double faceroexvis = fvis + mu0/froe;
				double mut0 = Viseddy2mut(fvis, froe, mu0); // 此处返回涡粘性量，应该在壁面速度更新处引入ODE求解
				mu0 += mut0;
				qk0 += mut0/(Re*Pr_tur*(gama-1)*pow(Ma,2));
				mutflux[level_n][bn] = faceroexvis*(dvisx*theface.keisa[0]+dvisy*theface.keisa[1]+dvisz*theface.keisa[2])/segma;
				mutflux[level_n][bn] *= a_mesh.m_level[level_n].m_face[bn].area;
#endif				
				tauxx = 2.0*mu0*(dux-(dux+dvy+dwz)/3.0)/Re;
				tauyy = 2.0*mu0*(dvy-(dux+dvy+dwz)/3.0)/Re;
				tauzz = 2.0*mu0*(dwz-(dux+dvy+dwz)/3.0)/Re;
				tauxy = mu0/Re*(duy+dvx);
				tauxz = mu0/Re*(duz+dwx);
				tauyz = mu0/Re*(dvz+dwy);
				// if (abs(tauzz) > 0.0000001)
				// {
				// 	printf("L%dF%d %f %f %f\n",level_n,bn,tauxx,tauyy,tauzz);
				// }
				// printf("L%dF%d (%f,%f,%f,%f,%f,%f)\n",level_n,bn,tauxx,tauyy,tauzz,tauxy,tauyz,tauxz);
				qx = qk0*dTx;
				qy = qk0*dTy;
				qz = qk0*dTz;

				facevel[0] = (vsflowptr[1][1][1]->u+vsflowptr[neibx][neiby][neibz]->u)*0.5;
				facevel[1] = (vsflowptr[1][1][1]->v+vsflowptr[neibx][neiby][neibz]->v)*0.5;
				facevel[2] = (vsflowptr[1][1][1]->w+vsflowptr[neibx][neiby][neibz]->w)*0.5;
				//printf("L%dF%dD%d (%f,%f,%f)\n",level_n,bn,theface.fnv,facevel[0],facevel[1],facevel[2]);
				// if (abs(dux) > 0 || abs(duy) > 0 || abs(duz) > 0)
				// {
				// 	//printf("L%dF%dD%d (%f,%f,%f)\n",level_n,bn,theface.fnv,dux,duy,duz);
				// 	//printf("L%dF%dD%d (%f,%f,%f)\n",level_n,bn,theface.fnv,facevel[0],facevel[1],facevel[2]);
				// }

				bbbx = facevel[0]*tauxx+facevel[1]*tauxy+facevel[2]*tauxz+qx;
				bbby = facevel[0]*tauxy+facevel[1]*tauyy+facevel[2]*tauyz+qy;
				bbbz = facevel[0]*tauxz+facevel[1]*tauyz+facevel[2]*tauzz+qz;

				vsflux[level_n][bn][0] = 0.0;
				vsflux[level_n][bn][1] = (tauxx*theface.keisa[0]+tauxy*theface.keisa[1]+tauxz*theface.keisa[2])*
					a_mesh.m_level[level_n].m_face[bn].area;
				vsflux[level_n][bn][2] = (tauxy*theface.keisa[0]+tauyy*theface.keisa[1]+tauyz*theface.keisa[2])*
					a_mesh.m_level[level_n].m_face[bn].area;
				vsflux[level_n][bn][3] = (tauxz*theface.keisa[0]+tauyz*theface.keisa[1]+tauzz*theface.keisa[2])*
					a_mesh.m_level[level_n].m_face[bn].area;
				vsflux[level_n][bn][4] = (bbbx*theface.keisa[0]+bbby*theface.keisa[1]+bbbz*theface.keisa[2])*
					a_mesh.m_level[level_n].m_face[bn].area;							
#ifdef DEBUG				
				if (bool flux_right = vsflux[level_n][bn].hasnan(level_n, bn, 
					a_mesh.Levelface(level_n, bn).fnv))
				{
					// printf("keisa is 0 (%f,%f,%f) 1 (%f,%f,%f) 2 (%f,%f,%f) face (%f,%f,%f)\n", 
					// 	fk0[0][0],fk0[0][1],fk0[0][2],
					// 	fk0[1][0],fk0[1][1],fk0[1][2],
					// 	fk0[2][0],fk0[2][1],fk0[2][2], theface.keisa[0], theface.keisa[1], theface.keisa[2]);
					for (int di = 0; di < neibe[faceptr[0]]; ++di)
					{
						for (int dj = 0; dj < neibe[faceptr[1]]; ++dj)
						{
							for (int dk = 0; dk < neibe[faceptr[2]]; ++dk)
							{
								Point rp(di+ddn[faceptr[0]], dj+ddn[faceptr[1]], dk+ddn[faceptr[2]]);								
								int aneib0 = a_mesh.m_level[level_n].m_box[abox].neib[rp[0]][rp[1]][rp[2]];
								printf("L%dB%d (%d,%d,%d) box center is (%f,%f,%f)\n", level_n, aneib0, rp[0],rp[1],rp[2],
									a_mesh.bc(level_n, aneib0)[0],a_mesh.bc(level_n, aneib0)[1],a_mesh.bc(level_n, aneib0)[2]);
								vsflowptr[rp[0]][rp[1]][rp[2]]->showdata("vsflow data");
								if (a_mesh.isghostcell(level_n,aneib0))
								{
									printf("Ghost cell L%dB%d at neib (%d, %d, %d)\n",level_n,aneib0,rp[0],rp[1],rp[2]);
								}
							}
						}
					}					
					MPI_Abort(MPI_COMM_WORLD, 552);
				}
#endif	
				faceptr[theface.fnv] = 0;
#if DIM == 2
				}
#endif									
			}
#ifdef TEMPORAL_REFINE
			}
#endif			
			//MPI_Win_unlock_all(vsflux[level_n].arraywin());
			//MPI_Win_fence(0,vsflux[level_n].arraywin());
		}
		MPI_Barrier(MPI_COMM_WORLD);
		//BalanceFaceFlux(a_mesh.mdyface, vsflux);
#ifdef SHOWTIME		
		double endtime = MPI_Wtime();
		step_vs_time = endtime -starttime;
		total_vs_time += step_vs_time;
#endif	
	}

#ifdef TURBULENCE
	void NS_Solver::ComptMutSource()
	{
#ifdef SHOWTIME
		double starttime = MPI_Wtime();
#endif				
		double dvisk0, dvisk1, dvisk2;
		double dvisx, dvisy, dvisz;
		double droek0, droek1, droek2;
		double droex, droey, droez;
		int xleft, xright, yleft, yright, zleft, zright;
		double duy, duz, dvx, dvz, dwx, dwy;
		double duk0, duk1, duk2, dvk0, dvk1, dvk2, dwk0, dwk1, dwk2;
		double bigx, fv1, fv2, bigs, r_ratio, g_ratio, ft2, fw;
		Pointxyz omiga;
		double vis0, vis2d, gsix, cw3six, rsix, total_face_flux, segdis;
		double dirsign[2] = {-1.0, 1.0};
		for (int level_n = 0; level_n < a_mesh.cur_level_num; ++level_n)
		{
#ifdef TEMPORAL_REFINE
			if (ts%marching_step[level_n] == marching_left_step[level_n])
			{
#endif			
			int start0 = a_mesh.m_level[level_n].m_box.ps();
			int end0 = a_mesh.m_level[level_n].m_box.pe();
			for (int bn = start0; bn < end0; ++bn)
			{
				if (!a_mesh.m_level[level_n].m_box[bn].solid)
				{					
					int b0 = bn-start0;
					xleft = a_mesh.m_level[level_n].m_box[bn].neib[0][1][1];
					xright = a_mesh.m_level[level_n].m_box[bn].neib[2][1][1];
					yleft = a_mesh.m_level[level_n].m_box[bn].neib[1][0][1];
					yright = a_mesh.m_level[level_n].m_box[bn].neib[1][2][1];
					zleft = a_mesh.m_level[level_n].m_box[bn].neib[1][1][0];
					zright = a_mesh.m_level[level_n].m_box[bn].neib[1][1][2];
					dvisk0 = (a_mesh.m_level[level_n].m_data[xright].viseddy - 
								a_mesh.m_level[level_n].m_data[xleft].viseddy)/2.0;
					dvisk1 = (a_mesh.m_level[level_n].m_data[yright].viseddy - 
								a_mesh.m_level[level_n].m_data[yleft].viseddy)/2.0;
					dvisk2 = (a_mesh.m_level[level_n].m_data[zright].viseddy - 
								a_mesh.m_level[level_n].m_data[zleft].viseddy)/2.0;
					VarReflect(dvisx, a_mesh.m_level[level_n].m_geom[bn].keisa, dvisk0, dvisk1, dvisk2, 0);
					VarReflect(dvisy, a_mesh.m_level[level_n].m_geom[bn].keisa, dvisk0, dvisk1, dvisk2, 1);
					VarReflect(dvisz, a_mesh.m_level[level_n].m_geom[bn].keisa, dvisk0, dvisk1, dvisk2, 2);
					duk0 = (a_mesh.m_level[level_n].m_data[xright].u - 
								a_mesh.m_level[level_n].m_data[xleft].u)/2.0;
					duk1 = (a_mesh.m_level[level_n].m_data[yright].u - 
								a_mesh.m_level[level_n].m_data[yleft].u)/2.0;
					duk2 = (a_mesh.m_level[level_n].m_data[zright].u - 
								a_mesh.m_level[level_n].m_data[zleft].u)/2.0;
					dvk0 = (a_mesh.m_level[level_n].m_data[xright].v - 
								a_mesh.m_level[level_n].m_data[xleft].v)/2.0;
					dvk1 = (a_mesh.m_level[level_n].m_data[yright].v - 
								a_mesh.m_level[level_n].m_data[yleft].v)/2.0;
					dvk2 = (a_mesh.m_level[level_n].m_data[zright].v - 
								a_mesh.m_level[level_n].m_data[zleft].v)/2.0;
					dwk0 = (a_mesh.m_level[level_n].m_data[xright].w - 
								a_mesh.m_level[level_n].m_data[xleft].w)/2.0;
					dwk1 = (a_mesh.m_level[level_n].m_data[yright].w - 
								a_mesh.m_level[level_n].m_data[yleft].w)/2.0;
					dwk2 = (a_mesh.m_level[level_n].m_data[zright].w - 
								a_mesh.m_level[level_n].m_data[zleft].w)/2.0;
					//VarReflect(dux, a_mesh.m_level[level_n].m_geom[bn].keisa, duk0, duk1, duk2, 0);
					VarReflect(duy, a_mesh.m_level[level_n].m_geom[bn].keisa, duk0, duk1, duk2, 1);
					VarReflect(duz, a_mesh.m_level[level_n].m_geom[bn].keisa, duk0, duk1, duk2, 2);
					VarReflect(dvx, a_mesh.m_level[level_n].m_geom[bn].keisa, dvk0, dvk1, dvk2, 0);
					//VarReflect(dvy, a_mesh.m_level[level_n].m_geom[bn].keisa, dvk0, dvk1, dvk2, 1);
					VarReflect(dvz, a_mesh.m_level[level_n].m_geom[bn].keisa, dvk0, dvk1, dvk2, 2);
					VarReflect(dwx, a_mesh.m_level[level_n].m_geom[bn].keisa, dwk0, dwk1, dwk2, 0);
					VarReflect(dwy, a_mesh.m_level[level_n].m_geom[bn].keisa, dwk0, dwk1, dwk2, 1);
					//VarReflect(dwz, a_mesh.m_level[level_n].m_geom[bn].keisa, dwk0, dwk1, dwk2, 2);
					omiga[0] = 0.5*(dwy - dvz);
					omiga[1] = 0.5*(duz - dwx);
					omiga[2] = 0.5*(dvx - duy);
					double mu0 = sqrt(pow(a_mesh.m_level[level_n].m_data[bn].T, 3))*
						((1.0+S_over_T_ref)/(a_mesh.m_level[level_n].m_data[bn].T+S_over_T_ref));
					vis0 = mu0/a_mesh.m_level[level_n].m_data[bn].roe;
					bigx = a_mesh.m_level[level_n].m_data[bn].viseddy/vis0;
					double bigxthree = bigx*bigx*bigx;
					fv1 = bigxthree/(bigxthree + cv1three);
					fv2 = 1.0 - bigx/(1.0+bigx*fv1);
					if (abs(a_mesh.m_level[level_n].m_box[bn].pair.signdis) < 0.00000001)
					{
#ifdef DEBUG						
						if (level_n != a_mesh.cur_level_num-1 || a_mesh.infectbox[bn] == -1)
						{
							printf("The zero dis cell is not a ib cell!!!\n");
							MPI_Abort(MPI_COMM_WORLD, 557);
						}
#endif						
						segdis = 1.0;
						vis2d = 1.0;
					}
					else
					{
						segdis = (kappa*a_mesh.m_level[level_n].m_box[bn].pair.signdis);
						vis2d = a_mesh.m_level[level_n].m_data[bn].viseddy/a_mesh.m_level[level_n].m_box[bn].pair.signdis;
					}
					segdis *= segdis;
					double omiga_length = sqrt(2.0)*omiga.length();
					bigs = omiga_length + a_mesh.m_level[level_n].m_data[bn].viseddy*fv2/segdis/Re;
					bigs = max(0.3*omiga_length, bigs);
					ft2 = ct3*exp(-ct4*bigx*bigx);
					r_ratio = min(a_mesh.m_level[level_n].m_data[bn].viseddy/(bigs*segdis)/Re, 10.0);
					rsix = pow(r_ratio, 6);
					g_ratio = r_ratio+cw2*(rsix - r_ratio);
					gsix = pow(g_ratio, 6);
					cw3six = pow(cw3, 6);
					fw = g_ratio*pow((1.0+cw3six)/(gsix+cw3six), 1.0/6.0);
					Pointxyz rel_vel = frame_rotate_vel.cross(a_mesh.m_level[level_n].m_geom[bn].boxcenter - frame_rotate_origin)
									 + frame_trans_vel;
					mutsource[level_n][b0] = -((a_mesh.m_level[level_n].m_data[bn].u-rel_vel[0])*dvisx+
										   	   (a_mesh.m_level[level_n].m_data[bn].v-rel_vel[1])*dvisy+
										   	   (a_mesh.m_level[level_n].m_data[bn].w-rel_vel[2])*dvisz);
#ifdef DEBUG
					if (isnan(mutsource[level_n][b0]))
					{
						PRINTFinLEVEL("Box (%d,%d,%d) mut source step 0 is nan!!!", level_n,
							a_mesh.m_level[level_n].m_box[bn].ix(),
							a_mesh.m_level[level_n].m_box[bn].iy(),
							a_mesh.m_level[level_n].m_box[bn].iz());
						MPI_Abort(MPI_COMM_WORLD, 1028);
					}
#endif				
					mutsource[level_n][b0] += cb1*(1.0-ft2)*bigs*a_mesh.m_level[level_n].m_data[bn].viseddy;
#ifdef DEBUG
					if (isnan(mutsource[level_n][b0]))
					{
						PRINTFinLEVEL("Box (%d,%d,%d) mut source step 1 is nan!!!", level_n,
							a_mesh.m_level[level_n].m_box[bn].ix(),
							a_mesh.m_level[level_n].m_box[bn].iy(),
							a_mesh.m_level[level_n].m_box[bn].iz());
						MPI_Abort(MPI_COMM_WORLD, 1028);
					}
#endif					
					mutsource[level_n][b0] -= (cw1*fw-cb1*ft2/(kappa*kappa))*(vis2d*vis2d)/Re;
#ifdef DEBUG
					if (isnan(mutsource[level_n][b0]))
					{
						PRINTFinLEVEL("Box (%d,%d,%d) mut source step 2 is nan!!! vis2d is %f viseddy %f dis %f", level_n,
							a_mesh.m_level[level_n].m_box[bn].ix(),
							a_mesh.m_level[level_n].m_box[bn].iy(),
							a_mesh.m_level[level_n].m_box[bn].iz(),
							vis2d, a_mesh.m_level[level_n].m_data[bn].viseddy,
							a_mesh.m_level[level_n].m_box[bn].pair.signdis);
						MPI_Abort(MPI_COMM_WORLD, 1028);
					}
#endif				
					mutsource[level_n][b0] += cb2*(dvisx*dvisx+dvisy*dvisy+dvisz*dvisz)/segma/Re;
#ifdef DEBUG
					if (isnan(mutsource[level_n][b0]))
					{
						PRINTFinLEVEL("Box (%d,%d,%d) mut source step 3 is nan!!!", level_n,
							a_mesh.m_level[level_n].m_box[bn].ix(),
							a_mesh.m_level[level_n].m_box[bn].iy(),
							a_mesh.m_level[level_n].m_box[bn].iz());
						MPI_Abort(MPI_COMM_WORLD, 1028);
					}
#endif				
					total_face_flux = 0.0;
					for (int di = 0; di < DIM; ++di)
					{
						for (int f0 = 0; f0 < 2; ++f0)
						{
							int faceindex = a_mesh.m_level[level_n].m_box[bn].faces[di][f0];
							total_face_flux += dirsign[f0]*mutflux[level_n][faceindex];
						}
					}
					mutsource[level_n][b0] += total_face_flux/a_mesh.m_level[level_n].m_geom[bn].v/Re;								
#ifdef DEBUG
					if (isnan(mutsource[level_n][b0]))
					{
						if (level_n == a_mesh.MyCurNum() - 1)
						{
							if (a_mesh.infectbox[bn] > -1 || a_mesh.m_level[level_n].m_box[bn].solid)
							{
								mutsource[level_n][b0] = 0.0;
							}
							else
							{
								PRINTFinLEVEL("Box (%d,%d,%d) signdis %f mut source is nan tag 0!!!", level_n,
									a_mesh.m_level[level_n].m_box[bn].ix(),
									a_mesh.m_level[level_n].m_box[bn].iy(),
									a_mesh.m_level[level_n].m_box[bn].iz(),
									a_mesh.m_level[level_n].m_box[bn].pair.signdis);
								MPI_Abort(MPI_COMM_WORLD, 1028);
							}
						}
						else
						{
							PRINTFinLEVEL("Box (%d,%d,%d) mut source is nan tag 1!!!", level_n,
								a_mesh.m_level[level_n].m_box[bn].ix(),
								a_mesh.m_level[level_n].m_box[bn].iy(),
								a_mesh.m_level[level_n].m_box[bn].iz());
							MPI_Abort(MPI_COMM_WORLD, 1028);
						}				
					}
#endif
				}					
			}
#ifdef TEMPORAL_REFINE
			}
#endif			
		}
#ifdef SHOWTIME
		double endtime = MPI_Wtime();
		step_tur_time += endtime - starttime;
#endif				
	}
#endif	

#ifndef LOCAL_TIME_STEPPING
	void NS_Solver::UpdateConsVars(const double & dt0)
	{
		int level0;
		if (bodynum > 0)
		{
			level0 = a_mesh.cur_level_num - 1;
		}
		else
		{
			level0 = a_mesh.cur_level_num;
		}
		//printf("Body number %d level0 %d\n", Body::Bodynum(), level0);
		for (int i = 0; i < level0; ++i)
		{
			double level_ts = dt0*step_num[i];
#ifdef TEMPORAL_REFINE
			if (ts%marching_step[i] == marching_left_step[i])
			{		
#endif			
			int start0 = a_mesh.LevelBoxStart(i);
			int end0 = a_mesh.LevelBoxEnd(i);
			for (int bn = start0; bn < end0; ++bn)
			{
				int b0 = bn - start0;
				for (int fv0 = 0; fv0 < 5; ++fv0)
				{
					cv_new[i][b0][fv0] = cv_old[i][b0][fv0]+source[0][i][b0][fv0]*level_ts/a_mesh.m_level[i].m_geom[bn].v;
					cv_new[i][b0][fv0] += source[1][i][b0][fv0]*level_ts;			
				}
#ifdef TURBULENCE
				a_mesh.m_level[i].m_data[bn].viseddy = c_vis_old[i][b0]+mutsource[i][b0]*level_ts;
				//a_mesh.m_level[i].m_data[bn].viseddy = min(c_vis_old[i][b0]+mutsource[i][b0]*level_ts, max_viseddy);
				a_mesh.m_level[i].m_data[bn].viseddy = max(a_mesh.m_level[i].m_data[bn].viseddy, min_viseddy);
#ifdef DEBUG
				if (isnan(a_mesh.m_level[i].m_data[bn].viseddy))
				{
					PRINTFinLEVEL("Box (%d,%d,%d) viseddy is NAN when updating new primary variable!!!", i,
						a_mesh.m_level[i].m_box[bn].ix(),
						a_mesh.m_level[i].m_box[bn].iy(),
						a_mesh.m_level[i].m_box[bn].iz());
					MPI_Abort(MPI_COMM_WORLD, 400);
				}			
#endif
#endif
				if (cv_new[i][b0][0] < 0.0)
				{
					printf("Level %d Box (%d,%d,%d) new density is %f old is %f!!!\n",
						i, a_mesh.m_level[i].m_box[bn].ix(),a_mesh.m_level[i].m_box[bn].iy(),a_mesh.m_level[i].m_box[bn].iz(),
						cv_new[i][b0][0], cv_old[i][b0][0]);
					MPI_Abort(MPI_COMM_WORLD, 523);
				}						
#ifdef DEBUG
				if (cv_new[i][b0][0] < 0.0)
				{
					PRINTFinLEVEL("Box %d (%d,%d,%d) new density is %f old is %f cell volume is %16.10f",i,
						bn, a_mesh.m_level[i].m_box[bn].ix(),a_mesh.m_level[i].m_box[bn].iy(),a_mesh.m_level[i].m_box[bn].iz(),
						cv_new[i][b0][0], cv_old[i][b0][0], a_mesh.m_level[i].m_geom[bn].v);
					for (Point_iterator p(0,3); p.end(); ++p)
					{
						int an0 = a_mesh.m_level[i].m_box[bn].neib[p.i][p.j][p.k];
						if (an0 > -1 && p.k == 1)
						{
							string tag0;
							int * ap = &p.i;
							for (int ap0 = 0; ap0 < 3; ++ap0)
							{
								char cd0[5];
								sprintf(cd0, "%d ", ap[ap0]);
								tag0 += cd0;
							}
							a_mesh.m_level[i].m_data[an0].showdata(tag0);
						}
					}
					for (int fi = 0; fi < DIM; ++fi)
					{
						for (int fj = 0; fj < 2; ++fj)
						{
							string tag0;
							char cd0[32];
							int face0 = a_mesh.m_level[i].m_box[bn].faces[fi][fj];
							sprintf(cd0, "D%dF%dFI%d", fi, fj, face0);
							tag0 += cd0;
							faceflux[i][face0].showdata(tag0);
						}
					}
				}		
				Assert(!cv_new[i][b0].hasnan(i, bn, 0), "new cv check after update", 164);
				Assert(cv_new[i][b0][0] > 0.0, "new cv must has a positive density!!!", 165);
#endif				
			}
#ifdef TEMPORAL_REFINE
			}
#endif			
		}
		if (bodynum > 0)
		{
			double level_ts = dt0*step_num[level0];
			int start0 = a_mesh.LevelBoxStart(level0);
			int end0 = a_mesh.LevelBoxEnd(level0);
			for (int bn = start0; bn < end0; ++bn)
			{
				if (!a_mesh.m_level[level0].m_box[bn].solid && a_mesh.infectbox[bn] == -1)					
				{
					int b0 = bn - start0;
					for (int fv0 = 0; fv0 < 5; ++fv0)
					{
						cv_new[level0][b0][fv0] = cv_old[level0][b0][fv0]+
							source[0][level0][b0][fv0]*level_ts/a_mesh.m_level[level0].m_geom[bn].v;
						cv_new[level0][b0][fv0] += source[1][level0][b0][fv0]*level_ts;
					}
#ifdef TURBULENCE
					a_mesh.m_level[level0].m_data[bn].viseddy = c_vis_old[level0][b0]+mutsource[level0][b0]*level_ts;
					//a_mesh.m_level[level0].m_data[bn].viseddy = min(c_vis_old[level0][b0]+mutsource[level0][b0]*level_ts, max_viseddy);
					a_mesh.m_level[level0].m_data[bn].viseddy = max(a_mesh.m_level[level0].m_data[bn].viseddy, min_viseddy_wall);
#ifdef DEBUG
					if (isnan(a_mesh.m_level[level0].m_data[bn].viseddy))
					{
						PRINTFinLEVEL("Box (%d,%d,%d) viseddy is nan when updating new primary variable!!!", level0,
							a_mesh.m_level[level0].m_box[bn].ix(),
							a_mesh.m_level[level0].m_box[bn].iy(),
							a_mesh.m_level[level0].m_box[bn].iz());
						MPI_Abort(MPI_COMM_WORLD, 468);
					}			
#endif					
#endif
					if (cv_new[level0][b0][0] < 0.0)
					{
						printf("Level %d Box (%d,%d,%d) new density is %f old is %f!!!\n",
							level0, 
							a_mesh.m_level[level0].m_box[bn].ix(),
							a_mesh.m_level[level0].m_box[bn].iy(),
							a_mesh.m_level[level0].m_box[bn].iz(),
							cv_new[level0][b0][0], cv_old[level0][b0][0]);
						MPI_Abort(MPI_COMM_WORLD, 611);
					}								
					Assert(!cv_new[level0][b0].hasnan(level0, bn, 0), "new cv check after update", 179);
					Assert(cv_new[level0][b0][0] > 0.0, "new cv must has a positive density!!!", 181);
				}
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
#else

	void NS_Solver::UpdateConsVars(const double & dt0)
	{
#ifdef TEMPORAL_REFINE
		printf("We can not have TEMPORAL_REFINE and LOCAL_TIME_STEPPING at the same time!!!");
		MPI_Abort(MPI_COMM_WORLD, 617);
#endif
#ifndef TURBULENCE
		printf("Why LOCAL_TIME_STEPPING is defined for laminar flow???");
		MPI_Abort(MPI_COMM_WORLD, 635);
#endif					
		int level0;
		if (bodynum > 0)
		{
			level0 = a_mesh.cur_level_num - 1;
		}
		else
		{
			level0 = a_mesh.cur_level_num;
		}
		//printf("Body number %d level0 %d\n", Body::Bodynum(), level0);
		for (int i = 0; i < level0; ++i)
		{		
			int start0 = a_mesh.LevelBoxStart(i);
			int end0 = a_mesh.LevelBoxEnd(i);
			for (int bn = start0; bn < end0; ++bn)
			{
				int b0 = bn - start0;
				for (int fv0 = 0; fv0 < 5; ++fv0)
				{				
					cv_new[i][b0][fv0] = cv_old[i][b0][fv0]+source[0][i][b0][fv0]*tp[subt]*local_dt[i][bn]/a_mesh.m_level[i].m_geom[bn].v;
					cv_new[i][b0][fv0] += source[1][i][b0][fv0]*tp[subt]*local_dt[i][bn];												
				}
				a_mesh.m_level[i].m_data[bn].viseddy = c_vis_old[i][b0]+mutsource[i][b0]*tp[subt]*local_dt[i][bn];
				//a_mesh.m_level[i].m_data[bn].viseddy = min(c_vis_old[i][b0]+mutsource[i][b0]*level_ts, max_viseddy);
				a_mesh.m_level[i].m_data[bn].viseddy = max(a_mesh.m_level[i].m_data[bn].viseddy, min_viseddy);
				// if (cv_new[i][b0][0] < 0.0 || isnan(a_mesh.m_level[i].m_data[bn].viseddy))
				// {
				// 	cv_new[i][b0] = cv_old[i][b0];
				// 	a_mesh.m_level[i].m_data[bn].viseddy = c_vis_old[i][b0];
				// }										
			}		
		}
		if (bodynum > 0)
		{
			int start0 = a_mesh.LevelBoxStart(level0);
			int end0 = a_mesh.LevelBoxEnd(level0);
			for (int bn = start0; bn < end0; ++bn)
			{
				if (!a_mesh.m_level[level0].m_box[bn].solid && a_mesh.infectbox[bn] == -1)					
				{
					int b0 = bn - start0;
					for (int fv0 = 0; fv0 < 5; ++fv0)
					{
						cv_new[level0][b0][fv0] = cv_old[level0][b0][fv0]+
							source[0][level0][b0][fv0]*tp[subt]*local_dt[level0][bn]/a_mesh.m_level[level0].m_geom[bn].v;
						cv_new[level0][b0][fv0] += source[1][level0][b0][fv0]*tp[subt]*local_dt[level0][bn];
					}
					a_mesh.m_level[level0].m_data[bn].viseddy = c_vis_old[level0][b0]+mutsource[level0][b0]*tp[subt]*local_dt[level0][bn];
					//a_mesh.m_level[level0].m_data[bn].viseddy = min(c_vis_old[level0][b0]+mutsource[level0][b0]*level_ts, max_viseddy);
					a_mesh.m_level[level0].m_data[bn].viseddy = max(a_mesh.m_level[level0].m_data[bn].viseddy, min_viseddy_wall);
					// if (cv_new[level0][b0][0] < 0.0 || isnan(a_mesh.m_level[level0].m_data[bn].viseddy))
					// {
					// 	cv_new[level0][b0] = cv_old[level0][b0];
					// 	a_mesh.m_level[level0].m_data[bn].viseddy = c_vis_old[level0][b0];
					// }																											
					Assert(!cv_new[level0][b0].hasnan(level0, bn, 0), "new cv check after update", 179);
					Assert(cv_new[level0][b0][0] > 0.0, "new cv must has a positive density!!!", 181);
				}
			}
		}
		//MPI_Allreduce(MPI_IN_PLACE, &restoreflag, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	}
#endif


	void NS_Solver::RestoreConsVar()
	{
		for (int i = 0; i < a_mesh.MyCurNum(); ++i)
		{			
			int start0 = a_mesh.LevelBoxStart(i);
			int end0 = a_mesh.LevelBoxEnd(i);
			for (int bn = start0; bn < end0; ++bn)
			{
				int b0 = bn - start0;
				cv_new[i][b0] = cv_old[i][b0];
#ifdef TURBULENCE				
				a_mesh.m_level[i].m_data[bn].viseddy = c_vis_old[i][b0];
#endif				
				Get_Cons_Vars(a_mesh.m_level[i].m_data[bn], cv_old[i][b0]);						
			}						
		}
		MPI_Barrier(share_comm);
	}	
