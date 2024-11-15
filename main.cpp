#include <iostream>
#include "Point.H"
#include "Domain.H"
#include "Mesh.H"
#include "Pltoutput.H"
#include "AMR.H"
#include "AMRmpi.H"
#include "AMRSpaceTime.H"
#include "NS_Solver.H"
#include "TwoLevelOp.H"
#include "InputFunc.H"
double f_osc0 = .20;
void TimeOperation(vector<Body> & sphere, vector<Body> & localbody, Pltoutput & pltfile, AMR & amr, NS_Solver & aslove);

void PreMesh(Domain & mydomain);

void testFunc(vector<Body> & abody);

int main(int argc, char ** argv)
{
	MPI_Init(&argc, &argv);
	MPI_CreateNode();
	GiveAFlag("InitSpaceParams...", 5);
	InitSpaceParams();

	GiveAFlag("Try to make a Domain!!!", 5);
	Domain p_domain(lowpt, highpt, periodic);

	PreMesh(p_domain);
	GiveAFlag("Finish PreMesh!!!", 5);
	Mesh mesh(p_domain, max_mesh_level);
	GiveAFlag("Finish make a mesh!!!", 5);
	
	vector<Body> sphere;
	vector<Body> localbody;
	Body::SetBodyParams(sphere, mesh, localbody);
	GiveAFlag("Sphere generated!!!", 5);

	AMR amr(mesh);

	GiveAFlag("Finish embeded the wall!!!", 5);

	Pltoutput pltfile(mesh);
	pltfile.Design_User_File();

	//pltfile.Output_Grid_Nobody("amr_grid-init-nobody");
	//pltfile.Output_Soln_Nobody("amr_soln-init-nobody");
	//pltfile.Output_Grid_Onlybody("amr_grid-init-localbody", localbody);
	//pltfile.Output_Soln_Onlybody("amr_soln-init-localbody", localbody);
	
	amr.SetupNewAMR(sphere);

	//pltfile.Output_Grid_Nobody("amr_cell-nobody");
	//pltfile.Output_GridParameters("amr_cell-params");

	GiveAFlag("AMR Generated!!!", 5);
	Body::ImmerseWallBoundary(mesh, sphere);

	GiveAFlag("Finish the first adaptive!!!", 5);
	
	NS_Solver flowsolver(mesh, p_domain, (int)sphere.size());
	amr.AllocateDataMem();
#ifndef INPUT_SOLN
#ifndef INITFLOW_USER	
	flowsolver.InitFlow();
#else
	flowsolver.InitFlow_User();
#endif		
#else
	InputSoln(mesh, flowsolver.BCinfo(), amr, pltfile);
	//amr.PerformAdaptive(sphere);
#endif
	GiveAFlag("Finish initialize flow...", 5);
	double start_time, end_time;
	double step_start_time, step_end_time;
	MPI_Barrier(MPI_COMM_WORLD);
	pltfile.Output_Grid("amr_grid-init", sphere);
	pltfile.Output_Soln("amr_soln-init", sphere);
	
	GiveAFlag("Finish initial grid output....\n", 5);
#ifdef SHOWTIME	
	start_time = MPI_Wtime();
#endif
	while (ts < ts_sum || t < t_sum)
	{		
		step_start_time = MPI_Wtime();
		amr.InterfaceExchange();	
		GiveAFlag("Step begin >>> start to obtain IBConditions...", 5);
		Body::IBConditions_Image(mesh, sphere);
		if (additional_output)
		{
			//pltfile.Output_Grid("amr_grid-periodic", sphere);
			//pltfile.Output_Soln("amr_soln-periodic", sphere);
			additional_output = false;
		}
		GiveAFlag("IBConditions has been finished...", 5);
		flowsolver.TimeMarching(amr);
		++ts;
		t += dt;
		step_end_time = MPI_Wtime();
		step_total_time = step_end_time - step_start_time;
#ifdef SHOWTIME
		
#ifdef TURBULENCE		
		if (nrank == 0) printf("step_total_time: %f vs: %f inv: %f ib: %f(move %f) mut: %f interface: %f exchange: %f\n", 
			step_total_time, step_vs_time, step_inv_time, step_ib_time, step_move_time, step_tur_time, step_interface_time, step_exchange_time);
#else
		if (nrank == 0) printf("step_total_time: %f vs: %f inv: %f ib: %f(move %f) interface: %f exchange: %f\n", 
			step_total_time, step_vs_time, step_inv_time, step_ib_time, step_move_time, step_interface_time, step_exchange_time);
#endif				
#else
		if (nrank == 0) printf("step_total_time: %f\n", step_total_time);		
#endif
		TimeOperation(sphere, localbody, pltfile, amr, flowsolver);		
		
	}
#ifdef SHOWTIME
	end_time = MPI_Wtime();
	PrintTime(start_time, end_time, "Total");
#endif	
	MPI_Barrier(share_comm);
	MPI_FreeDataType();
	MPI_Finalize();
	return 0;
}


// int main(int argc, char ** argv)
// {
// 	MPI_Init(&argc, &argv);
// 	MPI_CreateNode();
// 	GiveAFlag("InitSpaceParams...", 5);
// 	vector<Body> abody;
// 	InitSpaceParams();
// 	//testFunc(abody);
// 	DataArray<FlowVariables> f0;
// 	int i = 15000;
// 	while(true)
// 	{
// 		if (i > 8000) 
// 		{
// 			i = 5000;
// 			f0.setnum_nocopy(i*10000, 0);
// 		}
// 		else
// 		{
// 			i = 15000;
// 			f0.setnum_nocopy(i*10000, 0);
// 		}
// 	}
// 	MPI_Finalize();
// 	return 0;
// }

// void testFunc(vector<Body> & abody)
// {
// 	abody.resize(5);
// 	// abody.reserve(5);
// 	// Body a0 = Body(0);
// 	// // a0.set_pointnum(100);
// 	// // a0.set_patchnum(200);
// 	// // DataArray<FlowVariables> testarray;
// 	// // int i = 15000;
// 	// // // bool cycle = true;
// 	// // // while (cycle)
// 	// // // {
// 	// // // 	if (i < 10000) i=15000;
// 	// // // 	else i=5000;
// 	// // 	int arraysize = i*10000;
// 	// // 	testarray.setnum_nocopy(arraysize, 0);
// 	// // // }
// 	// int i0[5] = {100,200,300,400,500};
// 	// for (int i = 0; i < 5; ++i)
// 	// {
// 	// 	//Body a0 = ;
// 	// 	abody.push_back(Body(i));
// 	// 	abody[i].set_bodyindex(i+5);
// 	// 	abody[i].set_pointnum(i0[i]+10);
// 	// 	abody[i].set_patchnum(i0[i]+20);
// 	// }
// 	//abody.resize(3, a0);
// }

void TimeOperation(vector<Body> & sphere, vector<Body> & localbody, Pltoutput & pltfile, AMR & amr, NS_Solver & aslove)
{
	if (ts%dts_time_var == 0)
	{
		pltfile.Print_Var_as_Func_of_Time(sphere, aslove);
		// for (int i = 0; i < sphere.size(); ++i)
		// {
		// 	if (sphere[i].draglift(0) > 0.6)
		// 	{
		// 		pltfile.Output_Soln("amr_soln_emergency", sphere);
		// 	}
		// }
	}
	if (ts%dts_space_var == 0)
	{
		pltfile.Print_Var_as_Func_of_Space(sphere);
	}
	amr.ADPRule(sphere);
	if (adpflag)
	{
		// amr.InterfaceExchange();
		// pltfile.Output_Grid("before_adp_grid", sphere);
		// pltfile.Output_Soln("before_adp_soln", sphere);
		amr.PerformAdaptive(sphere);
		// pltfile.Output_Grid("after_adp_grid", sphere);
		// pltfile.Output_Soln("after_adp_soln", sphere);
	}
	if (ts%dts_plt == 0 && t > t_sum - 15.0)
	{
		pltfile.Output_Soln_Nobody("amr_soln");
        pltfile.Output_Soln_Onlybody("Body_soln",localbody);
		if (gridisnew) 
		{
			pltfile.Output_Grid_Nobody("amr_grid");
			pltfile.Output_Grid_Onlybody("amr_grid_body", sphere);
			gridisnew = false;
		}
	}
	adpflag = false;
}

void PreMesh(Domain & mydomain)
{
#ifdef PASSAGE_ANGLE
	if (nrank == 0)
	{
		printf("Start the PreMesh...\n");
	}
	vector<int> nodestart;
	vector<int> nodeend;
	Domain pre_domain(lowpt, highpt, periodic);
	Mesh mesh(pre_domain, max_mesh_level);
	vector<Body> sphere;
	vector<Body> localbody;
	Body::SetBodyParams(sphere, mesh, localbody);
	AMR amr(mesh);
	amr.SetupNewAMR(sphere);
	mesh.CounterSliceCellNumber(2, nodestart, nodeend);
	GiveAFlag("Finish CounterSliceCellNumber!!!", 5);
	mydomain.Set_Z_Bound(nodestart, nodeend);
	GiveAFlag("Finish Set_Z_Bound!!!", 5);
	MPI_Barrier(MPI_COMM_WORLD);
	if (nrank == 0)
	{
		printf("Finish the PreMesh!!!\n");
	}
#endif	
}
