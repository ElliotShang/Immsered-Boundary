#include "AMR.H"

void AMR::InterfaceExchange()
{
#ifdef SHOWTIME
	double starttime = MPI_Wtime();
#endif	
	for (int i = 0; i < m_mesh.cur_level_num-1; ++i)
	{
#ifdef TEMPORAL_REFINE
		if (ts%average_step[i+1] == average_left_step[i+1])
		{
#endif		
		InterfaceRestriction(i);
#ifdef TEMPORAL_REFINE
		}
#endif		
	}
	m_mesh.DataExchange_alllevels(ts+1);
	for (int i = 0; i < m_mesh.cur_level_num-1; ++i)
	{
#ifdef TEMPORAL_REFINE
		if (ts%prolongation_step[i+1] == prolongation_left_step[i+1])
		{
#endif		
		InterfaceProlongation(i);
#ifdef TEMPORAL_REFINE
		}
#endif		
	}
#ifdef SHOWTIME
	double endtime = MPI_Wtime();
	step_interface_time = endtime - starttime;
#endif	
}

void AMR::InterfaceProlongation(const int & nlevel)
{
	int bs = m_level[nlevel].f_pro_ghost.ps();
	int be = m_level[nlevel].f_pro_ghost.pe();
	int gbs = m_level[nlevel].g_pro.ps();
	int gbe = m_level[nlevel].g_pro.pe();
	if (!m_level[nlevel+1].twodflag)
	{
		for (int i = bs; i < be; ++i)
		{
			PairProlongation(m_mesh, nlevel, 
				m_level[nlevel].f_pro_ghost[i].ci,
				m_level[nlevel].f_pro_ghost[i].fi, tlop);
		}
		for (int i = gbs; i < gbe; ++i)
		{
			CtoFPair & cfp = m_level[nlevel].level_pro_ghost[m_level[nlevel].g_pro[i]];
			PairProlongation(m_mesh, nlevel, 
				cfp.ci,
				cfp.fi, tlop);
		}
	}
	else
	{
		for (int i = bs; i < be; ++i)
		{
			PairProlongation_2d(m_mesh, nlevel, 
				m_level[nlevel].f_pro_ghost[i].ci,
				m_level[nlevel].f_pro_ghost[i].fi, tlop);
		}
		for (int i = gbs; i < gbe; ++i)
		{
			CtoFPair & cfp = m_level[nlevel].level_pro_ghost[m_level[nlevel].g_pro[i]];
			PairProlongation_2d(m_mesh, nlevel, 
				cfp.ci,
				cfp.fi, tlop);
		}
	}
	MPI_Barrier(share_comm);
}

void AMR::InterfaceRestriction(const int & nlevel)
{
	if (!m_level[nlevel+1].twodflag)
	{
		for (int ig = 0; ig < ighost; ++ig)
		{
			int bs = m_level[nlevel].f_res_ghost[ig].ps();
			int be = m_level[nlevel].f_res_ghost[ig].pe();
			for (int i = bs; i < be; ++i)
			{
				PairRestriction(m_mesh, nlevel, 
					m_level[nlevel].f_res_ghost[ig][i].ci, 
					m_level[nlevel].f_res_ghost[ig][i].fi, tlop);
			}	
		}
	}
	else
	{
		for (int ig = 0; ig < ighost; ++ig)
		{
			int bs = m_level[nlevel].f_res_ghost[ig].ps();
			int be = m_level[nlevel].f_res_ghost[ig].pe();
			for (int i = bs; i < be; ++i)
			{
				PairRestriction_2d(m_mesh, nlevel, 
					m_level[nlevel].f_res_ghost[ig][i].ci, 
					m_level[nlevel].f_res_ghost[ig][i].fi, tlop);
			}	
		}
	}
	//MPI_Win_fence(0, m_mesh.m_level[nlevel].m_data.arraywin());
	MPI_Barrier(share_comm);
}
