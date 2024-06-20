#include "Filt.H"
#include "AMRSpaceTime.H"
#include "Body.H"

double filt_coef_a[FORDER+1] = {0.0, 0.0, 0.0, 0.0, 0.0};
double filt_coef_b[FORDER+1] = {1.0, 0.0, 0.0, 0.0, 0.0};

void FiltPatchPressure(Body & abody)
{
	DataArray<Surfpatch> & patch = abody.bodypatch();
	int pps = patch.ps();
	int ppe = patch.pe();
	if (ts == dts_time_var)
	{
		for (int i = pps; i < ppe; ++i)
		{
			for (int j = 0; j < FORDER; ++j)
			{
				patch[i].filt_pre_input[j] = 0.0;
				patch[i].filt_pre_output[j] = 0.0;
			}
		}
	}
	for (int i = pps; i < ppe; ++i)
	{
		double filt_pre0 = filt_coef_b[0]*patch[i].hgc.fv.p;
		for (int j = 0; j < FORDER; ++j)
		{
			filt_pre0 += filt_coef_b[j+1]*patch[i].filt_pre_input[j]-filt_coef_a[j+1]*patch[i].filt_pre_output[j];
		}
		for (int j = FORDER-1; j > 0; --j)
		{
			patch[i].filt_pre_input[j] = patch[i].filt_pre_input[j-1];
			patch[i].filt_pre_output[j] = patch[i].filt_pre_output[j-1];
		}
		patch[i].filt_pre_input[0] = patch[i].hgc.fv.p;
		patch[i].filt_pre_output[0] = filt_pre0;
		patch[i].hgc.fv.p = filt_pre0;
	}

	MPI_Barrier(share_comm);
}

void FiltSinglePatchPressure(Surfpatch & apatch)
{
	if (ts == dts_time_var)
	{
		for (int j = 0; j < FORDER; ++j)
		{
			apatch.filt_pre_input[j] = 0.0;
			apatch.filt_pre_output[j] = 0.0;
		}
	}
	double filt_pre0 = filt_coef_b[0]*apatch.hgc.fv.p;
	for (int j = 0; j < FORDER; ++j)
	{
		filt_pre0 += filt_coef_b[j+1]*apatch.filt_pre_input[j] - filt_coef_a[j+1]*apatch.filt_pre_output[j];
	}
	for (int j = FORDER-1; j > 0; --j)
	{
		apatch.filt_pre_input[j] = apatch.filt_pre_input[j-1];
		apatch.filt_pre_output[j] = apatch.filt_pre_output[j-1];
	}
	apatch.filt_pre_input[0] = apatch.hgc.fv.p;
	apatch.filt_pre_output[0] = filt_pre0;
	apatch.hgc.fv.p = filt_pre0;
}