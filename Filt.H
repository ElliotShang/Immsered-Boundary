#ifndef _FILT_H_
#define _FILT_H_

#include <iostream>
#include <fstream>
#include <cstdio>

#include "Box.H"
#include "Body.H"

extern double filt_coef_a[FORDER+1];
extern double filt_coef_b[FORDER+1];

void Get_Filt_Params();

void Renew_Filt_Params();

void FiltSinglePatchPressure(Surfpatch & apatch);

void FiltPatchPressure(Body & abody);

#endif
