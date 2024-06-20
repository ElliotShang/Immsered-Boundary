#ifndef _INPUTFUNC_H_
#define _INPUTFUNC_H_

#include "BCValues.H"
#include "AMR.H"
#include "Pltoutput.H"

extern const int init_jump;
extern const int filenodenum;
extern const int filelevelnum;
extern const int filebodynum;
extern const char soln_file_name[30];

void InputSoln(Mesh & amesh, BCValues & meshbc, AMR & myamr, Pltoutput & pltfile);

void ReadLevelSolution(Mesh & amesh, Pltoutput & pltfile);

void JumpIfstreamLine(ifstream & as, const int & jl);

int GetANumAfterString(char * a, const string & astring, const char & endchar);

void JumpIfstreamString(ifstream & as, const int & js);

#endif