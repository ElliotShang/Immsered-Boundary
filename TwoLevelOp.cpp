#include "TwoLevelOp.H"

// FlowVariables * TwoLevelOp::neibflow[3][3][3] = {{{NULL,NULL,NULL}, {NULL,NULL,NULL}, 
// 													{NULL,NULL,NULL}},
// 												 {{NULL,NULL,NULL}, {NULL,NULL,NULL}, {NULL,NULL,NULL}},
// 												 {{NULL,NULL,NULL}, {NULL,NULL,NULL}, {NULL,NULL,NULL}}};
// FlowVec TwoLevelOp::neibflowvec[3][3][3] = {{{FlowVec(), FlowVec(), FlowVec()}, 
// 											 {FlowVec(), FlowVec(), FlowVec()}, 
// 											 {FlowVec(), FlowVec(), FlowVec()}},
// 											{{FlowVec(), FlowVec(), FlowVec()}, {FlowVec(), FlowVec(), FlowVec()}, 
// 											{FlowVec(), FlowVec(), FlowVec()}},
// 											{{FlowVec(), FlowVec(), FlowVec()}, {FlowVec(), FlowVec(), FlowVec()}, 
// 											{FlowVec(), FlowVec(), FlowVec()}}};
// Boxson<FlowVariables *> TwoLevelOp::finedataptr = Boxson<FlowVariables *>();
// FlowVec TwoLevelOp::con_gradx = FlowVec(); 
// FlowVec TwoLevelOp::con_grady = FlowVec(); 
// FlowVec TwoLevelOp::con_gradz = FlowVec();

// Face * TwoLevelOp::finemcptr[3][2][4] = {{{NULL,NULL,NULL,NULL},{NULL,NULL,NULL,NULL}},
// 																				 {{NULL,NULL,NULL,NULL},{NULL,NULL,NULL,NULL}},
// 																				 {{NULL,NULL,NULL,NULL},{NULL,NULL,NULL,NULL}}};
// Boxson<double> TwoLevelOp::finevolume = Boxson<double>();

// Point_iterator TwoLevelOp::cneib[7] = {Point_iterator(1,1,0), Point_iterator(1,1,2), 
// 																			 Point_iterator(1,0,1), Point_iterator(1,2,1),
// 																			 Point_iterator(0,1,1), Point_iterator(2,1,1),
// 																			 Point_iterator(1,1,1)};

// CellMCoef * TwoLevelOp::neibmcoef[3][3][3] = {{{NULL,NULL,NULL}, {NULL,NULL,NULL}, 
// 													{NULL,NULL,NULL}},
// 												 {{NULL,NULL,NULL}, {NULL,NULL,NULL}, {NULL,NULL,NULL}},
// 												 {{NULL,NULL,NULL}, {NULL,NULL,NULL}, {NULL,NULL,NULL}}};

// // vector<int> TwoLevelOp::neibrefx(2+bcnx);
// // vector<int> TwoLevelOp::neibrefy(2+bcny);
// // vector<int> TwoLevelOp::neibrefz(2+bcnz);
// // vector<int> TwoLevelOp::neibptx(2+bcnx);
// // vector<int> TwoLevelOp::neibpty(2+bcny);
// // vector<int> TwoLevelOp::neibptz(2+bcnz);
// //int TwoLevelOp::off_ftoc[4] = {bcnx, 0, bcnx, 0};
// double TwoLevelOp::fincrem[2] = {-0.25, 0.25};
// int TwoLevelOp::off_ffneibx[3] = {-bcnx, 0, 0};
// int TwoLevelOp::off_ffneiby[3] = {-bcny, 0, 0};
// int TwoLevelOp::off_ffniebz[3] = {-bcnz, 0, 0};

// vector<FlowVariables> TwoLevelOp::tempfv(10);
// vector<CellMCoef>     TwoLevelOp::tempmc(10);
int neibson[4] = {1,0,1,0};
