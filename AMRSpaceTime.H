#ifndef _AMRSpaceTime_H
#define _AMRSpaceTime_H

#include <vector>

extern double t;
extern double dt;
extern int ts;
extern const int ts_sum;

extern const Point lowpt, highpt;
extern bool periodic[3];
extern const int max_mesh_level;
extern const int init_level_num;
extern const int first_distance_renew_level;
extern vector<vector<double> > dh;
extern double total_vs_time;
extern double total_inv_time;
extern double total_ib_time;
extern double total_force_time;
extern double step_inv_time;
extern double step_vs_time;
extern double step_ib_time;
extern double step_total_time;
extern double step_tur_time;
extern double step_interface_time;
extern double step_exchange_time;
extern double step_move_time;
extern const double CFL;
extern double CFL_dt;
extern const double dt_ratio;
extern int dt_var_step;
extern const int dts_time_var;
extern const int dts_space_var;
extern const int dts_plt;
extern vector<int> marching_step;
extern vector<int> average_step;
extern vector<int> prolongation_step;
extern vector<int> marching_left_step;
extern vector<int> average_left_step;
extern vector<int> prolongation_left_step;
extern vector<int> level_refine_ratio;
extern vector<Point> level_grid_ratio;
extern vector<double> step_num;
extern double dmlength[3];
extern double critic_ylength;
extern double t_sum;

void InitSpaceParams();

void ComptLevelgridsize(double * dh0);

void SetLevelOperationIndex();

void PeriodicLength(Pointxyz & dxyz);

double ComptPointAngle_Rotate_X(Pointxyz & p1, Pointxyz & p2);

#ifdef PASSAGE_ANGLE
void PeriodicAnnulaLength(Pointxyz & p1, Pointxyz & p2, Pointxyz & p1_new);
#endif

#endif