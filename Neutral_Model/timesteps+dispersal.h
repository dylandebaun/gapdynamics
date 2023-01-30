//
//  timesteps+dispersal.h
//  gapdynamics
//
//  Created by Dylan DeBaun on 1/21/20.
//  Copyright © 2020 Dylan DeBaun. All rights reserved.
//

#ifndef timesteps_dispersal_h
#define timesteps_dispersal_h

#include <stdio.h>
double *discrete_normalized_cummulative(double f(double, double*),
                                        double *parameters,
                                        double cutoff, int max_disp_dist,
                                        int *point_loc_cutoff);
void min_dist_bet_indiv(int J, double *LC_xcrd, double *LC_ycrd, double *threshold);
void run_actual(char *landscapefile, char *begoutfile, int *gap_xcrd, int *gap_ycrd, int *gap_sizex, int *gap_sizey, int gap_number, int **speciescount, int curr_year, int *gapgroup, int groupgapnum, int *groupgapsize, int num_species);
//void run_actual(char *landscapefile, char *begoutfile, int *gap_xcrd, int *gap_ycrd, int *gap_sizex, int *gap_sizey, int gap_number, int **speciescount, int curr_year, int *gapgroup, int groupgapnum, int *groupgapsize, int num_species, double **LC_xcrd, double **LC_ycrd, int **parent_spp_ID, int *numind);
//void run_actualnongap(char *landscapefile, char *begoutfile, int *gap_xcrd, int *gap_ycrd, int *gap_sizex, int *gap_sizey, int gap_number, int **speciescount, int curr_year, int *gapgroup, int groupgapnum, int *groupgapsize, int num_species, double **LC_xcrd, double **LC_ycrd, int **parent_spp_ID, int *numind);
void run_actualnongap(char *landscapefile, char *begoutfile, int **nongap_xcrd, int **nongap_ycrd, double **speciescount, int curr_year, int *gapgroup, int groupgapnum, int *groupgapsize, int num_species);
//void select_randloc_frm_cum_disp(int *ind1_coords, int sizex, int sizey, double *disc_norm_cum, int loc_cutoff, int *ind2_coords, int adjustmentx, int adjustmenty, int *i);

//double time_evolve(int ***LANDSCAPE, int size_x, int size_y, int *p_num_indiv, int *p_num_slots_inuse, int *p_curr_num_slots, int inc_num_slots, double *disc_norm_cum, int loc_kern_cutoff, int turnover_time, int do_age, short **p_year_born, int maturity_age, double curr_year);

//int choose_event_indiv(int ***LANDSCAPE, int sizex, int sizey);

//int choose_parent_indiv_wagestruct(int ***LANDSCAPE, int num_slots_inuse, double curr_year, int maturity_age, int sizex, int sizey);

//void sort_or_reallocate(short ***p_INDIV, int num_indiv, int *p_num_slots_inuse, int *p_curr_num_slots, int inc_num_slots, int do_age, short **p_year_born);

//void gen_offspring(int *event_ind, int ***LANDSCAPE, int size_x, int size_y, int *p_num_indiv, int *p_num_slots_inuse, double *disc_norm_cum, int loc_kern_cutoff, int do_age, short *year_born, double curr_year);

//void Kill(int *event_ind, int ***LANDSCAPE, int *p_num_indiv);

#endif /* timesteps_dispersal_h */
