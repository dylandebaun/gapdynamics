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

void select_randloc_frm_cum_disp(int *ind1_coords, int size_x, int size_y, double *disc_norm_cum, int loc_cutoff, int *ind2_coords);

void min_dist_bet_indiv(int J, double *LC_xcrd, double *LC_ycrd, double *threshold);

void run_actual(char *landscapefile, char *begoutfile, int *gap_xcrd, int *gap_ycrd, int *gap_sizex, int *gap_sizey, int gap_number, int curr_year, int *gapgroup, int groupgapnum, int *groupgapsize, int num_species);

void run_actualnongap(char *landscapefile, char *begoutfile, int **nongap_xcrd, int **nongap_ycrd, int curr_year, int *gapgroup, int groupgapnum, int *groupgapsize, int num_species);


#endif /* timesteps_dispersal_h */
