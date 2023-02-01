//
//  initializelandscape.h
//  gapdynamics
//
//  Created by Dylan DeBaun on 1/20/20.
//  Copyright © 2020 Dylan DeBaun. All rights reserved.
//

#ifndef initializelandscape_h
#define initializelandscape_h

#include <stdio.h>

void get_input_info(char *infile,
                    float *p_incr_percent,
                    int *p_size_x, int *p_size_y,
                    double *p_kernel_cutoff, int *p_max_disp_dist, int *p_gauss, double *p_sigma, int *p_fat_tail, double *p_u, double *p_p, int *p_global, int *p_near_neigh, double **p_disc_norm_cum, int *p_loc_kern_cutoff,
                    int *p_turnover_time, int *p_num_intev_to_meas, int *p_time_btwn_meas,
                    int *p_output_initial_spec_dist, int *p_output_final_spec_dist, int *p_output_interm_spec_dist, int *p_num_interm_outputs, double **p_when_to_output, int **p_outputted_interm, int *p_output_lowest_x, int *p_output_highest_x, int *p_output_lowest_y, int *p_output_highest_y);

void read_in_gaps(char *gap_file, int **gap_xcrdactual, int **gap_ycrdactual, int **gap_size_xactual, int **gap_size_yactual, int *numbergaps, int sizex, int sizey);

void read_in_nongaps(char *gap_file, int **gap_xcrdactual, int **gap_ycrdactual);

void read_in_landscape(char *landscapefile, double **LC_xcrd, double **LC_ycrd, int *num_indiv,int **parent_spp_ID, int *num_species);

void nongapsrandom(int *gapgroup, int *gapgroupsize, int gapgroupnum, int *gap_xcrd, int *gap_ycrd,  int numbersmallgaps, int **nongap_xcrdactual, int **nongap_ycrdactual,int size_x, int size_y);

#endif /* initializelandscape_h */
