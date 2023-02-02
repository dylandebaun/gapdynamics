//
//  stats_analysis.h
//  gapdynamics
//
//  Created by Dylan DeBaun on 3/25/20.
//  Copyright © 2020 Dylan DeBaun. All rights reserved.
//

#ifndef stats_analysis_h
#define stats_analysis_h

#include <stdio.h>

void print_landscape(char *begoutfile, double *LC_xcrd, double *LC_ycrd, int curr_year, int *species_ID, int J, int s, int *parentID);

void print_gaps(char *begoutfile, int *gap_xcrd, int *gap_ycrd, int *gap_sizex, int *gap_sizey, int gap_number);

void print_species_key(char *begoutfile, char **species_name, int species_number);

void ascendingarray(int *samples);

void print_SAD_per_simulation_groupgaps(char *begoutfile, int *gap_xcrd, int *gap_ycrd, int *gap_sizex, int *gap_sizey, int gap_number, int simnumber, double *LC_xcrd, double *LC_ycrd, int numindiv, int *spID, int num_species, double **speciescount, double **speciescountnongap, int curr_year, int *gapgroup, int groupgapnum, int *groupgapsize, int **nongap_xcrd1, int **nongap_ycrd1, int numindivstart);

#endif /* stats_analysis_h */
