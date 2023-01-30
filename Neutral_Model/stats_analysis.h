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

#endif /* stats_analysis_h */

void print_landscape(char *begoutfile, double *LC_xcrd, double *LC_ycrd, int curr_year, int *species_ID, int J, int s, int *parentID);

void print_gaps(char *begoutfile, int *gap_xcrd, int *gap_ycrd, int *gap_sizex, int *gap_sizey, int gap_number);

void print_species_key(char *begoutfile, char **species_name, int species_number);

void ascendingarray(int *samples);

void print_SAD_in_gaps(char *begoutfile, int *gap_xcrd, int *gap_ycrd, int *gap_sizex, int *gap_sizey, int gap_number, int simnumber, double *LC_xcrd, double *LC_ycrd, int numindiv, int *spID, int num_species, int ***speciescount, int curr_year);

void calc_SAD_in_gaps(char *begoutfile, int *gap_xcrd, int *gap_ycrd, int *gap_sizex, int *gap_sizey, int gap_number, double *LC_xcrd, double *LC_ycrd, int numindiv, int *spID, int num_species, int ***speciescount, int curr_year);

void print_SAD_per_simulation(char *begoutfile, int *gap_xcrd, int *gap_ycrd, int *gap_sizex, int *gap_sizey, int gap_number, int simnumber, double *LC_xcrd, double *LC_ycrd, int numindiv, int *spID, int num_species, int *speciescount, int curr_year);

void printgroupgapinfo(char *begoutfile, int *groupgap, int *groupgapsize, int groupgapnum , int numgaps);

//void print_SAD_per_simulation_groupgaps(char *begoutfile, char *simtype, int *gap_xcrd, int *gap_ycrd, int *gap_sizex, int *gap_sizey, int gap_number, int simnumber, double *LC_xcrd, double *LC_ycrd, int numindiv, int *spID, int num_species, int **speciescount, int curr_year, int *gapgroup, int groupgapnum, int *groupgapsize);

//void print_recruits(char *begoutfile, char *simtype, int *gap_xcrd, int *gap_ycrd, int *gap_sizex, int *gap_sizey, int gap_number, int simnumber, double *LC_xcrd, double *LC_ycrd, int numindiv, int *spID, int num_species, int **speciescount, int curr_year, int *gapgroup, int groupgapnum, int *groupgapsize, int numindivstart);

void print_SAD_per_simulation_groupgaps(char *begoutfile, int *gap_xcrd, int *gap_ycrd, int *gap_sizex, int *gap_sizey, int gap_number, int simnumber, double *LC_xcrd, double *LC_ycrd, int numindiv, int *spID, int num_species, double **speciescount, double **speciescountnongap, double **speciescountnongaprecruits, double **speciescountrecruit, int curr_year, int *gapgroup, int groupgapnum, int *groupgapsize, int **nongap_xcrd1, int **nongap_ycrd1, int numindivstart);

void print_SAD_per_simulation_nongapmeta(char *begoutfile, int *gap_xcrd, int *gap_ycrd, int *gap_sizex, int *gap_sizey, int gap_number, int simnumber, double *LC_xcrd, double *LC_ycrd, int numindiv, int *spID, int num_species, int curr_year, int numindivstart, int *speciescountmeta, int *speciescountmetarecruits);

void print_stats(char *begoutfile, int *gap_xcrd, int *gap_ycrd, int *gap_sizex, int *gap_sizey, int gap_number, int simnumber, double *LC_xcrd, double *LC_ycrd, int numindiv, int *spID, int num_species, double **speciescount, double **speciescountnongap, double **speciescountnongaprecruits, double **speciescountrecruit, int *gapgroup, int groupgapnum, int *groupgapsize, int **nongap_xcrd1, int **nongap_ycrd1, int numindivstart, int gen, double ***largearray, double ***smallarray, double ***allarray, double *wd, double *factor, int size, double **array, int mod);

void print_statsfull(char *begoutfile, int simnumber, double *LC_xcrd, double *LC_ycrd, int numindiv, int *spID, int num_species, int **speciescount, int gen, double ***allarray, double *wd, double *factor, double **array);
