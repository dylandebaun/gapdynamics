//
//  august.c
//  gapdynamics
//
//  Created by Dylan DeBaun on 8/6/20.
//  Copyright © 2020 Dylan DeBaun. All rights reserved.
//

#include <stdio.h>
# include <stdlib.h>
# include <time.h>
# include <string.h>
# include "general.h"
#include "initializelandscape.h"
#include "time_keeper.h"
#include <math.h>
#include "stats_analysis.h"
#include "timesteps+dispersal.h"


int main(int argc, char *argv[])
{
    int gapcode=1;
    int mod = 100;
    
    
    int ngen = 31594;//174753; //0.75*J; //change this each time rf=157861
    //number of simulations to run
    //int tenyr = 10375;//78584;//number events by yr 10
    int numdeaths = 15742;//109578;
    int numrecruits = 15852;//92832;
   // int numdeaths10 = 4819;
   // int numrecruits10 = tenyr-numdeaths10;
 //   int numrecruitsrest =numrecruits-numrecruits10;
 //   int numdeathsrest = numdeaths -numdeaths10;
    
    //number of simulations to run
    /*
    int tenyr = 19605+23067+19807+11657 + (1012*2);//number events by yr 10
    int numdeaths = 109351;
    int numrecruits = 93568;
    int ngen = numdeaths +numrecruits; //0.75*J; //change this each time rf=157861
    int numdeaths10 = 19605+19807+1012;
    int numrecruits10 = tenyr-numdeaths10;
    int numrecruitsrest =numrecruits-numrecruits10;
    int numdeathsrest = numdeaths -numdeaths10;
*/
    
    int start;
    int end;
    //INITIALIZING PARAMETERS+ARRAY ALLOCATION
    int size_x = 1000, size_y=500;
    
    int i,j,k;
    double *LC_xcrd, *LC_ycrd;//array of size 1xJ of the coordinate values for all individuals
    double *LC_xcrd_start, *LC_ycrd_start; //keep track of start landscape per simulation
    int *parent_spp_IDstart;
    int *parent_spp_ID; //array of size 1xJ of the species id values for all individuals
    int *age, *age_start; //matrix with age to choose adult
    int ***species_count; //for calculating avg SAD
    double **speciescountpersim; //for calculating SAD per simulation
    //GAP PARAMETERS
    int *gap_xcrd, *gap_ycrd,**nongap_xcrd, **nongap_ycrd, *gap_ycrdmaster,*gap_xcrdmaster;
    int *gap_size_x, *gap_size_y; //rectangular dimensions of gap
    int number_gaps; //total number of gaps
    int *gap_size; //tells if gap is 0=small or 1=large
    
    //ARRAY ALLOCATION
    LC_xcrd = calloc(236000, sizeof(double));
    for(i=0;i<236000; i++){
        LC_xcrd[i]=0;
    }
    LC_ycrd = calloc(236000, sizeof(double));
    for(i=0;i<236000; i++){
        LC_ycrd[i]=0;
    }
    LC_xcrd_start = calloc(236000, sizeof(double));
    for(i=0;i<236000; i++){
        LC_xcrd_start[i]=0;
    }
    LC_ycrd_start = calloc(236000, sizeof(double));
    for(i=0;i<236000; i++){
        LC_ycrd_start[i]=0;
    }
    parent_spp_IDstart = calloc(236000, sizeof(int));
    for(i=0;i<236000; i++){
        parent_spp_IDstart[i]=0;
    }
    age = calloc(236000, sizeof(int));
    for(i=0;i<236000; i++){
        age[i]=0;
    }
    age_start = calloc(236000, sizeof(int));
    for(i=0;i<236000; i++){
        age_start[i]=0;
    }
    parent_spp_ID = calloc(236000, sizeof(int));
    for(i=0;i<236000; i++){
        parent_spp_ID[i]=0;
    }
    species_count = calloc(1000, sizeof(int **));
    for (i = 0; i<1000; ++i) {
        species_count[i] = calloc(1000, sizeof(int *));
        for (j=0; j<1000; ++j) {
            species_count[i][j] = calloc(1000, sizeof(int));
            for(k = 0; k<1000; k++){
                species_count[i][j][k] = 0;
            }
        }
    }
    speciescountpersim = calloc(500, sizeof(double *));
    for(i=0;i<500; i++){
        speciescountpersim[i] = calloc(1000, sizeof(double ));
        for (j = 0; j<1000; ++j) {
            speciescountpersim[i][j] = 0;
        }
    }
    gap_xcrd = calloc(50000, sizeof(int));
    for(i=0;i<50000; i++){
        gap_xcrd[i]=0;
    }
    gap_ycrd = calloc(50000, sizeof(int));
    for(i=0;i<50000; i++){
        gap_ycrd[i]=0;
    }
    gap_xcrdmaster = calloc(50000, sizeof(int));
    for(i=0;i<50000; i++){
        gap_xcrdmaster[i]=0;
    }
    gap_ycrdmaster = calloc(50000, sizeof(int));
    for(i=0;i<50000; i++){
        gap_ycrdmaster[i]=0;
    }
    nongap_xcrd = calloc(20, sizeof(int *));
    for(i=0;i<20; i++){
      //  nongap_xcrd[i]=0;
        nongap_xcrd[i] = calloc(1000, sizeof(int ));
        for (j = 0; j<1000; ++j) {
            nongap_xcrd[i][j] = 0;
        }
    }
    nongap_ycrd = calloc(20, sizeof(int *));
    for(i=0;i<20; i++){
        //  nongap_xcrd[i]=0;
        nongap_ycrd[i] = calloc(1000, sizeof(int ));
        for (j = 0; j<1000; ++j) {
            nongap_ycrd[i][j] = 0;
        }
    }
    gap_size_x = calloc(50000, sizeof(int));
    for(i=0;i<50000; i++){
        gap_size_x[i]=0;
    }
    gap_size_y = calloc(50000, sizeof(int));
    for(i=0;i<50000; i++){
        gap_size_y[i]=0;
    }
    gap_size = calloc(1000, sizeof(int));
    for(i=0;i<1000; i++){
        gap_size[i]= 2;
    }
    int *groupgap;
    int *groupgapsize;
    groupgap = calloc(700, sizeof(int));
    for(i=0;i<700; i++){
        groupgap[i] = -1;
    }
    groupgapsize = calloc(700, sizeof(int));
    for(i=0;i<700; i++){
        groupgapsize[i] = -1;
    }
    int *groupnongap;
    int *groupnongapsize;
    groupnongap = calloc(700, sizeof(int));
    for(i=0;i<700; i++){
        groupnongap[i] = -1;
    }
    groupnongapsize = calloc(700, sizeof(int));
    for(i=0;i<700; i++){
        groupnongapsize[i] = -1;
    }
    double **speciescountnongap, **speciescountrecruit,**speciescountnongaprecruits;
    //, *speciescountmeta, *speciescountmetarecruits;
    speciescountnongap = calloc(500, sizeof(double *));
    for(i=0;i<500; i++){
        speciescountnongap[i] = calloc(1000, sizeof(double ));
        for (j = 0; j<1000; ++j) {
            speciescountnongap[i][j] = 0;
        }
    }
    speciescountnongaprecruits = calloc(500, sizeof(double *));
    for(i=0;i<500; i++){
        speciescountnongaprecruits[i] = calloc(1000, sizeof(double ));
        for (j = 0; j<1000; ++j) {
            speciescountnongaprecruits[i][j] = 0;
        }
    }
    speciescountrecruit = calloc(500, sizeof(double *));
    for(i=0;i<500; i++){
        speciescountrecruit[i] = calloc(1000, sizeof(double ));
        for (j = 0; j<1000; ++j) {
            speciescountrecruit[i][j] = 0;
        }
    }
    /*speciescountmeta = calloc(500, sizeof(int *));
    for(i=0;i<500; i++){
        speciescountmeta[i] = 0;
    }
    speciescountmetarecruits = calloc(500, sizeof(int *));
    for(i=0;i<500; i++){
        speciescountmetarecruits[i] = 0;
    }*/
    
    //INPUT PARAMETERS
    int maturity_age;
    int do_age;
    float increase_percent;
    double kernel_cutoff[2];
    int max_disp_dist[2];
    int gauss[2], fat_tail[2], global[2], near_neigh[2];
    double sigma[2], u[2], p[2];
    double *disc_norm_cum[2];
    int loc_kern_cutoff[2], turnover_time, num_intev_to_meas, time_btwn_meas;
    int output_initial,output_final, output_interm,num_interm, *outputted_interm, lowestx,highestx,lowesty, highesty;
    double *when_output;
    
    //argv1= input file, argv2 = start landscape name, argv3=gap file name, argv4= outfile name
    char inputfile[MAX_FLNAME_SIZE];
    strcpy(inputfile, argv[1]);
    printf("%s\n", inputfile);
    char startlandscape[MAX_FLNAME_SIZE];
    strcpy(startlandscape, argv[2]);
    printf("%s\n", startlandscape);
    char gapfile[MAX_FLNAME_SIZE];
    strcpy(gapfile, argv[3]);
    printf("%s\n", gapfile);
    char outfile[MAX_FLNAME_SIZE];
    strcpy(outfile, argv[4]);
    printf("%s\n", outfile);
    char start1[MAX_FLNAME_SIZE];
    strcpy(start1, argv[5]);
    printf("%s\n", start1);
    char end1[MAX_FLNAME_SIZE];
    strcpy(end1, argv[6]);
    printf("%s\n", end1);
    
    start = atoi(start1);
    end = atoi(end1);
   printf("start value %f ,end value %f.\n", start,end); 
    
    srand(time(NULL));
    //allocate arrays
    // LC_xcrd_start = calloc(500000, sizeof(int));
    
    get_input_info(inputfile, &increase_percent, &size_x, &size_y, kernel_cutoff, max_disp_dist, gauss, sigma, fat_tail, u, p, global, near_neigh, disc_norm_cum, loc_kern_cutoff, &turnover_time, &num_intev_to_meas, &time_btwn_meas, &do_age, &maturity_age, &output_initial, &output_final, &output_interm, &num_interm, &when_output, &outputted_interm, &lowestx, &highestx, &lowesty, &highesty);
    
    //2. INITIALIZE LANDSCAPE WITH EMPIRICAL DATA
    //read in table with columns (species ID, x, y, dbh, dispersal type)
    
    //SIMULATION PARAMETERS
    int J = 0; //# indivs in the landscape
    int num_species;//total number of species present
    double threshd =1.3635; // R in Sedio and Ostling paper. Didn't have much effect for them when ranging between 0-5m.
    //double threshd = 4.0;
    printf("Threshold value %f m.\n", threshd);
    
    read_in_landscape(startlandscape, &LC_xcrd_start, &LC_ycrd_start, &age_start, &maturity_age, &do_age, &J , &parent_spp_IDstart,  &num_species);
    
    double **largegaps, **smallgaps,**allgaps, **all;
    largegaps = calloc(200000, sizeof(double *));
    for(i=0;i<200000; i++){
        largegaps[i] = calloc(50, sizeof(double ));
        for (j = 0; j<50; ++j) {
            largegaps[i][j] = 0;
        }
    }
    smallgaps = calloc(200000, sizeof(double *));
    for(i=0;i<200000; i++){
        smallgaps[i] = calloc(50, sizeof(double ));
        for (j = 0; j<50; ++j) {
            smallgaps[i][j] = 0;
        }
    }
    allgaps = calloc(200000, sizeof(double *));
    for(i=0;i<200000; i++){
        allgaps[i] = calloc(50, sizeof(double ));
        for (j = 0; j<50; ++j) {
            allgaps[i][j] = 0;
        }
    }
    all = calloc(200000, sizeof(double *));
    for(i=0;i<200000; i++){
        all[i] = calloc(50, sizeof(double ));
        for (j = 0; j<50; ++j) {
            all[i][j] = 0;
        }
    }
    double **array;
    array = calloc(236000, sizeof(double *));
    for(i=0;i<236000; i++){
        array[i] = calloc(50, sizeof(double ));
        for (j = 0; j<50; ++j) {
            array[i][j] = -1;
        }
    }
    double **array1;
    array1 = calloc(236000, sizeof(double *));
    for(i=0;i<236000; i++){
        array1[i] = calloc(50, sizeof(double ));
        for (j = 0; j<50; ++j) {
            array1[i][j] = -1;
        }
    }
    int *deathsgap;
    deathsgap = calloc(400, sizeof(int *));
    for(i=0;i<400; i++){
        deathsgap[i] = 0;
    }
    int *deathsnongap;
    deathsnongap = calloc(400, sizeof(int *));
    for(i=0;i<400; i++){
        deathsnongap[i] = 0;
    }
    int *recgap;
    recgap = calloc(400, sizeof(int *));
    for(i=0;i<400; i++){
        recgap[i] = 0;
    }
    int *recnongap;
    recnongap = calloc(400, sizeof(int *));
    for(i=0;i<400; i++){
        recnongap[i] = 0;
    }
    
    
    //sets threshold to average minimum distance between individuals. only use when you change the start landscape because it takes too long to run
   // min_dist_bet_indiv(J,LC_xcrd_start,LC_ycrd_start, &threshd);
    
    //d. identify the gap locations and mark them as small/large
    int number_gapsmaster;
    read_in_gaps(gapfile, &gap_xcrd, &gap_ycrd, &gap_size_x, &gap_size_y, &number_gaps, size_x,size_y);
   // read_in_gaps("mastergaps.txt",&gap_xcrdmaster,&gap_ycrdmaster,&gap_size_x, &gap_size_y, &number_gapsmaster, size_x,size_y);
    
    //print informatiion on gaps and on species ID key
    //print_gaps("85", gap_xcrd, gap_ycrd, gap_size_x, gap_size_y, number_gaps);
    int groupgapnum = 0;
    //print_gaps("85", gap_xcrd, gap_ycrd, gap_size_x, gap_size_y, number_gaps);
    groupgaps(&groupgap, &groupgapsize, &groupgapnum, gap_xcrd, gap_ycrd,number_gaps);
    //group_gaps(&groupgap, &groupgapsize, &groupgapnum, number_gaps);
   // printgroupgapinfo("group", groupgap, groupgapsize, groupgapnum,number_gaps);
    //  print_species_key(outfile, species_name, num_species);
    num_species = 328;
    
    //nongaps(groupgap, groupgapsize, groupgapnum, gap_xcrd, gap_ycrd, number_gaps, &nongap_xcrd, &nongap_ycrd, size_x, size_y);
    
   char nongapgroup[MAX_FLNAME_SIZE];
   // for(i=gapcode; i<gapcode+2;i++){
    for(i=0; i<9;i++){
        /*
         //FOR CREATING SETS OF NONGAPS
        nongapsrandom(groupgap, groupgapsize, groupgapnum, gap_xcrd, gap_ycrd, number_gaps, &nongap_xcrd[1], &nongap_ycrd[1], size_x, size_y);
        sprintf(nongapgroup, "%d85non", i);
        print_gaps(nongapgroup, nongap_xcrd[1], nongap_ycrd[1], groupgapsize, groupgapsize, groupgapnum);
        
        //FOR CREATING SETS OF NONGAPS
        nongapsrandom(groupgap, groupgapsize, groupgapnum, gap_xcrdmaster, gap_ycrdmaster, number_gapsmaster, &nongap_xcrd[1], &nongap_ycrd[1], size_x, size_y);
        sprintf(nongapgroup, "%d85non", i);
        print_gaps(nongapgroup, nongap_xcrd[1], nongap_ycrd[1], groupgapsize, groupgapsize, groupgapnum);
         */
        
         //FOR READING IN NONGAPS
        sprintf(nongapgroup, "%d85nongaplocations.txt", i);
        read_in_nongaps(nongapgroup,&nongap_xcrd[i], &nongap_ycrd[i]);//changed
       /* for(j=0;j<200;j++){
            printf("%d\n",nongap_xcrd[i][j]);
        }*/
    }
    
    //ACTUAL
    //get the actual from the gaps for each year

    int numind = 0;
   // run_actual("start.txt", "/Users/dylandebaun/Desktop/statscluster/ng110cmactual/ng110cm", gap_xcrd, gap_ycrd, gap_size_x, gap_size_y, number_gaps, speciescountpersim, 0, groupgap, groupgapnum, groupgapsize, num_species);
  /*  run_actual("actualyr0.txt", "/Users/dylandebaun/Desktop/statscluster/ng110cmactual/ng110cm", gap_xcrd, gap_ycrd, gap_size_x, gap_size_y, number_gaps, speciescountpersim, 0, groupgap, groupgapnum, groupgapsize, num_species);
        run_actual("10cm85individuals5.txt", "/Users/dylandebaun/Desktop/statscluster/ng110cmactual/ng110cm", gap_xcrd, gap_ycrd, gap_size_x, gap_size_y, number_gaps, speciescountpersim, 5, groupgap, groupgapnum, groupgapsize, num_species);
       run_actual("10cm85individuals10.txt", "/Users/dylandebaun/Desktop/statscluster/ng110cmactual/ng110cm", gap_xcrd, gap_ycrd, gap_size_x, gap_size_y, number_gaps, speciescountpersim, 10, groupgap, groupgapnum, groupgapsize, num_species);
    run_actual("10cm85individuals15.txt", "/Users/dylandebaun/Desktop/statscluster/ng110cmactual/ng110cm", gap_xcrd, gap_ycrd, gap_size_x, gap_size_y, number_gaps, speciescountpersim, 15, groupgap, groupgapnum, groupgapsize, num_species);
        run_actual("10cm85individuals20.txt", "/Users/dylandebaun/Desktop/statscluster/ng110cmactual/ng110cm", gap_xcrd, gap_ycrd, gap_size_x, gap_size_y, number_gaps, speciescountpersim, 20, groupgap, groupgapnum, groupgapsize, num_species);
        run_actual("10cm85individuals25.txt", "/Users/dylandebaun/Desktop/statscluster/ng110cmactual/ng110cm", gap_xcrd, gap_ycrd, gap_size_x, gap_size_y, number_gaps, speciescountpersim, 25, groupgap, groupgapnum, groupgapsize, num_species);
        run_actual("10cm85individuals30.txt", "/Users/dylandebaun/Desktop/statscluster/ng110cmactual/ng110cm", gap_xcrd, gap_ycrd, gap_size_x, gap_size_y, number_gaps, speciescountpersim, 30, groupgap, groupgapnum, groupgapsize, num_species);
   */ /*run_actual("recruitsyr10.txt", "/Users/dylandebaun/Desktop/statscluster/ng110cmactual/recruitsng110cm", gap_xcrd, gap_ycrd, gap_size_x, gap_size_y, number_gaps, speciescountpersim, 10, groupgap, groupgapnum, groupgapsize, num_species);
    run_actualnongap("deathsyr10.txt", "/Users/dylandebaun/Desktop/statscluster/ng110cmactual/deathsnongapng110cm", nongap_xcrd, nongap_ycrd, speciescountpersim, 10, groupgap, groupgapnum, groupgapsize, num_species);*/
   // run_actualnongap("actualyr0.txt", "/Users/dylandebaun/Desktop/statscluster/ng110cmactual/nongapng110cm", nongap_xcrd, nongap_ycrd, speciescountpersim, 0, groupgap, groupgapnum, groupgapsize, num_species);
    //RF
    //gap
    /*
    run_actual("rf85individuals30.txt", "/Users/dylandebaun/Desktop/statscluster/ng1rfactual/ng1rf", gap_xcrd, gap_ycrd, gap_size_x, gap_size_y, number_gaps, speciescountpersim, 30, groupgap, groupgapnum, groupgapsize, num_species);
    run_actual("rf85individuals25.txt", "/Users/dylandebaun/Desktop/statscluster/ng1rfactual/ng1rf", gap_xcrd, gap_ycrd, gap_size_x, gap_size_y, number_gaps, speciescountpersim, 25, groupgap, groupgapnum, groupgapsize, num_species);
    run_actual("rf85individuals20.txt", "/Users/dylandebaun/Desktop/statscluster/ng2rfactual/ng2rf", gap_xcrd, gap_ycrd, gap_size_x, gap_size_y, number_gaps, speciescountpersim, 20, groupgap, groupgapnum, groupgapsize, num_species);
    run_actual("rf85individuals15.txt", "/Users/dylandebaun/Desktop/statscluster/ng2rfactual/ng2rf", gap_xcrd, gap_ycrd, gap_size_x, gap_size_y, number_gaps, speciescountpersim, 15, groupgap, groupgapnum, groupgapsize, num_species);
    run_actual("rf85individuals10.txt", "/Users/dylandebaun/Desktop/statscluster/ng2rfactual/ng2rf", gap_xcrd, gap_ycrd, gap_size_x, gap_size_y, number_gaps, speciescountpersim, 10, groupgap, groupgapnum, groupgapsize, num_species);
    run_actual("rf85individuals5.txt", "/Users/dylandebaun/Desktop/statscluster/ng2rfactual/ng2rf", gap_xcrd, gap_ycrd, gap_size_x, gap_size_y, number_gaps, speciescountpersim, 5, groupgap, groupgapnum, groupgapsize, num_species);
    //gaprecruit
    run_actual("rf85recruits30.txt", "/Users/dylandebaun/Desktop/statscluster/ng2rfactual/recruitng2rf", gap_xcrd, gap_ycrd, gap_size_x, gap_size_y, number_gaps, speciescountpersim, 30, groupgap, groupgapnum, groupgapsize, num_species);
    run_actual("rf85recruits25.txt", "/Users/dylandebaun/Desktop/statscluster/ng2rfactual/recruitng2rf", gap_xcrd, gap_ycrd, gap_size_x, gap_size_y, number_gaps, speciescountpersim, 25, groupgap, groupgapnum, groupgapsize, num_species);
    run_actual("rf85recruits20.txt", "/Users/dylandebaun/Desktop/statscluster/ng2rfactual/recruitng2rf", gap_xcrd, gap_ycrd, gap_size_x, gap_size_y, number_gaps, speciescountpersim, 20, groupgap, groupgapnum, groupgapsize, num_species);
    run_actual("rf85recruits15.txt", "/Users/dylandebaun/Desktop/statscluster/ng2rfactual/recruitng2rf", gap_xcrd, gap_ycrd, gap_size_x, gap_size_y, number_gaps, speciescountpersim, 15, groupgap, groupgapnum, groupgapsize, num_species);
    run_actual("rf85recruits10.txt", "/Users/dylandebaun/Desktop/statscluster/ng2rfactual/recruitng2rf", gap_xcrd, gap_ycrd, gap_size_x, gap_size_y, number_gaps, speciescountpersim, 10, groupgap, groupgapnum, groupgapsize, num_species);
    run_actual("rf85recruits5.txt", "/Users/dylandebaun/Desktop/statscluster/ng2rfactual/recruitng2rf", gap_xcrd, gap_ycrd, gap_size_x, gap_size_y, number_gaps, speciescountpersim, 5, groupgap, groupgapnum, groupgapsize, num_species);
    */
    //10CM
    //gap
   /* run_actual("10cm85individuals30.txt", "/Users/dylandebaun/Desktop/statscluster/ng110cmactual/ng110cm", gap_xcrd, gap_ycrd, gap_size_x, gap_size_y, number_gaps, speciescountpersim, 30, groupgap, groupgapnum, groupgapsize, num_species);
    run_actual("10cm85individuals25.txt", "/Users/dylandebaun/Desktop/statscluster/ng110cmactual/ng110cm", gap_xcrd, gap_ycrd, gap_size_x, gap_size_y, number_gaps, speciescountpersim, 25, groupgap, groupgapnum, groupgapsize, num_species);
    run_actual("10cm85individuals20.txt", "/Users/dylandebaun/Desktop/statscluster/ng110cmactual/ng110cm", gap_xcrd, gap_ycrd, gap_size_x, gap_size_y, number_gaps, speciescountpersim, 20, groupgap, groupgapnum, groupgapsize, num_species);
    run_actual("10cm85individuals15.txt", "/Users/dylandebaun/Desktop/statscluster/ng110cmactual/ng110cm", gap_xcrd, gap_ycrd, gap_size_x, gap_size_y, number_gaps, speciescountpersim, 15, groupgap, groupgapnum, groupgapsize, num_species);
    run_actual("10cm85individuals10.txt", "/Users/dylandebaun/Desktop/statscluster/ng110cmactual/ng110cm", gap_xcrd, gap_ycrd, gap_size_x, gap_size_y, number_gaps, speciescountpersim, 10, groupgap, groupgapnum, groupgapsize, num_species);
    
    run_actual("10cm85individuals5.txt", "/Users/dylandebaun/Desktop/statscluster/ng110cmactual/ng110cm", gap_xcrd, gap_ycrd, gap_size_x, gap_size_y, number_gaps, speciescountpersim, 5, groupgap, groupgapnum, groupgapsize, num_species);
    
    //gaprecruit
    run_actual("10cm85recruits30.txt", "/Users/dylandebaun/Desktop/statscluster/ng110cmactual/recruitng110cm", gap_xcrd, gap_ycrd, gap_size_x, gap_size_y, number_gaps, speciescountpersim, 30, groupgap, groupgapnum, groupgapsize, num_species);
    run_actual("10cm85recruits25.txt", "/Users/dylandebaun/Desktop/statscluster/ng110cmactual/recruitng110cm", gap_xcrd, gap_ycrd, gap_size_x, gap_size_y, number_gaps, speciescountpersim, 25, groupgap, groupgapnum, groupgapsize, num_species);
    run_actual("10cm85recruits20.txt", "/Users/dylandebaun/Desktop/statscluster/ng110cmactual/recruitng110cm", gap_xcrd, gap_ycrd, gap_size_x, gap_size_y, number_gaps, speciescountpersim, 20, groupgap, groupgapnum, groupgapsize, num_species);
    run_actual("10cm85recruits15.txt", "/Users/dylandebaun/Desktop/statscluster/ng110cmactual/recruitng110cm", gap_xcrd, gap_ycrd, gap_size_x, gap_size_y, number_gaps, speciescountpersim, 15, groupgap, groupgapnum, groupgapsize, num_species);
    run_actual("10cm85recruits10.txt", "/Users/dylandebaun/Desktop/statscluster/ng110cmactual/recruitng110cm", gap_xcrd, gap_ycrd, gap_size_x, gap_size_y, number_gaps, speciescountpersim, 10, groupgap, groupgapnum, groupgapsize, num_species);
    run_actual("10cm85recruits5.txt", "/Users/dylandebaun/Desktop/statscluster/ng110cmactual/recruitng110cm", gap_xcrd, gap_ycrd, gap_size_x, gap_size_y, number_gaps, speciescountpersim, 5, groupgap, groupgapnum, groupgapsize, num_species);
    */
    //10CM
    //nongap
   /*run_actualnongap("10cm85individuals30.txt", "/Users/dylandebaun/Desktop/statscluster/ng110cmactual/nongapng110cm", nongap_xcrd, nongap_ycrd, speciescountpersim, 30, groupgap, groupgapnum, groupgapsize, num_species);
    run_actualnongap("10cm85individuals25.txt", "/Users/dylandebaun/Desktop/statscluster/ng110cmactual/nongapng110cm", nongap_xcrd, nongap_ycrd, speciescountpersim, 25, groupgap, groupgapnum, groupgapsize, num_species);
    run_actualnongap("10cm85individuals20.txt", "/Users/dylandebaun/Desktop/statscluster/ng110cmactual/nongapng110cm", nongap_xcrd, nongap_ycrd, speciescountpersim, 20, groupgap, groupgapnum, groupgapsize, num_species);
    run_actualnongap("10cm85individuals15.txt", "/Users/dylandebaun/Desktop/statscluster/ng110cmactual/nongapng110cm", nongap_xcrd, nongap_ycrd, speciescountpersim, 15, groupgap, groupgapnum, groupgapsize, num_species);
    run_actualnongap("10cm85individuals10.txt", "/Users/dylandebaun/Desktop/statscluster/ng110cmactual/nongapng110cm", nongap_xcrd, nongap_ycrd,speciescountpersim, 10, groupgap, groupgapnum, groupgapsize, num_species);
    run_actualnongap("10cm85individuals5.txt", "/Users/dylandebaun/Desktop/statscluster/ng110cmactual/nongapng110cm", nongap_xcrd, nongap_ycrd, speciescountpersim, 5, groupgap, groupgapnum, groupgapsize, num_species);*/
    /*//nongaprecruit
    run_actualnongap("10cm85recruits30.txt", "/Users/dylandebaun/Desktop/statscluster/ng110cmactual/nongaprecruitng110cm", nongap_xcrd, nongap_ycrd, speciescountpersim, 30, groupgap, groupgapnum, groupgapsize, num_species);
    run_actualnongap("10cm85recruits25.txt", "/Users/dylandebaun/Desktop/statscluster/ng110cmactual/nongaprecruitng110cm", nongap_xcrd, nongap_ycrd, speciescountpersim, 25, groupgap, groupgapnum, groupgapsize, num_species);
    run_actualnongap("10cm85recruits20.txt", "/Users/dylandebaun/Desktop/statscluster/ng110cmactual/nongaprecruitng110cm", nongap_xcrd, nongap_ycrd, speciescountpersim, 20, groupgap, groupgapnum, groupgapsize, num_species);
    run_actualnongap("10cm85recruits15.txt", "/Users/dylandebaun/Desktop/statscluster/ng110cmactual/nongaprecruitng110cm", nongap_xcrd, nongap_ycrd, speciescountpersim, 15, groupgap, groupgapnum, groupgapsize, num_species);
    run_actualnongap("10cm85recruits10.txt", "/Users/dylandebaun/Desktop/statscluster/ng110cmactual/nongaprecruitng110cm", nongap_xcrd, nongap_ycrd, speciescountpersim, 10, groupgap, groupgapnum, groupgapsize, num_species);
    run_actualnongap("10cm85recruits5.txt", "/Users/dylandebaun/Desktop/statscluster/ng110cmactual/nongaprecruitng110cm", nongap_xcrd, nongap_ycrd, speciescountpersim, 5, groupgap, groupgapnum, groupgapsize, num_species);
 */
    printf("done with actual\n");
    //simtype,recruittype,method,"expectedingap",(gaplocations$gapnum[g]-1), "at",year,"yrs.txt"
    //1. read in data for trait values
    FILE *ifp;
    int status;
    int numspeciestemp;
    double *wd, *factor;
    double wood, fac;
    int increment =0;
    
    wd = calloc(236000, sizeof(double));
    for(i=0;i<236000; i++){
        wd[i]=0;
    }
    factor = calloc(236000, sizeof(double));
    for(i=0;i<236000; i++){
        factor[i]=0;
    }
    
    numspeciestemp = 0;
    
    ifp = fopen("wdfactor.txt", "r");
    status=fscanf(ifp, "%lf %lf \n", &fac, &wood);
    i=0;
    
    //fflush(ifp);
    while (status != EOF) {
        wd[i] = wood;
        factor[i] = fac;
        status=fscanf(ifp, "%lf %lf \n", &fac, &wood);
        i++;
       // fflush(ifp);
    }

    //run_actualstats
    
    //read in actual landscape
    
    //calculate in gap stuff
    
    //3. SIMULATION
    /*
     This section holds actual community dynamics.
     1) Random individual is chosen to experience death event.
     2) Must choose whether replacement is from immigration or local community.
     3) Replacement individual is chosen. Either randomly from LC (local replacement), or from regional pool (immigration).
     4) Determine whether or not offspring successful establishes itself. If not, go back to (3).
     */
    
    // DYNAMICS VARIABLES.
    int s, count;
    int rbinom_check =0; //add in 11.19
    int x;
    int my_gen; // what time step is it.
    int death_choice; // who dies.
    int repl_choice; // who reproduces.
    int repl_choicecoords[2];
    //   double imm_choice; // who immigrates into LC.
    double toss;
    
    double min_dist = 99999; // setting arbitrarily high.
    double cur_dist; // current distance of closest neighbor.
    
    int offspringcoords[2];
    double offspring_X, offspring_Y; // coordinates of new recruit within plot.
    int offspring_ID; // spp to which recruit belongs.
    int b; // test for successful birth.
    int curr_year;
    int *eventscape, *eventscapestart;
    
    //int hypbardeath = 0, hypbarbirth = 0;
    
    
    //initialize the eventscape array which tells whehter there is a death or recruitment event. making sure we have accurate number of d/r in each five year interval
    //are the appropriate number of deaths and recruitments
        int Jstart = J;
    int d; //event number for death/recruits
    int numdeaths5 = 2650, fiveyr = 5818, numdeaths10 = 2320, tenyr = 10677, numdeaths15 = 2525, fifteenyr = 15476, numdeaths20 = 2747, twentyyr = 20615, numdeaths25 = 2640, twentyfiveyr = 25846, numdeathsrest = 2860;
   int numdeaths101 =0, numdeaths10to30=0, numrecruits101=0, numrecruits10to30=0;
    eventscapestart = calloc(ngen, sizeof(int));
    for(i=0;i<numdeaths5; i++){
        eventscapestart[i] = 0;
    }
    for(i=(numdeaths5); i<fiveyr; i++){
        eventscapestart[i] = 1;
    }
    for(i=fiveyr; i<(fiveyr+numdeaths10); i++){
        eventscapestart[i] = 0;
    }
    for(i=(fiveyr+numdeaths10); i<tenyr; i++){
        eventscapestart[i] = 1;
    }
    for(i=tenyr; i<(tenyr+numdeaths15); i++){
        eventscapestart[i] = 0;
    }
    for(i=(tenyr+numdeaths15); i<fifteenyr; i++){
        eventscapestart[i] = 1;
    }
    for(i=fifteenyr; i<(fifteenyr+numdeaths20); i++){
        eventscapestart[i] = 0;
    }
    for(i=(fifteenyr+numdeaths20); i<twentyyr; i++){
        eventscapestart[i] = 1;
    }
    for(i=twentyyr; i<(twentyyr+numdeaths25); i++){
        eventscapestart[i] = 0;
    }
    for(i=(twentyyr+numdeaths25); i<twentyfiveyr; i++){
        eventscapestart[i] = 1;
    }
    for(i=twentyfiveyr; i<(twentyfiveyr+numdeathsrest); i++){
        eventscapestart[i] = 0;
        numdeaths10to30++;
    }
    for(i=(twentyfiveyr+numdeathsrest); i<ngen; i++){
        eventscapestart[i] = 1;
        numrecruits10to30++;
    }
    
    int r =numrecruits10to30+numrecruits101; //should equal numrecruits
    i = numdeaths10to30 + numdeaths101; //should equal numdeaths
    int n = r+i; //should equal ngen
    
    
    eventscape = calloc(ngen, sizeof(int));
    for(i=0;i<ngen; i++){
        eventscape[i] = 0;
    }
    int q=0;
    for(i = 0;i <ngen;i++){
        if(i <tenyr){
            if(eventscape[i] == 0){
                q++;
            }
        }
    }
    int countem,a;
    int deaths = 0, recruits =0;
    int o,count1 =0,count2=0;
    //127 to 250
    char outfile1[MAX_FLNAME_SIZE];
    FILE *ofp;
    char *metacommunity;
    for(s=start;s<end;s++){
        count2=0;
        //a. Reset landscape to starting configuration
        for(x=0; x< (Jstart+numrecruits); x++){
            LC_xcrd[x]=LC_xcrd_start[x];
            LC_ycrd[x]=LC_ycrd_start[x];
            // age[x] = age_start[x];
            parent_spp_ID[x] = parent_spp_IDstart[x];//added this in; might have been a major mistake about earlier trials?
            //printf("%f,%f\n", LC_xcrd[x],LC_ycrd[x]);
        }
        
        //check that we will have correct number of deaths and recruits
        for(x=0; x<ngen;x++){
            eventscape[x]=eventscapestart[x];
            if(eventscape[x] == 0){deaths++;}else{recruits++;}
        }
        for(x=0;x<400;x++){
            deathsnongap[x] = 0;
            recnongap[x] =0;
            deathsgap[x] = 0;
            recgap[x] =0;
        }
        curr_year=0;
        J = Jstart;
        
        //large JANUARY
      //  print_stats(outfile, gap_xcrd,  gap_ycrd, gap_size_x, gap_size_y, number_gaps, s, LC_xcrd, LC_ycrd, J, parent_spp_ID,num_species,speciescountpersim, speciescountnongap, speciescountnongaprecruits, speciescountrecruit, groupgap, groupgapnum, groupgapsize, nongap_xcrd, nongap_ycrd, Jstart+1, 0, &largegaps,&smallgaps,&allgaps, wd, factor, 1,array1,mod);
        //small
      //  print_stats(outfile, gap_xcrd,  gap_ycrd, gap_size_x, gap_size_y, number_gaps, s, LC_xcrd, LC_ycrd, J, parent_spp_ID,num_species,speciescountpersim, speciescountnongap, speciescountnongaprecruits, speciescountrecruit, groupgap, groupgapnum, groupgapsize, nongap_xcrd, nongap_ycrd, Jstart+1, 0, &largegaps,&smallgaps,&allgaps, wd, factor, 2,array1,mod);
        //  printf("small\n");
        //all
     //   print_stats(outfile, gap_xcrd,  gap_ycrd, gap_size_x, gap_size_y, number_gaps, s, LC_xcrd, LC_ycrd, J, parent_spp_ID,num_species,speciescountpersim, speciescountnongap, speciescountnongaprecruits, speciescountrecruit, groupgap, groupgapnum, groupgapsize, nongap_xcrd, nongap_ycrd, Jstart+1, 0, &largegaps,&smallgaps,&allgaps, wd, factor, 3,array1,mod);
        // printf("all\n");
         
       //print_statsfull(outfile, s, LC_xcrd, LC_ycrd, J, parent_spp_ID, num_species, speciescountpersim, 0, &all, wd, factor, array1);
        
        //one generation = J deaths
        //ngen = number of deaths
        for(my_gen = 1; my_gen <= ngen; my_gen++){ //change <=ngen
            d=0;
            //PICK IF DEATH EVENT OR RECRUIT EVENT
            if(my_gen <= fiveyr){
                d = get_rand_integ_intvl(1,fiveyr);
                while(eventscape[d-1] == -1){
                    d = get_rand_integ_intvl(1,fiveyr);
                }
            }else if(my_gen <= tenyr){
                d = get_rand_integ_intvl(fiveyr+1,tenyr);
                while(eventscape[d-1] == -1){
                    d = get_rand_integ_intvl(fiveyr+1,tenyr);
                }
            }else if(my_gen <= fifteenyr){
                d = get_rand_integ_intvl(tenyr+1,fifteenyr);
                while(eventscape[d-1] == -1){
                    d = get_rand_integ_intvl(tenyr+1,fifteenyr);
                }
            }else if(my_gen <= twentyyr){
                d = get_rand_integ_intvl(fifteenyr+1,twentyyr);
                while(eventscape[d-1] == -1){
                    d = get_rand_integ_intvl(fifteenyr+1,twentyyr);
                }
            }else if(my_gen <= twentyfiveyr){
                twentyyr= get_rand_integ_intvl(twentyyr+1,twentyfiveyr);
                while(eventscape[d-1] == -1){
                    d = get_rand_integ_intvl(twentyyr+1,twentyfiveyr);
                }
            }else{
                d = get_rand_integ_intvl(twentyfiveyr+1,ngen);
                while(eventscape[d-1] == -1){
                    d = get_rand_integ_intvl(twentyfiveyr+1,ngen);
                }
            }
            if(eventscape[d-1] == 0){
                count2++;
                // For death events, choose an index between 1 & J, inclusive. That member of LC will die.
                toss = get_rand_unit(); // U(0,1).
                toss *= J; // scale by LC size.
                death_choice = (int) ceil(toss); // convert to index.
                //make sure we are picking an alive individual
                countem = 0;
                while(LC_xcrd[ death_choice ] == -1){
                    toss = get_rand_unit(); // U(0,1).
                    toss *= J; // scale by LC size.
                    death_choice = (int) ceil(toss); // convert to index.
                    countem++;
                    /*if(countem == 1000){
                        for(x=0; x< (J); x++){
                            printf("%f,%f\n", LC_xcrd[x],LC_ycrd[x]);
                        }
                    }*/
                }
              //  ingapornongap(LC_xcrd[ death_choice ], LC_ycrd[ death_choice ], parent_spp_ID[ death_choice ], gap_xcrd, gap_ycrd, number_gaps, &deathsgap, nongap_xcrd[0], nongap_ycrd[0], groupgapnum, groupgapsize, &deathsnongap);
                for(x=0;x<400;x++){
                    if(deathsgap[x] > 0){
                        o = deathsgap[x];
                    }
                }
                LC_xcrd[ death_choice ] = -1;
                LC_ycrd[ death_choice ] = -1;
                parent_spp_ID[ death_choice ] = -1;
                
                
                count1= 0;
                for(j=0; j<J; j++){
                    o = LC_xcrd[j];
                    //printf("%f\n",LC_xcrd[j]);
                    if(o == -1){
                        count1++;
                    }
                }
                //printf("J=%d,deaths = %d, numdeathevents = %d\n",J,count1,count2);
            }else{
                // Do replacements here.
                b = 1;
                while(b){
                    count =0;
                   // printf("%d\n",count);
                    // ******************* Come back to this point if recruitment is unsuccessful. ******************************
                    //reset min dist every time we are trying to create a new individual
                    min_dist = 99999;
                    // local replacement
                    
                    // Choose individual to reproduce.
                    toss = get_rand_unit(); // U(0,1).
                    toss *= J; // scale by LC size.
                    repl_choice = (int) ceil(toss); // convert to index.
                    
                    while(LC_xcrd[ repl_choice ] == -1){
                        toss = get_rand_unit(); // U(0,1).
                        toss *= J; // scale by LC size.
                        repl_choice = (int) ceil(toss); // convert to index.
                    }
                    
                    repl_choicecoords[0]= LC_xcrd[ repl_choice];
                    repl_choicecoords[1]= LC_ycrd[ repl_choice];
                    
                    select_randloc_frm_cum_disp(repl_choicecoords, size_x, size_y, disc_norm_cum[0], loc_kern_cutoff[0], offspringcoords);
                    
                    offspring_X = offspringcoords[0];
                    offspring_Y = offspringcoords[1];
                    
                                        /*
                    // Eliminate edge effects by making plot taurus. moved to within select_randloc_frm_cum_disp function
                    if( offspring_X < 0 ){ // i.e. x-coordinate is a negative number.
                        offspring_X = size_x + offspring_X;
                    }
                    if( offspring_X >= size_x ){
                        offspring_X = offspring_X - size_x;
                    }
                    if( offspring_Y < 0 ){ // i.e. y-coordinate is a positive number.
                        offspring_Y = size_y + offspring_Y;
                    }
                    if( offspring_Y >= size_y ){
                        offspring_Y = offspring_Y - size_y;
                    }
                    */
                    
                    
                    // Assign spp identity. Temporary assignment as recuitment might not be successful.
                    offspring_ID = parent_spp_ID[ repl_choice ];
                    // printf( "Offspring is from spp. %d.\n", offspring_ID );
                    // Determine success of recruit.
                    // For each individual in the plot:
                    // Is that individual at least tphreshold distance away from the potential recruit? If not, stop & go back to reproduction section above.
                    // Calculate the density of enemies around the recruit due to the plot individual.
                    // Calculate the impact on the recruit due to the overlap in shared enemies between the recruit and the plot individual.
                    // Update the running total of this weighted enemy-mediated factor (big H, Eq. 4 in paper).
                   // count = 0;
                    for( k = 0; k < J; k++ ){
                        
                        cur_dist = sqrt( pow( (LC_xcrd[k] - offspring_X), 2 ) + pow( (LC_ycrd[k] - offspring_Y), 2 ) );
                        
                        if( cur_dist < threshd){//} && age[k] >= 10){ // if someone is too close, just stop immediately and go back to choosing replacement.
                            break;
                        }
                        if( cur_dist < min_dist ){ // Need to know after this loop whether everyone was far enough away.
                            // printf( " ...updating min. distance... \n" );
                            min_dist = cur_dist;
                        }

                    }
                    // So long as all individuals are outside threshold distance, we can then flip a coin to decide success (Eq. 5 in paper).
                    if( k >= J & min_dist > threshd ){ // Was able to run through the entire population before breaking and everyone was far enough away.
                        //add in 11.19
                        toss = get_rand_unit(); // U(0,1).

                       // if( toss >= 1/min_dist){
                            //end of add in
                        LC_xcrd[ J+1 ] = offspring_X;
                        LC_ycrd[ J+1 ] = offspring_Y;
                        parent_spp_ID[ J+1 ] = offspring_ID;
                        J++;
                       // ingapornongap(offspring_X, offspring_Y, offspring_ID, gap_xcrd, gap_ycrd, number_gaps, &recgap, nongap_xcrd[0], nongap_ycrd[0], groupgapnum, groupgapsize, &recnongap);
                        //age[death_choice] = 0;
                        
                        if(my_gen%1000 == 0){
                            printf( "Successful establishment @ generation: %d. min_dist is %f. \n", my_gen, min_dist );
                        }
                        b = 0; // exit while loop & proceed to next generation.
                       // printf("pass");
                      //  }
                        
                    }
                    else{
                        count++;
                    }
                    
            }
               // count++;
            }
            eventscape[d-1] = -1; //put at the end
            //Print out @ 10,20,30,40 year intervals.
            
            if(my_gen%mod == 0){
                /*if(my_gen == 5800){
                    printf("pause");
                }*/
                //large JANUARY
             //   print_stats(outfile, gap_xcrd,  gap_ycrd, gap_size_x, gap_size_y, number_gaps, s, LC_xcrd, LC_ycrd, J, parent_spp_ID,num_species,speciescountpersim, speciescountnongap, speciescountnongaprecruits, speciescountrecruit, groupgap, groupgapnum, groupgapsize, nongap_xcrd, nongap_ycrd, Jstart+1, my_gen, &largegaps,&smallgaps,&allgaps, wd, factor, 1,array1,mod);
                //printf("large\n");
                //small
              //  print_stats(outfile, gap_xcrd,  gap_ycrd, gap_size_x, gap_size_y, number_gaps, s, LC_xcrd, LC_ycrd, J, parent_spp_ID,num_species,speciescountpersim, speciescountnongap, speciescountnongaprecruits, speciescountrecruit, groupgap, groupgapnum, groupgapsize, nongap_xcrd, nongap_ycrd, Jstart+1, my_gen, &largegaps,&smallgaps,&allgaps, wd, factor, 2,array1,mod);
              //  printf("small\n");
            
                //all
             //   print_stats(outfile, gap_xcrd,  gap_ycrd, gap_size_x, gap_size_y, number_gaps, s, LC_xcrd, LC_ycrd, J, parent_spp_ID,num_species,speciescountpersim, speciescountnongap, speciescountnongaprecruits, speciescountrecruit, groupgap, groupgapnum, groupgapsize, nongap_xcrd, nongap_ycrd, Jstart+1, my_gen, &largegaps,&smallgaps,&allgaps, wd, factor, 3,array1,mod);
               // printf("all\n");
                
                //print_statsfull(outfile, s, LC_xcrd, LC_ycrd, J, parent_spp_ID, num_species, speciescountpersim, my_gen, &all, wd, factor, array1);
            }
printf("%f\n", my_gen);
            
            if(my_gen == tenyr){//(ngen*10/35)){ rf= 44637
                curr_year = 10;
                /*sprintf(outfile1, "%ssimulation%d%sdistributionsatyr%d.txt", outfile,s, "deathrec",curr_year);
                ofp = fopen(outfile1, "w");
                for(i =0; i <num_species; i++){
                    fprintf(ofp, "%d\t%d\t%d\t%d\t%d\n", i, deathsgap[i],deathsnongap[i], recgap[i],recnongap[i]);
                }
                fclose(ofp);*/
               print_SAD_per_simulation_groupgaps(outfile, gap_xcrd, gap_ycrd, gap_size_x, gap_size_y, number_gaps, s, LC_xcrd, LC_ycrd, J, parent_spp_ID,num_species,speciescountpersim, speciescountnongap,speciescountnongaprecruits,speciescountrecruit,curr_year, groupgap, groupgapnum, groupgapsize, nongap_xcrd, nongap_ycrd, Jstart+1);
               /* if(s<100){
                print_SAD_per_simulation_nongapmeta(outfile, gap_xcrd, gap_ycrd, gap_size_x, gap_size_y, number_gaps, s, LC_xcrd, LC_ycrd, J, parent_spp_ID,num_species, curr_year, Jstart, speciescountmeta,speciescountmetarecruits);
                }*/
                //add in nongap output here
                
            }
            if(my_gen == ngen){
                curr_year = 30;
               print_SAD_per_simulation_groupgaps(outfile, gap_xcrd, gap_ycrd, gap_size_x, gap_size_y, number_gaps, s, LC_xcrd, LC_ycrd, J, parent_spp_ID,num_species,speciescountpersim, speciescountnongap,speciescountnongaprecruits,speciescountrecruit,curr_year, groupgap, groupgapnum, groupgapsize, nongap_xcrd, nongap_ycrd, Jstart+1);
             
            }
            if(my_gen == 5818){
                curr_year = 5;
                print_SAD_per_simulation_groupgaps(outfile, gap_xcrd, gap_ycrd, gap_size_x, gap_size_y, number_gaps, s, LC_xcrd, LC_ycrd, J, parent_spp_ID,num_species,speciescountpersim, speciescountnongap,speciescountnongaprecruits,speciescountrecruit,curr_year, groupgap, groupgapnum, groupgapsize, nongap_xcrd, nongap_ycrd, Jstart+1);
                
            }
            if(my_gen == 15400){
                curr_year = 15;
                print_SAD_per_simulation_groupgaps(outfile, gap_xcrd, gap_ycrd, gap_size_x, gap_size_y, number_gaps, s, LC_xcrd, LC_ycrd, J, parent_spp_ID,num_species,speciescountpersim, speciescountnongap,speciescountnongaprecruits,speciescountrecruit,curr_year, groupgap, groupgapnum, groupgapsize, nongap_xcrd, nongap_ycrd, Jstart+1);
                
            }
            if(my_gen == 20600){
                curr_year = 20;
                print_SAD_per_simulation_groupgaps(outfile, gap_xcrd, gap_ycrd, gap_size_x, gap_size_y, number_gaps, s, LC_xcrd, LC_ycrd, J, parent_spp_ID,num_species,speciescountpersim, speciescountnongap,speciescountnongaprecruits,speciescountrecruit,curr_year, groupgap, groupgapnum, groupgapsize, nongap_xcrd, nongap_ycrd, Jstart+1);
            }
            if(my_gen == 25800){
                curr_year = 25;
                print_SAD_per_simulation_groupgaps(outfile, gap_xcrd, gap_ycrd, gap_size_x, gap_size_y, number_gaps, s, LC_xcrd, LC_ycrd, J, parent_spp_ID,num_species,speciescountpersim, speciescountnongap,speciescountnongaprecruits,speciescountrecruit,curr_year, groupgap, groupgapnum, groupgapsize, nongap_xcrd, nongap_ycrd, Jstart+1);
            }
        }
        printf( "Done with ngen loop. \n" );
       /* int z,a;
        char outfile1[MAX_FLNAME_SIZE];
        FILE *ofp;
        char *metacommunity;
        
       double k; //JANUARY
        for(a =0; a<3;a++){
            if(a == 1){
            metacommunity = "large";
            array = largegaps;
            }else if(a==2){
                metacommunity = "small";
                array = smallgaps;
            }else{
                metacommunity = "all";
                array = allgaps;
            }
            sprintf(outfile1, "%ssimulation%d%sdistributions.txt", outfile,s, metacommunity);
            ofp = fopen(outfile1, "w");
            //upon getting all q gaps in the group gap tallied, write the group gap to a file
            k = ngen/mod;
            for(i =0; i <ngen/mod; i++){
                //if(i%10 == 0){
                    //for(z=0;z<num_species; z++){
                    fprintf(ofp, "%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", i, array[i][0],array[i][1], array[i][2], array[i][3], array[i][4], array[i][5], array[i][6], array[i][7], array[i][8], array[i][9], array[i][10],array[i][11],  array[i][12], array[i][13], array[i][14], array[i][15], array[i][16], array[i][17], array[i][18], array[i][19], array[i][20],array[i][21],  array[i][22], array[i][23]);
                    //}
                //}
            }
            fclose(ofp);
        }
        */
        /*
            metacommunity = "all";
            array = all;
            sprintf(outfile1, "%ssimulation%d%sFULLPLOT.txt", outfile,s, metacommunity);
            ofp = fopen(outfile1, "w");
            //upon getting all q gaps in the group gap tallied, write the group gap to a file
            for(i =0; i <ngen/mod; i++){
                    //for(z=0;z<num_species; z++){
                //if(i%100 == 0){
                fprintf(ofp, "%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", i, array[i][0],array[i][1], array[i][2], array[i][3], array[i][4], array[i][5]);
                    //}
               // }
            }
            fclose(ofp);*/
        
    }
    printf("done with simulation\n");
    
}



