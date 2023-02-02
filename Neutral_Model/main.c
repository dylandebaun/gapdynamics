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
    int ngen = 31594; // deaths+recruits
    int numdeaths = 15742; 
    int numrecruits = 15852; 
    
    int start; //simulation nums
    int end;
    //INITIALIZING PARAMETERS+ARRAY ALLOCATION
    int size_x = 1000, size_y=500;
    
    int i,j,k;
    double *LC_xcrd, *LC_ycrd;//array of size 1xJ of the coordinate values for all individuals
    double *LC_xcrd_start, *LC_ycrd_start; //keep track of start landscape per simulation
    int *parent_spp_IDstart;
    int *parent_spp_ID; //array of size 1xJ of the species id values for all individuals
    double **speciescountpersim; //for calculating composition of species per simulation
    
    //GAP PARAMETERS
    int *gap_xcrd, *gap_ycrd,**nongap_xcrd, **nongap_ycrd;
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
    parent_spp_ID = calloc(236000, sizeof(int));
    for(i=0;i<236000; i++){
        parent_spp_ID[i]=0;
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
    double **speciescountnongap;
    //, *speciescountmeta, *speciescountmetarecruits;
    speciescountnongap = calloc(500, sizeof(double *));
    for(i=0;i<500; i++){
        speciescountnongap[i] = calloc(1000, sizeof(double ));
        for (j = 0; j<1000; ++j) {
            speciescountnongap[i][j] = 0;
        }
    }
    
    //INPUT PARAMETERS
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
    
    get_input_info(inputfile, &increase_percent, &size_x, &size_y, kernel_cutoff, max_disp_dist, gauss, sigma, fat_tail, u, p, global, near_neigh, disc_norm_cum, loc_kern_cutoff, &turnover_time, &num_intev_to_meas, &time_btwn_meas, &output_initial, &output_final, &output_interm, &num_interm, &when_output, &outputted_interm, &lowestx, &highestx, &lowesty, &highesty);
    
    //2. INITIALIZE LANDSCAPE WITH EMPIRICAL DATA
    //read in table with columns (species ID, x, y, dbh, dispersal type)
    
    //SIMULATION PARAMETERS
    int J = 0; //# indivs in the landscape
    int num_species =328;//total number of species present
    double threshd =1.3635; // R in Sedio and Ostling paper. Didn't have much effect for them when ranging between 0-5m.
    //double threshd = 4.0;
    printf("Threshold value %f m.\n", threshd);
    
    read_in_landscape(startlandscape, &LC_xcrd_start, &LC_ycrd_start, &J , &parent_spp_IDstart,  &num_species);
    
    //sets threshold to average minimum distance between individuals. only use when you change the start landscape because it takes too long to run
    min_dist_bet_indiv(J,LC_xcrd_start,LC_ycrd_start, &threshd);
    
    //identify the gap locations and mark them as small/large
    read_in_gaps(gapfile, &gap_xcrd, &gap_ycrd, &gap_size_x, &gap_size_y, &number_gaps, size_x,size_y);

	//Determine the number of 5x5 gaps that are touching to get gap 'size' this was used in the present paper to identify non-gap regions of similar size
    int groupgapnum = 0;
    groupgaps(&groupgap, &groupgapsize, &groupgapnum, gap_xcrd, gap_ycrd,number_gaps);
    
   char nongapgroup[MAX_FLNAME_SIZE];

    for(i=0; i<9;i++){
    	/*
    	//FOR CREATING SETS OF NONGAPS. NOTE: must mechanically alter these results after running for no overlap.
        nongapsrandom(groupgap, groupgapsize, groupgapnum, gap_xcrd, gap_ycrd, number_gaps, &nongap_xcrd[1], &nongap_ycrd[1], size_x, size_y);
        sprintf(nongapgroup, "%d85non", i);
        print_gaps(nongapgroup, nongap_xcrd[1], nongap_ycrd[1], groupgapsize, groupgapsize, groupgapnum);
         */
         //FOR READING IN NONGAPS
        sprintf(nongapgroup, "NonGap_Locations/%d_1985_nongaplocations.txt", i);
        read_in_nongaps(nongapgroup,&nongap_xcrd[i], &nongap_ycrd[i]);
    }
    
    //ACTUAL- get the observed compositions for each year
    //gap
    run_actual("Observed_Community_Values/10cm85individuals30.txt", outfile, gap_xcrd, gap_ycrd, gap_size_x, gap_size_y, number_gaps, 30, groupgap, groupgapnum, groupgapsize, num_species);
    run_actual("Observed_Community_Values/10cm85individuals25.txt", outfile, gap_xcrd, gap_ycrd, gap_size_x, gap_size_y, number_gaps,  25, groupgap, groupgapnum, groupgapsize, num_species);
    run_actual("Observed_Community_Values/10cm85individuals20.txt", outfile, gap_xcrd, gap_ycrd, gap_size_x, gap_size_y, number_gaps, 20, groupgap, groupgapnum, groupgapsize, num_species);
    run_actual("Observed_Community_Values/10cm85individuals15.txt", outfile, gap_xcrd, gap_ycrd, gap_size_x, gap_size_y, number_gaps, 15, groupgap, groupgapnum, groupgapsize, num_species);
    run_actual("Observed_Community_Values/10cm85individuals10.txt", outfile, gap_xcrd, gap_ycrd, gap_size_x, gap_size_y, number_gaps, 10, groupgap, groupgapnum, groupgapsize, num_species);
	run_actual("Observed_Community_Values/10cm85individuals5.txt", outfile, gap_xcrd, gap_ycrd, gap_size_x, gap_size_y, number_gaps, 5, groupgap, groupgapnum, groupgapsize, num_species);
    //nongap
    run_actualnongap("Observed_Community_Values/10cm85individuals30.txt", outfile, nongap_xcrd, nongap_ycrd,  30, groupgap, groupgapnum, groupgapsize, num_species);
    run_actualnongap("Observed_Community_Values/10cm85individuals25.txt", outfile, nongap_xcrd, nongap_ycrd,  25, groupgap, groupgapnum, groupgapsize, num_species);
    run_actualnongap("Observed_Community_Values/10cm85individuals20.txt", outfile, nongap_xcrd, nongap_ycrd,  20, groupgap, groupgapnum, groupgapsize, num_species);
    run_actualnongap("Observed_Community_Values/10cm85individuals15.txt", outfile, nongap_xcrd, nongap_ycrd,  15, groupgap, groupgapnum, groupgapsize, num_species);
    run_actualnongap("Observed_Community_Values/10cm85individuals10.txt", outfile, nongap_xcrd, nongap_ycrd, 10, groupgap, groupgapnum, groupgapsize, num_species);
    run_actualnongap("Observed_Community_Values/10cm85individuals5.txt", outfile, nongap_xcrd, nongap_ycrd,  5, groupgap, groupgapnum, groupgapsize, num_species);
    printf("done with actual\n");

    FILE *ifp;
    int status;
    
    //SIMULATION
    /*
     This section holds actual community dynamics.
     1) Random individual is chosen to experience death event.
     2) Must choose whether replacement is from immigration or local community.
     3) Replacement individual is chosen. Either randomly from LC (local replacement), or from regional pool (immigration).
     4) Determine whether or not offspring successful establishes itself. If not, go back to (3).
     */
    
    // DYNAMICS VARIABLES.
    int s, count; //simulation and count
    int x; //iterative variables
    int my_gen; // what time step is it.
    int death_choice; // who dies.
    int repl_choice; // who reproduces.
    int repl_choicecoords[2];
    double toss; //coin toss value
    
    double min_dist = 99999; // min distance to a neighbor, setting arbitrarily high.
    double cur_dist; // current distance of closest neighbor.
    
    int offspringcoords[2];
    double offspring_X, offspring_Y; // coordinates of new recruit within plot.
    int offspring_ID; // spp to which recruit belongs.
    int b; // test for successful birth.
    int curr_year;
    int *eventscape, *eventscapestart; //1xngen array of events, death =0, recruit = 1
    
    //int hypbardeath = 0, hypbarbirth = 0;
    
    
    //initialize the eventscape array which tells whehter there is a death or recruitment event. making sure we have accurate number of d/r in each five year interval
    //are the appropriate number of deaths and recruitments
    int Jstart = J;
    int d; //event number for death/recruits
    
    //number of events within each five year interval and total numbers of events at each 5 year interval
    int numdeaths5 = 2650, fiveyr = 5818, numdeaths10 = 2320, tenyr = 10677, numdeaths15 = 2525, fifteenyr = 15476, numdeaths20 = 2747, twentyyr = 20615, numdeaths25 = 2640, twentyfiveyr = 25846, numdeathsrest = 2860;
   	int numdeaths101 =0, numdeaths10to30=0, numrecruits101=0, numrecruits10to30=0;
    
    //initialize eventscape by 5 year partitions
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
    
    eventscape = calloc(ngen, sizeof(int));
    for(i=0;i<ngen; i++){
        eventscape[i] = 0;
    }
   /* int q=0;
    for(i = 0;i <ngen;i++){
        if(i <tenyr){
            if(eventscape[i] == 0){
                q++;
            }
        }
    }*/ //for checking that the eventscape is set correctly
    
    int a;
    int deaths = 0, recruits =0;//for checking number of death and recruit events that are running

    char outfile1[MAX_FLNAME_SIZE];
    FILE *ofp;
    
    //START SIMULATIONS
    for(s=start;s<end;s++){
        //Reset landscape to starting configuration
        for(x=0; x< (Jstart+numrecruits); x++){
            LC_xcrd[x]=LC_xcrd_start[x];
            LC_ycrd[x]=LC_ycrd_start[x];
            parent_spp_ID[x] = parent_spp_IDstart[x];//added this in; might have been a major mistake about earlier trials?
            //printf("%f,%f\n", LC_xcrd[x],LC_ycrd[x]);
        }
        
        //check that we will have correct number of deaths and recruits
        for(x=0; x<ngen;x++){
            eventscape[x]=eventscapestart[x];
            if(eventscape[x] == 0){deaths++;}else{recruits++;}
        }
    
        curr_year=0;
        J = Jstart;
        
        for(my_gen = 1; my_gen <= ngen; my_gen++){ //run through total number of death+recruits
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
                d = get_rand_integ_intvl(twentyyr+1,twentyfiveyr);
                if(my_gen == twentyfiveyr-1){
                  /*  for(i = twentyyr+1; i <= twentyfiveyr; i++){
                        if(eventscape[i-1] != -1){
                            printf( "%d \n", i);
                        }
                    }
                    printf( "done one \n");
                    for(i = 0; i <= twentyfiveyr+1; i++){
                        if(eventscape[i-1] != -1){
                            printf( "%d \n", i);
                        }
                    }
                    printf( "done two \n");*/
                }
                while(eventscape[d-1] == -1){
                    d = get_rand_integ_intvl(twentyyr+1,twentyfiveyr);
                }

            }else{
                d = get_rand_integ_intvl(twentyfiveyr+1,ngen);
                while(eventscape[d-1] == -1){
                    d = get_rand_integ_intvl(twentyfiveyr+1,ngen);
                }
            }
            //IF DEATH
            if(eventscape[d-1] == 0){
                // For death events, choose an index between 1 & J, inclusive. That member of LC will die.
                toss = get_rand_unit(); // U(0,1).
                toss *= J; // scale by LC size.
                death_choice = (int) ceil(toss); // convert to index.
                //make sure we are picking an alive individual
                while(LC_xcrd[ death_choice ] == -1){
                    toss = get_rand_unit(); // U(0,1).
                    toss *= J; // scale by LC size.
                    death_choice = (int) ceil(toss); // convert to index.
                }
                LC_xcrd[ death_choice ] = -1;
                LC_ycrd[ death_choice ] = -1;
                parent_spp_ID[ death_choice ] = -1;
            }else{
                // Do replacements here.
                b = 1;
				count =0;//count the number of times we have unsuccessful recruitment
                while(b){
                   // printf("%d\n",count); 
                    
                    // ******************* Come back to this point if recruitment is unsuccessful. ******************************
                    min_dist = 99999; //reset min dist every time we are trying to create a new individual
                    
                    // Choose individual to reproduce.
                    toss = get_rand_unit(); // U(0,1).
                    toss *= J; // scale by LC size.
                    repl_choice = (int) ceil(toss); // convert to index.
                    //make sure individual is not dead
                    while(LC_xcrd[ repl_choice ] == -1){
                        toss = get_rand_unit(); // U(0,1).
                        toss *= J; // scale by LC size.
                        repl_choice = (int) ceil(toss); // convert to index.
                    }
                    
                    repl_choicecoords[0]= LC_xcrd[ repl_choice];
                    repl_choicecoords[1]= LC_ycrd[ repl_choice];
                    
                    //pick a dispersal point from this chosen parent
                    select_randloc_frm_cum_disp(repl_choicecoords, size_x, size_y, disc_norm_cum[0], loc_kern_cutoff[0], offspringcoords);
                    
                    //set that to the offspring location
                    offspring_X = offspringcoords[0];
                    offspring_Y = offspringcoords[1];
                    
                    // Assign spp identity. Temporary assignment as recuitment might not be successful.
                    offspring_ID = parent_spp_ID[ repl_choice ];
                    // printf( "Offspring is from spp. %d.\n", offspring_ID );
                    
                    // Determine success of recruit.
                    // For each individual in the plot: Is that individual at least threshold distance away from the potential recruit? If not, stop & go back to reproduction section above.
                    for( k = 0; k < J; k++ ){
                        
                        cur_dist = sqrt( pow( (LC_xcrd[k] - offspring_X), 2 ) + pow( (LC_ycrd[k] - offspring_Y), 2 ) );
                        
                        if( cur_dist < threshd){// if someone is too close, just stop immediately and go back to choosing replacement.
                            break;
                        }
                        if( cur_dist < min_dist ){ // Need to know after this loop whether everyone was far enough away.
                            // printf( " ...updating min. distance... \n" );
                            min_dist = cur_dist;
                        }

                    }
                    
                    if( k >= J & min_dist > threshd ){ // Was able to run through the entire population before breaking loop and everyone was far enough away.
                        //assign offspring to LC
                        LC_xcrd[ J+1 ] = offspring_X;
                        LC_ycrd[ J+1 ] = offspring_Y;
                        parent_spp_ID[ J+1 ] = offspring_ID;
                        J++; //add individual
                        if(my_gen%1000 == 0){
                            printf( "Successful establishment @ generation: %d. min_dist is %f. \n", my_gen, min_dist );
                        }
                        b = 0; // exit while loop & proceed to next generation.
                       // printf("pass");
                        
                    }
                    else{
                        count++;
                    }
                    
            }
            }
            eventscape[d-1] = -1; //we won't use that event again
			//printf("%f\n", my_gen);
            
            //output the compositions at each 5 year interval
            if(my_gen == tenyr){
                curr_year = 10;
               print_SAD_per_simulation_groupgaps(outfile, gap_xcrd, gap_ycrd, gap_size_x, gap_size_y, number_gaps, s, LC_xcrd, LC_ycrd, J, parent_spp_ID,num_species,speciescountpersim, speciescountnongap,curr_year, groupgap, groupgapnum, groupgapsize, nongap_xcrd, nongap_ycrd, Jstart+1);
                
            }
            if(my_gen == ngen){
                curr_year = 30;
               print_SAD_per_simulation_groupgaps(outfile, gap_xcrd, gap_ycrd, gap_size_x, gap_size_y, number_gaps, s, LC_xcrd, LC_ycrd, J, parent_spp_ID,num_species,speciescountpersim, speciescountnongap,curr_year, groupgap, groupgapnum, groupgapsize, nongap_xcrd, nongap_ycrd, Jstart+1);
             
            }
            if(my_gen == 5818){
                curr_year = 5;
                print_SAD_per_simulation_groupgaps(outfile, gap_xcrd, gap_ycrd, gap_size_x, gap_size_y, number_gaps, s, LC_xcrd, LC_ycrd, J, parent_spp_ID,num_species,speciescountpersim, speciescountnongap,curr_year, groupgap, groupgapnum, groupgapsize, nongap_xcrd, nongap_ycrd, Jstart+1);
                
            }
            if(my_gen == 15400){
                curr_year = 15;
                print_SAD_per_simulation_groupgaps(outfile, gap_xcrd, gap_ycrd, gap_size_x, gap_size_y, number_gaps, s, LC_xcrd, LC_ycrd, J, parent_spp_ID,num_species,speciescountpersim, speciescountnongap,curr_year, groupgap, groupgapnum, groupgapsize, nongap_xcrd, nongap_ycrd, Jstart+1);
                
            }
            if(my_gen == 20600){
                curr_year = 20;
                print_SAD_per_simulation_groupgaps(outfile, gap_xcrd, gap_ycrd, gap_size_x, gap_size_y, number_gaps, s, LC_xcrd, LC_ycrd, J, parent_spp_ID,num_species,speciescountpersim, speciescountnongap,curr_year, groupgap, groupgapnum, groupgapsize, nongap_xcrd, nongap_ycrd, Jstart+1);
            }
            if(my_gen == 25800){
                curr_year = 25;
                print_SAD_per_simulation_groupgaps(outfile, gap_xcrd, gap_ycrd, gap_size_x, gap_size_y, number_gaps, s, LC_xcrd, LC_ycrd, J, parent_spp_ID,num_species,speciescountpersim, speciescountnongap,curr_year, groupgap, groupgapnum, groupgapsize, nongap_xcrd, nongap_ycrd, Jstart+1);
            }
        }
        printf( "Done with ngen loop. \n" );
        
    }
    printf("done with simulation\n");
    
}



