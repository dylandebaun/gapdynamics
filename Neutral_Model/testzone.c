//
//  testzone.c
//  gapdynamics
//
//  Created by Dylan DeBaun on 3/10/20.
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
    int start;
    int end;
    //INITIALIZING PARAMETERS+ARRAY ALLOCATION
    int size_x = 1000, size_y=500;
    
    int i,j,k;
    double *LC_xcrd, *LC_ycrd;//array of size 1xJ of the coordinate values for all individuals
    double *LC_xcrd_start, *LC_ycrd_start; //keep track of start landscape per simulation
    int *parent_spp_ID; //array of size 1xJ of the species id values for all individuals
    int *age, *age_start; //matrix with age to choose adult
    int ***species_count; //for calculating avg SAD
    int **speciescountpersim; //for calculating SAD per simulation
    //GAP PARAMETERS
    int *gap_xcrd, *gap_ycrd;
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
    speciescountpersim = calloc(1000, sizeof(int *));
    for (i = 0; i<1000; ++i) {
        speciescountpersim[i] = 0;
    }
    gap_xcrd = calloc(50000, sizeof(int));
    for(i=0;i<50000; i++){
        gap_xcrd[i]=0;
    }
    gap_ycrd = calloc(50000, sizeof(int));
    for(i=0;i<50000; i++){
        gap_ycrd[i]=0;
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
    speciescountpersim = calloc(500, sizeof(int *));
    for(i=0;i<500; i++){
        speciescountpersim[i] = calloc(1000, sizeof(int ));
        for (j = 0; j<1000; ++j) {
            speciescountpersim[i][j] = 0;
        }
    }
 
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
    
    
    srand(time(NULL));
    //allocate arrays
   // LC_xcrd_start = calloc(500000, sizeof(int));
    
    // get_input_info(inputfile, &increase_percent, &size_x, &size_y, kernel_cutoff, max_disp_dist, gauss, sigma, fat_tail, u, p, global, near_neigh, disc_norm_cum, loc_kern_cutoff, &turnover_time, &num_intev_to_meas, &time_btwn_meas, &do_age, &maturity_age, &output_initial, &output_final, &output_interm, &num_interm, &when_output, &outputted_interm, &lowestx, &highestx, &lowesty, &highesty);
    
    //2. INITIALIZE LANDSCAPE WITH EMPIRICAL DATA
    //read in table with columns (species ID, x, y, dbh, dispersal type)
    
    //SIMULATION PARAMETERS
    int J = 0; //# indivs in the landscape
    int num_species;//total number of species present
    double threshd; // R in Sedio and Ostling paper. Didn't have much effect for them when ranging between 0-5m.
    
   // read_in_landscape(startlandscape, &LC_xcrd_start, &LC_ycrd_start, &age_start, &maturity_age, &do_age, &J , &parent_spp_ID,  &num_species);

    //sets threshold to average minimum distance between individuals. only use when you change the start landscape because it takes too long to run
    min_dist_bet_indiv(J,LC_xcrd_start,LC_ycrd_start, &threshd);
    
    //d. identify the gap locations and mark them as small/large
    read_in_gaps(gapfile, &gap_xcrd, &gap_ycrd, &gap_size_x, &gap_size_y, &number_gaps, size_x,size_y);

    //Qualifying gaps as large or small- Not in Use currently
    /*
    num_species  = 328;
    int large = 0, small = 0;
    for(i=0; i < number_gaps; i++){
        if(gap_size_y[i]*gap_size_x[i] >= 100){ //here is where the threshold for small v large is set
            gap_size[i] = 1; //is a large gap
            large++;
        }
        else{
            gap_size[i] = 0; //is a small gap
            small++;
        }
    }
     */
    
    //print informatiion on gaps and on species ID key
    //print_gaps("85", gap_xcrd, gap_ycrd, gap_size_x, gap_size_y, number_gaps);
    int groupgapnum = 0;
    group_gaps(&groupgap, &groupgapsize, &groupgapnum, number_gaps);
   //printgroupgapinfo("group", groupgap, groupgapsize, groupgapnum);
  //  print_species_key(outfile, species_name, num_species);
    num_species = 328;
    
    //ACTUAL
    //get the actual from the gaps for each year
    run_actual("1start10cm85.txt", "10cm85", gap_xcrd, gap_ycrd, gap_size_x, gap_size_y, number_gaps, speciescountpersim, 0, groupgap, groupgapnum, groupgapsize, num_species);
    
    //read in actual landscape
    
    //calculate in gap stuff
    
    printf("done with actual\n");
    
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
    //int hypbardeath = 0, hypbarbirth = 0;
    int ngen = 28490; //0.75*J; //change this each time rf=157861
    //number of simulations to run
    for(s=start;s<end;s++){
        
        //a. Reset landscape to starting configuration
        for(x=0; x< J; x++){
            LC_xcrd[x]=LC_xcrd_start[x];
            LC_ycrd[x]=LC_ycrd_start[x];
           // age[x] = age_start[x];
        }
        curr_year=0;
        
        //one generation = J deaths
        //ngen = number of deaths
        for(my_gen = 1; my_gen <= ngen; my_gen++){
            // For death events, choose an index between 1 & J, inclusive. That member of LC will die.
            toss = get_rand_unit(); // U(0,1).
            toss *= J; // scale by LC size.
            death_choice = (int) ceil(toss); // convert to index.
            
            
            // Do replacements here.
            b = 1;
            while(b){
                // ******************* Come back to this point if recruitment is unsuccessful. ******************************
                //reset min dist every time we are trying to create a new individual
                min_dist = 99999;
                // local replacement
             
                    // Choose individual to reproduce.
                    toss = get_rand_unit(); // U(0,1).
                    toss *= J; // scale by LC size.
                    repl_choice = (int) ceil(toss); // convert to index.
                    
                    
                    //while(age[repl_choice]<10){
                        toss = get_rand_unit(); // U(0,1).
                        toss *= J; // scale by LC size.
                        repl_choice = (int) ceil(toss); // convert to index.

                    repl_choicecoords[0]= LC_xcrd[ repl_choice];
                    repl_choicecoords[1]= LC_ycrd[ repl_choice];

                    select_randloc_frm_cum_disp(repl_choicecoords, size_x, size_y, disc_norm_cum[0], loc_kern_cutoff[0], offspringcoords);
                    
                    offspring_X = offspringcoords[0];
                    offspring_Y = offspringcoords[1];
                    
                    // Eliminate edge effects by making plot spherical.
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
                
                    // Assign spp identity. Temporary assignment as recuitment might not be successful.
                    offspring_ID = parent_spp_ID[ repl_choice ];
                    // printf( "Offspring is from spp. %d.\n", offspring_ID );
                
                // Determine success of recruit.
                // For each individual in the plot:
                // Is that individual at least threshold distance away from the potential recruit? If not, stop & go back to reproduction section above.
                // Calculate the density of enemies around the recruit due to the plot individual.
                // Calculate the impact on the recruit due to the overlap in shared enemies between the recruit and the plot individual.
                // Update the running total of this weighted enemy-mediated factor (big H, Eq. 4 in paper).
                count = 0;
                for(k = 0; k < J; k++ ){

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

                    LC_xcrd[ death_choice ] = offspring_X;
                    LC_ycrd[ death_choice ] = offspring_Y;
                    parent_spp_ID[ death_choice ] = offspring_ID;
                    //age[death_choice] = 0;
                    
                    if(my_gen%1000 == 0){
                    printf( "Successful establishment @ generation: %d. min_dist is %f. \n", my_gen, min_dist );
                    }
                    b = 0; // exit while loop & proceed to next generation.
                
                }
                else{
                    count++;
                }
            
            }
            
            //Print out @ 10,20,30,40 year intervals.
            if(my_gen == 5300){//(ngen*10/35)){ rf= 44637
                curr_year = 10;
                //calc_SAD_in_gaps(outfile, gap_xcrd, gap_ycrd, gap_size_x, gap_size_y, number_gaps, LC_xcrd, LC_ycrd, J, parent_spp_ID,num_species, species_count, curr_year);
                //Reset the species abundances per sim
                //print_SAD_per_simulation(outfile, gap_xcrd, gap_ycrd, gap_size_x, gap_size_y, number_gaps, s, LC_xcrd, LC_ycrd, J, parent_spp_ID,num_species,speciescountpersim, curr_year);
                print_SAD_per_simulation_groupgaps(outfile, gap_xcrd, gap_ycrd, gap_size_x, gap_size_y, number_gaps, s, LC_xcrd, LC_ycrd, J, parent_spp_ID,num_species,speciescountpersim, curr_year, groupgap, groupgapnum, groupgapsize);

            }
            
           /* if(my_gen == 44789){//(ngen*10/35)){
                curr_year = 20;
                //calc_SAD_in_gaps(outfile, gap_xcrd, gap_ycrd, gap_size_x, gap_size_y, number_gaps, LC_xcrd, LC_ycrd, J, parent_spp_ID,num_species, species_count, curr_year);
                //Reset the species abundances per sim
                /*
                for(i=0; i< 1000; i++){
                    speciescountpersim[i] = 0;
                }
                //print_SAD_per_simulation(outfile, gap_xcrd, gap_ycrd, gap_size_x, gap_size_y, number_gaps, s, LC_xcrd, LC_ycrd, J, parent_spp_ID,num_species,speciescountpersim, curr_year);
                print_SAD_per_simulation_groupgaps(outfile, gap_xcrd, gap_ycrd, gap_size_x, gap_size_y, number_gaps, s, LC_xcrd, LC_ycrd, J, parent_spp_ID,num_species,speciescountpersim, curr_year, groupgap, groupgapnum, groupgapsize);
                
            }*/
            
            if(my_gen == ngen){
                curr_year = 30;
                //print_landscape(outfile, LC_xcrd, LC_ycrd, curr_year, parent_spp_ID, J, s, parent_spp_ID);
                //calc_SAD_in_gaps(outfile, gap_xcrd, gap_ycrd, gap_size_x, gap_size_y, number_gaps, LC_xcrd, LC_ycrd, J, parent_spp_ID,num_species, species_count, curr_year);
                //Reset the species abundances per sim

                //print_SAD_per_simulation(outfile, gap_xcrd, gap_ycrd, gap_size_x, gap_size_y, number_gaps, s, LC_xcrd, LC_ycrd, J, parent_spp_ID,num_species,speciescountpersim, curr_year);
               // printf("%d ", species_count[1][1]);
                print_SAD_per_simulation_groupgaps(outfile, gap_xcrd, gap_ycrd, gap_size_x, gap_size_y, number_gaps, s, LC_xcrd, LC_ycrd, J, parent_spp_ID,num_species,speciescountpersim, curr_year, groupgap, groupgapnum, groupgapsize);
                
                
            }
        }
        printf( "Done with ngen loop. \n" );

    }
    //printf("hybar: birth=%d, death=%d",hypbarbirth, hypbardeath);
   /* curr_year = 10;
    print_SAD_in_gaps(outfile, gap_xcrd, gap_ycrd, gap_size_x, gap_size_y, number_gaps, s, LC_xcrd, LC_ycrd, J, parent_spp_ID,num_species,species_count, curr_year);
    curr_year = 20;
    print_SAD_in_gaps(outfile, gap_xcrd, gap_ycrd, gap_size_x, gap_size_y, number_gaps, s, LC_xcrd, LC_ycrd, J, parent_spp_ID,num_species,species_count, curr_year);
    //curr_year = 30;
   // print_SAD_in_gaps(outfile, gap_xcrd, gap_ycrd, gap_size_x, gap_size_y, number_gaps, s, LC_xcrd, LC_ycrd, J, parent_spp_ID,num_species,species_count, curr_year);
    curr_year = 35;
    print_SAD_in_gaps(outfile, gap_xcrd, gap_ycrd, gap_size_x, gap_size_y, number_gaps, s, LC_xcrd, LC_ycrd, J, parent_spp_ID,num_species,species_count, curr_year);
    */
    printf("done with simulation\n");

}


/*y=0;
 for(k = 84; k <= 86; k++){
 gapgroupt[k] = p;
 y++;
 }
 gapgroupsizet[p] = y;
 p++;
 
 
 gapgroupnumt = p;
 for(i = 0; i < realgapnum; i++){
 if(gapgroupt[i] == -1){
 gapgroupt[i] = gapgroupnumt;
 gapgroupsizet[gapgroupnumt] = 1;
 gapgroupnumt++;
 
 }
 }
 
 *gapgroupnum = gapgroupnumt;
 *gapgroup = gapgroupt;
 *gapgroupsize = gapgroupsizet;
 gapgroupt[4] = 0;
 gapgroupt[5] = 0;
 gapgroupt[6] = 0;
 gapgroupsizet[0] = 3;
 
 gapgroupt[119] = 1;
 gapgroupt[120] = 1;
 gapgroupt[121] = 1;
 gapgroupsizet[1] = 3;
 
 gapgroupt[1] = 2;
 gapgroupt[2] = 2;
 gapgroupsizet[2] = 2;
 
 gapgroupt[9] = 3;
 gapgroupt[10] = 3;
 gapgroupsizet[3] = 2;
 
 gapgroupt[96] = 4;
 gapgroupt[97] = 4;
 gapgroupsizet[4] = 2;
 
 gapgroupt[100] = 5;
 gapgroupt[102] = 5;
 gapgroupsizet[5] = 2;
 
 gapgroupt[103] = 6;
 gapgroupt[104] = 6;
 gapgroupt[105] = 6;
 gapgroupt[106] = 6;
 gapgroupt[107] = 6;
 gapgroupt[108] = 6;
 gapgroupt[101] = 6;
 gapgroupt[110] = 6;
 gapgroupt[111] = 6;
 gapgroupt[112] = 6;
 gapgroupt[113] = 6;
 gapgroupt[114] = 6;
 gapgroupt[115] = 6;
 gapgroupt[122] = 6;
 gapgroupt[123] = 6;
 gapgroupsizet[6] = 15;
 
 gapgroupt[124] = 7;
 gapgroupt[125] = 7;
 gapgroupsizet[7] = 2;
 
 gapgroupt[135] = 8;
 gapgroupt[134] = 8;
 gapgroupsizet[8] = 2;
 
 gapgroupt[80] = 9;
 gapgroupt[81] = 9;
 gapgroupt[92] = 9;
 gapgroupt[93] = 9;
 gapgroupt[94] = 9;
 gapgroupsizet[9] = 5;
 
 int k;
 int num =0;
 for(k = 11; k <= 16; k++){
 gapgroupt[k] = 10;
 num++;
 }
 gapgroupsizet[10] = num;
 num = 0;
 
 for(k = 17; k <= 19; k++){
 gapgroupt[k] = 11;
 num++;
 }
 gapgroupsizet[11] = num;
 num = 0;
 
 gapgroupt[78] = 12;
 gapgroupt[79] = 12;
 gapgroupt[82] = 12;
 gapgroupsizet[12] = 3;
 
 int p = 13;
 gapgroupt[40] = p;
 //   gapgroupt[41] = p;
 gapgroupt[67] = p;
 gapgroupt[69] = p;
 gapgroupsizet[p] = 3;
 p++;
 
 int y;
 y=0;
 gapgroupt[70] = p;
 y++;
 gapgroupt[71] = p;
 y++;
 gapgroupsizet[p] = y;
 p++;
 
 y=0;
 for(k = 42; k <= 48; k++){
 gapgroupt[k] = p;
 y++;
 }
 gapgroupt[51] = p;
 y++;
 for(k = 73; k <= 76; k++){
 gapgroupt[k] = p;
 y++;
 }
 gapgroupsizet[p] = y;
 p++;
 
 y=0;
 // gapgroupt[34] = p;
 for(k = 36; k <= 38; k++){
 gapgroupt[k] = p;
 y++;
 }
 for(k = 57; k <= 66; k++){
 gapgroupt[k] = p;
 y++;
 }
 gapgroupsizet[p] = y;
 p++;
 
 y=0;
 for(k = 137; k <= 147; k++){
 gapgroupt[k] = p;
 y++;
 }
 for(k = 129; k <= 130; k++){
 gapgroupt[k] = p;
 y++;
 }
 gapgroupsizet[p] = y;
 p++;
 
 y=0;
 gapgroupt[249] = p;
 y++;
 gapgroupt[250] = p;
 y++;
 gapgroupsizet[p] = y;
 p++;
 
 y=0;
 gapgroupt[266] = p;
 y++;
 gapgroupt[267] = p;
 y++;
 gapgroupsizet[p] = y;
 p++;
 
 y=0;
 for(k = 258; k <= 260; k++){
 gapgroupt[k] = p;
 y++;
 }
 for(k = 244; k <= 245; k++){
 gapgroupt[k] = p;
 y++;
 }
 gapgroupt[242] = p;
 y++;
 gapgroupsizet[p] = y;
 p++;
 
 y=0;
 for(k = 230; k <= 231; k++){
 gapgroupt[k] = p;
 y++;
 }
 for(k = 251; k <= 252; k++){
 gapgroupt[k] = p;
 y++;
 }
 gapgroupsizet[p] = y;
 p++;
 
 y=0;
 for(k = 261; k <= 263; k++){
 gapgroupt[k] = p;
 y++;
 }
 gapgroupsizet[p] = y;
 p++;
 
 y=0;
 for(k = 207; k <= 209; k++){
 gapgroupt[k] = p;
 y++;
 }
 gapgroupsizet[p] = y;
 p++;
 
 y=0;
 for(k = 233; k <= 234; k++){
 gapgroupt[k] = p;
 y++;
 }
 gapgroupsizet[p] = y;
 p++;
 
 y=0;
 for(k = 205; k <= 206; k++){
 gapgroupt[k] = p;
 y++;
 }
 for(k = 216; k <= 217; k++){
 gapgroupt[k] = p;
 y++;
 }
 gapgroupsizet[p] = y;
 p++;
 
 y=0;
 for(k = 264; k <= 265; k++){
 gapgroupt[k] = p;
 y++;
 }
 gapgroupsizet[p] = y;
 p++;
 
 y=0;
 for(k = 149; k <= 158; k++){
 gapgroupt[k] = p;
 y++;
 }
 for(k = 165; k <= 180; k++){
 gapgroupt[k] = p;
 y++;
 }
 gapgroupsizet[p] = y;
 p++;
 
 y=0;
 for(k = 254; k <= 255; k++){
 gapgroupt[k] = p;
 y++;
 }
 gapgroupsizet[p] = y;
 p++;
 
 y=0;
 for(k = 271; k <= 274; k++){
 gapgroupt[k] = p;
 y++;
 }
 gapgroupsizet[p] = y;
 p++;
 
 y=0;
 for(k = 211; k <= 212; k++){
 gapgroupt[k] = p;
 y++;
 }
 gapgroupsizet[p] = y;
 p++;
 
 y=0;
 gapgroupt[275] = p;
 y++;
 gapgroupt[281] = p;
 y++;
 gapgroupsizet[p] = y;
 p++;
 
 y=0;
 gapgroupt[198] = p;
 y++;
 gapgroupt[199] = p;
 y++;
 gapgroupt[201] = p;
 y++;
 gapgroupsizet[p] = y;
 p++;
 
 y=0;
 for(k = 238; k <= 240; k++){
 gapgroupt[k] = p;
 y++;
 }
 gapgroupt[236] = p;
 y++;
 gapgroupsizet[p] = y;
 p++;
 
 y=0;
 for(k = 282; k <= 287; k++){
 gapgroupt[k] = p;
 y++;
 }
 for(k = 276; k <= 277; k++){
 gapgroupt[k] = p;
 y++;
 }
 gapgroupsizet[p] = y;
 p++;
 
 y=0;
 for(k = 222; k <= 225; k++){
 gapgroupt[k] = p;
 y++;
 }
 gapgroupsizet[p] = y;
 p++;
 
 y=0;
 gapgroupt[401] = p;
 y++;
 gapgroupt[421] = p;
 y++;
 gapgroupt[422] = p;
 y++;
 gapgroupsizet[p] = y;
 p++;
 
 y=0;
 for(k = 159; k <= 162; k++){
 gapgroupt[k] = p;
 y++;
 }
 gapgroupt[185] = p;
 y++;
 gapgroupsizet[p] = y;
 p++;
 
 y=0;
 for(k = 296; k <= 300; k++){
 gapgroupt[k] = p;
 y++;
 }
 gapgroupsizet[p] = y;
 p++;
 
 y=0;
 gapgroupt[330] = p;
 y++;
 gapgroupt[333] = p;
 y++;
 gapgroupsizet[p] = y;
 p++;
 
 y=0;
 gapgroupt[331] = p;
 y++;
 gapgroupt[332] = p;
 y++;
 gapgroupsizet[p] = y;
 p++;
 
 y=0;
 for(k = 301; k <= 303; k++){
 gapgroupt[k] = p;
 y++;
 }
 for(k = 311; k <= 322; k++){
 gapgroupt[k] = p;
 y++;
 }
 gapgroupt[327] = p;
 y++;
 gapgroupt[328] = p;
 y++;
 gapgroupsizet[p] = y;
 p++;
 
 y=0;
 for(k = 402; k <= 403; k++){
 gapgroupt[k] = p;
 y++;
 }
 gapgroupsizet[p] = y;
 p++;
 
 y=0;
 for(k = 323; k <= 324; k++){
 gapgroupt[k] = p;
 y++;
 }
 gapgroupsizet[p] = y;
 p++;
 
 y=0;
 for(k = 355; k <= 356; k++){
 gapgroupt[k] = p;
 y++;
 }
 gapgroupt[366] = p;
 y++;
 gapgroupsizet[p] = y;
 p++;
 
 y=0;
 for(k = 375; k <= 376; k++){
 gapgroupt[k] = p;
 y++;
 }
 gapgroupsizet[p] = y;
 p++;
 
 y=0;
 for(k = 377; k <= 379; k++){
 gapgroupt[k] = p;
 y++;
 }
 gapgroupsizet[p] = y;
 p++;
 
 y=0;
 for(k = 405; k <= 414; k++){
 gapgroupt[k] = p;
 y++;
 }
 gapgroupsizet[p] = y;
 p++;
 
 y=0;
 for(k = 340; k <= 341; k++){
 gapgroupt[k] = p;
 y++;
 }
 gapgroupsizet[p] = y;
 p++;
 
 y=0;
 for(k = 305; k <= 306; k++){
 gapgroupt[k] = p;
 y++;
 }
 gapgroupsizet[p] = y;
 p++;
 
 y=0;
 for(k = 350; k <= 351; k++){
 gapgroupt[k] = p;
 y++;
 }
 gapgroupsizet[p] = y;
 p++;
 
 y=0;
 for(k = 359; k <= 364; k++){
 gapgroupt[k] = p;
 y++;
 }
 gapgroupsizet[p] = y;
 p++;
 y=0;
 for(k = 307; k <= 308; k++){
 gapgroupt[k] = p;
 y++;
 }
 gapgroupsizet[p] = y;
 p++;
 
 y=0;
 for(k = 343; k <= 344; k++){
 gapgroupt[k] = p;
 y++;
 }
 gapgroupt[348] = p;
 y++;
 gapgroupsizet[p] = y;
 p++;
 
 gapgroupt[339] = p;
 y++;
 gapgroupt[346] = p;
 y++;
 gapgroupsizet[p] = y;
 p++;
 
 gapgroupt[417] = p;
 y++;
 gapgroupt[420] = p;
 y++;
 gapgroupsizet[p] = y;
 p++;
 
 y=0;
 for(k = 390; k <= 392; k++){
 gapgroupt[k] = p;
 y++;
 }
 for(k = 387; k <= 389; k++){
 gapgroupt[k] = p;
 y++;
 }
 gapgroupt[384] = p;
 y++;
 gapgroupt[385] = p;
 y++;
 gapgroupt[386] = p;
 y++;
 gapgroupt[371] = p;
 y++;
 for(k = 393; k <= 394; k++){
 gapgroupt[k] = p;
 y++;
 }
 gapgroupsizet[p] = y;
 p++;
 
 y=0;
 for(k = 535; k <= 536; k++){
 gapgroupt[k] = p;
 y++;
 }
 gapgroupsizet[p] = y;
 p++;
 
 y=0;
 for(k = 527; k <= 528; k++){
 gapgroupt[k] = p;
 y++;
 }
 gapgroupsizet[p] = y;
 p++;
 
 y=0;
 gapgroupt[530] = p;
 y++;
 gapgroupt[538] = p;
 y++;
 gapgroupsizet[p] = y;
 p++;
 
 y=0;
 for(k = 500; k <= 501; k++){
 gapgroupt[k] = p;
 y++;
 }
 gapgroupsizet[p] = y;
 p++;
 
 y=0;
 for(k = 458; k <= 459; k++){
 gapgroupt[k] = p;
 y++;
 }
 gapgroupt[460] = p;
 y++;
 gapgroupt[481] = p;
 y++;
 gapgroupsizet[p] = y;
 p++;
 
 y=0;
 for(k = 485; k <= 487; k++){
 gapgroupt[k] = p;
 y++;
 }
 gapgroupsizet[p] = y;
 p++;
 
 y=0;
 for(k = 490; k <= 491; k++){
 gapgroupt[k] = p;
 y++;
 }
 gapgroupt[488] = p;
 y++;
 gapgroupsizet[p] = y;
 p++;
 
 y=0;
 for(k = 461; k <= 462; k++){
 gapgroupt[k] = p;
 y++;
 }
 gapgroupt[489] = p;
 y++;
 gapgroupsizet[p] = y;
 p++;
 
 y=0;
 for(k = 463; k <= 465; k++){
 gapgroupt[k] = p;
 y++;
 }
 gapgroupsizet[p] = y;
 p++;
 
 y=0;
 for(k = 516; k <= 517; k++){
 gapgroupt[k] = p;
 y++;
 }
 gapgroupt[522] = p;
 y++;
 gapgroupsizet[p] = y;
 p++;
 
 y=0;
 for(k = 504; k <= 508; k++){
 gapgroupt[k] = p;
 y++;
 }
 for(k = 519; k <= 521; k++){
 gapgroupt[k] = p;
 y++;
 }
 for(k = 523; k <= 525; k++){
 gapgroupt[k] = p;
 y++;
 }
 gapgroupsizet[p] = y;
 p++;
 
 y=0;
 for(k = 541; k <= 542; k++){
 gapgroupt[k] = p;
 y++;
 }
 gapgroupsizet[p] = y;
 p++;
 
 y=0;
 for(k = 466; k <= 474; k++){
 gapgroupt[k] = p;
 y++;
 }
 gapgroupsizet[p] = y;
 p++;
 
 y=0;
 for(k = 449; k <= 450; k++){
 gapgroupt[k] = p;
 y++;
 }
 gapgroupsizet[p] = y;
 p++;
 
 y=0;
 for(k = 451; k <= 453; k++){
 gapgroupt[k] = p;
 y++;
 }
 gapgroupsizet[p] = y;
 p++;
 
 y=0;
 for(k = 493; k <= 495; k++){
 gapgroupt[k] = p;
 y++;
 }
 gapgroupt[475] = p;
 y++;
 gapgroupsizet[p] = y;
 p++;
 
 y=0;
 for(k = 439; k <= 440; k++){
 gapgroupt[k] = p;
 y++;
 }
 gapgroupsizet[p] = y;
 p++;
 
 y=0;
 for(k = 497; k <= 498; k++){
 gapgroupt[k] = p;
 y++;
 }
 gapgroupsizet[p] = y;
 p++;
 
 y=0;
 for(k = 455; k <= 457; k++){
 gapgroupt[k] = p;
 y++;
 }
 gapgroupsizet[p] = y;
 p++;
 
 y=0;
 for(k = 478; k <= 480; k++){
 gapgroupt[k] = p;
 y++;
 }
 gapgroupsizet[p] = y;
 p++;
 
 y=0;
 for(k = 429; k <= 430; k++){
 gapgroupt[k] = p;
 y++;
 }
 gapgroupt[443] = p;
 y++;
 gapgroupt[444] = p;
 y++;
 gapgroupsizet[p] = y;
 p++;
 
 y=0;
 for(k = 591; k <= 594; k++){
 gapgroupt[k] = p;
 y++;
 }
 gapgroupsizet[p] = y;
 p++;
 
 y=0;
 for(k = 557; k <= 558; k++){
 gapgroupt[k] = p;
 y++;
 }
 gapgroupt[553] = p;
 y++;
 gapgroupsizet[p] = y;
 p++;
 
 y=0;
 for(k = 559; k <= 560; k++){
 gapgroupt[k] = p;
 y++;
 }
 gapgroupt[565] = p;
 y++;
 gapgroupsizet[p] = y;
 p++;
 
 y=0;
 gapgroupt[573] = p;
 y++;
 gapgroupt[576] = p;
 y++;
 gapgroupsizet[p] = y;
 p++;
 
 y=0;
 for(k = 574; k <= 575; k++){
 gapgroupt[k] = p;
 y++;
 }
 gapgroupt[577] = p;
 y++;
 gapgroupsizet[p] = y;
 p++;
 
 y=0;
 gapgroupt[587] = p;
 y++;
 gapgroupt[595] = p;
 y++;
 gapgroupsizet[p] = y;
 p++;
 
 y=0;
 gapgroupt[588] = p;
 y++;
 gapgroupt[597] = p;
 y++;
 gapgroupsizet[p] = y;
 p++;
 
 y=0;
 gapgroupt[598] = p;
 y++;
 gapgroupt[599] = p;
 y++;
 gapgroupsizet[p] = y;
 p++;
 
 y=0;
 gapgroupt[571] = p;
 y++;
 gapgroupt[579] = p;
 y++;
 gapgroupsizet[p] = y;
 p++;
 
 y=0;
 gapgroupt[608] = p;
 y++;
 gapgroupt[609] = p;
 y++;
 gapgroupsizet[p] = y;
 p++;
 
 y=0;
 for(k = 605; k <= 607; k++){
 gapgroupt[k] = p;
 y++;
 }
 gapgroupsizet[p] = y;
 p++;
 
 /*y=0;
 for(k = 84; k <= 86; k++){
 gapgroupt[k] = p;
 y++;
 }
 gapgroupsizet[p] = y;
 p++;
 
 
 gapgroupnumt = p;
 for(i = 0; i < realgapnum; i++){
 if(gapgroupt[i] == -1){
 gapgroupt[i] = gapgroupnumt;
 gapgroupsizet[gapgroupnumt] = 1;
 gapgroupnumt++;
 
 }
 }
 
 *gapgroupnum = gapgroupnumt;
 *gapgroup = gapgroupt;
 *gapgroupsize = gapgroupsizet;
 */


//from read in gaps
/*
 int k,j;
 int n = i;
 //now we need to concatenate if some gaps are bigger than 5x5
 for(k=0;k<n;k++){
 if(gap_xcrd[k] + gap_size_x[k]== (gap_xcrd[k+1])){
 gap_size_x[k] += 5;
 //must then sort the matrix to get rid of gap_xcrd[k+1] since it combined with gap_xcrd[k]
 for(j=(k+1);j<n-1;j++){
 gap_xcrd[j] = gap_xcrd[j+1];
 gap_ycrd[j] = gap_ycrd[j+1];
 gap_size_x[j] = gap_size_x[j+1];
 gap_size_y[j] = gap_size_y[j+1];
 }
 n = n-1;
 k = k-1; //this is to recheck that spot now that we expanded it
 }
 }
 for(k=0;k<n;k++){
 //make sure we are creating rectangular gaps instead of L shaped gaps
 if((gap_ycrd[k] + gap_size_y[k] ==gap_ycrd[k+1]) && (gap_size_x[k] = gap_size_x[k+1])){
 gap_size_y[k] += 5;
 //must then sort the matrix to get rid of gap_xcrd[k+1] since it combined with gap_xcrd[k]
 for(j=(k+1);j<i-1;j++){
 gap_xcrd[j] = gap_xcrd[j+1];
 gap_ycrd[j] = gap_ycrd[j+1];
 gap_size_x[j] = gap_size_x[j+1];
 gap_size_y[j] = gap_size_y[j+1];
 }
 n = n-1;
 k = k-1; //this is to recheck that spot now that we expanded it
 }
 }
 */
