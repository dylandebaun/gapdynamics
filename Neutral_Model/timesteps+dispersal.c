//
//  timesteps+dispersal.c
//  gapdynamics
//
//  Created by Dylan DeBaun on 1/21/20.
//  Copyright © 2020 Dylan DeBaun. All rights reserved.
//

#include "timesteps+dispersal.h"
#include "general.h"

# include <stdio.h>
# include <math.h>
# include <stdlib.h>
# include <string.h>
# include "general.h"

//HELP make sure all lanscapes changed to int from char

//dispersal kernel info

/* discrete_normalized_cummulative will return the discrete, normalized
 cummulative probability distribution for the dispersal kernel f(r,
 parameters) out to the point where f(r, parameters) is smaller than
 the cutoff.  It will also return the location of this cutoff point
 in loc_cutoff. */
double *discrete_normalized_cummulative(double f(double, double*),
                                        double *parameters,
                                        double cutoff, int max_disp_dist,
                                        int *point_loc_cutoff)
{
    
    int curr_arraysize=100, pastmax=0;
    double *disc_norm_cum;
    
    int i, j;
    
    disc_norm_cum = calloc(curr_arraysize, sizeof(double));
    printf("\nCalculating dispersal kernel...\n");
   // fflush(NULL);
    
    i=1;
    /* while we are either not past the maximum or the value to add is
     not less than the cutoff, keep filling array*/
    while (!pastmax || (1-disc_norm_cum[i-1] >= cutoff && i<= max_disp_dist)) {
        /* if haven't already passed the maximum, check if you just did. */
        if (!pastmax)
            if (f(i,parameters) <= f(i-1,parameters)) pastmax=1;
        /* fill if at point less than array size, otherwise reallocate */
        if (i==curr_arraysize) {
            curr_arraysize = curr_arraysize + 100;
            disc_norm_cum = realloc(disc_norm_cum,
                                    curr_arraysize*sizeof(double));
        }
        disc_norm_cum[i] = disc_norm_cum[i-1];
        if (i==0)                                                                   // wouldn't i never be 0 b/c it is initialized at 1 and only has incrementing parts
            for(j=0; j<=4; ++j)
                disc_norm_cum[i] = disc_norm_cum[i] + 0.1*f(i+j*0.1, parameters);//0.1 is discretization adding a dr
        else
            for(j=0; j<=9; ++j)
                disc_norm_cum[i] = disc_norm_cum[i] + 0.1*f(i-0.5+j*0.1, parameters);
        /*    printf("i=%d, disc_norm_cum=%f\n", i, disc_norm_cum[i]); */
        ++i;
    }
    
    /* return the location where the cutoff was met */
    *point_loc_cutoff = i-1;
    
    /* normalize the cummulative probabilities, making the largest
     probability be 1 */
    printf("Normalizing w/ loc_cutoff=%d, remaining prob=%f...\n", *point_loc_cutoff, 1-disc_norm_cum[*point_loc_cutoff]);
    for (i=0; i<=*point_loc_cutoff; ++i) {
        disc_norm_cum[i] = disc_norm_cum[i]/disc_norm_cum[*point_loc_cutoff];
        /*    printf("i=%d, disc_norm_cum=%f\n", i, disc_norm_cum[i]);*/
    }
    
    /* return the pointer to the discrete normalized cummulative
     probability array */
    return disc_norm_cum;
    
}


void min_dist_bet_indiv(int J, double *LC_xcrd, double *LC_ycrd, double *threshold){
    int k,j;
    double min_dist = 99999; // setting arbitrarily high.
    double cur_dist; // current distance of closest neighbor.
    double min_dist_total = 0;
    
    for( k = 0; k < J; k++ ){
        min_dist=99999;
        for(j=0;j<J;j++){
        cur_dist = sqrt( pow( (LC_xcrd[k] - LC_xcrd[j]), 2 ) + pow( (LC_ycrd[k] - LC_ycrd[j] ), 2 ) );
       
        if( cur_dist < min_dist && cur_dist != 0){ // Need to know after this loop whether everyone was far enough away.
            // printf( " ...updating min. distance... \n" );
            min_dist = cur_dist;
            
        }
        
    }
        min_dist_total += min_dist;
    }
    *threshold = min_dist_total/J;
    printf("average minimum distance = %f",*threshold);
}


/*readTopology(((((M1,M2),(Ig,((Ip,Ip),(Im,Io)))),(La,(Lm,Lm))),Maf););
readTopology((((((M1,M2)#H1),(Ig,((Ip,Ip),(Im,Io)))),((La,(Lm,Lm))),#H1),Maf););
*/

//ACTUAL SIM STUFF
void run_actual(char *landscapefile, char *begoutfile, int *gap_xcrd, int *gap_ycrd, int *gap_sizex, int *gap_sizey, int gap_number, int **speciescount, int curr_year, int *gapgroup, int groupgapnum, int *groupgapsize, int num_species){

//}, double **LC_xcrd, double **LC_ycrd, int **parent_spp_ID, int *numind){
    //READ IN
    FILE *ifp;
    int dbh, sp;
    float x1,y1;
    int i;
    int status;
    double *LC_xcrdtemp, *LC_ycrdtemp;
    int *parent_spp_IDtemp;
    
    LC_xcrdtemp = calloc(236000, sizeof(double));
    for(i=0;i<236000; i++){
        LC_xcrdtemp[i]=0;
    }
    LC_ycrdtemp = calloc(236000, sizeof(double));
    for(i=0;i<236000; i++){
        LC_ycrdtemp[i]=0;
    }
    
    parent_spp_IDtemp = calloc(236000, sizeof(int));
    for(i=0;i<236000; i++){
        parent_spp_IDtemp[i]=0;
    }

    int n =0;
    int numspeciestemp = 0;
    printf("Now reading in actual landscape from %s\n", landscapefile);

    ifp = fopen(landscapefile, "r");
    status=fscanf(ifp, "%d %f %f", &sp, &x1, &y1);
    numspeciestemp++;
    //fflush(ifp);
    
    while (status != EOF) {
        //sets individual's coordinates and age
        LC_xcrdtemp[n] = x1;
        LC_ycrdtemp[n] = y1;
        parent_spp_IDtemp[n] = sp;
        if(sp >= numspeciestemp){
            numspeciestemp = sp;
        }
        //Read in the next individual:
        n++;
        status=fscanf(ifp, "%d %f %f", &sp, &x1, &y1);
        //fflush(ifp);
        
    }
    fclose(ifp);
    //fflush(NULL);
    
    //CALCULATE INDIVS IN GAPS
        int z,j,k,t,q;
        char outfile[MAX_FLNAME_SIZE];
        FILE *ofp;
    int ingap,g1,g2,g3,g4,p,r,op=0;
    double lc1,lc2,stop;
    int sumindiv = 0;
    int **speciescount1,gxsize,gysize;
    double x,y,gx,gy;
    speciescount1 = calloc(500, sizeof(int *));
    for(i=0;i<500; i++){
        speciescount1[i] = calloc(1000, sizeof(int));
        for (j = 0; j<1000; ++j) {
            speciescount1[i][j] = 0;
        }
    }
    int count = 0;
    for(j=0; j<n; j++){
        x = LC_xcrdtemp[j];
        y = LC_ycrdtemp[j];
        if(x!= -1){
            count++;
            //printf("%f\t%f\n",LC_xcrdtemp[j],LC_ycrdtemp[j]);
        }
    }
    printf("%d\n",count);
    int b,thresholdlower,thresholdupper,countt;
    char *type;
        //for a grouped gap
    
    for(i=0;i<631;i++){
        printf("%d\t%d\n",gap_xcrd[i],gap_ycrd[i]);
    }
    for(b =0; b<3;b++){
        if(b==0){
            type = "all";
            thresholdlower = 0;
            thresholdupper = 100;
        }else if(b==1){
            type = "large";
            thresholdlower = 2;
            thresholdupper = 100;
        }else if(b==2){
            type = "small";
            thresholdlower = 0;
            thresholdupper = 2;
        }
        for(k=0; k< 1000; k++){
            speciescount1[0][k] = 0;
        }
        sumindiv = 0;
      //  countt=0;
        for(t = 0; t < groupgapnum; t++){
           // if(groupgapsize[t] > thresholdlower && groupgapsize[t] <= thresholdupper){
            //    op++;
            stop = 0;
            if(b==2){
                if(groupgapsize[t]<3){
                    stop = 1;
                }else{
                    op++;
                }
            }else if(b==1){
                if(groupgapsize[t]>2){
                    stop = 1;
                }else{
                    op++;
                }
            }
            if(stop == 0){
             //   countt++;
            //find out how many real gaps are in that grouped gap
           // for(q = 0; q< groupgapsize[t]; q++){
                //for the real gap in the grouped gap, count the species
                for(i = 0; i < gap_number; i++){
                   // ingap = gapgroup[i];
                    if(gapgroup[i] == t){
                        for(j=0; j<n; j++){
                            x = LC_xcrdtemp[j];
                            y = LC_ycrdtemp[j];
                            gx = gap_xcrd[i];
                            gy = gap_ycrd[i];
                            gxsize =gap_sizex[i];
                            gysize =gap_sizey[i];
                            if(gxsize != 5 || gysize != 5){
                                printf("STOP! gxsize = %d,gysize=%d",gxsize,gysize);
                            }
                            if((x > gx) && (y > gy) && (y <= (gy+gysize)) && (x <= (gx+gxsize))){
                                // printf("%d ", speciescount[i][spID[j]]);
                               // p = parent_spp_IDtemp[j];
                               // printf("%d\n",p);
                                speciescount1[0][parent_spp_IDtemp[j]]++;
                               // r=speciescount1[0][parent_spp_IDtemp[j]];
                                sumindiv++;
                                //if(speciescount[curr_year][i][spID[j]] > 0){
                                
                                // }
                                // printf("%d,%d\n", parent_spp_IDtemp[j],r);
                            }
                        }
                 //   }
                }
            }
           
            //printf("done\n");
           }
           // printf("%d\n", speciescount1[t][22]);
            //upon getting all q gaps in the group gap tallied, write the group gap to a file
        }
        printf("count here should be 554. is it: %d\n",sumindiv);
       // printf("%s: %d", type,countt);
        sprintf(outfile, "%sexpectedin%sat%dyrs.txt", begoutfile,type,curr_year);
        // printf("%d\n", speciescount1[t][22]);
        ofp = fopen(outfile, "w");
        // printf("%d\n", speciescount1[t][22]);
        for(z=0;z<num_species; z++){
            // printf("%d\n", speciescount1[t][22]);
            fprintf(ofp, "%d\t%d\n", z, speciescount1[0][z]);
            //printf("%d\n", speciescount1[t][z]);
        }
        fclose(ofp);
        printf("%s: %d", type,sumindiv);
    }
    
   /* *LC_xcrd= LC_xcrdtemp;
    *LC_ycrd= LC_ycrdtemp;
    *parent_spp_ID = parent_spp_IDtemp;
    *numind = n;*/
}


void run_actualnongap(char *landscapefile, char *begoutfile, int **nongap_xcrd, int **nongap_ycrd, double **speciescount, int curr_year, int *gapgroup, int groupgapnum, int *groupgapsize, int num_species){
   
    //}, double **LC_xcrd, double **LC_ycrd, int **parent_spp_ID, int *numind){
    //READ IN
    FILE *ifp;
    int dbh, sp;
    float x1,y1;
    int i;
    int status;
    double *LC_xcrdtemp, *LC_ycrdtemp;
    int *parent_spp_IDtemp;
    
    LC_xcrdtemp = calloc(236000, sizeof(double));
    for(i=0;i<236000; i++){
        LC_xcrdtemp[i]=0;
    }
    LC_ycrdtemp = calloc(236000, sizeof(double));
    for(i=0;i<236000; i++){
        LC_ycrdtemp[i]=0;
    }
    
    parent_spp_IDtemp = calloc(236000, sizeof(int));
    for(i=0;i<236000; i++){
        parent_spp_IDtemp[i]=0;
    }
    
    int n =0;
    int numspeciestemp = 0;
    printf("Now reading in actual landscape from %s\n", landscapefile);
    
    ifp = fopen(landscapefile, "r");
    status=fscanf(ifp, "%d %f %f", &sp, &x1, &y1);
    numspeciestemp++;
    //fflush(ifp);
    
    while (status != EOF) {
        //sets individual's coordinates and age
        LC_xcrdtemp[n] = x1;
        LC_ycrdtemp[n] = y1;
        parent_spp_IDtemp[n] = sp;
        if(sp >= numspeciestemp){
            numspeciestemp = sp;
        }
        //Read in the next individual:
        n++;
        status=fscanf(ifp, "%d %f %f", &sp, &x1, &y1);
        //fflush(ifp);
        
    }
    fclose(ifp);
    //fflush(NULL);
    
    //CALCULATE INDIVS IN GAPS
    int z,j,k,t,q;
    char outfile[MAX_FLNAME_SIZE];
    FILE *ofp;
    int ingap,g1,g2,g3,g4,p,r;
    double lc1,lc2;
    int sumindiv = 0;
    int x,y,gx,gy,gxsize,gysize;
    double **speciescount1;
    speciescount1 = calloc(500, sizeof(double *));
    for(i=0;i<500; i++){
        speciescount1[i] = calloc(1000, sizeof(double));
        for (j = 0; j<1000; ++j) {
            speciescount1[i][j] = 0;
        }
    }
    int count = 0;
    for(j=0; j<n; j++){
        x = LC_xcrdtemp[j];
        y = LC_ycrdtemp[j];
        if(x!= -1){
            count++;
            //printf("%f\t%f\n",LC_xcrdtemp[j],LC_ycrdtemp[j]);
        }
    }
 //   printf("%d\n",count);
    int b,thresholdlower,thresholdupper,countt;
    char *type;
    //for a grouped gap
    for(q=0;q < 10;q++){
    for(b =0; b<3;b++){
        if(b==0){
            type = "all";
            thresholdlower = 0;
            thresholdupper = 100;
        }else if(b==1){
            type = "large";
            thresholdlower = 2;
            thresholdupper = 100;
        }else if(b==2){
            type = "small";
            thresholdlower = 0;
            thresholdupper = 2;
        }
        for(k=0; k< 1000; k++){
            speciescount1[0][k] = 0;
        }
        sumindiv = 0;
        countt=0;
        for(t = 0; t < groupgapnum; t++){
            if(groupgapsize[t] > thresholdlower && groupgapsize[t] < thresholdupper){
                countt++;
                //find out how many real gaps are in that grouped gap
                // for(q = 0; q< groupgapsize[t]; q++){
                //for the real gap in the grouped gap, count the species
                //for(i = 0; i < gap_number; i++){
                   // ingap = gapgroup[i];
                    //if(gapgroup[i] == t){
                for(i =0;i < 1;i++){
                        for(j=0; j<n; j++){
                            x = LC_xcrdtemp[j];
                            y = LC_ycrdtemp[j];
                            gx = nongap_xcrd[q][t];
                            gy = nongap_ycrd[q][t];
                            gxsize =sqrt(groupgapsize[t]*25);
                            gysize =sqrt(groupgapsize[t]*25);
                            if(x > gx && y > gy && y <= gy+gysize && x <= gx+gxsize){
                                // printf("%d ", speciescount[i][spID[j]]);
                                p = parent_spp_IDtemp[j];
                               // printf("%d\n",p);
                                speciescount1[0][parent_spp_IDtemp[j]]++;
                                r=speciescount1[0][parent_spp_IDtemp[j]];
                                sumindiv++;
                            }
                        }
                }
             //   printf("done\n");
            }
            // printf("%d\n", speciescount1[t][22]);
            //upon getting all q gaps in the group gap tallied, write the group gap to a file
        }
        /*
        for(i =0;i<num_species;i++){
            speciescount1[0][i] /= 10;
        }*/
        
        //printf("%s: %d", type,countt);
        sprintf(outfile, "%s%dexpectedin%sat%dyrs.txt", begoutfile,q,type,curr_year);
        // printf("%d\n", speciescount1[t][22]);
        ofp = fopen(outfile, "w");
        // printf("%d\n", speciescount1[t][22]);
        for(z=0;z<num_species; z++){
            // printf("%d\n", speciescount1[t][22]);
            fprintf(ofp, "%d\t%f\n", z, speciescount1[0][z]);
            //printf("%d\n", speciescount1[t][z]);
        }
        fclose(ofp);
       // printf("%s: %d", type,sumindiv);
    }
    }
    
    /* *LC_xcrd= LC_xcrdtemp;
     *LC_ycrd= LC_ycrdtemp;
     *parent_spp_ID = parent_spp_IDtemp;
     *numind = n;*/
}










/*

void select_randloc_frm_cum_disp(int *ind1_coords, int sizex, int sizey, double *disc_norm_cum, int loc_cutoff, int *ind2_coords, int adjustmentx, int adjustmenty, int *i)
{
    
    double toss;
    int y = 0;
    toss = get_rand_unit();
    while((toss > disc_norm_cum[y]) && (y<loc_cutoff)) {
        y += 1;
    }
    choose_randloc_atdist(y, ind2_coords);
    ind2_coords[0] = ind2_coords[0] + ind1_coords[0];
    ind2_coords[1] = ind2_coords[1] + ind1_coords[1];
    
    if (ind2_coords[0] >= sizex || ind2_coords[0] < 0 || ind2_coords[1] >= sizey || ind2_coords[1] < 0) {
        ind2_coords[0] = 0;
     }
    
}*/

//run time step
/**
 * Modifies: time. elapses the time by a factor of probabilities
 * Effects: kills/generates individual
 *
double time_evolve(int ***LANDSCAPE, int size_x, int size_y, int *p_num_indiv, int *p_num_slots_inuse, int *p_curr_num_slots, int inc_num_slots, double *disc_norm_cum, int loc_kern_cutoff, int turnover_time, int do_age, short **p_year_born, int maturity_age, double curr_year) //HELP check reqs for all of these
{
    double toss;
    double dmsk2, dmsk3;
    double p_same;
    
    short *yb_TEMP;
    
    int *event_ind;
    
    //GENERATE OFFSPRING/KILL INDIV
    
     50-50 death/birth event*
    toss = get_rand_unit();
    //kills an individual
    if (*p_num_indiv == size_x*size_y){
        printf("Forest full.\n"); //HELP don't want this to happen.
    }
    if ( ((toss < 0.5) && (*p_num_indiv >1)) || (*p_num_indiv == size_x*size_y)) {
        *event_ind = choose_event_indiv(LANDSCAPE, size_x, size_y);
        Kill(&event_ind[2], LANDSCAPE, p_num_indiv);
    }
    //generate offspring by picking parent and picking location w. disp kernel from parent
    else {
        if (do_age) {
            *event_ind=choose_parent_indiv_wagestruct(LANDSCAPE, *p_num_slots_inuse, curr_year, maturity_age, size_x, size_y);
            yb_TEMP=*p_year_born; //HELP make sure this works
        }
        else {
            *event_ind = choose_event_indiv(LANDSCAPE, size_x, size_y); //HELP make sure this has the right number of pointers
        }
        
        if (*p_num_slots_inuse == *p_curr_num_slots-1) {

            printf("HELP! Should not have had to sort! \n");
        }
        
        gen_offspring(event_ind, LANDSCAPE, size_x, size_y, p_num_indiv, p_num_slots_inuse, disc_norm_cum, loc_kern_cutoff, do_age, *p_year_born, curr_year);
        
    }
    
    //MODIFYING TIME
     approximate the likelihood nothing happens to the species using the global dispersal approximation to that probability*
    dmsk2=(double) *p_num_indiv;
    dmsk3=(double) size_x*size_y;
    p_same = (1-dmsk2/dmsk3)*(1-dmsk2/dmsk3)+(dmsk2/dmsk3)*(dmsk2/dmsk3);
    
     must update the time according to the likelihood that the species distribution would have remained the same *
    return (p_same/(1-p_same))*(((double) turnover_time)/((double) size_x*size_y));
    
}

*
 * Effects: Chooses random individual for death/birth
 *
int choose_event_indiv(int ***LANDSCAPE, int sizex, int sizey)
{
    
    int *event_ind;
    
    event_ind[0]=get_rand_integ_intvl(1,sizex);
    event_ind[1]=get_rand_integ_intvl(1,sizey);
    
    while (LANDSCAPE[event_ind[0]][event_ind[1]][0] == 0){
        event_ind[0]=get_rand_integ_intvl(1,sizex);
        event_ind[1]=get_rand_integ_intvl(1,sizey);
    }
    return *event_ind;
    
}

void Kill(int *event_ind, int ***LANDSCAPE, int *p_num_indiv)
{
    
    *p_num_indiv -= 1;//HELP may be able to delete
    //removes individual from landscape
    LANDSCAPE[event_ind[0]][event_ind[1]][0] = 0;
    //resets age in case that was set
    LANDSCAPE[event_ind[0]][event_ind[1]][1] = 0;
    
}


*
 * Effects: picks parent that is older than x years old (i.e. not random)
 *
int choose_parent_indiv_wagestruct(int ***LANDSCAPE, int num_slots_inuse, double curr_year, int maturity_age, int sizex, int sizey)
{
    
    int event_ind[2], num_tries;
    
    *this function tries at most num_slots_inuse times to look for an adult*
    
    event_ind[0]=get_rand_integ_intvl(1,sizex);
    event_ind[1]=get_rand_integ_intvl(1,sizey);
    
    num_tries=1;
    //while we havent looked at all the individuals and the age of the last individual is under the maturity age:
    while ((num_tries <= num_slots_inuse) & //changed this to 'and' and unsure about it
           ((double) LANDSCAPE[event_ind[0]][event_ind[1]][1]<(double) maturity_age)) { //HELP check this over
        event_ind[0]=get_rand_integ_intvl(1,sizex);
        event_ind[1]=get_rand_integ_intvl(1,sizey);
               ++num_tries;
           }
    * if maxed out on trying to look for an adult, atleast give back a slot that is NOT empty *
    if (num_tries == num_slots_inuse) {
        printf("\nCOULDN'T FIND AND ADULT!!!\n\n");
        event_ind[0]=get_rand_integ_intvl(1,sizex);
        event_ind[1]=get_rand_integ_intvl(1,sizey);
        while (LANDSCAPE[event_ind[0]][event_ind[1]][1] == 0) {
            event_ind[0]=get_rand_integ_intvl(1,sizex);
            event_ind[1]=get_rand_integ_intvl(1,sizey);
        }
    }
    
    return *event_ind;
    
}

**
 * Requires: landscape, indiv array
 * Modifies: ""
 * Effects: Uses dispersal kernel to pick offspring location. If lands outside, does not pick a new parent because that would bias the center
 *
void gen_offspring(int *event_ind, int ***LANDSCAPE, int size_x, int size_y, int *p_num_indiv, int *p_num_slots_inuse, double *disc_norm_cum, int loc_kern_cutoff, int do_age, short *year_born, double curr_year)
{
    
    int ind2_coords[2];
    int i;
    
    *p_num_slots_inuse += 1;
    *p_num_indiv += 1;
    
    //disc_norm_cum[0] = 0 unless set explicitly in the input file. if set explicitly means that dispersal kernel is not used
    if (disc_norm_cum[0] == 2) {
        choose_random_location(size_x, size_y, ind2_coords);
    }
    else if (disc_norm_cum[0] == 3){
        choose_random_nearest_neighbor(event_ind, size_x, size_y, ind2_coords);
    }
    else{
        i = 0;
    select_randloc_frm_cum_disp(event_ind, size_x, size_y, disc_norm_cum, loc_kern_cutoff, ind2_coords); //HELP make sure event ind gets put into cum disp perfectly
         }
    * is  site of offspring already occupied? if yes, chooses a completely new random location from the parent *
    while (LANDSCAPE[ind2_coords[0]][ind2_coords[1]][0]== 1 && disc_norm_cum[0] != 3) {
        
        if (disc_norm_cum[0] == 2) {
            choose_random_location(size_x, size_y, ind2_coords);
        }
        else if (disc_norm_cum[0] == 3)
            choose_random_nearest_neighbor(event_ind, size_x, size_y, ind2_coords);
        else
            select_randloc_frm_cum_disp(event_ind, size_x, size_y, disc_norm_cum, loc_kern_cutoff, ind2_coords);
    }
    
    //Sets the landscape and indiv array to occupied.
    LANDSCAPE[ind2_coords[0]][ind2_coords[1]][0] = 1;
    
    if (do_age)
        LANDSCAPE[ind2_coords[0]][ind2_coords[1]][1] = 0;//HELP used to be:  = myround(curr_year);changed to 0 so it would add maturity age
    
}

*/

