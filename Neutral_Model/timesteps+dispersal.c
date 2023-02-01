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

//sets threshold to average minimum distance between individuals. only use when you change the start landscape because it takes too long to run
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

//Calculates the observed compositions in gaps based on the given coordinate file for the whole landscape
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
        //sets individual's coordinates
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
    for(b =0; b<1;b++){
        if(b==0){
            type = "all";
            thresholdlower = 0;
            thresholdupper = 100;
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

//Calculates the observed compositions in nongaps based on the given coordinate file for the whole landscape
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
        //sets individual's coordinates
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

    for(q=0;q < 9;q++){ //run through all nongaps
    for(b =0; b<1;b++){ 
        if(b==0){
            type = "all";
            thresholdlower = 0;
            thresholdupper = 100;
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

void select_randloc_frm_cum_disp(int *ind1_coords, int size_x, int size_y, double *disc_norm_cum, int loc_cutoff, int *ind2_coords)
{
    // small-dynamic:
    //         *ind1_coords: death location
    //         *ind2_coords: parent location
    int i;
    double toss;
    int offspring_X = 0, offspring_Y = 0;
    
    toss = get_rand_unit();
    // get_rand_unit returns a real random number in the interval [0,1]
   /* for(int u =0; u < 20; u++){
        printf("%f\n", disc_norm_cum[u]);
    }
    */
    i=0;
    
    while((toss > disc_norm_cum[i]) && (i<loc_cutoff)) {
        ++i;
    }
    
    choose_randloc_atdist(i, ind2_coords);
    ind2_coords[0] = ind2_coords[0] + ind1_coords[0];
    ind2_coords[1] = ind2_coords[1] + ind1_coords[1];
    
    
    int count = 0;
     //eliminate edge effects by rechoosing angle at which offspring should land
    while( ind2_coords[0] < 0 || ind2_coords[0] >= size_x || ind2_coords[1] < 0 || ind2_coords[1] >= size_y){
            choose_randloc_atdist(i, ind2_coords);
            ind2_coords[0] = ind2_coords[0] + ind1_coords[0];
            ind2_coords[1] = ind2_coords[1] + ind1_coords[1];
            count++;
            
            if(count > 5000){
                //if we can't land it in the plot, make the plot a tarus and put in the plot
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
            }
    }
}
