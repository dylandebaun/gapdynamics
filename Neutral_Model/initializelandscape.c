//
//  initializelandscape.c
//  gapdynamics
//
//  Created by Dylan DeBaun on 1/20/20.
//  Copyright © 2020 Dylan DeBaun. All rights reserved.
//

#include "initializelandscape.h"
#include "timesteps+dispersal.h"
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include "general.h"
# include <math.h>

void get_input_info(char *infile,
                    float *p_incr_percent,
                    int *p_size_x, int *p_size_y,

                    double *p_kernel_cutoff, int *p_max_disp_dist, int *p_gauss, double *p_sigma, int *p_fat_tail, double *p_u, double *p_p, int *p_global, int *p_near_neigh, double **p_disc_norm_cum, int *p_loc_kern_cutoff,
                    
                    
                    int *p_turnover_time, int *p_num_intev_to_meas, int *p_time_btwn_meas,

                    int *p_output_initial_spec_dist, int *p_output_final_spec_dist, int *p_output_interm_spec_dist, int *p_num_interm_outputs, double **p_when_to_output, int **p_outputted_interm, int *p_output_lowest_x, int *p_output_highest_x, int *p_output_lowest_y, int *p_output_highest_y)
{
    /* note that I've written the above variable names the same as in the program, even though they are actually pointers to those variables in this function */
    
    FILE *ifp;
    char startstr[MAX_STR_SIZE], tempstr[MAX_STR_SIZE];
    char c;
    float fltemp1, fltemp2;
    double parameters[2][2], *disc_norm_cum_TEMP[2];
    int status;
 
    printf("\n--------------------------------\n");
    printf("GETTING INPUT INFO FROM %s...\n", infile);
    ifp = fopen(infile, "r");
    status=fscanf(ifp, "%s", startstr);
    while (status != EOF) {
        if (strcmp(startstr, "PRINT_PROGRESS")==0) {
            fscanf(ifp, "%f", p_incr_percent);
            printf("Will print progress every %f percent of the run.\n", *p_incr_percent);
        }
        if (strcmp(startstr, "SIZE")==0) {
            fscanf(ifp, "%d%d", p_size_x, p_size_y);
            printf("The landscape will be %d by %d in size\n", *p_size_x, *p_size_y);
        }
        
        if (strcmp(startstr, "DISPERSAL")==0) {
            fscanf(ifp, "%s", tempstr);
            if (strcmp(tempstr, "gaussian")==0) {
                p_gauss[0] = TRUE;
                fscanf(ifp, "%f", &fltemp1);
                p_sigma[0] = (double) fltemp1;
                parameters[0][0] = p_sigma[0];
                printf("LARGE: Gaussian dispersal with sigma=%f\n", p_sigma[0]);
                fscanf(ifp, "%f", &fltemp1);
                p_kernel_cutoff[0]= (double) fltemp1;
                fscanf(ifp, "%d", &p_max_disp_dist[0]);
                printf("\nCutoff dispersal kernel at remaining probability %f or max_disp_dist %d\n", p_kernel_cutoff[0], p_max_disp_dist[0]);
                disc_norm_cum_TEMP[0] = discrete_normalized_cummulative(gaussian, parameters[0], p_kernel_cutoff[0], p_max_disp_dist[0], &p_loc_kern_cutoff[0]);
            }
            else if (strcmp(tempstr, "fat_tailed")==0) {
                p_fat_tail[0] = TRUE;
                fscanf(ifp, "%f%f", &fltemp1, &fltemp2);
                p_u[0] = (double) fltemp1;
                p_p[0] = (double) fltemp2;
                parameters[0][0] = p_u[0];
                parameters[0][1] = p_p[0];
                printf("Fat tailed dispersal with u=%f, p=%f\n", p_u[0], p_p[0]);
                fscanf(ifp, "%f", &fltemp1);
                p_kernel_cutoff[0] = (double) fltemp1;
                fscanf(ifp, "%d", &p_max_disp_dist[0]);
                printf("\nCutoff dispersal kernel at remaining probability %f or max_disp_dist %d\n", p_kernel_cutoff[0], p_max_disp_dist[0]);
                disc_norm_cum_TEMP[0] = discrete_normalized_cummulative(fat_tailed, parameters[0], p_kernel_cutoff[0], p_max_disp_dist[0], &p_loc_kern_cutoff[0]);
        
               // fflush(NULL);
            }
            else if (strcmp(tempstr, "global")==0) {
                p_global[0]=TRUE;
                printf("Global dispersal\n");
                disc_norm_cum_TEMP[0] = calloc(1, sizeof(double));
                disc_norm_cum_TEMP[0][0] = 2;
            }
            else if (strcmp(tempstr, "nearest_neighbor")==0) {
                p_near_neigh[0]=TRUE;
                printf("Nearest neighbor dispersal\n");
                disc_norm_cum_TEMP[0] = calloc(1, sizeof(double));
                disc_norm_cum_TEMP[0][0] = 3;
            }
            else {
                printf("\n\nNO DISPERSAL KERNEL WAS CHOSEN!!\n\n");
                exit(1);
            }
            p_disc_norm_cum[0] = disc_norm_cum_TEMP[0];
            
        }
        
        c = getc(ifp);
        while (c != '\n'){
            c=getc(ifp);
        }
        status=fscanf(ifp, "%s", startstr);
        if(status == EOF){
            fclose(ifp);
            printf("closed ifp \n");
        }
    }
    printf("FINISHED GETTING INPUT INFO!\n");
    printf("--------------------------------\n\n");
   // fflush(NULL);
    
}

void read_in_landscape(char *landscapefile, double **LC_xcrd, double **LC_ycrd, int *num_indiv,int **parent_spp_ID, int *num_species){
    FILE *ifp;
    int dbh, sp;
    float x1,y1;
    int i;
    int status;
    double *LC_xcrdtemp, *LC_ycrdtemp;
    int *parent_spp_IDtemp;
    int numspeciestemp;
    
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
    
    int n = 0; //for keeping track of number of species
    numspeciestemp = 0;
    
    printf("Now initializing the landscape from %s\n", landscapefile);
    //col1 = species, col 1=x coord, col2=y coord, col4=dbh (delete the rest of the columns)
    ifp = fopen(landscapefile, "r");
    status=fscanf(ifp, "%d %f %f %d", &sp, &x1, &y1,&dbh);

    //fflush(ifp);
    
    while (status != EOF) {
        //sets individual's coordinates and 
        LC_xcrdtemp[n] = x1;
        LC_ycrdtemp[n] = y1;
        parent_spp_IDtemp[n] = sp;
        if(sp >= numspeciestemp){
            numspeciestemp = sp;
        }
        //Read in the next individual:
        n++;
        status=fscanf(ifp, "%d %f %f %d", &sp, &x1, &y1,&dbh);
        if(status == EOF){
             fclose(ifp);
            printf("closed ifp \n");
        }
        //fflush(ifp);
        //fflush(NULL);
    }
    
   // fclose(ifp);
  //  fflush(NULL);
    *LC_xcrd = LC_xcrdtemp;
    *LC_ycrd = LC_ycrdtemp;
    *parent_spp_ID = parent_spp_IDtemp;
    *num_species = numspeciestemp+1;
    *num_indiv = n;
}

void read_in_gaps(char *gap_file, int **gap_xcrdactual, int **gap_ycrdactual, int **gap_size_xactual, int **gap_size_yactual, int *numbergaps, int sizex, int sizey){
    FILE *ifp;
    int status;
    int xcrd;
    int ycrd;
    int ht1, ht2, ht3,ht4;
    int i = 0;
    int z;
    int *gap_xcrd, *gap_ycrd, *gap_size_x, *gap_size_y;
    int *nongap_xcrd, *nongap_ycrd, *nongap_size_x, *nongap_size_y;
    
    gap_xcrd = calloc(22000, sizeof(int));
    for(z=0;z<22000; z++){
        gap_xcrd[z]=0;
    }
    gap_ycrd = calloc(22000, sizeof(int));
    for(z=0;z<22000; z++){
        gap_ycrd[z]=0;
    }
    gap_size_x = calloc(22000, sizeof(int));
    for(z=0;z<22000; z++){
        gap_size_x[z]=0;
    }
    gap_size_y = calloc(22000, sizeof(int));
    for(z=0;z<22000; z++){
        gap_size_y[z]=0;
    }
    
    ifp = fopen(gap_file, "r");

    status=fscanf(ifp, "%d %d %d %d %d %d\n", &xcrd, &ycrd, &ht1,&ht2,&ht3,&ht4);
    if(xcrd < sizex - 10 && xcrd > 10 && ycrd < sizey - 10 && ycrd > 10){
        if(ht1==0 && ht2==0 && ht3==0 && ht4==0){
        //then it is a gap! so we document the x and y coordinate + 5meters since it is a grid
            gap_xcrd[i]=xcrd;
            gap_ycrd[i]=ycrd;
            gap_size_x[i]=5;
            gap_size_y[i]=5;
            i++;
        }
    }
    while (status != EOF) {
        status=fscanf(ifp, "%d %d %d %d %d %d\n", &xcrd, &ycrd, &ht1,&ht2,&ht3,&ht4);
        if(xcrd < sizex - 10 && xcrd > 10 && ycrd < sizey - 10 && ycrd > 10){
            if(ht1==0 && ht2==0 && ht3==0 && ht4==0){
            //then it is a gap! so we document the x and y coordinate + 5meters since it is a grid
            gap_xcrd[i]=xcrd;
            gap_ycrd[i]=ycrd;
            gap_size_x[i]=5;
            gap_size_y[i]=5;
            i++;
            }
        }
    }
    
    *gap_xcrdactual = gap_xcrd;
    *numbergaps= i; //sets total number of gaps in the landscape
    *gap_ycrdactual = gap_ycrd;
    *gap_size_xactual = gap_size_x;
    *gap_size_yactual = gap_size_y;
    
}


void read_in_nongaps(char *gap_file, int **gap_xcrdactual, int **gap_ycrdactual){
    FILE *ifp;
    int status;
    int xcrd;
    int ycrd;
    int ht1, ht2;
    int i = 0;
    int z;
    int *gap_xcrd, *gap_ycrd;
    
    gap_xcrd = calloc(22000, sizeof(int));
    for(z=0;z<22000; z++){
        gap_xcrd[z]=0;
    }
    gap_ycrd = calloc(22000, sizeof(int));
    for(z=0;z<22000; z++){
        gap_ycrd[z]=0;
    }
    
    ifp = fopen(gap_file, "r");
    
    status=fscanf(ifp, "%d %d %d %d\n", &xcrd, &ycrd, &ht1,&ht2);
            //then it is a gap! so we document the x and y coordinate + 5meters since it is a grid
    gap_xcrd[i]=xcrd;
    gap_ycrd[i]=ycrd;
    i++;
    while (status != EOF) {
        status=fscanf(ifp, "%d %d %d %d\n", &xcrd, &ycrd, &ht1,&ht2);
        gap_xcrd[i]=xcrd;
        gap_ycrd[i]=ycrd;
        i++;
    }
    
    *gap_xcrdactual = gap_xcrd;
    *gap_ycrdactual = gap_ycrd;
}

//THIS FUNCTION REQUIRES INSPECTION. it doesn't do a perfect job, so one must inspect and correct the incorrect information manually (can be seen in below)
void groupgaps(int **gapgroup, int **gapgroupsize, int *gapgroupnum,int *gap_xcrd, int *gap_ycrd,  int numbersmallgaps){
    int *gapgroupt, *gapgroupsizet,*gapgrouptemp;
    int i,j,k;
    int gapgroupnumt;
    
    gapgroupt = calloc(700, sizeof(int));
    for(i=0;i<700; i++){
        gapgroupt[i] = -1;
    }
    gapgrouptemp = calloc(700, sizeof(int));
    for(i=0;i<700; i++){
        gapgrouptemp[i] = -1;
    }
    gapgroupsizet = calloc(700, sizeof(int));
    for(i=0;i<700; i++){
        gapgroupsizet[i] = -1;
    }
    int jx,jy,ix,iy,y;
    int stop =0; //checks if we already have a group or not for our 2 gaps
    int fullstop =0; //checks if our 1 gap is completely standalone
    //set the first group:
    gapgroupt[0] = 0;
    gapgroupnumt = 1;
    gapgroupsizet[0] = 1;
    int x;
    for(i=0; i< numbersmallgaps; i++){
        y=gapgroupt[134];
        fullstop = 0;
        
        //EDIT HERE: here is where the grouping is manually edited
        if(i == 134){
            gapgroupt[134] = 71; //remove gap 72
            gapgroupsize[71] += 1;
        } else if(i==171){
            gapgroupt[171] = 73;//was 74, remove gap 80
            gapgroupsize[73] += 1;
        }else if(i==262){
            gapgroupt[262] = 117;//was 119, remove gap 120
        }else if(i==263){
            gapgroupt[263] = 117;
            gapgroupsize[117] += 2;
        }else if(i==342){
            gapgroupt[342] = 155;//was 158,remove gap 159
            gapgroupsize[155] += 1;
            
        }else{
        for(j=0; j< numbersmallgaps; j++){
            stop = 0;
            if(i == j){
                //printf("");
            }else{
            jx= gap_xcrd[j];
            jy =gap_ycrd[j];
            ix= gap_xcrd[i];
            iy =gap_ycrd[i];
                x=distance(jx, jy, ix, iy);
            //if these two gaps are touching add them to the same group
            if(distance(jx, jy, ix, iy) <= 7.08){
                //to do that we must first make sure they aren't already part of a group
                for(k=0; k<gapgroupnumt; k++){
                    if(gapgroupt[j] == k || gapgroupt[i] == k){
                        //if they are, add the odd one out to the group and take note of that assignment
                        if(gapgroupt[j] == k && gapgroupt[i] == k){
                            gapgroupsizet[k] += 0;
                        }else{
                         gapgroupsizet[k] += 1;
                        }
                        gapgroupt[j] = k;
                        gapgroupt[i] = k;
                        stop = 1;
                        fullstop = 1;
                        break;
                    }
                }
                if(stop == 0){
                    //if they aren't add them both to their own new group
                    gapgroupt[j] = gapgroupnumt;
                    gapgroupt[i] = gapgroupnumt;
                    gapgroupsizet[gapgroupnumt] = 2;
                    gapgroupnumt++;
                    fullstop = 1;
                }
            }
            }
        }
            //if i was not touching any other gaps, add it to it's own group
            if(fullstop == 0){
                gapgroupt[i] = gapgroupnumt;
                gapgroupsizet[gapgroupnumt] = 1;
                gapgroupnumt++;
            }
    }
    }
    
    *gapgroupnum = gapgroupnumt;
    *gapgroup = gapgroupt;
    *gapgroupsize = gapgroupsizet;
}

//for choosing nongap locations
//NOTE: this code may still produce overlapping nongap communities. one must manually inspect nongaps for overlap
void nongapsrandom(int *gapgroup, int *gapgroupsize, int gapgroupnum, int *gap_xcrd, int *gap_ycrd,  int numbersmallgaps, int **nongap_xcrdactual, int **nongap_ycrdactual,int size_x, int size_y){
    int i,j,k;
    
    int *nongap_ycrd, *nongap_xcrd;
    int *listbefore,*listafter;
    
    int z,listend = gapgroupnum;
    
    nongap_xcrd = calloc(22000, sizeof(int));
    for(z=0;z<22000; z++){
        nongap_xcrd[z]=0;
    }
    nongap_ycrd = calloc(22000, sizeof(int));
    for(z=0;z<22000; z++){
        nongap_ycrd[z]=0;
    }
    listbefore = calloc(1000, sizeof(int));
    for(z=0;z<gapgroupnum; z++){
        listbefore[z]= z;
    }
    listafter = calloc(1000, sizeof(int));
    for(z=0;z<gapgroupnum; z++){
        listafter[z]=0;
    }
    
    int lowerx,lowery,size,jx,jy,kx,ky;
    int lowerleftdist,upperleftdist,upperrightdist,lowerrightdist;
    int RectAX1, RectAX2, RectAY1, RectAY2;
    int RectBX1, RectBX2, RectBY1, RectBY2;
    double x;
    int u =0;
    int stop;
    //for each full gap, find a square gap of equivalent area that doesn't overlap gaps or nongaps
    for(i=0;i <gapgroupnum; i++){
        u = gapgroupsize[i];
        size = gapgroupsize[i] * 25;
        stop = 1;
        //check that doesn't overlap gaps
        while(stop == 1){
            lowerx = get_rand_integ_intvl(10,990);
            lowery = get_rand_integ_intvl(10,490);
            RectAX1 = lowerx;
            RectAX2 = lowerx + sqrt(size);
            RectAY2 = lowery;
            RectAY1 = lowery + sqrt(size);
            while(RectAX2 >size_x || RectAY2 > size_y){
                lowerx = get_rand_integ_intvl(10,990);
                lowery = get_rand_integ_intvl(10,490);
                RectAX1 = lowerx;
                RectAX2 = lowerx + sqrt(size);
                RectAY2 = lowery;
                RectAY1 = lowery + sqrt(size);
            }
            for(j =0;j<numbersmallgaps;j++){
                if(i ==244 && j == 39){
                    printf("pause");
                }
                kx= gap_xcrd[j];
                ky =gap_ycrd[j];
                RectBX1 = kx;
                RectBX2 = kx + 5;
                RectBY2 = ky;
                RectBY1 = ky + 5;
                //rectangles overlap
                if((RectAX1 < RectBX2 && RectAX2 > RectBX1 &&
                    RectAY1 > RectBY2 && RectAY2 < RectBY1 )){
                    stop = 1;
                    break;
                }else{
                    stop = 0;
                }
            }
            //if didn't overlap with any small gaps, check if overlaps with other gaps
            if(stop == 0){
                for(j = 0; j< i; j++){
                    kx= nongap_xcrd[j];
                    ky =nongap_ycrd[j];
                    RectBX1 = kx;
                    RectBX2 = kx + sqrt(gapgroupsize[j] * 25);
                    RectBY2 = ky;
                    RectBY1 = ky + sqrt(gapgroupsize[j] * 25);
                    //rectangles  overlap
                    if((RectAX1 < RectBX2 && RectAX2 > RectBX1 &&
                       RectAY1 > RectBY2 && RectAY2 < RectBY1 )){
                        stop = 1;
                        break;
                    }
                }
            }
        }
        //successfully exited while loop
        nongap_xcrd[i] = RectAX1;
        nongap_ycrd[i] = RectAY2;
    }

    
    //now check for nongaps with 0,0
    *nongap_xcrdactual = nongap_xcrd;
    *nongap_ycrdactual = nongap_ycrd;
}
