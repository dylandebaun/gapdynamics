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
                    int *p_do_age, int *p_maturity_age,
                    /*
                    int *p_spec_gap_meas, int *p_num_gap, int ***p_gap_locs, int ***p_gap_sizes,*/
                    

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
        
 //0= SMALL, 1= LARGE
        if (strcmp(startstr, "DISPERSAL_SMALL")==0) {
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
            
        if (strcmp(startstr, "DISPERSAL_LARGE")==0) {
            fscanf(ifp, "%s", tempstr);
            if (strcmp(tempstr, "gaussian")==0) {
                p_gauss[1] = TRUE;
                fscanf(ifp, "%f", &fltemp1);
                p_sigma[1] = (double) fltemp1;
                parameters[1][0] = p_sigma[1];
                printf("LARGE: Gaussian dispersal with sigma=%f\n", p_sigma[1]);
                fscanf(ifp, "%f", &fltemp1);
                p_kernel_cutoff[1]= (double) fltemp1;
                fscanf(ifp, "%d", &p_max_disp_dist[1]);
                printf("\nCutoff dispersal kernel at remaining probability %f or max_disp_dist %d\n", p_kernel_cutoff[1], p_max_disp_dist[1]);
                disc_norm_cum_TEMP[1] = discrete_normalized_cummulative(gaussian, parameters[1], p_kernel_cutoff[1], p_max_disp_dist[1], &p_loc_kern_cutoff[1]);
            }
            
            else if (strcmp(tempstr, "fat_tailed")==0) {
                p_fat_tail[1] = TRUE;
                fscanf(ifp, "%f%f", &fltemp1, &fltemp2);
                p_u[1] = (double) fltemp1;
                p_p[1] = (double) fltemp2;
                parameters[1][0] = p_u[1];
                parameters[1][1] = p_p[1];
                printf("Fat tailed dispersal with u=%f, p=%f\n", p_u[1], p_p[1]);
                fscanf(ifp, "%f", &fltemp1);
                p_kernel_cutoff[1] = (double) fltemp1;
                fscanf(ifp, "%d", &p_max_disp_dist[1]);
                printf("\nCutoff dispersal kernel at remaining probability %f or max_disp_dist %d\n", p_kernel_cutoff[1], p_max_disp_dist[1]);
                disc_norm_cum_TEMP[1] = discrete_normalized_cummulative(fat_tailed, parameters[1], p_kernel_cutoff[1], p_max_disp_dist[1], &p_loc_kern_cutoff[1]);
                printf("Finished calculating dispersal kernel.\n");
                //fflush(NULL);
            }
            else if (strcmp(tempstr, "global")==0) {
                p_global[1]=TRUE;
                printf("Global dispersal\n");
                disc_norm_cum_TEMP[1] = calloc(1, sizeof(double));
                disc_norm_cum_TEMP[1][0] = 2;
            }
            else if (strcmp(tempstr, "nearest_neighbor")==0) {
                p_near_neigh[1]=TRUE;
                printf("Nearest neighbor dispersal\n");
                disc_norm_cum_TEMP[1] = calloc(1, sizeof(double));
                disc_norm_cum_TEMP[1][0] = 3;
            }
            else {
                printf("\n\nNO DISPERSAL KERNEL WAS CHOSEN!!\n\n");
                exit(1);
            }
            p_disc_norm_cum[1] = disc_norm_cum_TEMP[1];
            
        }
        
        if(strcmp(startstr, "AGE_STRUCTURE")==0) {
            *p_do_age=TRUE;
            fscanf(ifp, "%d", p_maturity_age);
            printf("\nWill not let individuals have offspring until they are older than %d yrs\n", *p_maturity_age);
        }
        
        c = getc(ifp);
        while (c != '\n'){
            c=getc(ifp);
        }
        status=fscanf(ifp, "%s", startstr);
        if(status == EOF){
            fclose(ifp);
            printf("closed ifp");
        }
    }
    printf("FINISHED GETTING INPUT INFO!\n");
    printf("--------------------------------\n\n");
   // fflush(NULL);
    
}

void read_in_landscape(char *landscapefile, double **LC_xcrd, double **LC_ycrd, int **age, int *p_maturity_age, int *p_do_age, int *num_indiv,int **parent_spp_ID, int *num_species){
    FILE *ifp;
    int dbh, sp;
    float x1,y1;
    int i;
    int status;
    double *LC_xcrdtemp, *LC_ycrdtemp;
    int *agetemp, *parent_spp_IDtemp;
    int numspeciestemp;
    
   LC_xcrdtemp = calloc(236000, sizeof(double));
    for(i=0;i<236000; i++){
        LC_xcrdtemp[i]=0;
    }
    LC_ycrdtemp = calloc(236000, sizeof(double));
    for(i=0;i<236000; i++){
        LC_ycrdtemp[i]=0;
    }
    
    agetemp = calloc(236000, sizeof(int));
    for(i=0;i<236000; i++){
        agetemp[i]=0;
    }
    
    parent_spp_IDtemp = calloc(236000, sizeof(int));
    for(i=0;i<236000; i++){
        parent_spp_IDtemp[i]=0;
    }
    
    int adult = *p_maturity_age + 1; //HELP Check this
    int n = 0; //for keeping track of number of species
    numspeciestemp = 0;
    
    printf("Now initializing the landscape from %s\n", landscapefile);
    //col1 = species, col 1=x coord, col2=y coord, col4=dbh (delete the rest of the columns)
    ifp = fopen(landscapefile, "r");
    status=fscanf(ifp, "%d %f %f %d", &sp, &x1, &y1,&dbh);

    //fflush(ifp);
    
    while (status != EOF) {
        //sets individual's coordinates and age
        LC_xcrdtemp[n] = x1;
        LC_ycrdtemp[n] = y1;
        parent_spp_IDtemp[n] = sp;
        if(sp >= numspeciestemp){
            numspeciestemp = sp;
        }
        if(p_do_age){
            agetemp[n] = adult;
            /*if(dbh>10){
                agetemp[n] = adult; //sets age
            }
            else if (dbh > 5){
                agetemp[n] = 5; //I think this choice will need to be better assigned
            }
            else{
                agetemp[n] = 1;
            }*/
        }
        //Read in the next individual:
        n++;
        if(n == 20717){
          printf("no");
        }
        status=fscanf(ifp, "%d %f %f %d", &sp, &x1, &y1,&dbh);
        if(status == EOF){
             fclose(ifp);
            printf("closed ifp");
        }
        //fflush(ifp);
        //fflush(NULL);
    }
    
   // fclose(ifp);
  //  fflush(NULL);
    *LC_xcrd = LC_xcrdtemp;
    *LC_ycrd = LC_ycrdtemp;
    *age = agetemp;
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

//THIS FUNCTION REQUIRES INSPECTION. it doesn't do a perfect job, sometimes misses a grouping, must inspect and set the incorrect groups below
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
                if(gapgroupt[134] != -1){
                    //printf("pause");
                }
                if(i == 154 && j == 134){
                    //printf("%d,%d\n", gapgroupt[i], gapgroupt[j]);
                }
                if(j == 134){
                   // printf("%d\n", gapgroupt[134]);
                }
                if(i == 133 && j == 134){
                   // printf("%d,%d\n", gapgroupt[i], gapgroupt[j]);
                }
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

void movegaps(int *listend, int numbersmallgaps, int *gapgroup, int *gap_xcrd,int *gap_ycrd, int size_x, int size_y, int *listbefore, int **listafter, int incx, int incy, int **nongap_xcrdreal,int **nongap_ycrdreal){
    int *listA, z;
    listA = calloc(1000, sizeof(int));
    for(z=0;z<1000; z++){
        listA[z]=0;
    }
    int jx, jy,kx,ky, i, t,j,k;
    int test1,test2,test3,test4;
    int stop;
    int *nongap_xcrd,*nongap_ycrd;
    int b=0;

    nongap_xcrd = calloc(22000, sizeof(int));
    for(z=0;z<22000; z++){
        nongap_xcrd[z] = 0;
    }
    nongap_ycrd = calloc(22000, sizeof(int));
    for(z=0;z<22000; z++){
        nongap_ycrd[z] = 0;
    }
    nongap_xcrd = *nongap_xcrdreal;
    nongap_ycrd = *nongap_ycrdreal;
    int x;
    for(t = 0; t < *listend; t++){
        i = listbefore[t];
        if(i == 53){
           // printf("pause");
        }
        stop = 0;
        //for all gaps in groupgap i
        for(j=0; j< numbersmallgaps; j++){
            if(gapgroup[j] == i){
                if(j == 134){
                   // printf("pause");
                }
                // printf("%d\n",1);
                //if the gap +x lands in another gap
                for(k =0; k< numbersmallgaps;k++){
                    if(k == 110){
                       // printf("pause");
                    }
                    //also make sure it lands 10 m away from other gaps
                    jx= gap_xcrd[j];
                    jy =gap_ycrd[j];
                    kx= gap_xcrd[k];
                    ky =gap_ycrd[k];
                    x = distance(jx + incx, jy + incy, kx, ky);
                    if((distance(jx + incx, jy + incy, kx, ky) < 6) || jx + incx > size_x - 10 || jx + incx < 10 || jy + incy > size_y - 10 || jy + incy < 10){
                        test1 =gap_xcrd[j];
                        test2 =gap_xcrd[k];
                        test3 =gap_ycrd[j];
                        test4 =  gap_ycrd[k];
                        stop = 1;
                        listA[b] = i;
                        b++;
                        break;
                    }
                }
            }
            //if one of the gaps fell outside, break. no need to run through rest
            if(stop == 1){
                break;
            }
        }
        //if we never broke, assign that value to all nongaps
        if(stop == 0){
            //for all gaps in groupgap i we are going to give them this assignment
            for(j=0; j< numbersmallgaps; j++){
                if(gapgroup[j] == i){
                    test1 =gap_xcrd[j];
                    test2 =gap_ycrd[j];
                    nongap_xcrd[j] = gap_xcrd[j] + incx;
                    nongap_ycrd[j] = gap_ycrd[j] + incy;
                    test3 =nongap_xcrd[j];
                }
            }
        }
    }
    *nongap_xcrdreal = nongap_xcrd;
    *nongap_ycrdreal = nongap_ycrd;
    *listafter = listA;
    *listend = b;
}

void nongaps(int *gapgroup, int *gapgroupsize, int gapgroupnum, int *gap_xcrd, int *gap_ycrd,  int numbersmallgaps, int **nongap_xcrdactual, int **nongap_ycrdactual,int size_x, int size_y){
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
   // printf("%d\n",listbefore[2]);
    movegaps(&listend, numbersmallgaps, gapgroup, gap_xcrd, gap_ycrd, size_x, size_y, listbefore, &listafter, 0, -10, &nongap_xcrd, &nongap_ycrd);
    movegaps(&listend, numbersmallgaps, gapgroup, gap_xcrd, gap_ycrd, size_x, size_y, listafter, &listafter, 10, 0, &nongap_xcrd, &nongap_ycrd);
    movegaps(&listend, numbersmallgaps, gapgroup, gap_xcrd, gap_ycrd, size_x, size_y, listafter, &listafter, 0, 10, &nongap_xcrd, &nongap_ycrd);
    movegaps(&listend, numbersmallgaps, gapgroup, gap_xcrd, gap_ycrd, size_x, size_y, listafter, &listafter, -10, 0, &nongap_xcrd, &nongap_ycrd);
    movegaps(&listend, numbersmallgaps, gapgroup, gap_xcrd, gap_ycrd, size_x, size_y, listafter, &listafter, -10, -10, &nongap_xcrd, &nongap_ycrd);
    movegaps(&listend, numbersmallgaps, gapgroup, gap_xcrd, gap_ycrd, size_x, size_y, listafter, &listafter, 10, 10, &nongap_xcrd, &nongap_ycrd);
    movegaps(&listend, numbersmallgaps, gapgroup, gap_xcrd, gap_ycrd, size_x, size_y, listafter, &listafter, 20, 0, &nongap_xcrd, &nongap_ycrd);
    movegaps(&listend, numbersmallgaps, gapgroup, gap_xcrd, gap_ycrd, size_x, size_y, listafter, &listafter, 0, 20, &nongap_xcrd, &nongap_ycrd);
    movegaps(&listend, numbersmallgaps, gapgroup, gap_xcrd, gap_ycrd, size_x, size_y, listafter, &listafter, -20, 0, &nongap_xcrd, &nongap_ycrd);
    movegaps(&listend, numbersmallgaps, gapgroup, gap_xcrd, gap_ycrd, size_x, size_y, listafter, &listafter, 0, -20, &nongap_xcrd, &nongap_ycrd);
    movegaps(&listend, numbersmallgaps, gapgroup, gap_xcrd, gap_ycrd, size_x, size_y, listafter, &listafter, 20, 20, &nongap_xcrd, &nongap_ycrd);
    movegaps(&listend, numbersmallgaps, gapgroup, gap_xcrd, gap_ycrd, size_x, size_y, listafter, &listafter, -20, -20, &nongap_xcrd, &nongap_ycrd);
    //printf("%d\n",listafter[0]);
    movegaps(&listend, numbersmallgaps, gapgroup, gap_xcrd, gap_ycrd, size_x, size_y, listafter, &listafter, 50, 0, &nongap_xcrd, &nongap_ycrd);
    movegaps(&listend, numbersmallgaps, gapgroup, gap_xcrd, gap_ycrd, size_x, size_y, listafter, &listafter, 0, 50, &nongap_xcrd, &nongap_ycrd);
    movegaps(&listend, numbersmallgaps, gapgroup, gap_xcrd, gap_ycrd, size_x, size_y, listafter, &listafter, 60, 0, &nongap_xcrd, &nongap_ycrd);

    //now check for nongaps with 0,0
    *nongap_xcrdactual = nongap_xcrd;
    *nongap_ycrdactual = nongap_ycrd;
}

void movegapsrandom(int *listend, int numbersmallgaps, int *gapgroup, int *gap_xcrd,int *gap_ycrd, int size_x, int size_y, int *listbefore, int **listafter, int **nongap_xcrdreal,int **nongap_ycrdreal){
    int *listA, z;
    listA = calloc(1000, sizeof(int));
    for(z=0;z<1000; z++){
        listA[z]=0;
    }
    int jx, jy,kx,ky, i, t,j,k;
    int test1,test2,test3,test4;
    int stop;
    int *nongap_xcrd,*nongap_ycrd;
    int b=0;
    int incx = 0, incy =0;
    int numgaps = 0;
    
    nongap_xcrd = calloc(22000, sizeof(int));
    for(z=0;z<22000; z++){
        nongap_xcrd[z] = 0;
    }
    nongap_ycrd = calloc(22000, sizeof(int));
    for(z=0;z<22000; z++){
        nongap_ycrd[z] = 0;
    }
    nongap_xcrd = *nongap_xcrdreal;
    nongap_ycrd = *nongap_ycrdreal;
    int x;
    for(t = 0; t < *listend; t++){
        i = listbefore[t];
        incx = get_rand_integ_intvl(10,490);
        incy = get_rand_integ_intvl(10,490);
        numgaps = 0;
        if(i == 53){
            // printf("pause");
        }
        stop = 0;
        //for all gaps in groupgap i
        for(j=0; j< numbersmallgaps; j++){
            if(gapgroup[j] == i){
                numgaps++;
                if(j == 134){
                    // printf("pause");
                }
                // printf("%d\n",1);
                //if the gap +x lands in another gap
                for(k =0; k< numbersmallgaps;k++){
                    if(k == 110){
                        // printf("pause");
                    }
                    //also make sure it lands 10 m away from other gaps
                    jx= gap_xcrd[j];
                    jy =gap_ycrd[j];
                    kx= gap_xcrd[k];
                    ky =gap_ycrd[k];
                    x = distance(jx + incx, jy + incy, kx, ky);
                    if(numgaps == 1){
                        while(jx + incx > size_x || jx + incx < 10 || jy + incy > size_y || jy + incy < 10 || x < 6){
                            incx = get_rand_integ_intvl(10,490);
                            incy = get_rand_integ_intvl(10,490);
                        }
                    }
                    if((distance(jx + incx, jy + incy, kx, ky) < 6) || jx + incx > size_x || jx + incx < 10 || jy + incy > size_y || jy + incy < 10){
                        test1 =gap_xcrd[j];
                        test2 =gap_xcrd[k];
                        test3 =gap_ycrd[j];
                        test4 =  gap_ycrd[k];
                        stop = 1;
                        listA[b] = i;
                        b++;
                        break;
                    }
                }
            }
            //if one of the gaps fell outside, break. no need to run through rest
            if(stop == 1){
                break;
            }
        }
        //if we never broke, assign that value to all nongaps
        if(stop == 0){
            //for all gaps in groupgap i we are going to give them this assignment
            for(j=0; j< numbersmallgaps; j++){
                if(gapgroup[j] == i){
                    test1 =gap_xcrd[j];
                    test2 =gap_ycrd[j];
                    nongap_xcrd[j] = gap_xcrd[j] + incx;
                    nongap_ycrd[j] = gap_ycrd[j] + incy;
                    test3 =nongap_xcrd[j];
                }
            }
        }
    }
    *nongap_xcrdreal = nongap_xcrd;
    *nongap_ycrdreal = nongap_ycrd;
    *listafter = listA;
    *listend = b;
}

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
    // printf("%d\n",listbefore[2]);
    /*movegapsrandom(&listend, numbersmallgaps, gapgroup, gap_xcrd, gap_ycrd, size_x, size_y, listbefore, &listafter,  &nongap_xcrd, &nongap_ycrd);
    while(listend != 0){
        listbefore = listafter;
        movegapsrandom(&listend, numbersmallgaps, gapgroup, gap_xcrd, gap_ycrd, size_x, size_y, listbefore, &listafter,  &nongap_xcrd, &nongap_ycrd);
    }*/
    
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








/*
 listB= calloc(1000, sizeof(int));
 for(z=0;z<1000; z++){
 listB[z]=0;
 }
 listC= calloc(1000, sizeof(int));
 for(z=0;z<1000; z++){
 listC[z]=0;
 }
 listD= calloc(1000, sizeof(int));
 for(z=0;z<1000; z++){
 listD[z]=0;
 }
 listE= calloc(1000, sizeof(int));
 for(z=0;z<1000; z++){
 listE[z]=0;
 }
 listF= calloc(1000, sizeof(int));
 for(z=0;z<1000; z++){
 listF[z]=0;
 }
 for(i = 0; i < gapgroupnum; i++){
 /*   if(i == 296){
 printf("%d\n",1);
 }
stop = 0;
//for all gaps in groupgap i
for(j=0; j< numbersmallgaps; j++){
    if(gapgroup[j] == i){
        // printf("%d\n",1);
        //if the gap +x lands in another gap
        for(k =0; k< numbersmallgaps;k++){
            //also make sure it lands 10 m away from other gaps
            if(k == 311 && j==309){
                printf("pause");
            }
            jx= gap_xcrd[j];
            jy =gap_ycrd[j];
            kx= gap_xcrd[k];
            ky =gap_ycrd[k];
            if((distance(jx + 10, jy, kx, ky) < 6) || jx + 10 > size_x - 10){
                test1 =gap_xcrd[j];
                test2 =gap_xcrd[k];
                test3 =gap_ycrd[j];
                test4 =  gap_ycrd[k];
                stop = 1;
                
                listA[b] = i;
                b++;
                break;
            }
        }
    }
    //if one of the gaps fell outside, break. no need to run through rest
    if(stop == 1){
        break;
    }
}
//if we never broke, assign that value to all nongaps
/*if(i >= 64){
 printf("%d\n",nongap_xcrd[115]);
 }
if(stop == 0){
    //for all gaps in groupgap i we are going to give them this assignment
    for(j=0; j< numbersmallgaps; j++){
        if(gapgroup[j] == i){
            test1 =gap_xcrd[j];
            test2 =gap_ycrd[j];
            nongap_xcrd[j] = gap_xcrd[j] + 10;
            nongap_ycrd[j] = gap_ycrd[j];
            test3 =nongap_xcrd[j];
        }
    }
}
}
//printf("break\n");
/*for(j=0; j< numbersmallgaps; j++){
 printf("%d\n",nongap_xcrd[j]);
 }
//now let's look at the things that got rejected
for(t = 0; t < (b); t++){
    i = listA[t];
    stop = 0;
    //for all gaps in groupgap i
    for(j=0; j< numbersmallgaps; j++){
        if(gapgroup[j] == i){
            //if the gap +5 lands in another gap
            for(k =0; k< numbersmallgaps;k++){
                jx= gap_xcrd[j];
                jy =gap_ycrd[j];
                kx= gap_xcrd[k];
                ky =gap_ycrd[k];
                if((distance(jx, jy + 10, kx, ky) < 6) || jy+10 > size_y -10){
                    // if(gap_xcrd[j] == gap_xcrd[k] && gap_ycrd[j] + 20 == gap_ycrd[k] &&  gap_ycrd[j] + 20 < size_y -20){
                    listB[o] = i;
                    o++;
                    stop = 1;
                    break;
                }
            }
        }
        //if one of the gaps fell outside, break. no need to run through rest
        if(stop == 1){
            break;
        }
    }
    //if we never broke, assign that value to all nongaps
    if(stop == 0){
        //for all gaps in groupgap i we are going to give them this assignment
        for(j=0; j< numbersmallgaps; j++){
            if(gapgroup[j] == i){
                nongap_xcrd[j] = gap_xcrd[j];
                nongap_ycrd[j] = gap_ycrd[j] + 10;
            }
        }
    }
}
//now let's look at the things that got rejected...again

for(t = 0; t < (o); t++){
    i = listB[t];
    stop = 0;
    //for all gaps in groupgap i
    for(j=0; j< numbersmallgaps; j++){
        if(gapgroup[j] == i){
            //if the gap +5 lands in another gap
            for(k =0; k< numbersmallgaps;k++){
                jx= gap_xcrd[j];
                jy =gap_ycrd[j];
                kx= gap_xcrd[k];
                ky =gap_ycrd[k];
                if((distance(jx - 40, jy, kx, ky) < 6) || jx-40 < 10){//if(gap_xcrd[j] == gap_xcrd[k]-20 && gap_ycrd[j] == gap_ycrd[k] &&  gap_xcrd[j] - 20 > 20){
                    listC[h] = i;
                    h++;
                    stop = 1;
                    break;
                }
            }
        }
        //if one of the gaps fell outside, break. no need to run through rest
        if(stop == 1){
            break;
        }
    }
    //if we never broke, assign that value to all nongaps
    if(stop == 0){
        //for all gaps in groupgap i we are going to give them this assignment
        for(j=0; j< numbersmallgaps; j++){
            if(gapgroup[j] == i){
                nongap_xcrd[j] = gap_xcrd[j] - 40;
                nongap_ycrd[j] = gap_ycrd[j];
            }
        }
    }
}
int l=0;
for(t = 0; t < (h); t++){
    i = listC[t];
    stop = 0;
    //for all gaps in groupgap i
    for(j=0; j< numbersmallgaps; j++){
        if(gapgroup[j] == i){
            //if the gap +5 lands in another gap
            for(k =0; k< numbersmallgaps;k++){
                jx= gap_xcrd[j];
                jy =gap_ycrd[j];
                kx= gap_xcrd[k];
                ky =gap_ycrd[k];
                if((distance(jx, jy - 40, kx, ky) < 6) || jy-40 < 10){
                    listD[m] = i;
                    m++;
                    stop = 1;
                    break;
                }
            }
        }
        //if one of the gaps fell outside, break. no need to run through rest
        if(stop == 1){
            break;
        }
    }
    //if we never broke, assign that value to all nongaps
    if(stop == 0){
        //for all gaps in groupgap i we are going to give them this assignment
        for(j=0; j< numbersmallgaps; j++){
            if(gapgroup[j] == i){
                nongap_xcrd[j] = gap_xcrd[j] ;
                nongap_ycrd[j] = gap_ycrd[j]-40;
            }
        }
    }
}
for(t = 0; t < m; t++){
    i = listD[t];
    stop = 0;
    //for all gaps in groupgap i
    for(j=0; j< numbersmallgaps; j++){
        if(gapgroup[j] == i){
            //if the gap +5 lands in another gap
            for(k =0; k< numbersmallgaps;k++){
                jx= gap_xcrd[j];
                jy =gap_ycrd[j];
                kx= gap_xcrd[k];
                ky =gap_ycrd[k];
                if((distance(jx + 40, jy, kx, ky) < 6) || jx+40 > size_x-10){
                    // if(gap_xcrd[j] + 10 == gap_xcrd[k] && gap_ycrd[j]+10 == gap_ycrd[k] &&  gap_xcrd[j] + 10 <size_x+ 10 && gap_ycrd[j] + 10 <size_y+ 10){
                    listE[m] = i;
                    l++;
                    stop = 1;
                    break;
                }
            }
        }
        //if one of the gaps fell outside, break. no need to run through rest
        if(stop == 1){
            break;
        }
    }
    //if we never broke, assign that value to all nongaps
    if(stop == 0){
        //for all gaps in groupgap i we are going to give them this assignment
        for(j=0; j< numbersmallgaps; j++){
            if(gapgroup[j] == i){
                nongap_xcrd[j] = gap_xcrd[j]+20;
                nongap_ycrd[j] = gap_ycrd[j]+20;
            }
        }
    }
}
for(t = 0; t < l; t++){
    i = listE[t];
    stop = 0;
    //for all gaps in groupgap i
    for(j=0; j< numbersmallgaps; j++){
        if(gapgroup[j] == i){
            //if the gap +5 lands in another gap
            for(k =0; k< numbersmallgaps;k++){
                jx= gap_xcrd[j];
                jy =gap_ycrd[j];
                kx= gap_xcrd[k];
                ky =gap_ycrd[k];
                if((distance(jx - 20, jy -20, kx, ky) < 6) || jx -20 <10|| jy - 20 < 10){
                    //listE[m] = i;
                    c++;
                    stop = 1;
                    break;
                }
            }
        }
        //if one of the gaps fell outside, break. no need to run through rest
        if(stop == 1){
            break;
        }
    }
    //if we never broke, assign that value to all nongaps
    if(stop == 0){
        //for all gaps in groupgap i we are going to give them this assignment
        for(j=0; j< numbersmallgaps; j++){
            if(gapgroup[j] == i){
                nongap_xcrd[j] = gap_xcrd[j]-20;
                nongap_ycrd[j] = gap_ycrd[j]-20;
            }
        }
    }
}

for(j=0; j< numbersmallgaps; j++){
 if(gapgroup[j] == 71){
 nongap_xcrd[j] = gap_xcrd[j];
 nongap_ycrd[j] = gap_ycrd[j]-55;
 }
 }
 for(j=0; j< numbersmallgaps; j++){
 if(gapgroup[j] == 117){
 nongap_xcrd[j] = gap_xcrd[j]+50;
 nongap_ycrd[j] = gap_ycrd[j];
 }
 }
 for(j=0; j< numbersmallgaps; j++){
 if(gapgroup[j] == 299){
 nongap_xcrd[j] = gap_xcrd[j];
 nongap_ycrd[j] = gap_ycrd[j]-45;
 }
 }*/
/*void group_gaps(int **gapgroup, int **gapgroupsize, int *gapgroupnum, int realgapnum){
    int *gapgroupt, *gapgroupsizet,*gapgrouptemp;
    int i;
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
    int p = 0;
    int k;
    int y;
    
    y=0;
    gapgroupt[20] = p;
    y++;
    gapgroupsizet[p] = y;
    p++;
    
    y=0;
    gapgroupt[27] = p;
    y++;
    gapgroupsizet[p] = y;
    p++;
    
    y=0;
    gapgroupt[21] = p;
    y++;
    gapgroupsizet[p] = y;
    p++;
    
    y=0;
    gapgroupt[50] = p;
    y++;
    gapgroupsizet[p] = y;
    p++;
    
    y=0;
    gapgroupt[95] = p;
    y++;
    gapgroupsizet[p] = y;
    p++;
    
    y=0;
    gapgroupt[163] = p;
    y++;
    gapgroupsizet[p] = y;
    p++;
    
    y=0;
    gapgroupt[215] = p;
    y++;
    gapgroupsizet[p] = y;
    p++;
    
    y=0;
    gapgroupt[243] = p;
    y++;
    gapgroupsizet[p] = y;
    p++;
    
    y=0;
    gapgroupt[336] = p;
    y++;
    gapgroupsizet[p] = y;
    p++;
    
    y=0;
    gapgroupt[400] = p;
    y++;
    gapgroupsizet[p] = y;
    p++;
    
    y=0;
    gapgroupt[431] = p;
    y++;
    gapgroupsizet[p] = y;
    p++;
    
    y=0;
    gapgroupt[596] = p;
    y++;
    gapgroupsizet[p] = y;
    p++;
    
    y=0;
    gapgroupt[603] = p;
    y++;
    gapgroupsizet[p] = y;
    p++;
    
    y=0;
    gapgroupt[610] = p;
    y++;
    gapgroupsizet[p] = y;
    p++;
    
    y=0;
    gapgroupt[611] = p;
    y++;
    gapgroupsizet[p] = y;
    p++;
    
    y=0;
    gapgroupt[612] = p;
    y++;
    gapgroupsizet[p] = y;
    p++;
    
    y=0;
    gapgroupt[586] = p;
    y++;
    gapgroupsizet[p] = y;
    p++;
    
    y=0;
    gapgroupt[551] = p;
    y++;
    gapgroupsizet[p] = y;
    p++;
    
    y=0;
    gapgroupt[550] = p;
    y++;
    gapgroupsizet[p] = y;
    p++;
    
    y=0;
    gapgroupt[511] = p;
    y++;
    gapgroupsizet[p] = y;
    p++;
    
    y=0;
    gapgroupt[510] = p;
    y++;
    gapgroupsizet[p] = y;
    p++;
    
    y=0;
    gapgroupt[428] = p;
    y++;
    gapgroupsizet[p] = y;
    p++;
    
    y=0;
    gapgroupt[204] = p;
    y++;
    gapgroupsizet[p] = y;
    p++;
    
    y=0;
    gapgroupt[226] = p;
    y++;
    gapgroupsizet[p] = y;
    p++;
    
    y=0;
    for(k = 7; k <= 8; k++){
        gapgroupt[k] = p;
        y++;
    }
    gapgroupsizet[p] = y;
    p++;
    
    y=0;
    for(k = 1; k <= 2; k++){
        gapgroupt[k] = p;
        y++;
    }
    gapgroupsizet[p] = y;
    p++;
    
    y=0;
    for(k = 9; k <= 10; k++){
        gapgroupt[k] = p;
        y++;
    }
    gapgroupsizet[p] = y;
    p++;
    
    y=0;
    for(k = 11; k <= 16; k++){
        gapgroupt[k] = p;
        y++;
    }
    gapgroupsizet[p] = y;
    p++;
    
    y=0;
    for(k = 17; k <= 19; k++){
        gapgroupt[k] = p;
        y++;
    }
    gapgroupsizet[p] = y;
    p++;
    
    y=0;
    for(k = 128; k <= 130; k++){
        gapgroupt[k] = p;
        y++;
    }
    for(k = 137; k <= 147; k++){
        gapgroupt[k] = p;
        y++;
    }
    gapgroupsizet[p] = y;
    p++;
    
    y=0;
    for(k = 159; k <= 162; k++){
        gapgroupt[k] = p;
        y++;
    }
    for(k = 184; k <= 186; k++){
        gapgroupt[k] = p;
        y++;
    }
    gapgroupsizet[p] = y;
    p++;
    
    y=0;
    for(k = 244; k <= 245; k++){
        gapgroupt[k] = p;
        y++;
    }
    for(k = 258; k <= 260; k++){
        gapgroupt[k] = p;
        y++;
    }
    gapgroupt[242] = p;
    gapgroupsizet[p] = y;
    p++;
    
    
    y=0;
    for(k = 590; k <= 594; k++){
        gapgroupt[k] = p;
        y++;
    }
    gapgroupsizet[p] = y;
    p++;
    
    
    y=0;
    for(k = 598; k <= 602; k++){
        gapgroupt[k] = p;
        y++;
    }
    gapgroupsizet[p] = y;
    p++;
    
    
    y=0;
    for(k = 605; k <= 607; k++){
        gapgroupt[k] = p;
        y++;
    }
    gapgroupsizet[p] = y;
    p++;
    
    
    y=0;
    for(k = 581; k <= 585; k++){
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
    for(k = 350; k <= 351; k++){
        gapgroupt[k] = p;
        y++;
    }
    for(k = 358; k <= 364; k++){
        gapgroupt[k] = p;
        y++;
    }
    gapgroupsizet[p] = y;
    p++;
    
    y=0;
    for(k = 383; k <= 394; k++){
        gapgroupt[k] = p;
        y++;
    }
    gapgroupt[371] = p;
    y++;
    gapgroupsizet[p] = y;
    p++;
    
    y=0;
    for(k = 309; k <= 310; k++){
        gapgroupt[k] = p;
        y++;
    }
    gapgroupsizet[p] = y;
    p++;
    
    y=0;
    for(k = 276; k <= 277; k++){
        gapgroupt[k] = p;
        y++;
    }
    for(k = 282; k <= 287; k++){
        gapgroupt[k] = p;
        y++;
    }
    gapgroupsizet[p] = y;
    p++;
    
    y=0;
    for(k = 535; k <= 537; k++){
        gapgroupt[k] = p;
        y++;
    }
    gapgroupsizet[p] = y;
    p++;
    
    y=0;
    for(k = 421; k <= 422; k++){
        gapgroupt[k] = p;
        y++;
    }
    gapgroupt[401] = p;
    gapgroupsizet[p] = y;
    p++;
    
    
    
    y=0;
    for(k = 22; k <= 23; k++){
        gapgroupt[k] = p;
        y++;
    }
    gapgroupsizet[p] = y;
    p++;
    
    y=0;
    for(k = 4; k <= 6; k++){
        gapgroupt[k] = p;
        y++;
    }
    gapgroupsizet[p] = y;
    p++;
    
    y=0;
    for(k = 42; k <= 49; k++){
        gapgroupt[k] = p;
        y++;
    }
    gapgroupt[51] = p;
    y++;
    for(k = 72; k <= 76; k++){
        gapgroupt[k] = p;
        y++;
    }
    gapgroupsizet[p] = y;
    p++;
    
    y=0;
    for(k = 70; k <= 71; k++){
        gapgroupt[k] = p;
        y++;
    }
    gapgroupsizet[p] = y;
    p++;
    
    y=0;
    for(k = 119; k <= 121; k++){
        gapgroupt[k] = p;
        y++;
    }
    gapgroupsizet[p] = y;
    p++;
    
    y=0;
    for(k = 96; k <= 97; k++){
        gapgroupt[k] = p;
        y++;
    }
    gapgroupsizet[p] = y;
    p++;
    
    y=0;
    for(k = 78; k <= 79; k++){
        gapgroupt[k] = p;
        y++;
    }
    gapgroupt[82] = p;
    y++;
    gapgroupsizet[p] = y;
    p++;
    
    y=0;
    for(k = 80; k <= 81; k++){
        gapgroupt[k] = p;
        y++;
    }
    for(k = 92; k <= 94; k++){
        gapgroupt[k] = p;
        y++;
    }
    gapgroupsizet[p] = y;
    p++;
    
    y=0;
    for(k = 99; k <= 108; k++){
        gapgroupt[k] = p;
        y++;
    }
    for(k = 110; k <= 115; k++){
        gapgroupt[k] = p;
        y++;
    }
    for(k = 122; k <= 123; k++){
        gapgroupt[k] = p;
        y++;
    }
    gapgroupsizet[p] = y;
    p++;
    
    y=0;
    for(k = 124; k <= 125; k++){
        gapgroupt[k] = p;
        y++;
    }
    gapgroupsizet[p] = y;
    p++;
    
    y=0;
    for(k = 40; k <= 41; k++){
        gapgroupt[k] = p;
        y++;
    }
    gapgroupt[67] = p;
    y++;
    gapgroupt[69] = p;
    y++;
    gapgroupsizet[p] = y;
    p++;
    
    y=0;
    for(k = 34; k <= 39; k++){
        gapgroupt[k] = p;
        y++;
    }
    for(k = 57; k <= 66; k++){
        gapgroupt[k] = p;
        y++;
    }
    gapgroupt[68] = p;
    y++;
    for(k = 84; k <= 86; k++){
        gapgroupt[k] = p;
        y++;
    }
    gapgroupsizet[p] = y;
    p++;
    
    y=0;
    gapgroupt[116] = p;
    y++;
    gapgroupt[126] = p;
    y++;
    gapgroupsizet[p] = y;
    p++;
    
    y=0;
    for(k = 134; k <= 135; k++){
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
    for(k = 165; k <= 181; k++){
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
    for(k = 261; k <= 263; k++){
        gapgroupt[k] = p;
        y++;
    }
    gapgroupsizet[p] = y;
    p++;
    
    y=0;
    for(k = 249; k <= 250; k++){
        gapgroupt[k] = p;
        y++;
    }
    gapgroupsizet[p] = y;
    p++;
    
    y=0;
    for(k = 266; k <= 267; k++){
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
    for(k = 323; k <= 324; k++){
        gapgroupt[k] = p;
        y++;
    }
    gapgroupsizet[p] = y;
    p++;
    
    y=0;
    for(k = 330; k <= 333; k++){
        gapgroupt[k] = p;
        y++;
    }
    gapgroupsizet[p] = y;
    p++;
    
    y=0;
    for(k = 296; k <= 303; k++){
        gapgroupt[k] = p;
        y++;
    }
    for(k = 311; k <= 322; k++){
        gapgroupt[k] = p;
        y++;
    }
    for(k = 327; k <= 328; k++){
        gapgroupt[k] = p;
        y++;
    }
    gapgroupsizet[p] = y;
    p++;
    
    y=0;
    gapgroupt[354] = p;
    y++;
    gapgroupt[365] = p;
    y++;
    gapgroupsizet[p] = y;
    p++;
    
    y=0;
    gapgroupt[339] = p;
    y++;
    gapgroupt[346] = p;
    y++;
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
    for(k = 340; k <= 341; k++){
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
    gapgroupt[369] = p;
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
    for(k = 402; k <= 403; k++){
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
    for(k = 485; k <= 487; k++){
        gapgroupt[k] = p;
        y++;
    }
    gapgroupt[483] = p;
    y++;
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
    for(k = 437; k <= 438; k++){
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
    for(k = 466; k <= 474; k++){
        gapgroupt[k] = p;
        y++;
    }
    gapgroupsizet[p] = y;
    p++;
    
    y=0;
    for(k = 434; k <= 436; k++){
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
    for(k = 417; k <= 420; k++){
        gapgroupt[k] = p;
        y++;
    }
    gapgroupt[426] = p;
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
    for(k = 503; k <= 508; k++){
        gapgroupt[k] = p;
        y++;
    }
    for(k = 516; k <= 526; k++){
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
    gapgroupt[549] = p;
    y++;
    gapgroupsizet[p] = y;
    p++;
    
    y=0;
    for(k = 562; k <= 563; k++){
        gapgroupt[k] = p;
        y++;
    }
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
    for(k = 500; k <= 501; k++){
        gapgroupt[k] = p;
        y++;
    }
    gapgroupsizet[p] = y;
    p++;
    
    y=0;
    for(k = 512; k <= 513; k++){
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
    for(k = 527; k <= 528; k++){
        gapgroupt[k] = p;
        y++;
    }
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
    gapgroupt[573] = p;
    y++;
    gapgroupt[576] = p;
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
    for(k = 608; k <= 609; k++){
        gapgroupt[k] = p;
        y++;
    }
    gapgroupt[589] = p;
    y++;
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
    for(k = 211; k <= 212; k++){
        gapgroupt[k] = p;
        y++;
    }
    gapgroupt[200] = p;
    y++;
    gapgroupsizet[p] = y;
    p++;
    
    y=0;
    for(k = 198; k <= 199; k++){
        gapgroupt[k] = p;
        y++;
    }
    gapgroupt[201] = p;
    y++;
    gapgroupsizet[p] = y;
    p++;
    
    y=0;
    for(k = 202; k <= 203; k++){
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
    for(k = 238; k <= 240; k++){
        gapgroupt[k] = p;
        y++;
    }
    gapgroupt[236] = p;
    y++;
    gapgroupsizet[p] = y;
    p++;
    
    y=0;
    gapgroupt[235] = p;
    y++;
    gapgroupt[237] = p;
    y++;
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
    for(k = 279; k <= 280; k++){
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
    for(k = 305; k <= 306; k++){
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
    for(k = 348; k <= 349; k++){
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
    for(k = 373; k <= 374; k++){
        gapgroupt[k] = p;
        y++;
    }
    gapgroupsizet[p] = y;
    p++;
    
    
    y=0;
    for(k = 458; k <= 460; k++){
        gapgroupt[k] = p;
        y++;
    }
    gapgroupt[481] = p;
    y++;
    gapgroupsizet[p] = y;
    p++;
    
    gapgroupnumt = p;
    for( i = 1; i < realgapnum+1; i++){
        gapgrouptemp[i-1] = gapgroupt[i];
    }
    for(i = 0; i < realgapnum; i++){
        if(gapgrouptemp[i] == -1){
            gapgrouptemp[i] = gapgroupnumt;
            gapgroupsizet[gapgroupnumt] = 1;
            gapgroupnumt++;
            
        }
    }
    printf("%d/n", gapgrouptemp[630]);
    *gapgroupnum = gapgroupnumt;
    *gapgroup = gapgrouptemp;
    *gapgroupsize = gapgroupsizet;
    
}
*/
/*for(i =0; i < numbersmallgaps; i++){
    nongap_xcrd[i] += 5;
    
    for(j =0; j < numbersmallgaps; j++){
        if(nongap_xcrd[i] == nongap_xcrd[j] && gapgroup[nongap_xcrd[i]] != gapgroup[nongap_xcrd[j]]){
            nongap_ycrd[i] += 5;
        }
    }
}*/
//put extra stuff in the testzone
