//
//  stats_analysis.c
//  gapdynamics
//
//  Created by Dylan DeBaun on 3/25/20.
//  Copyright © 2020 Dylan DeBaun. All rights reserved.
//

#include "stats_analysis.h"
#include <stdio.h>
# include <stdlib.h>
# include <time.h>
# include <string.h>
# include "general.h"
#include "initializelandscape.h"
#include "time_keeper.h"
#include <math.h>

void print_landscape(char *begoutfile, double *LC_xcrd, double *LC_ycrd, int curr_year, int *species_ID, int J, int s, int *parentID){
    int i;
    char outfile[MAX_FLNAME_SIZE];
    FILE *ofp;
    
    printf("\nPrinting out landscape configuration at %d yrs...\n", curr_year);
    sprintf(outfile, "%slandscape%dat%dyrs.txt", begoutfile, s, curr_year);
    ofp = fopen(outfile, "w");
    
    for (i = 0; i < J; ++i){
        fprintf(ofp, "%f\t%f\t%d\n", LC_xcrd[i], LC_ycrd[i], parentID[i]);
    }
    fclose(ofp);
}

void print_gaps(char *begoutfile, int *gap_xcrd, int *gap_ycrd, int *gap_sizex, int *gap_sizey, int gap_number){
    int i;
    char outfile[MAX_FLNAME_SIZE];
    FILE *ofp;
    printf("\nPrinting out gap location data...\n");
    sprintf(outfile, "%sgaplocations.txt", begoutfile);
    ofp = fopen(outfile, "w");
    
    for (i = 0; i < gap_number; ++i){
        //printf("%d\t%d\n",gap_xcrd[115], gap_ycrd[115]);
        fprintf(ofp, "%d\t%d\t%d\t%d\n", gap_xcrd[i],gap_ycrd[i], gap_sizex[i], gap_sizey[i]);
    }
    fclose(ofp);
}

void print_species_key(char *begoutfile, char **species_name, int species_number){
    int i;
    char outfile[MAX_FLNAME_SIZE];
    FILE *ofp;
    
    printf("\nPrinting out species key...\n");
    sprintf(outfile, "%sspecieskey.txt", begoutfile);
    ofp = fopen(outfile, "w");
    printf("%s\n", species_name[301]);
    for (i = 0; i < species_number; ++i){
        fprintf(ofp, "%d\t%s\n", i, species_name[i]);
    }
    fclose(ofp);
}

void ascendingarray(int *samples){
    int n,i,j;
    n=301;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            if (samples[j] < samples[i])
            {
                int tmp = samples[i];
                samples[i] = samples[j];
                samples[j] = tmp;
            }
        }
    }
}

void print_SAD_per_simulation_groupgaps(char *begoutfile, int *gap_xcrd, int *gap_ycrd, int *gap_sizex, int *gap_sizey, int gap_number, int simnumber, double *LC_xcrd, double *LC_ycrd, int numindiv, int *spID, int num_species, double **speciescount, double **speciescountnongap, int curr_year, int *gapgroup, int groupgapnum, int *groupgapsize, int **nongap_xcrd1, int **nongap_ycrd1, int numindivstart){
    int i,z,j,k,t,q;
    char outfile[MAX_FLNAME_SIZE];
    FILE *ofp;
    int *nongap_xcrd,*nongap_ycrd;
    nongap_xcrd = calloc(1000, sizeof(int));
    for(i=0;i<1000; i++){
        nongap_xcrd[i]=0;
    }
    nongap_ycrd = calloc(1000, sizeof(int));
    for(i=0;i<1000; i++){
        nongap_ycrd[i]=0;
    }
    /*int **speciescountnongap, **speciescountrecruit,**speciescountnongaprecruits;
    speciescountnongap = calloc(500, sizeof(int *));
    for(i=0;i<500; i++){
        speciescountnongap[i] = calloc(1000, sizeof(int ));
        for (j = 0; j<1000; ++j) {
            speciescountnongap[i][j] = 0;
        }
    }
    speciescountnongaprecruits = calloc(500, sizeof(int *));
    for(i=0;i<500; i++){
        speciescountnongaprecruits[i] = calloc(1000, sizeof(int ));
        for (j = 0; j<1000; ++j) {
            speciescountnongaprecruits[i][j] = 0;
        }
    }
    speciescountrecruit = calloc(500, sizeof(int *));
    for(i=0;i<500; i++){
        speciescountrecruit[i] = calloc(1000, sizeof(int ));
        for (j = 0; j<1000; ++j) {
            speciescountrecruit[i][j] = 0;
        }
    }*/
    int ingap,p,n,r,a,stop;
    int op= 0;
    double b,c,d,e,f1,g1,h1,i1,j1,k1,l1;;
    double x,y, gx,gy,gxsize,gysize;
    n= numindiv;
    //for a grouped gap
    int count = 0, count1= 0;
    for(j=0; j<n; j++){
        x = LC_xcrd[j];
        y = LC_ycrd[j];
        if(x!= -1){
            count++;
        }
    }
    char *type;
    for(a =0;a<3;a++){
        if(a==0){
            type = "large";
        }else if(a==1){
            type= "small";
        }else{
            type="all";
        }
        for(k=0; k< 1000; k++){
            speciescount[0][k] = 0;
            for(i =0;i <10;i++){
            speciescountnongap[i][k] = 0;
            }
        }
        for(t = 0; t < groupgapnum; t++){
            stop = 0;
            if(a==0){
                if(groupgapsize[t]<3){
                    stop = 1;
                }else{
                    op++;
                }
            }else if(a==1){
                if(groupgapsize[t]>2){
                    stop = 1;
                }else{
                    op++;
                }
            }
            //find out how many real gaps are in that grouped gap
            if(stop == 0){
                for(i = 0; i < gap_number; i++){
                    ingap = gapgroup[i];
                    if(gapgroup[i] == t){
                        for(j=0; j<n; j++){
                            x = LC_xcrd[j];
                            y = LC_ycrd[j];
                            gx = gap_xcrd[i];
                            gy = gap_ycrd[i];
                            gxsize =gap_sizex[i];
                            gysize =gap_sizey[i];
                            if(x!=-1){
                                if(x > gx && y > gy && y <= gy+gysize && x <= gx+gxsize){
                                    // printf("%d ", speciescount[i][spID[j]]);
                                    p = spID[j];
                                    speciescount[0][spID[j]]++;
                                    r=speciescount[0][spID[j]];
                                    count1++;
                                }
                            }
                        }
                    }
                    //}
                }
                //FOR NON GAPS
                //for(i = 0; i < gap_number; i++){
                //for(q = 0; q < 9; i++){
                //  ingap = gapgroup[i];
                //if(gapgroup[i] == t){
                for(i =0;i<10;i++){
                    nongap_xcrd[t] = nongap_xcrd1[i][t];
                    nongap_ycrd[t] = nongap_ycrd1[i][t];
                for(j=0; j<n; j++){
                    x = LC_xcrd[j];
                    y = LC_ycrd[j];
                    gx = nongap_xcrd[t];
                    gy = nongap_ycrd[t];
                    gxsize =sqrt(groupgapsize[t]*25);
                    gysize =sqrt(groupgapsize[t]*25);
                    if(LC_xcrd[j] != -1 && LC_xcrd[j] > nongap_xcrd[t] && LC_ycrd[j] > nongap_ycrd[t] && LC_ycrd[j] <= (nongap_ycrd[t]+gysize) && LC_xcrd[j] <= (nongap_xcrd[t]+gxsize)){
                        // printf("%d ", speciescount[i][spID[j]]);
                        p = spID[j];
                       // speciescountnongap[0][spID[j]]++;
                        speciescountnongap[i][spID[j]]++;
                        r=speciescountnongap[0][spID[j]];
                    }
                        //}
                    // }
                    //}
                }
                }
                //FOR RECRUITS
                for(i = 0; i < gap_number; i++){
                    ingap = gapgroup[i];
                    if(gapgroup[i] == t){
                        for(j=numindivstart; j<n; j++){
                            x = LC_xcrd[j];
                            y = LC_ycrd[j];
                            gx = gap_xcrd[i];
                            gy = gap_ycrd[i];
                            gxsize =gap_sizex[i];
                            gysize =gap_sizey[i];
                            if(LC_xcrd[j] != -1 && LC_xcrd[j] > gap_xcrd[i] && LC_ycrd[j] > gap_ycrd[i] && LC_ycrd[j] <= (gap_ycrd[i]+gap_sizey[i]) && LC_xcrd[j] <= (gap_xcrd[i]+gap_sizex[i])){
                                // printf("%d ", speciescount[i][spID[j]]);
                                p = spID[j];
                            }
                        }
                    }
                    //}
                }
                //FOR NON GAPS recruits
                for(i =0;i<1;i++){
                    nongap_xcrd[t] = nongap_xcrd1[i][t];
                    nongap_ycrd[t] = nongap_ycrd1[i][t];
                for(j=numindivstart; j<n; j++){
                    x = LC_xcrd[j];
                    y = LC_ycrd[j];
                    gx = nongap_xcrd[t];
                    gy = nongap_ycrd[t];
                    gxsize =sqrt(groupgapsize[t]*25);
                    gysize =sqrt(groupgapsize[t]*25);
                    if(LC_xcrd[j] != -1 && LC_xcrd[j] > nongap_xcrd[t] && LC_ycrd[j] > nongap_ycrd[t] && LC_ycrd[j] <= (nongap_ycrd[t]+gysize) && LC_xcrd[j] <= (nongap_xcrd[t]+gxsize)){
                        // printf("%d ", speciescount[i][spID[j]]);
                        p = spID[j];
                    }
                    //}
                }
                }

            }
        }
        
        //Average 10 nongap sets and send to int (floor) for downstream analysis
       /* for(i =0; i < num_species; i++){
            printf("%f\n",speciescountnongap[0][i]);
            speciescountnongaprecruits[0][i] /= 10;
            speciescountnongap[0][i] /= 10;
            printf("%f\n",speciescountnongap[0][i]);
            //speciescountnongap[0][i]= floor(speciescountnongap[0][i]);
            //printf("%d\n",speciescountnongap[0][i]);
        }*/
            
                //upon getting all q gaps in the group gap tallied, write the group gap to a file
                sprintf(outfile, "%ssimulation%d%sspeciescountyr%d.txt", begoutfile,simnumber,type,curr_year);
                ofp = fopen(outfile, "w");
                for(z=0;z<num_species; z++){
                    e=speciescount[0][z];
                    b=speciescountnongap[0][z];
                    //c=speciescountrecruit[0][z];
                    //d=speciescountnongaprecruits[0][z];
                    c=speciescountnongap[1][z];
                    d=speciescountnongap[2][z];
                    f1=speciescountnongap[3][z];
                    g1=speciescountnongap[4][z];
                    h1=speciescountnongap[5][z];
                    i1=speciescountnongap[6][z];
                    j1=speciescountnongap[7][z];
                    k1=speciescountnongap[8][z];
                    l1=speciescountnongap[9][z];
                    fprintf(ofp, "%d\t%f\t%d\t%f\t%d\t%f\t%d\t%f\t%d\t%f\t%d\t%f\t%d\t%f\t%d\t%f\t%d\t%f\t%d\t%f\t%d\t%f\n", z, e, z, b, z, c, z, d, z, f1, z, g1, z, h1, z, i1,z, j1, z, k1, z, l1);
                }
                fclose(ofp);
    }
    /*printf("%d\n",count1);
    char outfile2[MAX_FLNAME_SIZE];
    FILE *ofp2;
    sprintf(outfile2, "%ssimulation%dfulllandscapeyr%d.txt", begoutfile,simnumber,curr_year);
    count =0;
    ofp2 = fopen(outfile2, "w");
        for(j=0; j<n; j++){
            x = LC_xcrd[j];
            y = LC_ycrd[j];
            if(x!= -1){
                count++;
            fprintf(ofp2, "%f\t%f\t%d\n",x,y, spID[j]);
            }
        }
    printf("%d\n",count);
    fclose(ofp2);*/
        
}




