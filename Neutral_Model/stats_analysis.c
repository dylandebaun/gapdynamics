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

void print_SAD_in_gaps(char *begoutfile, int *gap_xcrd, int *gap_ycrd, int *gap_sizex, int *gap_sizey, int gap_number, int simnumber, double *LC_xcrd, double *LC_ycrd, int numindiv, int *spID, int num_species, int ***speciescount, int curr_year){
    int i,z;
    char outfile[MAX_FLNAME_SIZE];
    FILE *ofp;
    double spcount;
    
    for(i = 0; i<gap_number; i++){
        sprintf(outfile, "%sgap%dspeciescountyr%d.txt", begoutfile, i,curr_year);
        ofp = fopen(outfile, "w");
        
        for(z=0;z<num_species; z++){
            if(speciescount[curr_year][i][z] > 0){
               // printf("%d\n",speciescount[curr_year][i][z]);
            }
            spcount = (double) (speciescount[curr_year][i][z])/(simnumber);
            //printf("%f\n",speciescount[curr_year][i][z]/simnumber);
            fprintf(ofp, "%d\t%f\n", z, spcount);
        }
        //should print rank abundance as well
        
        fclose(ofp);
        //printf("\n");
    }
    printf("\nDone printing out gap individuals...\n");
}

 
void calc_SAD_in_gaps(char *begoutfile, int *gap_xcrd, int *gap_ycrd, int *gap_sizex, int *gap_sizey, int gap_number, double *LC_xcrd, double *LC_ycrd, int numindiv, int *spID, int num_species, int ***speciescount, int curr_year){
        
    int i,j;

    for(i = 0; i<gap_number; i++){
        for(j=0; j<numindiv; j++){
            if(LC_xcrd[j] > gap_xcrd[i] && LC_ycrd[j] > gap_ycrd[i] && LC_ycrd[j] <= (gap_ycrd[i]+gap_sizey[i]) && LC_xcrd[j] <= (gap_xcrd[i]+gap_sizex[i])){
               // printf("%d ", speciescount[i][spID[j]]);
                speciescount[curr_year][i][spID[j]]++;
                //if(speciescount[curr_year][i][spID[j]] > 0){
                
               // }
               // printf("%d\n", speciescount[i][spID[j]]);
            }
        }
    }
    
}
void print_SAD_per_simulation(char *begoutfile, int *gap_xcrd, int *gap_ycrd, int *gap_sizex, int *gap_sizey, int gap_number, int simnumber, double *LC_xcrd, double *LC_ycrd, int numindiv, int *spID, int num_species, int *speciescount, int curr_year){
    int i,z,j,k,a;
    char outfile[MAX_FLNAME_SIZE];
    FILE *ofp;
    
    for(a =0;a <3;a++){
    for(i = 0; i<gap_number; i++){
        for(k=0; k< 1000; k++){
            speciescount[k] = 0;
        }
        for(j=0; j<numindiv; j++){
            if(LC_xcrd[j] > gap_xcrd[i] && LC_ycrd[j] > gap_ycrd[i] && LC_ycrd[j] <= (gap_ycrd[i]+gap_sizey[i]) && LC_xcrd[j] <= (gap_xcrd[i]+gap_sizex[i])){
                // printf("%d ", speciescount[i][spID[j]]);
                speciescount[spID[j]]++;
                //if(speciescount[curr_year][i][spID[j]] > 0){
                
                // }
                // printf("%d\n", speciescount[i][spID[j]]);
            }
        }
        sprintf(outfile, "%ssimulation%dgap%dspeciescountyr%d.txt", begoutfile,simnumber, i,curr_year);
        ofp = fopen(outfile, "w");
     
        for(z=0;z<num_species; z++){
            fprintf(ofp, "%d\t%d\n", z, speciescount[z]);
        }
        //should print rank abundance as well
        
        fclose(ofp);
        //printf("\n");
    }
    }
    printf("\nDone printing out gap individuals...\n");
}

void printgroupgapinfo(char *begoutfile, int *groupgap, int *groupgapsize, int groupgapnum, int numgaps){
    int i;
    char outfile[MAX_FLNAME_SIZE];
    FILE *ofp;
    
    printf("\nPrinting out gap location data...\n");
    sprintf(outfile, "%sgaplocations.txt", begoutfile);
    ofp = fopen(outfile, "w");
    
   /* for (i = 0; i < groupgapnum; ++i){
        fprintf(ofp, "%d\t%d\n", i, groupgapsize[i]);
    }*/
    for(i =0; i < numgaps; i++){
        fprintf(ofp, "%d\t%d\t%d\n", i, groupgap[i], groupgapsize[groupgap[i]]);
    }
    fclose(ofp);
    
    sprintf(outfile, "groupgaps.txt");
    ofp = fopen(outfile, "w");
    
    /* for (i = 0; i < groupgapnum; ++i){
     fprintf(ofp, "%d\t%d\n", i, groupgapsize[i]);
     }*/
    for(i =0; i < groupgapnum; i++){
        fprintf(ofp, "%d\t%d\n", i, groupgapsize[i]);
    }
    fclose(ofp);
}

void print_SAD_per_simulation_groupgaps(char *begoutfile, int *gap_xcrd, int *gap_ycrd, int *gap_sizex, int *gap_sizey, int gap_number, int simnumber, double *LC_xcrd, double *LC_ycrd, int numindiv, int *spID, int num_species, double **speciescount, double **speciescountnongap, double **speciescountnongaprecruits, double **speciescountrecruit, int curr_year, int *gapgroup, int groupgapnum, int *groupgapsize, int **nongap_xcrd1, int **nongap_ycrd1, int numindivstart){
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
            speciescountrecruit[i][k] = 0;
            }
            speciescountnongaprecruits[0][k] = 0;
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
                                speciescountrecruit[0][spID[j]]++;
                                r=speciescountrecruit[0][spID[j]];
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
                        speciescountnongaprecruits[0][spID[j]]++;
                        r=speciescountnongaprecruits[0][spID[j]];
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

void print_SAD_per_simulation_nongapmeta(char *begoutfile, int *gap_xcrd, int *gap_ycrd, int *gap_sizex, int *gap_sizey, int gap_number, int simnumber, double *LC_xcrd, double *LC_ycrd, int numindiv, int *spID, int num_species, int curr_year, int numindivstart, int *speciescountmeta, int *speciescountmetarecruits){
    int i,z,j,k;
    char outfile[MAX_FLNAME_SIZE];
    FILE *ofp;
    int *speciescount,*speciescountrecruit;
    int n = numindiv;
    //for a grouped gap
    for(k=0; k< 1000; k++){
        speciescountmeta[k] = 0;
        speciescountmetarecruits[k] = 0;
    }
    int stop = 0,count =0;
    //scan for individuals not in any gap
    for(j=0; j<n; j++){
        stop = 0;
        for(i = 0; i < gap_number; i++){
            if(LC_xcrd[j] != -1 && LC_xcrd[j] > gap_xcrd[i] && LC_ycrd[j] > gap_ycrd[i] && LC_ycrd[j] <= (gap_ycrd[i]+gap_sizey[i]) && LC_xcrd[j] <= (gap_xcrd[i]+gap_sizex[i])){
                stop = 1;
                break;
            }
        }
        if(stop == 0 && LC_xcrd[j] != -1){
            count++;
            speciescountmeta[spID[j]]++;
            if(j >= numindivstart){
                speciescountmetarecruits[spID[j]]++;
            }
        }
    }
    printf("%d\n",count);
        
        //upon getting all q gaps in the group gap tallied, write the group gap to a file
        sprintf(outfile, "%ssimulation%dmetanongapspeciescountyr%d.txt", begoutfile,simnumber,curr_year);
        ofp = fopen(outfile, "w");
        for(z=0;z<num_species; z++){
            fprintf(ofp, "%d\t%d\t%d\t%d\n", z, speciescountmeta[z], z, speciescountmetarecruits[z]);
        }
        fclose(ofp);
}

void print_stats(char *begoutfile, int *gap_xcrd, int *gap_ycrd, int *gap_sizex, int *gap_sizey, int gap_number, int simnumber, double *LC_xcrd, double *LC_ycrd, int numindiv, int *spID, int num_species, double **speciescount, double **speciescountnongap, double **speciescountnongaprecruits, double **speciescountrecruit, int *gapgroup, int groupgapnum, int *groupgapsize, int **nongap_xcrd1, int **nongap_ycrd1, int numindivstart, int gen, double ***largearray, double ***smallarray, double ***allarray, double *wd, double *factor, int size, double **array, int mod){
    
    int i,z,j,k,t,q,a;
    char outfile[MAX_FLNAME_SIZE];
    FILE *ofp;
    char *metacommunity;
    int ingap,p,n,r,b,c,d;
    double x,y, gx,gy,gxsize,gysize;
    n= numindiv;
    int richgaps = 0,richnongaps = 0,richgapsrec = 0,richnongapsrec = 0;
    int densgaps = 0,densnongaps = 0,densgapsrec = 0,densnongapsrec = 0;
    double Hgaps = 0,Hnongaps = 0,Hgapsrec = 0,Hnongapsrec = 0;
    double Egaps = 0,Enongaps = 0,Egapsrec = 0,Enongapsrec = 0;
    int sumofspgaps = 0,sumofspnongaps = 0,sumofspgapsrec = 0,sumofspnongapsrec = 0;
    double avgwdgaps = 0,avgwdnongaps = 0,avgwdgapsrec = 0,avgwdnongapsrec = 0;
    double avgfacgaps = 0,avgfacnongaps = 0,avgfacgapsrec = 0,avgfacnongapsrec = 0;
    int cwdgaps = 0,cwdnongaps = 0,cwdgapsrec = 0,cwdnongapsrec = 0;
    int cfacgaps = 0,cfacnongaps = 0,cfacgapsrec = 0,cfacnongapsrec = 0;
    int *nongap_xcrd,*nongap_ycrd;
    int doit;
    nongap_xcrd = calloc(1000, sizeof(int));
    for(i=0;i<1000; i++){
        nongap_xcrd[i]=0;
    }
    nongap_ycrd = calloc(1000, sizeof(int));
    for(i=0;i<1000; i++){
        nongap_ycrd[i]=0;
    }
    //check
   /* for(i = 0; i < 400; i++){
        printf("%f %f \n",wd[i],factor[i]);
    }*/
    /*
    for(i =0; i <20;i++){
        printf("%lf\n", array[gen][i]);
    }*/
    int count =0,sizer;
    for(a=size; a<size+1;a++){
        for(k=0; k< 1000; k++){
            speciescount[0][k] = 0;
            speciescountnongap[0][k] = 0;
            speciescountrecruit[0][k] = 0;
            speciescountnongaprecruits[0][k] = 0;
        }
        if(a == 1){//LARGE
            array = *largearray;
        }else if(a == 2){
            array = *smallarray;
        }else{
            array = *allarray;
        }
        
        
        for(t = 0; t < groupgapnum; t++){
            doit = 1;
            if(a == 1){//LARGE
                if(groupgapsize[t] < 3){
                    doit = 0;
                }
            }else if(a == 2){
                if(groupgapsize[t] > 2){
                    doit = 0;
                }
            }
            if(doit != 0){
            for(i = 0; i < gap_number; i++){
                ingap = gapgroup[i];
                if(gapgroup[i] == t){
                    for(j=0; j<numindiv; j++){
                        x = LC_xcrd[j];
                        y = LC_ycrd[j];
                        gx = gap_xcrd[i];
                        gy = gap_ycrd[i];
                        gxsize =gap_sizex[i];
                        gysize =gap_sizey[i];
                        if(gxsize != 5 || gysize != 5){
                            printf("error gapsize\n");
                        }
                        if(x!=-1){
                          //  if(x >= gx && y >= gy && y <= (gy+gysize) && x <= (gx+gxsize)){
                            if(LC_xcrd[j] != -1 && LC_xcrd[j] > gap_xcrd[i] && LC_ycrd[j] > gap_ycrd[i] && LC_ycrd[j] <= (gap_ycrd[i]+5) && LC_xcrd[j] <= (gap_xcrd[i]+5)){
                                p = spID[j];
                                speciescount[0][spID[j]]++;
                                r=speciescount[0][spID[j]];
                                count++;
                            }
                        }
                    }
                }
            }
            //FOR NON GAPS
            /*for(i = 0; i < gap_number; i++){
                ingap = gapgroup[i];
                if(gapgroup[i] == t){
                    for(j=0; j<n; j++){
                        x = LC_xcrd[j];
                        y = LC_ycrd[j];
                        gx = gap_xcrd[i];
                        gy = gap_ycrd[i];
                        gxsize =gap_sizex[i];
                        gysize =gap_sizey[i];
                        if(LC_xcrd[j] != -1 && LC_xcrd[j] > nongap_xcrd[i] && LC_ycrd[j] > nongap_ycrd[i] && LC_ycrd[j] <= (nongap_ycrd[i]+gap_sizey[i]) && LC_xcrd[j] <= (nongap_xcrd[i]+gap_sizex[i])){
                            // printf("%d ", speciescount[i][spID[j]]);
                            p = spID[j];
                            speciescountnongap[0][spID[j]]++;
                            r=speciescountnongap[0][spID[j]];
                            //if(speciescount[curr_year][i][spID[j]] > 0){
                            //printf("%d,%d\n",p,r);
                            
                            // }
                            // printf("%d\n", speciescount[i][spID[j]]);
                        }
                    }
                }
                //}
            }*/
                
            //get species counts in all nongaps
            for(i =0;i<1;i++){
                nongap_xcrd[t] = nongap_xcrd1[i][t];
                nongap_ycrd[t] = nongap_ycrd1[i][t];
            for(j=0; j<n; j++){
                x = LC_xcrd[j];
                y = LC_ycrd[j];
                gx = nongap_xcrd[t];
                gy = nongap_ycrd[t];
                sizer = groupgapsize[t];
                gxsize =(sqrt(groupgapsize[t]*25));
                gysize =(sqrt(groupgapsize[t]*25));
                if(LC_xcrd[j] != -1 && LC_xcrd[j] > nongap_xcrd[t] && LC_ycrd[j] > nongap_ycrd[t] && LC_ycrd[j] <= (nongap_ycrd[t]+gxsize) && LC_xcrd[j] <= (nongap_xcrd[t]+gysize)){
                    // printf("%d ", speciescount[i][spID[j]]);
                    p = spID[j];
                    speciescountnongap[0][spID[j]]++;
                    r=speciescountnongap[0][spID[j]];
                }
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
                            speciescountrecruit[0][spID[j]]++;
                            r=speciescountrecruit[0][spID[j]];
                            //if(speciescount[curr_year][i][spID[j]] > 0){
                            //printf("%d,%d\n",p,r);
                            // }
                            // printf("%d\n", speciescount[i][spID[j]]);
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
                gxsize =(sqrt(groupgapsize[t]*25));
                gysize =(sqrt(groupgapsize[t]*25));
                if(LC_xcrd[j] != -1 && LC_xcrd[j] > nongap_xcrd[t] && LC_ycrd[j] > nongap_ycrd[t] && LC_ycrd[j] <= (nongap_ycrd[t]+gxsize) && LC_xcrd[j] <= (nongap_xcrd[t]+gysize)){
                            // printf("%d ", speciescount[i][spID[j]]);
                            p = spID[j];
                            speciescountnongaprecruits[0][spID[j]]++;
                            r=speciescountnongaprecruits[0][spID[j]];
                            //printf("%d,%d\n",p,r);
                            //if(speciescount[curr_year][i][spID[j]] > 0){
                            
                            // }
                            // printf("%d\n", speciescount[i][spID[j]]);
                        }
            }
                }
            
            }
        }
        //Average 10 nongap sets
       /* for(i =0; i < num_species; i++){
            speciescountnongaprecruits[0][i] /= 10;
            speciescountnongap[0][i] /= 10;
           // printf("%d\n",speciescountnongaprecruits[0][i]);
        }*/
        //2. calculate H,richness
        richgaps = 0;
        richnongaps = 0;
        richgapsrec = 0;
        richnongapsrec = 0;
        densgaps = 0;
       densnongaps = 0;
        densgapsrec = 0;
        densnongapsrec = 0;
        Hgaps = 0;
        Hnongaps = 0;
        Hgapsrec = 0;
        Hnongapsrec = 0;
        Egaps = 0;
        Enongaps = 0;
        Egapsrec = 0;
        Enongapsrec = 0;
        sumofspgaps = 0;
        sumofspnongaps = 0;
        sumofspgapsrec = 0;
        sumofspnongapsrec = 0;
        avgwdgaps = 0;
        avgwdnongaps = 0;
        avgwdgapsrec = 0;
        avgwdnongapsrec = 0;
        avgfacgaps = 0;
        avgfacnongaps = 0;
        avgfacgapsrec = 0;
        avgfacnongapsrec = 0;
        cwdgaps = 0;
        cwdnongaps = 0;
        cwdgapsrec = 0;
        cwdnongapsrec = 0;
        cfacgaps = 0;
        cfacnongaps = 0;
        cfacgapsrec = 0;
        cfacnongapsrec = 0;
        x=0;

        for(k=0; k < num_species; k++){
            sumofspgaps+=speciescount[0][k];
            sumofspnongaps+=speciescountnongap[0][k];
            sumofspgapsrec +=speciescountrecruit[0][k];
            sumofspnongapsrec+=speciescountnongaprecruits[0][k];
        }
        for(k=0; k < num_species; k++){
            if(speciescount[0][k] > 0.5){
                richgaps++;
                b =speciescount[0][k];
                x =(double) b/sumofspgaps;
                Hgaps += x*log(x);
                densgaps += speciescount[0][k] ;
            }
            if(speciescountnongap[0][k] > 0.5){
                richnongaps++;
                x =(double) speciescountnongap[0][k]/sumofspnongaps;
                Hnongaps += x*log(x);
                densnongaps += speciescountnongap[0][k] ;
            }
            if(speciescountrecruit[0][k] > 0.5){
                richgapsrec++;
                x =(double) speciescountrecruit[0][k]/sumofspgapsrec;
                Hgapsrec += (double) x*log(x);
                densgapsrec += speciescountrecruit[0][k] ;
            }
            if(speciescountnongaprecruits[0][k] > 0.5){
                richnongapsrec++;
                x = (double) speciescountnongaprecruits[0][k]/sumofspnongapsrec;
                Hnongapsrec += (double) x*log(x);
                densnongapsrec += speciescountnongaprecruits[0][k] ;
            }
            if(wd[k] != -10){
                if(speciescount[0][k] != 0){
                    b =speciescount[0][k];
                    c = wd[k];
                    avgwdgaps += (speciescount[0][k] * wd[k]);
                    cwdgaps += speciescount[0][k];
                }
                if(speciescountnongap[0][k] != 0){
                    avgwdnongaps += (speciescountnongap[0][k] * wd[k]);
                    cwdnongaps += speciescountnongap[0][k];
                }
                if(speciescountrecruit[0][k] != 0){
                avgwdgapsrec += (speciescountrecruit[0][k] * wd[k]);
                cwdgapsrec += speciescountrecruit[0][k];
                }
                if(speciescountnongaprecruits[0][k] != 0){
                avgwdnongapsrec += (speciescountnongaprecruits[0][k] * wd[k]);
                cwdnongapsrec += speciescountnongaprecruits[0][k];
                }
            }
            if(factor[k] != -10){
                if(speciescount[0][k] != 0){
                avgfacgaps += (speciescount[0][k] * factor[k]);
                cfacgaps += speciescount[0][k];
                }
                 if(speciescountnongap[0][k] != 0){
                avgfacnongaps += (speciescountnongap[0][k] * factor[k]);
                cfacnongaps += speciescountnongap[0][k];
                 }
                if(speciescountrecruit[0][k] != 0){
                avgfacgapsrec += (speciescountrecruit[0][k] * factor[k]);
                cfacgapsrec += speciescountrecruit[0][k];
                }
                if(speciescountnongaprecruits[0][k] != 0){
                avgfacnongapsrec += (speciescountnongaprecruits[0][k] * factor[k]);
                cfacnongapsrec += speciescountnongaprecruits[0][k];
                }
            }
        }
        Hgaps *= -1;
        Hnongapsrec *= -1;
        Hgapsrec *= -1;
        Hnongaps *= -1;
        
        //3. calculate evenness
        Egaps = Hgaps/(log(sumofspgaps));
        Enongaps = Hnongaps/(log(sumofspnongaps));
        Egapsrec = Hgapsrec/(log(sumofspgapsrec));
        Enongapsrec = Hnongapsrec/(log(sumofspnongapsrec));
        //4. calculate trait values
        avgfacgaps /= cfacgaps;
        avgfacnongaps /= cfacnongaps;
        if(cfacgapsrec != 0){avgfacgapsrec /= cfacgapsrec;}
        if(cfacgapsrec != 0){avgfacnongapsrec /= cfacnongapsrec;}
        
        avgwdgaps /= cwdgaps;
        avgwdnongaps /= cwdnongaps;
        if(cwdgapsrec != 0){avgwdgapsrec /= cwdgapsrec;}
        if(cwdnongapsrec != 0){avgwdnongapsrec /= cwdnongapsrec;}
        
        z =0;
        array[gen/mod][z] = richgaps;
        array[gen/mod][z+1] = Egaps;
        array[gen/mod][z+2] = Hgaps;
        array[gen/mod][z+3] = avgwdgaps;
        array[gen/mod][z+4] = avgfacgaps;
        array[gen/mod][z+5] = densgaps;
        z=6;
        array[gen/mod][z] = richnongaps;
        array[gen/mod][z+1] = Enongaps;
        array[gen/mod][z+2] = Hnongaps;
        array[gen/mod][z+3] = avgwdnongaps;
        array[gen/mod][z+4] = avgfacnongaps;
        array[gen/mod][z+5] = densnongaps;
        z=12;
        array[gen/mod][z] = richgapsrec;
        array[gen/mod][z+1] = Egapsrec;
        array[gen/mod][z+2] = Hgapsrec;
        array[gen/mod][z+3] = avgwdgapsrec;
        array[gen/mod][z+4] = avgfacgapsrec;
        array[gen/mod][z+5] = densgapsrec;
        z=18;
        array[gen/mod][z] = richnongapsrec;
        array[gen/mod][z+1] = Enongapsrec;
        array[gen/mod][z+2] = Hnongapsrec;
        array[gen/mod][z+3] = avgwdnongapsrec;
        array[gen/mod][z+4] = avgfacnongapsrec;
        array[gen/mod][z+5] = densnongapsrec;
        
        
        if(a == 1){
            //free(*largearray);
            *largearray = array;
        }else if(a==2){
            *smallarray = array;
            
        }else{
            *allarray = array;
        }
    }
}

void ingapornongap(double LC_xcrdchoice, double LC_ycrdchoice, int parent_spp_IDchoice, int *gap_xcrd, int *gap_ycrd, int gap_number, int **speciescountgap, int *nongap_xcrd, int *nongap_ycrd, int groupgapnum, int *groupgapsize, int **speciescountnongap){
    
    int i,p,r,t,sizer;
    double x,y,gx,gy;
    int b = 0;
    int *array;
    array = calloc(400, sizeof(int *));
    for(i=0;i<400; i++){
        array[i] = 0;
    }
    
    x = LC_xcrdchoice;
    y = LC_ycrdchoice;
    p = parent_spp_IDchoice;
    for(i = 0; i < gap_number; i++){
        gx = gap_xcrd[i];
        gy = gap_ycrd[i];
        if(x > gap_xcrd[i] && y > gap_ycrd[i] && y <= (gap_ycrd[i]+5) && x <= (gap_xcrd[i]+5)){
            b = 1;
            array = *speciescountgap;
            array[p]++;
            r=array[p];
            break;
        }
    }
    double gxsize,gysize;
    if(b != 1){
    for(t = 0; t < groupgapnum; t++){
    //get species counts in all nongaps
            gx = nongap_xcrd[t];
            gy = nongap_ycrd[t];
            sizer = groupgapsize[t];
            gxsize =(sqrt(groupgapsize[t]*25));
            gysize =(sqrt(groupgapsize[t]*25));
            if(x > nongap_xcrd[t] && y > nongap_ycrd[t] && y <= (nongap_ycrd[t]+gxsize) && x <= (nongap_xcrd[t]+gysize)){
                b =2;
                array = *speciescountnongap;
                array[p]++;
                r=array[p];
                break;
            }
    }
    }
    if(b == 1){
        *speciescountgap= array;
    }else if(b==2){
        *speciescountnongap= array;
    }
    //r=*speciescountnongap[p];
}

void print_statsfull(char *begoutfile, int simnumber, double *LC_xcrd, double *LC_ycrd, int numindiv, int *spID, int num_species, int **speciescount, int gen, double ***allarray, double *wd, double *factor, double **array){
    
    int i,z,j,k,t,q,a;
    char outfile[MAX_FLNAME_SIZE];
    FILE *ofp;
    char *metacommunity;
    int ingap,p,n,r,b,c,d;
    double x,y, gx,gy,gxsize,gysize;
    n= numindiv;
    int richgaps = 0,densgaps=0;
    double Hgaps = 0;
    double Egaps = 0;
    int sumofspgaps = 0;
    double avgwdgaps = 0;
    double avgfacgaps = 0;
    int cwdgaps = 0;
    int cfacgaps = 0;
    
    int count =0;
        for(k=0; k< 1000; k++){
            speciescount[0][k] = 0;
        }
        array = *allarray;

        for(j=0; j<numindiv; j++){
            x = LC_xcrd[j];
            y = LC_ycrd[j];
            if(x!=-1){
                p = spID[j];
                speciescount[0][spID[j]]++;
                r=speciescount[0][spID[j]];
                count++;
            }
        }
    
        
        //2. calculate H,richness
        richgaps = 0;
        Hgaps = 0;
        Egaps = 0;
        sumofspgaps = 0;
        avgwdgaps = 0;
        avgfacgaps = 0;
        cwdgaps = 0;
        cfacgaps = 0;
        x=0;
        densgaps=0;
        
        for(k=0; k < num_species; k++){
            sumofspgaps+=speciescount[0][k];
        }
        for(k=0; k < num_species; k++){
            if(speciescount[0][k] != 0){
                richgaps++;
                b =speciescount[0][k];
                x =(double) b/sumofspgaps;
                Hgaps += x*log(x);
                densgaps += b;
            }
            if(wd[k] != -10){
                if(speciescount[0][k] != 0){
                    b =speciescount[0][k];
                    c = wd[k];
                    avgwdgaps += (speciescount[0][k] * wd[k]);
                    cwdgaps += speciescount[0][k];
                }
            }
            if(factor[k] != -10){
                if(speciescount[0][k] != 0){
                    avgfacgaps += (speciescount[0][k] * factor[k]);
                    cfacgaps += speciescount[0][k];
                }
            }
        }
        Hgaps *= -1;
        
        //3. calculate evenness
        Egaps = Hgaps/(log(sumofspgaps));

        //4. calculate trait values
        avgfacgaps /= cfacgaps;
        
        avgwdgaps /= cwdgaps;
    
        z =0;
        array[gen/100][z] = richgaps;
        array[gen/100][z+1] = Egaps;
        array[gen/100][z+2] = Hgaps;
        array[gen/100][z+3] = avgwdgaps;
        array[gen/100][z+4] = avgfacgaps;
        array[gen/100][z+5] = densgaps;
        

        *allarray = array;

   /* if(simnumber == 0){
        sprintf(outfile, "%ssimulation%dtimestamp%d.txt", begoutfile, simnumber, gen/10);
        ofp = fopen(outfile, "w");
        //upon getting all q gaps in the group gap tallied, write the group gap to a file
        for(i =0; i <numindiv; i++){
            if(LC_xcrd[i]!= -1){
                x = LC_xcrd[i];
                y = LC_ycrd[i];
            //for(z=0;z<num_species; z++){
                fprintf(ofp, "%lf\t%lf\n", x,y);
            //}
            }
        }
        fclose(ofp);
    }*/
}

/*
void print_statsactual(char *begoutfile, int *gap_xcrd, int *gap_ycrd, int *gap_sizex, int *gap_sizey, int gap_number, int simnumber, double *LC_xcrd, double *LC_ycrd, int numindiv, int *spID, int num_species, int **speciescount, int **speciescountnongap, int **speciescountnongaprecruits, int **speciescountrecruit, int *gapgroup, int groupgapnum, int *groupgapsize, int *nongap_xcrd, int *nongap_ycrd, int numindivstart, int gen, double ***largearray, double ***smallarray, double ***allarray, double *wd, double *factor, int size, double **array){
    
    int i,z,j,k,t,q,a;
    int ingap,p,n,r,b,c,d;
    double x,y, gx,gy,gxsize,gysize;
    n= numindiv;
    int richgaps = 0,richnongaps = 0,richgapsrec = 0,richnongapsrec = 0;
    double Hgaps = 0,Hnongaps = 0,Hgapsrec = 0,Hnongapsrec = 0;
    double Egaps = 0,Enongaps = 0,Egapsrec = 0,Enongapsrec = 0;
    int sumofspgaps = 0,sumofspnongaps = 0,sumofspgapsrec = 0,sumofspnongapsrec = 0;
    double avgwdgaps = 0,avgwdnongaps = 0,avgwdgapsrec = 0,avgwdnongapsrec = 0;
    double avgfacgaps = 0,avgfacnongaps = 0,avgfacgapsrec = 0,avgfacnongapsrec = 0;
    int cwdgaps = 0,cwdnongaps = 0,cwdgapsrec = 0,cwdnongapsrec = 0;
    int cfacgaps = 0,cfacnongaps = 0,cfacgapsrec = 0,cfacnongapsrec = 0;
    int doit;
    
   /* for(i =0; i <20;i++){
        printf("%lf\n", array[gen][i]);
    }
    int count =0;
    for(a=size; a<size+1;a++){
        for(k=0; k< 1000; k++){
            speciescount[0][k] = 0;
            speciescountnongap[0][k] = 0;
            speciescountrecruit[0][k] = 0;
            speciescountnongaprecruits[0][k] = 0;
        }
        if(a == 1){//LARGE
            for(i=0;i<150000; i++){
                array = *largearray;
            }
        }else if(a == 2){
            array = *smallarray;
        }else{
            array = *allarray;
        }
        for(t = 0; t < groupgapnum; t++){
            doit = 1;
            if(a == 1){//LARGE
                if(groupgapsize[t] < 3){
                    doit = 0;
                }
            }else if(a == 2){
                if(groupgapsize[t] > 2){
                    doit = 0;
                }
            }
            if(doit != 0){
                for(i = 0; i < gap_number; i++){
                    ingap = gapgroup[i];
                    if(gapgroup[i] == t){
                        for(j=0; j<numindiv; j++){
                            x = LC_xcrd[j];
                            y = LC_ycrd[j];
                            gx = gap_xcrd[i];
                            gy = gap_ycrd[i];
                            gxsize =gap_sizex[i];
                            gysize =gap_sizey[i];
                            if(x!=-1){
                                if(x > gx && y > gy && y <= gy+gysize && x <= gx+gxsize){
                                    p = spID[j];
                                    speciescount[0][spID[j]]++;
                                    r=speciescount[0][spID[j]];
                                    count++;
                                }
                            }
                        }
                    }
                }
                //FOR NON GAPS
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
                            if(LC_xcrd[j] != -1 && LC_xcrd[j] > nongap_xcrd[i] && LC_ycrd[j] > nongap_ycrd[i] && LC_ycrd[j] <= (nongap_ycrd[i]+gap_sizey[i]) && LC_xcrd[j] <= (nongap_xcrd[i]+gap_sizex[i])){
                                // printf("%d ", speciescount[i][spID[j]]);
                                p = spID[j];
                                speciescountnongap[0][spID[j]]++;
                                r=speciescountnongap[0][spID[j]];
                                //if(speciescount[curr_year][i][spID[j]] > 0){
                                //printf("%d,%d\n",p,r);
                                
                                // }
                                // printf("%d\n", speciescount[i][spID[j]]);
                            }
                        }
                    }
                    //}
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
                                speciescountrecruit[0][spID[j]]++;
                                r=speciescountrecruit[0][spID[j]];
                                //if(speciescount[curr_year][i][spID[j]] > 0){
                                //printf("%d,%d\n",p,r);
                                // }
                                // printf("%d\n", speciescount[i][spID[j]]);
                            }
                        }
                    }
                    //}
                }
                //FOR NON GAPS recruits
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
                            if(LC_xcrd[j] != -1 && LC_xcrd[j] > nongap_xcrd[i] && LC_ycrd[j] > nongap_ycrd[i] && LC_ycrd[j] <= (nongap_ycrd[i]+gap_sizey[i]) && LC_xcrd[j] <= (nongap_xcrd[i]+gap_sizex[i])){
                                // printf("%d ", speciescount[i][spID[j]]);
                                p = spID[j];
                                speciescountnongaprecruits[0][spID[j]]++;
                                r=speciescountnongaprecruits[0][spID[j]];
                                //printf("%d,%d\n",p,r);
                                //if(speciescount[curr_year][i][spID[j]] > 0){
                                
                                // }
                                // printf("%d\n", speciescount[i][spID[j]]);
                            }
                        }
                    }
                    //}
                }
            }
        }
        //2. calculate H,richness
        richgaps = 0;
        richnongaps = 0;
        richgapsrec = 0;
        richnongapsrec = 0;
        Hgaps = 0;
        Hnongaps = 0;
        Hgapsrec = 0;
        Hnongapsrec = 0;
        Egaps = 0;
        Enongaps = 0;
        Egapsrec = 0;
        Enongapsrec = 0;
        sumofspgaps = 0;
        sumofspnongaps = 0;
        sumofspgapsrec = 0;
        sumofspnongapsrec = 0;
        avgwdgaps = 0;
        avgwdnongaps = 0;
        avgwdgapsrec = 0;
        avgwdnongapsrec = 0;
        avgfacgaps = 0;
        avgfacnongaps = 0;
        avgfacgapsrec = 0;
        avgfacnongapsrec = 0;
        cwdgaps = 0;
        cwdnongaps = 0;
        cwdgapsrec = 0;
        cwdnongapsrec = 0;
        cfacgaps = 0;
        cfacnongaps = 0;
        cfacgapsrec = 0;
        cfacnongapsrec = 0;
        x=0;
        
        for(k=0; k < num_species; k++){
            sumofspgaps+=speciescount[0][k];
            sumofspnongaps+=speciescountnongap[0][k];
            sumofspgapsrec +=speciescountrecruit[0][k];
            sumofspnongapsrec+=speciescountnongaprecruits[0][k];
        }
        for(k=0; k < num_species; k++){
            if(speciescount[0][k] != 0){
                richgaps++;
                b =speciescount[0][k];
                x =(double) b/sumofspgaps;
                Hgaps += x*log(x);
            }
            if(speciescountnongap[0][k] != 0){
                richnongaps++;
                x =(double) speciescountnongap[0][k]/sumofspnongaps;
                Hnongaps += x*log(x);
                
            }
            if(speciescountrecruit[0][k] != 0){
                richgapsrec++;
                x =speciescountrecruit[0][k]/sumofspgapsrec;
                Hgapsrec += (double) x*log(x);
                
            }
            if(speciescountnongaprecruits[0][k] != 0){
                richnongapsrec++;
                x =speciescountnongaprecruits[0][k]/sumofspnongapsrec;
                Hnongapsrec += (double) x*log(x);
            }
            if(wd[k] != -10){
                if(speciescount[0][k] != 0){
                    b =speciescount[0][k];
                    c = wd[k];
                    avgwdgaps += (speciescount[0][k] * wd[k]);
                    cwdgaps += speciescount[0][k];
                }
                if(speciescountnongap[0][k] != 0){
                    avgwdnongaps += (speciescountnongap[0][k] * wd[k]);
                    cwdnongaps += speciescountnongap[0][k];
                }
                if(speciescountrecruit[0][k] != 0){
                    avgwdgapsrec += (speciescountrecruit[0][k] * wd[k]);
                    cwdgapsrec += speciescountrecruit[0][k];
                }
                if(speciescountnongaprecruits[0][k] != 0){
                    avgwdnongapsrec += (speciescountnongaprecruits[0][k] * wd[k]);
                    cwdnongapsrec += speciescountnongaprecruits[0][k];
                }
            }
            if(factor[k] != -10){
                if(speciescount[0][k] != 0){
                    avgfacgaps += (speciescount[0][k] * factor[k]);
                    cfacgaps += speciescount[0][k];
                }
                if(speciescountnongap[0][k] != 0){
                    avgfacnongaps += (speciescountnongap[0][k] * factor[k]);
                    cfacnongaps += speciescountnongap[0][k];
                }
                if(speciescountrecruit[0][k] != 0){
                    avgfacgapsrec += (speciescountrecruit[0][k] * factor[k]);
                    cfacgapsrec += speciescountrecruit[0][k];
                }
                if(speciescountnongaprecruits[0][k] != 0){
                    avgfacnongapsrec += (speciescountnongaprecruits[0][k] * factor[k]);
                    cfacnongapsrec += speciescountnongaprecruits[0][k];
                }
            }
        }
        Hgaps *= -1;
        Hnongapsrec *= -1;
        Hgapsrec *= -1;
        Hnongaps *= -1;
        
        //3. calculate evenness
        Egaps = Hgaps/(log(sumofspgaps));
        Enongaps = Hnongaps/(log(sumofspnongaps));
        Egapsrec = Hgapsrec/(log(sumofspgapsrec));
        Enongapsrec = Hnongapsrec/(log(sumofspnongapsrec));
        //4. calculate trait values
        avgfacgaps /= cfacgaps;
        avgfacnongaps /= cfacnongaps;
        if(cfacgapsrec != 0){avgfacgapsrec /= cfacgapsrec;}
        if(cfacgapsrec != 0){avgfacnongapsrec /= cfacnongapsrec;}
        
        avgwdgaps /= cwdgaps;
        avgwdnongaps /= cwdnongaps;
        if(cwdgapsrec != 0){avgwdgapsrec /= cwdgapsrec;}
        if(cwdnongapsrec != 0){avgwdnongapsrec /= cwdnongapsrec;}
        
        z =0;
        array[gen][z] = richgaps;
        array[gen][z+1] = Egaps;
        array[gen][z+2] = Hgaps;
        array[gen][z+3] = avgwdgaps;
        array[gen][z+4] = avgfacgaps;
        z=5;
        array[gen][z] = richnongaps;
        array[gen][z+1] = Enongaps;
        array[gen][z+2] = Hnongaps;
        array[gen][z+3] = avgwdnongaps;
        array[gen][z+4] = avgfacnongaps;
        z=10;
        array[gen][z] = richgapsrec;
        array[gen][z+1] = Egapsrec;
        array[gen][z+2] = Hgapsrec;
        array[gen][z+3] = avgwdgapsrec;
        array[gen][z+4] = avgfacgapsrec;
        z=15;
        array[gen][z] = richnongapsrec;
        array[gen][z+1] = Enongapsrec;
        array[gen][z+2] = Hnongapsrec;
        array[gen][z+3] = avgwdnongapsrec;
        array[gen][z+4] = avgfacnongapsrec;
        
        
        if(a == 1){
            //free(*largearray);
            *largearray = array;
            /*for(i =0; i <20;i++){
             printf("%lf\n", array[gen][i]);
             }
             array = *largearray;
             for(i =0; i <20;i++){
             printf("%lf\n", array[gen][i]);
             }
        }else if(a==2){
            *smallarray = array;
            
        }else{
            *allarray = array;
        }
        
        /*for(i=0;i<200000; i++){
         for (j = 0; j<50; ++j) {
         array[i][j] = 0;
         }
         }
    }
}
*/
//was in print SAD group gaps
//for a grouped gap
/* for(t = 0; t < groupgapnum; t++){
 for(k=0; k< 1000; k++){
 speciescount[k] = 0;
 }
 //find out how many real gaps are in that grouped gap
 for(q = 0; q< groupgapsize[t]; q++){
 //for the real gap in the grouped gap, count the species
 for(i = 0; i < gap_number; i++){
 if(gapgroup[i] == t){
 for(j=0; j<numindiv; j++){
 if(LC_xcrd[j] > gap_xcrd[i] && LC_ycrd[j] > gap_ycrd[i] && LC_ycrd[j] <= (gap_ycrd[i]+gap_sizey[i]) && LC_xcrd[j] <= (gap_xcrd[i]+gap_sizex[i])){
 // printf("%d ", speciescount[i][spID[j]]);
 speciescount[spID[j]]++;
 //if(speciescount[curr_year][i][spID[j]] > 0){
 
 // }
 // printf("%d\n", speciescount[i][spID[j]]);
 }
 }
 }
 }
 }
 //upon getting all q gaps in the group gap tallied, write the group gap to a file
 sprintf(outfile, "%ssimulation%dgap%dspeciescountyr%d.txt", begoutfile,simnumber, t,curr_year);
 ofp = fopen(outfile, "w");
 for(z=0;z<num_species; z++){
 fprintf(ofp, "%d\t%d\n", z, speciescount[z]);
 }
 fclose(ofp);
 }*/

 /*
void print_smallvslarge_gaps(char *begoutfile, int *gap_xcrd, int *gap_ycrd, int *gap_sizex, int *gap_sizey, int gap_number, double *LC_xcrd, double *LC_ycrd, int numindiv, int *spID, int num_species, int **speciescount){
    
    
    
}

void print_SAD_in_actual(char *begoutfile, int *gap_xcrd, int *gap_ycrd, int *gap_sizex, int *gap_sizey, int gap_number, int simnumber, double *LC_xcrd, double *LC_ycrd, int numindiv, int *spID, int num_species){
    
    FILE *ifp;
    double *LC_xcrdtemp, *LC_ycrdtemp;
    int i;
    int status;
    
    LC_xcrdtemp = calloc(236000, sizeof(double));
    for(i=0;i<236000; i++){
        LC_xcrdtemp[i]=0;
    }
    LC_ycrdtemp = calloc(236000, sizeof(double));
    for(i=0;i<236000; i++){
        LC_ycrdtemp[i]=0;
    }
    
    ifp = fopen("individualson35yr50haplot.txt", "r");
    status=fscanf(ifp, "%s %f %f %d", speciesnameTEMP, &x1, &y1,&dbh);
    fflush(ifp);
    
    while (status != EOF) {
        status=fscanf(ifp, "%s %f %f %d", speciesnameTEMP, &x1, &y1,&dbh);
    }

}
  
  void print_recruits(char *begoutfile, char *simtype, int *gap_xcrd, int *gap_ycrd, int *gap_sizex, int *gap_sizey, int gap_number, int simnumber, double *LC_xcrd, double *LC_ycrd, int numindiv, int *spID, int num_species, int **speciescount, int curr_year, int *gapgroup, int groupgapnum, int *groupgapsize, int numindivstart){
  int i,z,j,k,t,q;
  char outfile[MAX_FLNAME_SIZE];
  FILE *ofp;
  
  int ingap,g1,g2,g3,g4,p,n,r;
  double lc1,lc2;
  n= numindiv;
  //for a grouped gap
  for(t = 0; t < groupgapnum; t++){
  for(k=0; k< 1000; k++){
  speciescount[t][k] = 0;
  }
  for(i = 0; i < gap_number; i++){
  ingap = gapgroup[i];
  if(gapgroup[i] == t){
  for(j=numindiv; j<n; j++){ //only change from print SAD is here
  lc1 = LC_xcrd[j];
  lc2 = LC_ycrd[j];
  g1 = gap_xcrd[i];
  g2 = gap_ycrd[i];
  g3 =gap_sizex[i];
  g4 =gap_sizey[i];
  if(LC_xcrd[j] != -1 && LC_xcrd[j] > gap_xcrd[i] && LC_ycrd[j] > gap_ycrd[i] && LC_ycrd[j] <= (gap_ycrd[i]+gap_sizey[i]) && LC_xcrd[j] <= (gap_xcrd[i]+gap_sizex[i])){
  p = spID[j];
  speciescount[t][spID[j]]++;
  r=speciescount[t][spID[j]];
  }
  }
  }
  }
  //upon getting all q gaps in the group gap tallied, write the group gap to a file
  sprintf(outfile, "%srecruitssimulation%d%sgap%dspeciescountyr%d.txt", begoutfile,simnumber, simtype,t,curr_year);
  ofp = fopen(outfile, "w");
  for(z=0;z<num_species; z++){
  fprintf(ofp, "%d\t%d\n", z, speciescount[t][z]);
  }
  fclose(ofp);
  }
  
  
  }
*/
