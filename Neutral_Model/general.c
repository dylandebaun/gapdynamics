/* standard header files */
# include <stdio.h>
# include <stdlib.h>
# include <math.h>

/* my header files */
# include "general.h"

/*function to input if want qsort to sort in descending order*/
int is_x_lower_than_y(void *vp, void *vq)
{
  int *point_x=vp, *point_y=vq;

  /*  printf("x=%d, y=%d\n", (int)*point_x, (int)*point_y);*/
  return *point_y-*point_x;


}

/* used to have rounddown and roundup--use floor(x) and ceil(x) isntead, they are math library functions!!*/

/* take care not to use this on doubles which are out of the INT_MAX range */
/* I've read that it is more efficient to use a #def command instead of a function here, so if wanting to find a way to speed things up, may want to do that instead.  I'm not sure how much it really matters.  Most of my lattice model functions don't use quite this anyway, which is were speed would matter. */
int myround(double x)
{
  //int intcast;

  if (fabs(((double) ceil(x)) - x) > 0.5) 
    return (int) floor(x);
  else
    return (int) ceil(x);
}

void modular(int x, int y, int *xmody)
{

  xmody[0] = floor(((double) x)/((double) y));
  xmody[1] = x - xmody[0]*y;
  
}


void modular_dbl(double x, double y, double *xmody)
{

  xmody[0] = x/y;
  xmody[1] = x - xmody[0]*y;

}


double abs_dbl(double x)
{
  if (x<0)
    return -x;
  else
    return x;
}


/* get_rand_unit returns a real random number in the interval [0,1],
   that is, including 0 and 1 as possiblities.  rand() returns a
   random integer between 0 and RAND_MAX.  Hence there is a discrete
   set of real numbers that this function can return--this set
   contains about 2 billion (2X10^9) numbers.  */
/* REMEMBER TO SEED THE RANDOM NUMBER GENERATOR when using this
   function! --srand(time(NULL)) */
double get_rand_unit()
{

  double toss;

  toss = (double) rand();
  toss = toss/((double) RAND_MAX);
  
  return toss;

}

/* get_rand_integ_intvl returns a random integer in the interval [x,y], that
   is, including x and y as possibilities. */
/* REMEMBER TO SEED THE RANDOM NUMBER GENERATOR when using this
   function! --srand(time(NULL)) */

int get_rand_integ_intvl(int x, int y)
{

  double toss;

  toss = get_rand_unit();
  /*  printf("toss=%f, y-x+1=%d, toss*(y-x+1)=%f, ceil(toss*(y-x+1))=%f, returning %d.\n", toss, y-x+1, toss*(y-x+1), ceil(toss*(y-x+1)), (int) ceil(toss*(y-x+1))-1 + x);*/

  /* technically should throw out zeros, because should choose a particular integer based on a region non-inclusive of the one before that, even when choosing the lowest value in the integer range */
  while (toss==0) toss=get_rand_unit();

  return (int) ceil(toss*(y-x+1))-1 + x;


}

/* get_rand_real_intvl returns a random real number in the interval
   [x,y], that is, including x and y as possibilities. */
/* REMEMBER TO SEED THE RANDOM NUMBER GENERATOR when using this
   function! --srand(time(NULL)) */
double get_rand_real_intvl(double x, double y)
{

  double toss;

  toss = rand();
  return ((toss/RAND_MAX)*(y-x) + x);

}

double area_ring(double r1, double r2)
{

  double area;

  area = 3.14*(pow(r2,2) - pow(r1,2));

  return area;

}


void descending_order(int *old_array, int *new_array, int array_size)
{
  int i, j, place;

  for (i=0; i<array_size; ++i)
    new_array[i] = 0;
  for (i=0; i<array_size; ++i) {
    place=0;
    while (old_array[i] < new_array[place] && place<array_size-1) ++place;
    for(j=array_size-1; j>=place; --j) {
      new_array[j+1] = new_array[j];
    }
    new_array[place] = old_array[i];
  }

}
						   
double deg_to_radians(double x)
{

  return (double)((x/360)*2*3.14);

}

double frac_circle_in_region(double x, double y, double r, double l_x, double l_y)
{
  
  int i;
  double frac=0;
  
  for (i=1; i<=360; ++i) 
    if ((x+r*cos(deg_to_radians(i)) < l_x) && (x+r*cos(deg_to_radians(i)) > 0) && (y+r*sin(deg_to_radians(i)) < l_y) && (y+r*sin(deg_to_radians(i)) > 0)) 
      frac = frac +1;
    
  frac = frac/360;

  return frac;

}

double frac_circle_in_region_faster(double x, double y, double r, double l_x, double l_y)
{

  int i;
  double x_1=r, x_2=r, y_1=r, y_2=r, theta_q[5], frac=0;

  if (r==0)
    return 1;
  else {
    if (x-r < 0) x_1 = x;
    if (x+r > l_x) x_2 = l_x-x;
    if (y-r < 0) y_1 = y;
    if (y+r > l_y) y_2 = l_y-y;

    theta_q[1] = max(3.14/2 - (acos(x_1/r) + acos(y_1/r)), 0);
    theta_q[2] = max(3.14/2 - (acos(y_1/r) + acos(x_2/r)), 0);
    theta_q[3] = max(3.14/2 - (acos(x_2/r) + acos(y_2/r)), 0);
    theta_q[4] = max(3.14/2 - (acos(y_2/r) + acos(x_1/r)), 0);
  
    for (i=1; i<=4; ++i)
      frac = frac + theta_q[i];

    frac = frac/(2*3.14);

    return frac;
  }

}

double max(double x, double y)
{
  if (x>y) return x;
  else return y;

}

int intmax3(int x, int y, int z)
{
  if (x >= y && x >= z) return x;
  if (y >= x && y >= z) return y;
  if (z >= x && z >= y) return z;
  else return 0;
    
}

int intmax4(int w, int x, int y, int z)
{
  if (w>= x && w>=y && w>=z) return w;
  if (x>= w && x>=y && x>=z) return x;
  if (y>= w && y>=w && y>=z) return y;
  if (z>= w && z>=x && z>=y) return w;
  else return 0;
}

int intmax_array(int *max_array, int num)
{

  int i, max=-RAND_MAX;

  for (i=1; i<=num; ++i) 
    if (max_array[i] > max) max=max_array[i];

  return max;

}

void printout_2d_dist(double **occurrence_list, int num_occurrences, double x_dist_range, double x_dist_bin_size, double y_dist_range, double y_dist_bin_size, char *beg_outfilename, char *x_name, char *y_name)
{

  int i, j, num_x_bins, num_y_bins, which_x_bin, which_y_bin, num_occ_in_range=0;
  double *x_dist, *y_dist, **xy_dist, x, y;
  char outfilename[MAX_FLNAME_SIZE];
  FILE *ofp;

  num_x_bins = (int) floor(x_dist_range/x_dist_bin_size);
  num_y_bins = (int) floor(y_dist_range/y_dist_bin_size);
  
  x_dist = calloc(num_x_bins+1, sizeof(double));
  y_dist = calloc(num_y_bins+1, sizeof(double));
  xy_dist = calloc(num_x_bins+1, sizeof(double *));
  for (i=1; i<=num_x_bins; ++i) 
    xy_dist[i] = calloc(num_y_bins+1, sizeof(double));

  for (i=1; i<=num_occurrences; ++i) {
    x=occurrence_list[i][1];
    y=occurrence_list[i][2];
    which_x_bin = (int) ceil(x/x_dist_bin_size);
    which_y_bin = (int) ceil(y/y_dist_bin_size);
    if (which_x_bin <= num_x_bins && which_y_bin <= num_y_bins) {
      x_dist[which_x_bin] = x_dist[which_x_bin] +1;
      y_dist[which_y_bin] = y_dist[which_y_bin] +1;
      xy_dist[which_x_bin][which_y_bin] = xy_dist[which_x_bin][which_y_bin] +1;
      num_occ_in_range = num_occ_in_range +1;
    }
  }
  printf("%d occurrences were out of the distribution calculation range.\n", num_occurrences-num_occ_in_range); 

  for (i=1; i<=num_x_bins; ++i)
    x_dist[i] = x_dist[i]/num_occurrences;

  for (i=1; i<=num_y_bins; ++i)
    y_dist[i] = y_dist[i]/num_occurrences;

  for (i=1; i<=num_x_bins; ++i)
    for (j=1; j<=num_y_bins; ++j)
      xy_dist[i][j] = xy_dist[i][j]/num_occurrences;

  sprintf(outfilename, "%s_%s_dist.dat", beg_outfilename, x_name);
  ofp = fopen(outfilename, "w");
  for (i=1; i<=num_x_bins; ++i)
    fprintf(ofp, "%f\t%f\n", ((double) i)*x_dist_bin_size, x_dist[i]);  
  fclose(ofp);

  sprintf(outfilename, "%s_%s_dist.dat", beg_outfilename, y_name);
  ofp = fopen(outfilename, "w");
  for (i=1; i<=num_y_bins; ++i)
    fprintf(ofp, "%f\t%f\n", ((double) i)*y_dist_bin_size, y_dist[i]);  
  fclose(ofp);
  
  sprintf(outfilename, "%s_%s%s_dist.dat", beg_outfilename, x_name, y_name);
  ofp = fopen(outfilename, "w");
  for (i=1; i<=num_x_bins; ++i)
    for (j=1; j<=num_y_bins; ++j)
      fprintf(ofp, "%f\t%f\t%f\n", ((double) i)*x_dist_bin_size, ((double) i)*y_dist_bin_size, xy_dist[i][j]);  
  fclose(ofp);

}

double sample_var(int num_samples, double *samples, double *point_mean, double *p_sum_of_samples, double *p_sum_of_sqrd_samples)
{

  int i;
  double sum_of_squares;
  
  *p_sum_of_samples =0;
  *p_sum_of_sqrd_samples =0;

  for (i=0; i<num_samples; ++i) {
    *p_sum_of_samples += samples[i];
    *p_sum_of_sqrd_samples += pow(samples[i],2);
  }
  
  *point_mean = *p_sum_of_samples/((double) num_samples);
  sum_of_squares = *p_sum_of_sqrd_samples-(1/((double) num_samples))*pow(*p_sum_of_samples,2);

  return sum_of_squares/((double) num_samples -1);

}

  /* compute means and confidence intervals of means, using hard-coded info from student's t-distribution at 95% (because 5% on each side) confidence see mathworld.wolfram.com/Studentst-Distribution.html and Jim Kirchner's notes and Zar to understand why take 95% one*/
void get_90perc_CL(int num_samples, double *samples, double *point_mean, double *point_lower_CL, double *point_upper_CL, double *p_sum_of_samples, double *p_sum_of_sqrd_sam) //was written *p_sum_of_sqrd_sam=0.0ples why? i deleted it b/c the variable can't be defined
{

  double variance;
   // int i;
  double t;

    variance = sample_var(num_samples, samples, point_mean, p_sum_of_samples, p_sum_of_sqrd_sam);

  if (num_samples==2) t=6.31371;
  if (num_samples==3) t=2.91999;
  if (num_samples==4) t=2.35336;
  if (num_samples==5) t=2.13185;
  if (num_samples>=6 && num_samples < 11) t=2.01505;
  if (num_samples>=11 && num_samples < 31) t=1.81246;
  if (num_samples>=31 && num_samples < 101) t=1.69726;
  if (num_samples>=101 && num_samples < 10000) t=1.66023;
  else t=1.64487;  /* for when num_samples>=10000 */
  
  *point_lower_CL = *point_mean - t*sqrt(variance/num_samples);
  *point_upper_CL = *point_mean + t*sqrt(variance/num_samples);
  
}
  
double get_median(int num_samples, double *samples)
{
  int i;
  double median, samples_copy[num_samples+1];

  for (i=1; i<=num_samples; ++i) 
    samples_copy[i] = samples[i];

  for (i=1; i<=num_samples/2; ++i) {
    delete_lowest(num_samples, samples_copy);
  }

  median=get_lowest(num_samples, samples_copy);

  return median;

}


/* computes the quantile boundaries of the distribution, which is different from the confidence interval of the mean (above), and in contrast to the std deviation can reflect sqewness of a distribution*/
void get_quantile_boundaries(int num_samples, double *samples, double interval_size, double *point_mean, double *point_lower_bound, double *point_upper_bound)
{

  int i, num_delete;
  double samples_copy[num_samples+1]; /*!!! Having trouble with this as dynamically allocated, not sure why, so hard coded for now. !!!!*/


  *point_mean=get_mean(num_samples, samples);
  /*   printf("The mean is %f.\n", *point_mean);*/

  /*  printf("Copying samples.\n");*/
  /*  samples_copy = calloc(num_samples+1, sizeof(double));*/
  /*  printf("Have allocated samples_copy\n");*/
  for (i=1; i<=num_samples; ++i) 
    samples_copy[i] = samples[i];

  num_delete = myround(((100-interval_size)/2)*0.01*num_samples);
  /*  printf("Deleting the %d lowest and %d highest samples.\n", num_delete, num_delete);*/

  for (i=1; i<=num_delete; ++i) {
    /*    printf("deleting the %dth lowest and highest\n", i);*/
    delete_lowest(num_samples, samples_copy);
    delete_highest(num_samples, samples_copy);
  }
  /*  printf("finsihed deleting the %d lowest and %d highest samples.\n", num_delete, num_delete);*/
  *point_lower_bound = get_lowest(num_samples, samples_copy);
  *point_upper_bound = get_highest(num_samples, samples_copy);
  
  /*  free(samples_copy);*/
  
}

double get_mean(int num_samples, double *samples)
{

  int i;
  double mean=0;

  for (i=1; i<=num_samples; ++i)
    mean = mean+samples[i];

  return mean/((double) num_samples);

}

/* works for positive samples.  deletes lowest by putting a -1 in the place of that sample */

void delete_lowest(int num_samples, double *samples)
{

  int i, j, pos_lowest;
  double lowest;

  i=1;
  while (samples[i] == -1 && i<=num_samples) ++i;
  if (i==num_samples && samples[i] == -1) {
    printf("The array is empty already.  Can't delete the minimum.\n");
    exit(1);
  }
  pos_lowest=i;
  lowest=samples[i];
  /*  printf("the starting point for looking for the lowest  is at %d, with value %f\n", i, lowest);*/

  for (j=i; j<=num_samples; ++j) 
    if (samples[j] > -1 && samples[j] < lowest) {
      lowest=samples[j];
      pos_lowest=j;
    }

  /*  printf("the position of the lowest is %d, with value %f\n", pos_lowest, lowest);*/

  /*  printf("Deleting the lowest, which is at position %d, value %f.\n", pos_lowest, samples[pos_lowest]);*/
  samples[pos_lowest] = -1;

}

double get_lowest(int num_samples, double *samples)
{

  int i, j, pos_lowest;
  double lowest;

  i=1;
  while (samples[i] == -1 && i<=num_samples) ++i;
  if (i==num_samples && samples[i] == -1) {
    printf("The array is empty already.  Can't delete the minimum.\n");
    exit(1);
  }
  pos_lowest=i;
  lowest=samples[i];

  for (j=i; j<=num_samples; ++j) 
    if (samples[j] > -1 && samples[j] < lowest) {
      lowest=samples[j];
      pos_lowest=j;
    }

  /*  printf("The lowest is at %d, with value %f.\n", pos_lowest, samples[pos_lowest]);*/
  return samples[pos_lowest];

}

void delete_highest(int num_samples, double *samples)
{

  int i, j, pos_highest=1;
  double highest=0;

  i=1;
  while (samples[i] == -1 && i<=num_samples) ++i;
  if (i==num_samples && samples[i] == -1) {
    printf("The array is empty already.  Can't delete the minimum.\n");
    exit(1);
  }
  pos_highest=i;
  highest=samples[i];

  for (j=i; j<=num_samples; ++j)
    if (samples[j] > -1 && samples[j] > highest) {
      highest=samples[j];
      pos_highest=j;
    }
  /*  printf("Deleting the highest, which is at position %d, value %f.\n", pos_highest, samples[pos_highest]);*/

  samples[pos_highest] = -1;

}

double get_highest(int num_samples, double *samples)
{

  int i, j, pos_highest;
  double highest;

  i=1;
  while (samples[i] == -1 && i<=num_samples) ++i;
  if (i==num_samples && samples[i] == -1) {
    printf("The array is empty already.  Can't delete the minimum.\n");
    exit(1);
  }
  pos_highest=i;
  highest=samples[i];

  for (j=i; j<=num_samples; ++j)
    if (samples[j] > -1 && samples[j] > highest) {
      highest=samples[j];
      pos_highest=j;
    }

  /*  printf("The highest is at %d, with value %f.\n", pos_highest, samples[pos_highest]);*/
  return samples[pos_highest];

}

double correlation_coefficient(int num_samples, double *X, double *Y)
{

  int i;
  double sumX, sumXsq, sumY, sumYsq, sumprod;

  sumX=0;
  sumXsq=0;
  sumY=0;
  sumYsq=0;
  sumprod=0;

  for (i=1; i<=num_samples; ++i) {
    sumX += X[i];
    sumXsq += pow(X[i],2);
    sumY += Y[i];
    sumYsq += pow(Y[i],2);
    sumprod += X[i]*Y[i];
  }

  return (sumprod - (sumX*sumY)/((double) num_samples))/sqrt((sumXsq - pow(sumX,2)/((double) num_samples))*(sumYsq - pow(sumY,2)/((double) num_samples)));
  
}
    

int conv_rad_to_deg(double radians)
{
  return myround((radians/PI)*180);
}

double conv_deg_to_rad(int degrees)
{
  return (((double) degrees)/((double) 180))*PI;
}

/* works for positive numbers */
int mymax(int *list, int size_list)
{
  int i, max_value=0;
  
  for (i=1; i<=size_list; ++i)
    if (list[i] > max_value) max_value=list[i];

  return max_value;

}

void convert_to_coords(int label, int size_y, int *coords) {
    
    modular(label-1,size_y, coords);
    
}

int convert_to_label(int size_y, int *coords) {
    
    int x, y;
    
    x=coords[0];
    y=coords[1];
    
    return x*size_y + y+1;
    
}


/* choose_random_location chooses at random one of the sites on the
 lattice of size size_x by size_y.  The coordinates of the randomly chosen
 site are returned in the ind1_coords array.  Only one random number
 generation is used by using the convert_to_coords function.  The
 coordinates go from 0 to size_x-1 and 0 to size_y-1.
 REMEMBER TO SEED THE RANDOM NUMBER GENERATOR when using this
 function! --srand(time(NULL)) */
void choose_random_location(int size_x, int size_y, int *coords)
{
    
    int chosen_location;
    
    chosen_location = get_rand_integ_intvl(1, size_x*size_y);
    convert_to_coords(chosen_location, size_y, coords);
    
}


/* choose_randloc_atdist chooses a random location on the lattice
 approximately a distance r from the origin.  It does this by
 picking an angle at random and rounding the coordinats xcos(angle)
 and xsin(angle).  This rounding picks the closest point on the
 lattice to these real coordinates, and returns the integer
 coordinates of that closest point in coords.*/

void choose_randloc_atdist(double r, int *coords)
{
    double angle;
    
    angle=get_rand_real_intvl(0,2*PI);
    coords[0] = myround(r*cos(angle));
    coords[1] = myround(r*sin(angle));
    
    /*  printf("angle chosen is %f radians\n", angle); */
    /*  printf("r*cos(angle)=%f, coords[0]=%d, r*sin(angle)=%f, coords[1]=%d\n", r*cos(angle), coords[0], r*sin(angle), coords[1]);*/
    
}


void choose_random_nearest_neighbor(int *ind1_coords, int size_x, int size_y, int *ind2_coords)
{
    
    int chosen_location;
    
    
    chosen_location = get_rand_integ_intvl(1, 4);
    if (chosen_location == 1) {
        ind2_coords[0] = ind1_coords[0] + 0;
        ind2_coords[1] = ind1_coords[1] + 1;
    }
    if (chosen_location == 2) {
        ind2_coords[0] = ind1_coords[0] + 1;
        ind2_coords[1] = ind1_coords[1] + 0;
    }
    if (chosen_location == 3) {
        ind2_coords[0] = ind1_coords[0] + 0;
        ind2_coords[1] = ind1_coords[1] + -1;
    }
    if (chosen_location == 4) {
        ind2_coords[0] = ind1_coords[0] + -1;
        ind2_coords[1] = ind1_coords[1] + 0;
    }
    while (ind2_coords[0] >= size_x || ind2_coords[0] < 0 ||
           ind2_coords[1] >= size_y || ind2_coords[1] < 0) {
        chosen_location = get_rand_integ_intvl(1, 4);
        chosen_location = get_rand_integ_intvl(1, 4);
        if (chosen_location == 1) {
            ind2_coords[0] = ind1_coords[0] + 0;
            ind2_coords[1] = ind1_coords[1] + 1;
        }
        if (chosen_location == 2) {
            ind2_coords[0] = ind1_coords[0] + 1;
            ind2_coords[1] = ind1_coords[1] + 0;
        }
        if (chosen_location == 3) {
            ind2_coords[0] = ind1_coords[0] + 0;
            ind2_coords[1] = ind1_coords[1] + -1;
        }
        if (chosen_location == 4) {
            ind2_coords[0] = ind1_coords[0] + -1;
            ind2_coords[1] = ind1_coords[1] + 0;
        }
    }
}


/*guassian(r, parameters) returns the gaussian probability for a seed
 to land in a ring between r and r+dr.  It is 2*pi*r times the one
 dimensional guassian probability distribution. */
double gaussian(double r, double *parameters)
{
    double sigma;
    
    sigma = parameters[0];
    
    return ((r/pow(sigma,2)) * (exp(-pow(r,2)/(2*pow(sigma,2))))); //this equation used is different from the equation cited in the article but seems closer to actual gaussian equation
    
}

double gaussianxy(double x, double y, double *parameters)
{
    double sigma;
    
    sigma = parameters[0];
    
    return (1/(2*PI*pow(sigma,2)))*(exp(-pow(distance(x,y,0,0),2)/(2*pow(sigma,2))));
}


/* fat_tailed(r, parameters) returns the probability for a seed to
 land in a ring between r and r+dr according to a fat-tailed
 dispersal kernel developped in Clark et al. 1999 Ecology
 80:1475. It is 2*pi*r times the one dimensional probability
 distribution. */
double fat_tailed(double r, double *parameters)
{
    double u, p;
    
    u = parameters[0];
    p = parameters[1];
    
    return ((2*p*r)/(u*pow(1 + (pow(r,2)/u),p+1)));
    
}

double fat_tailedxy(double x, double y, double *parameters)
{
    
    double u, p;
    
    u = parameters[0];
    p = parameters[1];
    
    
    return (p/(PI*u*pow(1 + (pow(distance(x,y,0,0),2)/u),p+1)));
    
}


// Example:
//In small-dynamic.c : select_randloc_frm_cum_disp(death_location, 200, 100, disc_norm_cum, loc_kern_cutoff, parent);

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



double distance(double x_1, double y_1, double x_2, double y_2)
{
    
    double x_dist, y_dist;
    
    x_dist = x_2-x_1;
    y_dist = y_2-y_1;
    
    /*  printf("x_dist=%f, pow(x_dist,2)=%f.\n", x_dist, pow(x_dist,2));
     printf("y_dist=%f, pow(y_dist,2)=%f.\n", y_dist, pow(y_dist,2));*/
    
    return sqrt(pow(x_dist,2) + pow(y_dist, 2));
    
}


