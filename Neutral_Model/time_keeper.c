/* Copied this from 'A Book on C' by Kelley and Pohl.  Use it to time
   programs.*/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "time_keeper.h"

# define MAXSTRING 100

typedef struct {
  clock_t begin_clock, save_clock;
  time_t begin_time, save_time;
} time_keeper;

static time_keeper tk;

void start_time(void)
{
  tk.begin_clock = tk.save_clock = clock();
  tk.begin_time = tk.save_time = time(NULL);
}

void prn_time(void)
{
  char s1[MAXSTRING], s2[MAXSTRING];
  int field_width, n1, n2;
  double clocks_per_second = (double) CLOCKS_PER_SEC, user_time, real_time;
  
  user_time = (clock() - tk.save_clock) / clocks_per_second;
  real_time = difftime(time(NULL), tk.save_time);
  tk.save_clock = clock();
  tk.save_time = time(NULL);

  /* print the values found, and do it neatly */

  n1 = sprintf(s1, "%.1f", user_time);
  n2 = sprintf(s2, "%.1f", real_time);
  field_width = (n1 > n2) ? n1 : n2;
  printf("SINCE LAST PRN_TIME:\n");
  printf("%s%*.1f%s%f%s%f%s\n%s%*.1f%s%f%s%f%s\n", 
	 "User time: ", field_width, user_time, " seconds, ", user_time/60, " minutes, ", user_time/3600, " hours ", 
	 "Real time: ", field_width, real_time, " seconds, ", real_time/60, " minutes, ", real_time/3600, " hours ");

}

void prn_total_time(void)
{
  char s1[MAXSTRING], s2[MAXSTRING];
  int field_width, n1, n2;
  double clocks_per_second = (double) CLOCKS_PER_SEC, user_time, real_time;
  
  user_time = (clock() - tk.begin_clock) / clocks_per_second;
  real_time = difftime(time(NULL), tk.begin_time);

  /* print the values found, and do it neatly */

  n1 = sprintf(s1, "%.1f", user_time);
  n2 = sprintf(s2, "%.1f", real_time);
  field_width = (n1 > n2) ? n1 : n2;
  printf("SINCE START:\n");
  printf("%s%*.1f%s%f%s%f%s\n%s%*.1f%s%f%s%f%s\n\n", 
	 "User time: ", field_width, user_time, " seconds, ", user_time/60, " minutes, ", user_time/3600, " hours ", 
	 "Real time: ", field_width, real_time, " seconds, ", real_time/60, " minutes, ", real_time/3600, " hours ");

}
