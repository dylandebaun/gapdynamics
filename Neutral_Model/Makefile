BASE    = .
MAIN    = .
EFILE   = gap
MYLIBS  = $(MAIN)
CC      = cc
CFLAGS  = -g -I $(MYLIBS)

# List the source files in the specified directory
OBJS = $(MAIN)/main.o \
	$(MAIN)/initializelandscape.o \
 	$(MAIN)/stats_analysis.o \
	$(MAIN)/timesteps+dispersal.o \
	$(MAIN)/time_keeper.o \
	$(MAIN)/general.o

$(EFILE): $(OBJS)
	$(CC) $(CFLAGS) -o $(EFILE) $(OBJS) -lm

clean:
	/bin/rm -f $(OBJS) $(EFILE)
