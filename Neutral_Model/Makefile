BASE    = .
MAIN    = .
EFILE   = gap
MYLIBS  = .
CC      = cc
CFLAGS  = -g -I $(MYLIBS)


OBJS = $(MAIN)/main.o $(MAIN)/initializelandscape.o $(MAIN)/stats_analysis.o $(MAIN)/timesteps+dispersal.o $(MAIN)/time_keeper.o $(MYLIBS)/general.o

$(EFILE): $(OBJS)
	$(CC) $(CFLAGS) -o $(EFILE) $(OBJS) -lm

clean:
	/bin/rm -f $(OBJS) $(EFILE)