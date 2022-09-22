# ----- Makefile FOR gsm CODE -----
#

OPTIMIZE =
OPT1    =

OPTIONS =  $(OPTIMIZE) \
	   $(OPT1)


EXEC = fkpt

all: $(EXEC)

# C compiler:
CC = gcc
# With icc (Intel):
#CC = icc

# Default CFLAGS:
#~ INCLUDES    =  -I ./ -I ./libs/
CFLAGS = -g -O3 $(OPTIONS)
#
# With OpenMP:
#CFLAGS = -g -O3 -fopenmp $(OPTIONS)
# For icc use instead:
#CFLAGS = -g -O3 -qopenmp $(OPTIONS)


#
# Nothing to do below
#
H_PATH = libs

S_PATH = src


OBJS	= $(S_PATH)/main.o $(S_PATH)/gsm_io.o \
	$(S_PATH)/write.o $(S_PATH)/startrun.o \
	$(S_PATH)/models.o $(S_PATH)/gsm_diffeqs.o \
	$(S_PATH)/kfunctions.o $(S_PATH)/global.o $(S_PATH)/rsd.o\
	$(H_PATH)/diffeqs.o $(H_PATH)/mathutil.o \
	$(H_PATH)/clib.o $(H_PATH)/getparam.o \
	$(H_PATH)/mathfns.o $(H_PATH)/numrec.o \
	$(H_PATH)/inout.o $(H_PATH)/quads.o

INCL	= $(S_PATH)/globaldefs.h $(S_PATH)/cmdline_defs.h \
    $(S_PATH)/protodefs.h $(S_PATH)/models.h $(S_PATH)/rsd.h \
    $(H_PATH)/diffeqs.h \
	$(H_PATH)/getparam.h $(H_PATH)/machines.h \
	$(H_PATH)/mathfns.h $(H_PATH)/stdinc.h \
	$(H_PATH)/strings.h $(H_PATH)/precision.h \
	$(H_PATH)/vectdefs.h $(H_PATH)/vectmath.h \
	$(H_PATH)/numrec.h \
	$(H_PATH)/inout.h $(H_PATH)/mathutil.h \
	$(H_PATH)/switchs.h $(H_PATH)/quads.h

#~ $(EXEC): $(OBJS) 
#~ 	($(CC) $(OBJS) $(LIBS) $(CFLAGS) -o $@ -lm; cp $(EXEC) ../)
$(EXEC): $(OBJS) 
	($(CC) $(OBJS) $(LIBS) $(CFLAGS) -o $@ -lm)

$(OBJS): $(INCL)

clean:
	(rm -f $(OBJS) $(EXEC); rm -fR mgpt.dSYM; rm $(EXEC))
#~ clean:
#~ 	(rm -f $(OBJS) $(EXEC); rm -fR mgpt.dSYM; rm ../$(EXEC))


.PHONY : all clean check install uninstall

