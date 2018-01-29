#-----------------------------------------------------------------------
# This Makefile is for building dependency generator "makdep" on a CRAY, 
# SGI or SUN.  The executable will be put in $(HOME)/bin.  If this directory
# does not exist or is not in $PATH, problems will occur when the CCM is built.
#------------------------------------------------------------------------

#include ../configure.lvt # from EMK for NU-WRF

EXEDIR  := .
EXENAME := makdep

#------------------------------------------------------------------------
# Linux
#------------------------------------------------------------------------

ifeq ("$(LIS_ARCH)","AIX")
CFLAGS := -O -qcpluscmt
else
CFLAGS := -O
#CFLAGS := $(CFLAGS_OPT) # from EMK for NU-WRF
endif

#------------------------------------------------------------------------
# Default rules
#------------------------------------------------------------------------

$(EXEDIR)/$(EXENAME):	main.o
	$(CC) -o $@ main.o 

main.o: main.c
	$(CC) -c $(CFLAGS) $<

RM := rm

clean:
	$(RM) -f *.o $(EXEDIR)/$(EXENAME)
