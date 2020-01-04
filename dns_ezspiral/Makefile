# ----------------------------------- #
# This is the Makefile for EZ-Spiral. #
# ----------------------------------- #

# NX, NY defined?
# -------------------
# If you wish to specify NX, NY (the number of grid points in each direction)
# at compile time, then uncomment the following and set the values to what you
# want.  They need not all be set.  Those that are set will be ignored in
# task.dat.

# NOPTION = -DNX=61 -DNY=61 

# -------------------
# Choose the compiler
# -------------------

# EZ-Spiral conforms to the ANSI C standard and *should* compile with any C
# compiler.  Use at least the -O2 compiler option, the others are optional
# options.  If you use warning options such as -Wall (which you should use if
# you modify the code), then you may get some warnings. These are probably OK.

# Compiler flags appropriate for the GNU C compiler
# -------------------------------------------------
CC = cc 
CFLAGS = -O2 -Wall -I/usr/X11R6/include

# ------------------------
# Set the linked libraries
# ------------------------

# OpenGL with X. You may need to change the line below depending on
# where libraries are found on your system.
# -----------------------------------------------------------------

LDFLAGS = -L/usr/X11R6/lib -lGL -lX11 -lXext -lm 

# otherwise for no graphics
# -------------------------

#LDFLAGS = -lm 

# -----------------------------------------------
# Probably you wont need to change anything below.
# ------------------------------------------------

OBJ=ezspiral.o ezstep.o ezgraphGL.o eztip.o

.c.o :
	$(CC) -c $(CFLAGS) $(NOPTION) $*.c

ezspiral: $(OBJ) 
	$(CC) $(CFLAGS) $(OBJ) -o ezspiral $(LDFLAGS)


# Dependencies
# ------------

ezspiral.o   :ezspiral.h 
ezstep.o     :ezspiral.h ezstep.h 
ezgraphGL.o  :ezspiral.h ezgraphGL.h 
eztip.o      :ezspiral.h 

clean:	tidy
	rm -f ezspiral *.o *~
	rm -f core *~ .*~ 	
	rm -f *.log *.dvi *.aux

tidy:
	rm -f ezspiral *.o 

tar: 
	tar -cvf ezspiral.tar README ezspiral.c ezstep.c ezgraphGL.c eztip.c \
	ezspiral.h ezstep.h ezgraphGL.h task.dat Makefile
