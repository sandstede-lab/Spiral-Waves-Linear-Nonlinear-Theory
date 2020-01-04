#
# This is the Makefile for EZ-Spiral.
#
# This Makefile should work with gcc (GNU C) on almost any machine.
# The -Wall option may give you warnings, but this is probably ok.
# If you don't have gcc you should try something like:
# CC=cc -O2 
# or whatever C compiler you want.  You may have to tell cc where 
# the Xlib libraries are, e.g:
#	$(CC) $(OBJ) -o ezspiral -lm -L/usr/X11/lib -lX11

#
# Variables
#
OBJ     = ezspiral.o ezstep.o eztip.o ezgraph.o
CC      = gcc -O
CFLAGS  = -I/usr/X11R6/include
LDFLAGS = -lm -L/usr/X11R6/lib -lX11 -lXpm
#
# Rules
#
ezspiral: $(OBJ) 
	$(CC) $(OBJ) -o ezspiral $(LDFLAGS)
#
superclean:	clean
	-@rm -f ezspiral *.o
#
clean:
	-@rm -f core
#
# Dependencies
#
ezspiral.o :ezspiral.h 
ezstep.o   :ezspiral.h ezstep.h 
ezgraph.o  :ezspiral.h ezgraph.h 
eztip.o    :ezspiral.h 
