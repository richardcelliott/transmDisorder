#===================================================
# Makefile
#===================================================


# statp program
PROGRAM = xExciton
PROGRAM_OBJS = \
	main.o 

PROGRAM_FILES = 

# custom library and its objects
MY_LIB = libxcrnk.a
LIB_OBJS = \

LIBSRC = \
	
# other lib and include
LIBS = $(MY_LIB) 
INCLUDES = -lm -larmadillo -lblas -llapack -DARMA_DONT_USE_WRAPPER -ffast-math -fno-math-errno -I./ 
#-DARMA_DONT_USE_WRAPPER -lblas -llapack -I./

# compiler
CXX = g++
CC = g++
LD = g++

ifeq ($(EX_DEBUG), 1)
EX_DEBUG_CFLAGS = -DEX_DEBUG
endif

#CXXFLAGS = -g -Wall 
CXXFLAGS = -O3 \
	-Wall \
	-std=c++14 #\
	-ffast-math
	#-larmadillo
	#-DARMA_DONT_USE_WRAPPER
	#-static \
	#-pipe \
	#-march=pentium4 \
	#-mfpmath=sse \
	#-fomit-frame-pointer \
	#-funroll-loops 

SUFFIXCXX = .cpp
SUFFIXC = .c

#===================================================
# Rules
#===================================================
ALL: $(MY_LIB) $(PROGRAM)

# build program
$(PROGRAM) : $(PROGRAM_OBJS)
	$(LD) -o $(PROGRAM) $(PROGRAM_OBJS) $(INCLUDES) $(LIBS)

$(MY_LIB): $(LIB_OBJS)
	ar rv $(MY_LIB) $(LIB_OBJS)
	ranlib $(MY_LIB)

# compile .o
%.o: %$(SUFFIXCXX)
	$(CXX) -c $(CXXFLAGS) $< -o $@

%.o: %$(SUFFIXC)
	$(CC) -c $(CFLAGS) $< -o $@

#%.o: %$(SUFFIXF)
#	$(FC) -c $(CFLAGS) $< -o $@

# clean
.PHONY: clean
clean :
	rm -f *.o $(MY_LIB) $(PROGRAM) core* 
#=============================================================
#	DEPENDENCIES
#=============================================================

#exciton_prop.o: exciton_prop.cpp \
	#exciton_prop.h \
	#consts.h 
main.o: main.cpp \
	consts.h 

