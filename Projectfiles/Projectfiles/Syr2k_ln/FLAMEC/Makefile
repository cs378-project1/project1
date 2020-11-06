#BLAS_LIB   := $(HOME)/lib/libopenblas_sandybridgep-r0.2.20.a
#FLAME_LIB  := $(HOME)/lib/libflame-x86_64-5.1.0-36.a
#FLAME_INC  := $(HOME)/include
BLAS_LIB   := $(HOME)/blis/lib/libblis.a
FLAME_LIB  := $(HOME)/libflame/lib/x86_64-unknown-linux-gnu/libflame.a
FLAME_INC  := $(HOME)/libflame/include

# indicate where the object files are to be created
CC         := gcc
LINKER     := $(CC)
CFLAGS     := -O0 -g -Wall -I$(FLAME_INC) -I$(CBLAS_INC) -m64 -msse3 
FFLAGS     := $(CFLAGS) 

# set the range of experiments to be performed
NREPEATS   := 3#       number of times each experiment is repeated.
NFIRST     := 100#     smallest size to be timed.
NLAST      := 1500#    largest size to be timed.
NINC       := 100#     increment between sizes.

LDFLAGS    := -m64 -lm -fopenmp

# indicate where the FLAME include files reside

# c files
C_FILES=$(wildcard *.c)
C_O_FILES=$(subst .c,.o,$(C_FILES))

# add files in which your implementations are coded
TEST_OBJS  := $(C_O_FILES)

# $%.o: %.f
#	$(FC) $(FFLAGS) -c $< -o $@

driver.x: $(TEST_OBJS) 
	$(LINKER) $(TEST_OBJS) $(FLAME_LIB) $(BLAS_LIB) $(LDFLAGS) -o driver.x

test:   driver.x
	echo "$(NREPEATS) $(NFIRST) $(NLAST) $(NINC)" | ./driver.x > output.m

$(C_O_FILES) : %.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f *.o *~ core *.x output.m
