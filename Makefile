export CC = mpicxx
export CCFLAGS = -D DIM=2 -O3 -D QUA_ELEMENT -D FORDER=4

#CCFLAGS += -D SHOWTIME
CCFLAGS += -D SHN=32
#CCFLAGS += -D TWO_SOLID_GHOSTS
#CCFLAGS += -D TURBULENCE
CCFLAGS += -D TEMPORAL_REFINE
#CCFLAGS += -D PRINTGHOST
#CCFLAGS += -D DEBUG
#CCFLAGS += -D TIMECOUNTING
#CCFLAGS += -D IMPORT_MESH
#CCFLAGS += -D PASSAGE_ANGLE=0.33069396

mainprog := main.cpp
cmd := main.exe

source := $(wildcard : *.cpp)
objects := $(source:.cpp=.o)
mainobject := $(mainprog:.cpp=.o)
MD := o

%.o:%.cpp
	$(CC) -c $(CCFLAGS) $<

$(cmd) : $(mainobject) $(objects)
	$(CC) $(CCFLAGS) -o $@ $^

$(MD) :
	mkdir $(MD)

all : $(MD) $(cmd) 
	mv *.o o

.PHONY: clean cleanfile
clean:
	rm -rf *.o *.exe o core.* *.err log* *.plt *.out machinefile.* *.dat

cleanfile:
	rm -rf o core.* *.err log* *.plt *.out machinefile.* *.dat
