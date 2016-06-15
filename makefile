SDIR = src
IDIR = include
ODIR = obj
INC = -I./include -I/usr/local/eigen

CXX = g++ -g -fopenmp -O3 -std=c++11

LFLAGS = -lgsl -lgslcblas -DHAVE_INLINE -lgomp 

_OBJS = Settings.o Mvngen.o Farm.o Grid.o Params.o Model.o Smc.o multiscale.o
OBJS = $(patsubst %,$(ODIR)/%,$(_OBJS))

$(ODIR)/%.o : $(SDIR)/%.cpp
	$(CXX) -c $(INC) -o $@ $< $(CFLAGS) $(LFLAGS)

disease : $(OBJS)
	$(CXX) $(OBJS) -o disease $(CLFAGS) $(LFLAGS)

clean :
	rm -f disease ./obj/*.o



