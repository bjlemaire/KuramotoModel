# Load the common configuration file
include ./config.mk

objs=rk4_solver.o
src=$(patsubst %.o,%.cc,$(objs))
execs=model1

all:
	$(MAKE) executables

executables: $(execs)

depend: $(src)
	$(cxx) -MM $(src) >Makefile.dep

-include Makefile.dep

%.o: %.cc
	$(cxx) $(cflags) -c $<

model1: Model1.cc rk4_solver.o
	$(cxx) $(cflags) -o $@ $< $(objs)


clean:
	rm -f $(execs) $(objs)

clean-all:
	rm -f $(execs) $(objs)

.PHONY: clean all executables depend
