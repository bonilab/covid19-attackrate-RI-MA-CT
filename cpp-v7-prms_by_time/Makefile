# driver: driver.o rkf.o
# 	lcc  driver.o rkf.o -o driver -lm
# makefile
# Tabs *must* be used for the indentations below;
# spaces cause make syntax errors.

.PHONY: clean

CC=g++
CFLAGS=-std=c++11 -O3 -flto
#CFLAGS=-fast -xO4 -xdepend -xarch=v8plusa -xprefetch -xvector -xunroll=8 -fsimple=2 -xsafe=mem
#LIBS=-lm -lgsl -lgslcblas
LIBS=-lm

odesim: $(wildcard *.cpp *.h)
	$(CC) $(CFLAGS) -c -o rkf.o rkf.cpp 							$(LIBS)
	$(CC) $(CFLAGS) -c -o derivs.o derivs.cpp 						$(LIBS)
	$(CC) $(CFLAGS) -c -o driver.o driver.cpp 						$(LIBS)
	$(CC) $(CFLAGS) -c -o parseargs.o parseargs.cpp 					$(LIBS)
	$(CC) $(CFLAGS) -c -o generate_trajectories.o generate_trajectories.cpp 		$(LIBS)
	$(CC) $(CFLAGS) -c -o time_varying_prms.o time_varying_prms.cpp				$(LIBS)
	$(CC) $(CFLAGS) -c -o time_varying_prms_nonage.o time_varying_prms_nonage.cpp				$(LIBS)
	$(CC) $(CFLAGS) -c -o time_varying_prms_age.o time_varying_prms_age.cpp 				$(LIBS)
	$(CC) $(CFLAGS) -c -o prms.o prms.cpp 						$(LIBS)
	$(CC) $(CFLAGS) -o $@ driver.o derivs.o rkf.o generate_trajectories.o time_varying_prms.o time_varying_prms_nonage.o time_varying_prms_age.o prms.o parseargs.o	$(LIBS)

test: $(wildcard *.cpp *.h)
	$(CC) $(CFLAGS) -c -o rkf.o rkf.cpp 							$(LIBS)
	$(CC) $(CFLAGS) -c -o derivs.o derivs.cpp 						$(LIBS)
	$(CC) $(CFLAGS) -c -o test.o test.cpp 						$(LIBS)
	$(CC) $(CFLAGS) -c -o parseargs.o parseargs.cpp 					$(LIBS)
	$(CC) $(CFLAGS) -c -o generate_trajectories.o generate_trajectories.cpp 		$(LIBS)
	$(CC) $(CFLAGS) -c -o time_varying_prms.o time_varying_prms.cpp				$(LIBS)
	$(CC) $(CFLAGS) -c -o time_varying_prms_nonage.o time_varying_prms_nonage.cpp				$(LIBS)
	$(CC) $(CFLAGS) -c -o time_varying_prms_age.o time_varying_prms_age.cpp				$(LIBS)
	$(CC) $(CFLAGS) -c -o prms.o prms.cpp 						$(LIBS)
	$(CC) $(CFLAGS) -o $@ test.o derivs.o rkf.o generate_trajectories.o time_varying_prms.o time_varying_prms_nonage.o time_varying_prms_age.o prms.o parseargs.o	$(LIBS)

clean:
	rm -f odesim test *.o *~

