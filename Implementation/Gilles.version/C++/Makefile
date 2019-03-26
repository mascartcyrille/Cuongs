UNAME = $(shell uname)
host = $(shell hostname)
#CXXFLAGS = -std=c++11 -g
CXXFLAGS = -std=c++11 -g -O2
ifneq (,$(findstring $(UNAME), "Darwin"))
armaPATH= $(HOME)/MyInstall/armadillo-9.200.7/usr/local
armalib= $(armaPATH)/lib
else
ifneq (,$(findstring $(host), math11.unice.fr math5.unice.fr))
armaPATH= /workdir/gscarella/myinstall/armadillo-9.200.7/usr/
armalib= $(armaPATH)/lib64
endif
endif

all: testV2.exe test11.exe computeGk1k2.exe test3.exe test_DataSpike.exe test_DataNeur.exe

computeGk1k2.exe: computeGk1k2_V2.o misc.o computeG_V2.o computeb_V2.o test_computeGk1k2.o
	g++ $(CXXFLAGS) $^ -o $@ -L$(armalib) -larmadillo

testV2.exe: testV2.o computeG_V2.o computeb_V2.o misc.o DataNeur.o
	g++ $(CXXFLAGS) $^ -o $@ -L$(armalib) -larmadillo

test11.exe: test11.o computeG_V2.o computeb_V2.o misc.o computeGk1k2_V2.o DataNeur.o run_short_test.o
	g++ $(CXXFLAGS) $^ -o $@ -L$(armalib) -larmadillo

test3.exe: test3.o computeG_V2.o computeb_V2.o misc.o computeGk1k2_V2.o DataNeur.o run_short_test.o
	g++ $(CXXFLAGS) $^ -o $@ -L$(armalib) -larmadillo

test_DataSpike.exe: test_DataSpike.o DataNeur.o misc.o
	g++ $(CXXFLAGS) $^ -o $@ -L$(armalib) -larmadillo

test_DataNeur.exe: test_DataNeur.o DataNeur.o misc.o
	g++ $(CXXFLAGS) $^ -o $@ -L$(armalib) -larmadillo

test_DataSpike.o: test_DataSpike.cpp DataSpike.hpp 
	g++ $(CXXFLAGS) -c -I$(armaPATH)/include -I. $< -o $@

%.o: %.cpp
	g++ $(CXXFLAGS) -c -I$(armaPATH)/include -I. $< -o $@

clean:
	- rm -f $(wildcard *~ *.o *.exe) 

