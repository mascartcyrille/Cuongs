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
else
ifeq ($(host),PcNice)
armaPATH= /usr/local/armadillo-9.200.7/usr/
armalib= $(armaPATH)/lib/x86_64-linux-gnu
else
ifeq (,$(findstring $(host), creos))
armaPATH= /users/scarella/my_install/armadillo-9.200.7/usr/
armalib= $(armaPATH)/lib/x86_64-linux-gnu
endif
endif
endif
endif


##all: testV2.exe test11.exe pw.exe computemu.exe test1.exe computeGpw.exe computeGk1k2.exe test_series.exe test3.exe test31.exe test_DataSpike.exe test_DataNeur.exe test_series_cas9.exe
#all: testV2.exe test11.exe pw.exe computemu.exe test1.exe computeGpw.exe computeGk1k2.exe test_series.exe test3.exe test_DataSpike.exe test_DataNeur.exe test_series_cas9.exe mytest3.exe test9_750.exe
all: testV2.exe test11.exe pw.exe computemu.exe test1.exe computeGpw.exe computeGk1k2.exe test_series.exe test3.exe test_DataSpike.exe test_DataNeur.exe test_series_cas9.exe test9_750.exe test9_125.exe test9_short_750.exe test_series_short_cas10.exe test10_1500.exe test3_full.exe test1_full.exe test11_full.exe test_series_short_cas10_2.exe test_series_short_cas11_2.exe testDN1.exe test_series_short_cas11_4.exe test_lassoshooting.exe

computemu.exe: computemu_V2.o misc.o PwcFunc.o computeb_V2.o test_computemu.o
	g++ $(CXXFLAGS) $^ -o $@ -L$(armalib) -larmadillo

computeGpw.exe: computeGpw_V2.o misc.o PwcFunc.o computemu_V2.o test_computeGpw.o computeG_V2.o computeb_V2.o
	g++ $(CXXFLAGS) $^ -o $@ -L$(armalib) -larmadillo

computeGk1k2.exe: computeGk1k2_V2.o misc.o PwcFunc.o computemu_V2.o computeG_V2.o computeb_V2.o test_computeGk1k2.o
	g++ $(CXXFLAGS) $^ -o $@ -L$(armalib) -larmadillo

testV2.exe: testV2.o computeG_V2.o computeb_V2.o misc.o DataNeur.o
	g++ $(CXXFLAGS) $^ -o $@ -L$(armalib) -larmadillo

testDN1.exe: testDN1.o computeG_V2.o computeb_V2.o misc.o computeGk1k2_V2.o DataNeur.o run_short_test.o
	g++ $(CXXFLAGS) $^ -o $@ -L$(armalib) -larmadillo

test11.exe: test11.o computeG_V2.o computeb_V2.o misc.o computeGk1k2_V2.o DataNeur.o run_short_test.o
	g++ $(CXXFLAGS) $^ -o $@ -L$(armalib) -larmadillo

test11_full.exe: test11_full.o computeG_V2.o computeb_V2.o misc.o computeGk1k2_V2.o DataNeur.o run_test.o run_short_test.o computemu_V2.o PwcFunc.o computeGpw_V2.o 
	g++ $(CXXFLAGS) $^ -o $@ -L$(armalib) -larmadillo

test3.exe: test3.o computeG_V2.o computeb_V2.o misc.o computeGk1k2_V2.o DataNeur.o run_short_test.o
	g++ $(CXXFLAGS) $^ -o $@ -L$(armalib) -larmadillo

test3_full.exe: test3_full.o computeG_V2.o computeb_V2.o misc.o computeGk1k2_V2.o DataNeur.o run_test.o run_short_test.o computemu_V2.o PwcFunc.o computeGpw_V2.o 
	g++ $(CXXFLAGS) $^ -o $@ -L$(armalib) -larmadillo

#test31.exe: test31.o computeG_V2.o computeb_V2.o misc.o computemu_V2.o PwcFunc.o computeGpw_V2.o computeGk1k2_V2.o DataNeur.o run_test.o run_short_test.o
#	g++ $(CXXFLAGS) $^ -o $@ -L$(armalib) -larmadillo

test_series.exe: test_series.o computeG_V2.o computeb_V2.o misc.o computemu_V2.o PwcFunc.o computeGpw_V2.o computeGk1k2_V2.o DataNeur.o run_test.o run_short_test.o
	g++ $(CXXFLAGS) $^ -o $@ -L$(armalib) -larmadillo

test_series_cas9.exe: test_series_cas9.o computeG_V2.o computeb_V2.o misc.o computemu_V2.o PwcFunc.o computeGpw_V2.o computeGk1k2_V2.o DataNeur.o run_test.o run_short_test.o
	g++ $(CXXFLAGS) $^ -o $@ -L$(armalib) -larmadillo

test1.exe: test1.o computeG_V2.o computeb_V2.o misc.o computeGk1k2_V2.o DataNeur.o run_short_test.o
	g++ $(CXXFLAGS) $^ -o $@ -L$(armalib) -larmadillo

test1_full.exe: test1_full.o computeG_V2.o computeb_V2.o misc.o computeGk1k2_V2.o DataNeur.o run_test.o run_short_test.o computemu_V2.o PwcFunc.o computeGpw_V2.o 
	g++ $(CXXFLAGS) $^ -o $@ -L$(armalib) -larmadillo

test_DataSpike.exe: test_DataSpike.o DataNeur.o misc.o
	g++ $(CXXFLAGS) $^ -o $@ -L$(armalib) -larmadillo

test_DataNeur.exe: test_DataNeur.o DataNeur.o misc.o
	g++ $(CXXFLAGS) $^ -o $@ -L$(armalib) -larmadillo

#mytest3.exe: computeG_V2.o computeb_V2.o misc.o computemu_V2.o PwcFunc.o computeGpw_V2.o computeGk1k2_V2.o DataNeur.o run_test.o run_short_test.o mytest3.o
#	g++ $(CXXFLAGS) $^ -o $@ -L$(armalib) -larmadillo

test9_750.exe: computeG_V2.o computeb_V2.o misc.o computemu_V2.o PwcFunc.o computeGpw_V2.o computeGk1k2_V2.o DataNeur.o run_test.o run_short_test.o test9_750.o
	g++ $(CXXFLAGS) $^ -o $@ -L$(armalib) -larmadillo

test10_1500.exe: computeG_V2.o computeb_V2.o misc.o computemu_V2.o PwcFunc.o computeGpw_V2.o computeGk1k2_V2.o DataNeur.o run_test.o run_short_test.o test10_1500.o
	g++ $(CXXFLAGS) $^ -o $@ -L$(armalib) -larmadillo

test9_short_750.exe: computeG_V2.o computeb_V2.o misc.o computeGk1k2_V2.o DataNeur.o run_short_test.o test9_short_750.o
	g++ $(CXXFLAGS) $^ -o $@ -L$(armalib) -larmadillo

test_series_short_cas10.exe: test_series_short_cas10.o computeG_V2.o computeb_V2.o misc.o computeGk1k2_V2.o DataNeur.o run_short_test.o
	g++ $(CXXFLAGS) $^ -o $@ -L$(armalib) -larmadillo

test_series_short_cas10_2.exe: test_series_short_cas10_2.o computeG_V2.o computeb_V2.o misc.o computeGk1k2_V2.o DataNeur.o run_short_test.o
	g++ $(CXXFLAGS) $^ -o $@ -L$(armalib) -larmadillo

test_series_short_cas11_2.exe: test_series_short_cas11_2.o computeG_V2.o computeb_V2.o misc.o computeGk1k2_V2.o DataNeur.o run_short_test.o
	g++ $(CXXFLAGS) $^ -o $@ -L$(armalib) -larmadillo

test_series_short_cas11_3.exe: test_series_short_cas11_3.o computeG_V2.o computeb_V2.o misc.o computeGk1k2_V2.o DataNeur.o run_short_test.o
	g++ $(CXXFLAGS) $^ -o $@ -L$(armalib) -larmadillo

test_series_short_cas11_4.exe: test_series_short_cas11_4.o computeG_V2.o computeb_V2.o misc.o computeGk1k2_V2.o DataNeur.o run_short_test.o
	g++ $(CXXFLAGS) $^ -o $@ -L$(armalib) -larmadillo

test9_125.exe: computeG_V2.o computeb_V2.o misc.o computemu_V2.o PwcFunc.o computeGpw_V2.o computeGk1k2_V2.o DataNeur.o run_test.o run_short_test.o test9_125.o
	g++ $(CXXFLAGS) $^ -o $@ -L$(armalib) -larmadillo

test_DataSpike.o: test_DataSpike.cpp DataSpike.hpp 
	g++ $(CXXFLAGS) -c -I$(armaPATH)/include -I. $< -o $@

test_lassoshooting.exe: test_lassoshooting.o lassoshooting.o
	g++ $(CXXFLAGS) $^ -o $@ -L$(armalib) -larmadillo

%.o: %.cpp
	g++ $(CXXFLAGS) -c -I$(armaPATH)/include -I. $< -o $@

clean:
	- rm -f $(wildcard *~ *.o *.exe) 

pw.exe: PwcFunc.o misc.o test_PwcFunc.o
	g++ $(CXXFLAGS) $^ -o $@ -L$(armalib) -larmadillo
