CC=g++
INCLUDES = -I /sw/include -I/sw/opt/boost-1_58/include/ -I /home/vskokov/lib/include -I /home/vskokov/lib/include  -I /sw/include -I/sw/opt/boost-1_58/include/ -I/usr/local/include 
LIBS= -L /sw/lib/ -L/home/vskokov/lib/lib -L  /home/vskokov/MC_DiJet/interp2d -L /home/vskokov/lib/lib  -L/usr/local/lib
CFLAGS= -std=c++11 -O2 -lgsl -lgslcblas -lfftw3



v2: 
	g++  $(INCLUDES) $(LIBS) -O3 -w -o v2.x v2_elliptic.cpp zheevc3.c zheevv3.c   -lgsl -lgslcblas -lfftw3 -lm

test: 
	g++  $(INCLUDES) $(LIBS) -w  -O3 -o v2.x v2_elliptic_tests.cpp zheevc3.c zheevv3.c   -lgsl -lgslcblas -lfftw3 -lm


HBT: 
	g++  $(INCLUDES) $(LIBS) -w  -O3 -o v2_fund.x v2_fund.cpp zheevc3.c zheevv3.c   -lgsl -lgslcblas -lfftw3 -lm



all05: 
	g++  $(INCLUDES) $(LIBS) -O3 -o v3p05.x v3_mp_par.cpp zheevc3.c zheevv3.c -fopenmp  -lboost_system -lgsl -lgslcblas -lfftw3 -lm

all: 
	g++  $(INCLUDES) $(LIBS) -O3 -o v3p.x v3_mp_par.cpp zheevc3.c zheevv3.c -fopenmp  -lboost_system -lgsl -lgslcblas -lfftw3 -lm
allq: 
	g++  $(INCLUDES) $(LIBS) -O3 -o v3.x v3_mv.cpp zheevc3.c zheevv3.c  -lgsl -lgslcblas -lfftw3
MV: 
	g++ -O3 -o mv.x JIMWLK_TMD_beta_LATEST_dipole_probability_mv.cpp zheevc3.c zheevv3.c  -lgsl -lgslcblas -lfftw3
Qs:
	g++ -O3 -o jimwlk_qs.x JIMWLK_Qs.cpp zheevc3.c zheevv3.c  -lgsl -lgslcblas -lfftw3
clean:
	rm -rf jimwlk.x
