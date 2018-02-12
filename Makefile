CC=g++
INCLUDES = -I /sw/include -I/sw/opt/boost-1_58/include/ -I /home/vskokov/lib/include -I /home/vskokov/lib/include  -I /sw/include -I/sw/opt/boost-1_58/include/ -I/usr/local/include 
LIBS= -L /sw/lib/ -L/home/vskokov/lib/lib -L  /home/vskokov/MC_DiJet/interp2d -L /home/vskokov/lib/lib  -L/usr/local/lib
CFLAGS= -std=c++11 -O2 -lgsl -lgslcblas -lfftw3



v3: 
	/usr/local/Cellar/llvm/5.0.0/bin/clang++ -I /sw/include -I/sw/opt/boost-1_58/include/ -I /home/vskokov/lib/include -I /home/vskokov/lib/include  -I /sw/include -I/sw/opt/boost-1_58/include/ -I/usr/local/include  -L /sw/lib/ -L/home/vskokov/lib/lib -L  /home/vskokov/MC_DiJet/interp2d -L /home/vskokov/lib/lib  -L/usr/local/lib -L /usr/local/Cellar/llvm/5.0.0/lib  -O3 -o v3p.x v3_mp_par.cpp zheevc3.c zheevv3.c -fopenmp  -lboost_system -lgsl -lgslcblas -lfftw3 -lm

v2: 
	g++  $(INCLUDES) $(LIBS) -O3 -w -o v2.x v2_elliptic.cpp zheevc3.c zheevv3.c   -lgsl -lgslcblas -lfftw3 -lm

test: 
	g++  $(INCLUDES) $(LIBS) -w  -O3 -o v2.x v2_elliptic_tests.cpp zheevc3.c zheevv3.c   -lgsl -lgslcblas -lfftw3 -lm


HBT: 
	g++  $(INCLUDES) $(LIBS) -w  -O3 -o v2_fund.x v2_fund.cpp zheevc3.c zheevv3.c   -lgsl -lgslcblas -lfftw3 -lm



all05: 
	g++  $(INCLUDES) $(LIBS) -O3 -o v3p05.x v3_mp_par.cpp zheevc3.c zheevv3.c -fopenmp  -lboost_system -lgsl -lgslcblas -lfftw3 -lm

allD: 
	g++  $(INCLUDES) $(LIBS) -O3 -o v3D.x v3_dilute.cpp zheevc3.c zheevv3.c    -lgsl -lgslcblas -lfftw3 -lm

allDq: 
	g++  $(INCLUDES) $(LIBS) -O3 -o v3D.x v3_dilute.cpp zheevc3.c zheevv3.c  -fopenmp -lboost_system -lgsl -lgslcblas -lfftw3 -lm

all: 
	g++  $(INCLUDES) $(LIBS) -O3 -o v3.x v3_mp_par.cpp zheevc3.c zheevv3.c -fopenmp  -lboost_system    -lgsl -lgslcblas -lfftw3 -lm


corr: 
	g++  $(INCLUDES) $(LIBS) -O3 -o v3_corr.x v3_corr.cpp zheevc3.c zheevv3.c -fopenmp  -lboost_system    -lgsl -lgslcblas -lfftw3 -lm

corrSP: 
	g++  $(INCLUDES) $(LIBS) -O3 -o corr.x v3_corr_SP.cpp zheevc3.c zheevv3.c   -lgsl -lgslcblas -lfftw3 -lm


corrSPNew: 
	g++  $(INCLUDES) $(LIBS) -O3 -o corrNew.x v3_corr_SP.cpp zheevc3.c zheevv3.c   -lgsl -lgslcblas -lfftw3 -lm


corrSPq: 
	g++  $(INCLUDES) $(LIBS) -O3 -o corrq.x v3_corr_SP.cpp zheevc3.c zheevv3.c   -lgsl -lgslcblas -lfftw3 -lm



corrq: 
	g++  $(INCLUDES) $(LIBS) -O3 -o v3_corrq.x v3_corr.cpp zheevc3.c zheevv3.c -fopenmp  -lboost_system    -lgsl -lgslcblas -lfftw3 -lm

allq: 
	g++  $(INCLUDES) $(LIBS) -O3 -o v3.x v3_mv.cpp zheevc3.c zheevv3.c  -lgsl -lgslcblas -lfftw3
MV: 
	g++ -O3 -o mv.x JIMWLK_TMD_beta_LATEST_dipole_probability_mv.cpp zheevc3.c zheevv3.c  -lgsl -lgslcblas -lfftw3
Qs:
	g++ -O3 -o jimwlk_qs.x JIMWLK_Qs.cpp zheevc3.c zheevv3.c  -lgsl -lgslcblas -lfftw3
clean:
	rm -rf jimwlk.x
