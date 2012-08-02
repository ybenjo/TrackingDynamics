# Makefile

CXX=g++-mp-4.7
CXX_FLAGS=-O3 -std=c++11

ttfmm: 
	${CXX} ${CXX_FLAGS} ./src/main.cc -o ttfmm
clean: 
	rm ttfmm 
