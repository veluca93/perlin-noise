CPY_TGT=cperlin$(shell python3-config --extension-suffix)

CXX=c++
CXXFLAGS=-O3 -Wall -shared -std=c++14 -fPIC -g
SRCS=perlin-py.cpp

all: ${CPY_TGT}

${CPY_TGT}: ${SRCS} perlin.hpp Makefile
	${CXX} ${CXXFLAGS} $(shell python3 -m pybind11 --includes) ${SRCS} -o $@


clean:
	rm -f ${CPY_TGT}
