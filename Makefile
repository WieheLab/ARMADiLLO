#ARMADiLLO_quick.make

cxx=/usr/bin/g++
flags=-ggdb -O3 -std=c++11

ifeq (,$(findstring $(shell uname -n),dcc-dhvi))
	links=-I /datacommons/dhvi/scripts/lib/boost/boost_1_70_system/ -L /datacommons/dhvi/scripts/lib/boost/boost_1_70_system/lib/ -lboost_filesystem -lboost_system -lboost_serialization -pthread
	libs=/datacommons/dhvi/scripts/lib/boost/boost_1_70_system/lib/libboost_serialization.a
else
	links=-L/usr/lib/x86_64-linux-gnu/ -lboost_filesystem -lboost_system -lboost_serialization -pthread
	libs=/usr/lib/x86_64-linux-gnu/libboost_serialization.a
endif


ARMADiLLO: ARMADiLLO_main.o HTML.o utilities.o
	${cxx} ${flags} ${links} -o ARMADiLLO ARMADiLLO_main.o HTML.o utilities.o ${libs}

ARMADiLLO_main.o: ARMADiLLO_main.cpp HTML.hpp utilities.hpp
	${cxx} ${flags} ${links} -c ARMADiLLO_main.cpp

utilities.o: utilities.cpp utilities.hpp
	${cxx} ${flags} ${links} -c utilities.cpp

HTML.o: HTML.cpp HTML.hpp
	${cxx} ${flags} ${links} -c HTML.cpp

clean:
	rm *.o ARMADiLLO

#END OF MAKEFILE
