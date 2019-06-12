#ARMADiLLO.make

cxx=/usr/bin/g++ -O3 -std=c++11 -ggdb -I /usr/local/lib/boost_1_59_0/

ARMADiLLO: ARMADiLLO_main.o HTML.o utilities.o
	${cxx} -o ARMADiLLO ARMADiLLO_main.o HTML.o utilities.o

ARMADiLLO_main.o: ARMADiLLO_main.cpp HTML.hpp utilities.hpp
	${cxx} -c ARMADiLLO_main.cpp

utilities.o: utilities.cpp utilities.hpp
	${cxx} -c utilities.cpp

HTML.o: HTML.cpp HTML.hpp
	${cxx} -c HTML.cpp

clean:
	rm *.o ARMADiLLO

#END OF MAKEFILE
