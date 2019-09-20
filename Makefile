#ARMADiLLO_quick.make

cxx=/usr/bin/g++ -ggdb -O3 -std=c++11 -L/usr/lib/x86_64-linux-gnu/ -lboost_filesystem -lboost_system -lboost_serialization -pthread

libs=/usr/lib/x86_64-linux-gnu/libboost_serialization.a


ARMADiLLO: ARMADiLLO_main.o HTML.o utilities.o
	${cxx} -o ARMADiLLO ARMADiLLO_main.o HTML.o utilities.o /usr/lib/x86_64-linux-gnu/libboost_serialization.a #/usr/lib/x86_64-linux-gnu/libboost_filesystem.a

ARMADiLLO_main.o: ARMADiLLO_main.cpp HTML.hpp utilities.hpp
	${cxx} -c ARMADiLLO_main.cpp #/usr/lib/x86_64-linux-gnu/libboost_serialization.a #/usr/lib/x86_64-linux-gnu/libboost_filesystem.a

utilities.o: utilities.cpp utilities.hpp
	${cxx} -c utilities.cpp

HTML.o: HTML.cpp HTML.hpp
	${cxx} -c HTML.cpp

clean:
	rm *.o ARMADiLLO

#END OF MAKEFILE
