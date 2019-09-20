#ARMADiLLO_quick.make


ifeq (,$(findstring dcc-dhvi,$(shell uname -n)))
	cxx=/usr/bin/g++ -ggdb -O3 -std=c++11 -I /datacommons/dhvi/scripts/lib/boost/boost_1_70_system	-lboost_filesystem -lboost_system -lboost_serialization -pthread
	libs=/datacommons/dhvi/scripts/lib/boost/boost_1_70_system/lib/libboost_serialization.a
else
	cxx=/usr/bin/g++ -ggdb -O3 -std=c++11 -L/usr/lib/x86_64-linux-gnu/ -lboost_filesystem -lboost_system -lboost_serialization -pthread
	libs=/usr/lib/x86_64-linux-gnu/libboost_serialization.a
endif


ARMADiLLO: ARMADiLLO_main.o HTML.o utilities.o
	${cxx} -o ARMADiLLO ARMADiLLO_main.o HTML.o utilities.o ${libs}

ARMADiLLO_main.o: ARMADiLLO_main.cpp HTML.hpp utilities.hpp
	${cxx} -c ARMADiLLO_main.cpp

utilities.o: utilities.cpp utilities.hpp
	${cxx} -c utilities.cpp

HTML.o: HTML.cpp HTML.hpp
	${cxx} -c HTML.cpp

clean:
	rm *.o ARMADiLLO

#END OF MAKEFILE
