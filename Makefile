#ARMADiLLO make

cxx=/usr/bin/g++
flags=-ggdb -O3 -std=c++11

ifneq (,$(findstring dhvi,$(shell uname -n)))#links and libs for duke cluster
	links=-I/datacommons/dhvi/scripts/lib/boost/boost_1_70_0/ -L/datacommons/dhvi/scripts/lib/boost/boost_1_70_0/stage/lib/ -lboost_filesystem -lboost_system -lboost_serialization -pthread
	#links=-I/datacommons/dhvi/scripts/lib/boost/boost_1_70_system/include/ -L/datacommons/dhvi/scripts/lib/boost/boost_1_70_system/lib/ -lboost_filesystem -lboost_system -lboost_serialization -pthread
	libs=/datacommons/dhvi/scripts/lib/boost/boost_1_70_0/stage/lib/libboost_serialization.a
	#libs=/datacommons/dhvi/scripts/lib/boost/boost_1_70_system/lib/libboost_serialization.a
else ifneq (,$(findstring Darwin,$(shell uname -s)))#links and libs for osx
	#links=-I/opt/local/include/  #-L/opt/local/include/ -lboost_filesystem -lboost_system -lboost_serialization -pthread
	links=-I/usr/local/lib
	libs=/usr/local/lib/libboost_serialization-mt.a
	#libs=/opt/local/lib/libboost_serialization-mt.dylib
else #links for linux
	links=-L/usr/lib/x86_64-linux-gnu/ -lboost_filesystem -lboost_system -lboost_serialization -pthread
	libs=/usr/lib/x86_64-linux-gnu/libboost_serialization.a
endif


ARMADiLLO: ARMADiLLO_main.o HTML.o utilities.o readInputFiles.o #linking of ARMADiLLO
	${cxx} ${flags} ${links} -o ARMADiLLO ARMADiLLO_main.o HTML.o utilities.o readInputFiles.o ${libs}

ARMADiLLO_main.o: ARMADiLLO_main.cpp ARMADiLLO_main.hpp HTML.hpp utilities.hpp nab.hpp readInputFiles.hpp
	${cxx} ${flags} ${links} -c ARMADiLLO_main.cpp

readInputFiles.o: readInputFiles.cpp readInputFiles.hpp utilities.hpp
	${cxx} ${flags} ${links} -c readInputFiles.cpp

utilities.o: utilities.cpp utilities.hpp
	${cxx} ${flags} ${links} -c utilities.cpp

HTML.o: HTML.cpp HTML.hpp
	${cxx} ${flags} ${links} -c HTML.cpp

clean: #cleaning out old compiled files
	rm -f *.o ARMADiLLO

#END OF MAKEFILE
