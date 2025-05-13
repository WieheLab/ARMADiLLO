#ARMADiLLO make

cxx=/usr/bin/g++
flags=-ggdb -O3 -std=c++11

ifneq (,$(findstring dcc-,$(shell uname -n)))#links and libs for duke cluster
	links=-I/datacommons/dhvi/scripts/lib/boost/boost_1_70_0/ -L/datacommons/dhvi/scripts/lib/boost/boost_1_70_0/stage/lib/ -lboost_filesystem -lboost_system -lboost_serialization -pthread
	libs= /datacommons/dhvi/scripts/lib/boost/boost_1_70_0/stage/lib/libboost_serialization.a
else ifneq (,$(findstring Darwin,$(shell uname -s)))#links and libs for osx
	links=-I/opt/homebrew/Cellar/boost/1.87.0/include/ 
	libs=-L/opt/homebrew/Cellar/boost/1.87.0/lib -lboost_filesystem -lboost_system -lboost_serialization -pthread /opt/homebrew/Cellar/boost/1.87.0/lib/libboost_serialization.a
else #links for linux
	links=-L/usr/lib/x86_64-linux-gnu/ 
	libs=-lboost_filesystem -lboost_system -lboost_serialization -pthread /usr/lib/x86_64-linux-gnu/libboost_serialization.a
endif

all: ARMADiLLO_main.o HTML.o utilities.o readInputFiles.o license.o #linking of ARMADiLLO
	${cxx} ${flags} ${links} -o ARMADiLLO ARMADiLLO_main.o HTML.o utilities.o readInputFiles.o license.o ${libs}

ARMADiLLO: ARMADiLLO_main.o HTML.o utilities.o readInputFiles.o license.o #linking of ARMADiLLO
	${cxx} ${flags} ${links} -o ARMADiLLO ARMADiLLO_main.o HTML.o utilities.o readInputFiles.o license.o ${libs}

ARMADiLLO_main.o: ARMADiLLO_main.cpp ARMADiLLO_main.hpp HTML.hpp utilities.hpp nab.hpp readInputFiles.hpp license.hpp
	${cxx} ${flags} ${links} -c ARMADiLLO_main.cpp

readInputFiles.o: readInputFiles.cpp readInputFiles.hpp utilities.hpp
	${cxx} ${flags} ${links} -c readInputFiles.cpp

utilities.o: utilities.cpp utilities.hpp
	${cxx} ${flags} ${links} -c utilities.cpp

HTML.o: HTML.cpp HTML.hpp
	${cxx} ${flags} ${links} -c HTML.cpp

generate_license: Generate_license.cpp license.hpp license.o
	@if [ -f Generate_license.cpp ]; then \
		echo "${cxx} ${flags} ${links} Generate_license.cpp license.o -o generate_license" ; \
		${cxx} ${flags} ${links} Generate_license.cpp license.o -o generate_license; \
	fi

FILES := $(filter-out license.o, $(wildcard license*))
license.o : $(FILES)
	@if [ -f license.cpp ]; then \
		echo "${cxx} ${flags} ${links} -c license.cpp"; \
		${cxx} ${flags} ${links} -c license.cpp; \
	else \
		touch license.o;\
	fi

clean:
	@if [ -f license.cpp ]; then \
		echo "rm -f *.o ARMADiLLO" ;\
		rm -f *.o ARMADiLLO ; \
	else \
		echo "rm HTML.o utilities.o readInputFiles.o ARMADiLLO_main.o ARMADiLLO" ;\
		rm HTML.o utilities.o readInputFiles.o ARMADiLLO_main.o ARMADiLLO ; \
	fi

#END OF MAKEFILE
