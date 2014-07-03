# Requires:
#   BOOST libraries
#   Gnu scientific libraries (GSL)
#   HDF5 library
#   Semi-Markov library

BOOST=/usr/local/boost_1_54_0mt
# Different Boost installations have different suffixes.
# If there is no suffix, use "BOOSTVARIANT=".
BOOSTVARIANT=-mt
SEMIMARKOV=/usr/local/include/semimarkov-0.1
HDF5=/usr/local/hdf5-1.8.11

CXX=g++
OPT=-g -O2
INCLUDES=-I$(SEMIMARKOV) -I. -I$(BOOST)/include -I$(HDF5)/include
LIBS=-L$(BOOST)/lib -L$(HDF5)/lib -lboost_unit_test_framework$(BOOSTVARIANT) \
	-lboost_log_setup$(BOOSTVARIANT) -lboost_log$(BOOSTVARIANT) \
	-lboost_chrono$(BOOSTVARIANT) -lboost_thread$(BOOSTVARIANT) \
	-lboost_date_time$(BOOSTVARIANT) -lboost_filesystem$(BOOSTVARIANT) \
	-lboost_program_options$(BOOSTVARIANT) -lboost_random$(BOOSTVARIANT) \
	-lboost_system$(BOOSTVARIANT) -lhdf5 -lhdf5_hl -lgsl -lgslcblas -lpthread


sirexp: sir_exp.o main.o seasonal.o hdf_file.o
	g++ $(OPT) -fPIC -o sirexp sir_exp.o main.o seasonal.o hdf_file.o $(LIBS)

sir_exp.o: sir_exp.cpp sir_exp.hpp
	g++ sir_exp.cpp -DHAVE_CONFIG_H -std=c++11 -fPIC $(INCLUDES) $(OPT) \
	-c -o sir_exp.o

seasonal.o: seasonal.cpp seasonal.hpp
	g++ seasonal.cpp -DHAVE_CONFIG_H -std=c++11 -fPIC $(INCLUDES) $(OPT) \
	-c -o seasonal.o

hdf_file.o: hdf_file.cpp hdf_file.hpp
	g++ hdf_file.cpp -DHAVE_CONFIG_H -std=c++11 -fPIC $(INCLUDES) $(OPT) \
	-c -o hdf_file.o

main.o: main.cpp
	g++ main.cpp -DHAVE_CONFIG_H -std=c++11 -fPIC $(INCLUDES) $(OPT) \
	-c -o main.o

clean:
	rm -f *.o sirexp
