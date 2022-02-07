#include <math.h>
#include <cstdlib>
#include <random>
#include <iostream>
#include <ctime>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <set>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <algorithm>
#include <boost/filesystem/operations.hpp>
//#include "boost/filesystem.hpp"
#include <sys/types.h>
#include <dirent.h>
//#include <boost/filesystem/fstream.hpp>
//#define BOOST_FILESYSTEM_NO_DEPRECATED
// include headers that implement a archive in simple text format
using namespace std;
using namespace boost;
using namespace boost::filesystem;

void read_SMUA_file(string, vector<vector<string> > &);

void read_PARTIScsv_file(string, vector<vector<string> > &);
void read_PARTISyaml_file(string, vector<vector<string> > &);

string cleanYAMLline(string);
void cleanSeqs(string &,string &);
void readIndividualFiles(string, string, vector<vector<string> > &);
void readIndividualFiles(string, string, string, vector<vector<string> > &);

