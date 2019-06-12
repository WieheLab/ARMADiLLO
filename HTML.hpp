///HTML.hpp
///Header file for HTML objects

#ifndef HTML
#define HTML

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <map>
using namespace std;

class Td
{
public:
  string hclass;
  string style;
  string id;
  string colspan;
  string value;

  //secondaries (not in constructor)
  string rowspan;
  string valign;
  string div_class;

  Td(){};
  Td(string, string, string, string, string);
  ~Td(){ };
  void print(string &);
};

class Tr
{
public:
  string hclass;
  string style;
  string id;
  vector<Td> cols;
  Tr(){};
  Tr(string, string, string, vector<Td>);
  ~Tr(){ };
  void print(string &);
};

class Table
{
public:
 
  string hclass;
  string style;
  string id;
  vector <Tr> rows;
  
  Table(){};
  Table(string, string, string, vector<Tr>);
  ~Table(){ };
  void print(string &);
  
};

#endif
