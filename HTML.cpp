#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <map>
#include "HTML.hpp"

using namespace std;

//Td
Td::Td(string _hclass, string _style, string _id, string _colspan, string _value)
{
  hclass=_hclass;
  style=_style;
  id=_id;
  colspan=_colspan;
  value=_value;
}

void Td::print(string &file_string)
{
  file_string+="  <td";
  if (!hclass.empty()){file_string+=" class=\""+hclass+"\"";}
  if (!style.empty()){file_string+=" style=\""+style+"\"";}
  if (!id.empty()){file_string+=" id=\""+id+"\"";}
  if (!colspan.empty()){file_string+=" colspan=\""+colspan+"\"";}
  if (!rowspan.empty()){file_string+=" rowspan=\""+rowspan+"\"";}
  if (!valign.empty()){file_string+=" valign=\""+valign+"\"";}
  if (!div_class.empty()){file_string+="><div class=\""+div_class+"\"";}
  file_string+=">"+value;
  if (!div_class.empty()){file_string+="</div>";}
  file_string+="</td>\n";
}

//Tr
Tr::Tr(string _hclass, string _style, string _id, vector<Td> _cols)
{
  hclass=_hclass;
  style=_style;
  id=_id;
  cols=_cols;
}

void Tr::print(string &file_string)
{
  file_string+=" <tr";
  if (!hclass.empty()){file_string+=" class=\""+hclass+"\"";}
  if (!style.empty()){file_string+=" style=\""+style+"\"";}
  if (!id.empty()){file_string+=" id=\""+id+"\"";}
  file_string+=">\n";
  for(int i=0; i<cols.size(); i++)
    {
      cols[i].print(file_string);
    }
  file_string+=" </tr>\n";
}

//Table
Table::Table(string _hclass, string _style, string _id, vector<Tr> _rows)
{
  hclass=_hclass;
  style=_style;
  id=_id;
  rows=_rows;
}
void Table::print(string &file_string)
{
  file_string+="<table";
  if (!hclass.empty()){file_string+=" class=\""+hclass+"\"";}
  if (!style.empty()){file_string+=" style=\""+style+"\"";}
  if (!id.empty()){file_string+=" id=\""+id+"\"";}
  file_string+=">";
  for(int i=0; i<rows.size(); i++)
    {
      rows[i].print(file_string);
    }
  file_string+="</table>\n"; 
}

bool writeAMA()
{

  ofstream amaFile("AMA.css");
  if(amaFile.is_open())
    {
      amaFile<< AMA_css;
      amaFile.close();
    }
  else
    return false;
  return true;
}

bool writeColor()
{
  ofstream colorFile("sequence_color.css");
    if(colorFile.is_open())
    {
      colorFile<< sequence_color_css;
      colorFile.close();
    }
    else
      return false;
  return true;
}

bool writeError(string filename, string seqName)
{

  string outString;
  
  outString+="<html xmlns='http://www.w3.org/1999/xhtml' xml:lang='en'>\n";
  outString+="<head>\n";
  outString+="<meta http-equiv='Content-Type' content='text/html; charset=utf-8' />\n";
  outString+="<title>Antibody Mutation Analysis</title>\n";
  outString+="<link rel='stylesheet' href='AMA.css' />\n";
  outString+="</head>\n";
  outString+="<body>\n";
  
  outString+="<br><br><div align=\"center\"><font size=6><b>\n";
  outString+= "Unidentified nucleotide in sequence: "+seqName+" <br>\n";
  outString+="ARMADiLLO cannot accept amino acid sequences or nucleotide sequences with ambuguity codes\n<br>\n";
  outString+="</font>\n</b>\n</div>\n</body>\n</html>";

  ofstream file_out;
  file_out.open(filename.c_str());
  file_out << outString;
  file_out.close();
  return true;
}
