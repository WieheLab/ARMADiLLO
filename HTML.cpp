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
