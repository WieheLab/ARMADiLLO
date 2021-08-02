#!/usr/bin/python

#libraries
import xlrd, re
import glob
import os
import sys
import PySimpleGUI as sg


def main():
   sg.theme('DarkTeal6')
   layout=[ [sg.Text("ARMADiLLO")],[sg.Text('enter somthing on Row 2'), sg.InputText()],[sg.FilesBrowse('Select')],[sg.Text('second text lin')],[sg.Button('Run'),sg.Button('OK'),sg.Button('Cancel')] ]
   window=sg.Window('ARMADiLLO GUI',layout)
   while True:
      event, values=window.read()
      if event == sg.WIN_CLOSED or event == 'Cancel':
         break
      print('you enter ',values[0])
   window.close()

if __name__=="__main__":
   main()
