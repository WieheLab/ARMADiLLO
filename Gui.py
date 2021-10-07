#!/usr/bin/python

#libraries
import xlrd, re
import glob
import os
import sys
import string
import PySimpleGUI as sg

def run_armadillo():

def main():
   sg.theme('DarkTeal6')
   
   layout=[ [sg.Text("ARMADiLLO")], [sg.Text('Filename')],
            [sg.Input(), sg.FileBrowse()],
            [sg.Text('enter somthing on Row 2'), sg.InputText()],[sg.FilesBrowse('Select')],[sg.Text('second text lin')],[sg.Button('Run'),sg.Button('OK'),sg.Button('Cancel')] ]
   window=sg.Window('ARMADiLLO GUI',layout)
   while True:
      event, values=window.read()
      if event == sg.WIN_CLOSED or event == 'Cancel':
         break
      print('you enter ',values[0])
   window.close()

if __name__=="__main__":
   main()

if [ -f ${UCA_FILE} ]
then
    $EXECUTABLE_LOCATION/$EXECUTABLE -seq ${FASTA_FILE} \
	-uca ${UCA_FILE} \
	-w $LINE_WRAP_LENGTH \
	-m $EXECUTABLE_LOCATION/Mutability.csv \
	-s $EXECUTABLE_LOCATION/Substitution.csv \
	-max_iter $MAX_ITER \
	-chain $CHAIN \
	-species $SPECIES 	
else
    $EXECUTABLE_LOCATION/$EXECUTABLE -SMUA ${FASTA_FILE/fasta/SimpleMarkedUAs.fasta} \
	-w $LINE_WRAP_LENGTH \
	-m $EXECUTABLE_LOCATION/Mutability.csv \
	-s $EXECUTABLE_LOCATION/Substitution.csv \
	-max_iter $MAX_ITER \
	-chain $CHAIN \
	-species $SPECIES
    
fi
