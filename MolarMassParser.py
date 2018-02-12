# -*- coding: utf-8 -*-
"""
This Python script is intended to scrape a given SDF file using rdkit, and 
utilize the "ExactMolWt" functionality of rdkit to get the molecular weight. 

If using a Unix Shell (or equivalent) and wish to call this .py with an input
file and desired output file names, uncomment appropriate lines below "file" 
and "newfile" in code, then run:
    >> python.exe <parser>.py <input filename> <output filename>
...with the "<>" not included, and the appropriate name of the python parser 
script included (NOTE: use EXTENSIONS for both filenames as well). 
Path is required in Windows for call to python.exe, unless added to PATH 
environment variable.

IF USING FOR LARGE DATASET (more than say 5 compounds), change "testing" value 
to FALSE, otherwise print statements will run rampant.
"""
#imports the sys module for input to be given in-line of unix shell
import sys
#Pandas is used in tandem with rdkit's PandasTools to load and SDF into a frame
import pandas as pd
from rdkit import Chem
from rdkit.Chem import PandasTools
from rdkit.Chem import Descriptors

#Trigger print statements with True setting (NOTE: change to a sys arg?)
testing = False #testing values (e.g., printing portions of dataframe)
helpful = True #helpful prints statement (input file, output file, etc...)

#Variables changed upon input type, controls later if/else statement for output
givenargs = False #2 args not given
onearg = False #1 arg not given
defofile = False #default ofile not in use

if testing or helpful == True:
    print("Number of arguments given: ", len(sys.argv), "\n")
##Requests the file to be accessed, takes arguments, or selects default file
#assumes that only 1 arg given is input, and requests output later on.
if len(sys.argv) == 2:
    file = str(sys.argv[1])
    #below response either a filename.csv, 'no', or enter to continue properly
    ofile = input("Output file name? (no OR enter to skip) \n") 
    print('\n')
    if ofile == 'no' or ofile == '':
        print("Default output file naming scheme in use. \n") #.sdf --> .csv
        defofile = True
    onearg = True
elif len(sys.argv) == 3:
	file = str(sys.argv[1])
	givenargs = True
elif len(sys.argv) > 3:
	print("Error, too many arguments. Only 2 allowed. \n")
	for x in sys.argv:
		print("Argument: ", x, "\n")
else:
	print("Using default file input. \n")
	file = 'IC_non-excl_89191.sdf'

if helpful == True:
    print("This is your input file name: ", file, "\n")
    
#Example testing variable in RDkit LoadSDF. Could just write it into
#the parameters as well.
getname = "Molecule"

if helpful == True:
    print("Now loading SDF file through RDKit... \n")
    
#Loads the SDF into a Pandas dataframe. smilesName and molColName are 
#the labels for SMILES and the rdkit molecule object number
frame = PandasTools.LoadSDF(file,
                            smilesName='SMILES',
                            molColName=getname,
                            includeFingerprints=True)

if testing == True:
	#Test print: prints index 0 (i.e. row 1) at column labeled "SMILES", giving
	#the desired SMILES label for the molecule, which can be used to calculate
	#molecular weight.
	print(frame.loc[0,'SMILES'])

#Storing the smiles and names in separate variables/arrays, with empty mw set
smiles = frame.loc[:,'SMILES']
#name = frame.loc[:,'name']

if testing == True:
	#Printing for visual testing
	print('\n')
	print("SMILES log", '\n', smiles)
    #commented out due to some datasets not having 'name' column, caused error
#	print('\n')
#	print("Compound Name log", '\n', name)

	##Example of printing the smiles value at a particular index (row)
	#print(smiles[1])
	print('\n')

#storage vector of molecular masses
molmass = []

#Loops through all molecules (in SMILES format) and calculates the mass, 
#storing it in either a dict or list (see above the index 'i').
for item in smiles:
    mm = Descriptors.ExactMolWt(Chem.MolFromSmiles(item))
    molmass.append(mm)

##uses pandas append to frame (does it?)
#frame2 = frame.append(molmass, ignore_index=True)
    
#Name desired for column header of molecular mass
mmtitle = 'molmass'

if helpful == True:
    print("Column title for Molecular Masses: ", mmtitle, "\n")
    
#appends molecular mass column to dataframe from RDkit
frame[mmtitle] = molmass

#sort dataframe by molecular mass
sframe = frame.sort_values(mmtitle)

##Prints the desired column of information from frame. Replace 'molmass' with 
##any other column name to get it printed out.
if testing == True:
    print(frame.molmass)

#obtains output file name if given, OR generates it from input file name 
#(either the defaulted one or given one)
if givenargs == True:
    ofile = str(sys.argv[2])
elif onearg == True and defofile == False:
    ofile = ofile #use given ofile from question prompt earlier
else:
    #for any file extension, this will replace them with
    #.csv for file naming below. -> occurs if defofile = True OR at full 
    #default settings.
    splitfile = file.split(sep='.') #cuts off extension
    ofile = str(splitfile[0]) + '_sorted.csv' #replaces with .csv

#Print the file name generated in if statement above, test check
if testing == True or helpful == True:
    print("This is your output file: ", ofile, "\n")
    print("Now writing to CSV file... \n")
#converts dataframe from SDF file to a CSV file
sframe.to_csv(ofile, na_rep='0')