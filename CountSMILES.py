# -*- coding: utf-8 -*-
"""
Created on Thu Jan 25 13:24:47 2018

@author: lsprague
"""
#Grab pandas package
import pandas as pd

#Select CSV file to peruse, needs SMILES column to search for atoms
samplesfile = "stackedtogether_samples.csv"
#samplesfile = "stackedtogether.csv"

#Read the requested CSV into a pd DataFrame
df = pd.read_csv(samplesfile)

#Due to datasets each having a different identifier for an ID number, these
#steps swap the resultant 'nan' values in 'idnumber' column with whatever the
#appropriate values are for the molecule's ID that ended up in a different 
#column label.
df.idnumber.fillna(df.IDNUMBER, inplace=True)
df.idnumber.fillna(df.MFCDNUMBER, inplace=True)
df.idnumber.fillna(df.CLNUMBER, inplace=True)

#Store the SMILES column and molmass column as separate lists
smilesList = df['SMILES'].tolist()
massframe = df['molmass'].tolist()
idlist = df['idnumber'].tolist()

#Printing for test reasons
print("The following print statements are test samples ONLY.\n")
print("Printing the example compound's smile ID: ",smilesList[1])
print("\nPrinting the counts for each atom requested in script: ")
countC = []


'''
FOR loops, below, for finding and counting particular atoms of interest.
To add a new atom, create an equivalent for loop, then add the new list
to the countlist at bottom, and the appropriate label in indexlist.
'''

#counts Carbons in each compound
for name in smilesList:
    num = smilesList[smilesList.index(name)].count('C')
    num2 = smilesList[smilesList.index(name)].count('c')
    fnum = num + num2
    countC.append(fnum)
    
print("C: ",countC[1])

countN = []

for name in smilesList:
    num = smilesList[smilesList.index(name)].count('N')
    num2 = smilesList[smilesList.index(name)].count('n')
    fnum = num + num2
    countN.append(fnum)

print("N: ",countN[1])
    
countO = []

for name in smilesList:
    num = smilesList[smilesList.index(name)].count('O')
    num2 = smilesList[smilesList.index(name)].count('o')
    fnum = num + num2
    countO.append(fnum)

print("O: ",countO[1])


countF = []

for name in smilesList:
    num = smilesList[smilesList.index(name)].count('F')
    num2 = smilesList[smilesList.index(name)].count('f')
    fnum = num + num2
    countF.append(fnum)

print("F: ",countF[1])

countP = []

for name in smilesList:
    num = smilesList[smilesList.index(name)].count('P')
    num2 = smilesList[smilesList.index(name)].count('p')
    fnum = num + num2
    countP.append(fnum)

print("P: ",countP[1])

countS = []

for name in smilesList:
    num = smilesList[smilesList.index(name)].count('S')
    num2 = smilesList[smilesList.index(name)].count('s')
    fnum = num + num2
    countS.append(fnum)

print("S: ",countS[1])

countCl = []

for name in smilesList:
    num = smilesList[smilesList.index(name)].count('Cl')
    num2 = smilesList[smilesList.index(name)].count('cl')
    fnum = num + num2
    countCl.append(fnum)

print("Cl: ",countCl[1])

countBr = []

for name in smilesList:
    num = smilesList[smilesList.index(name)].count('Br')
    num2 = smilesList[smilesList.index(name)].count('br')
    fnum = num + num2
    countBr.append(fnum)

print("Br: ",countBr[1])


countI = []

for name in smilesList:
    num = smilesList[smilesList.index(name)].count('I')
    num2 = smilesList[smilesList.index(name)].count('i')
    fnum = num + num2
    countI.append(fnum)

print("I: ",countI[1])

'''
End of FOR loops to search SMILES identifiers list.
Below is the generation of output dataframe to CSV.
'''

#List containing the lists above with atoms counted and stored
countlists = [countC,countN,countO,countF,countP,countS,countCl,countBr,countI,massframe,idlist]
#list for row/index notation, should match with appropriate countlist order
indexlist = ['C','N','O','F','P','S','Cl','Br','I','Mass','ID']
#Puts countlists into a dataframe using the indexlist as index
dfcounts = pd.DataFrame(countlists, index = (indexlist))
#Renames the columns by the SMILES IDs
dfcounts.columns = smilesList
#Names the appropriate axes 
dfcounts.rename_axis('Count').rename_axis('SMILES ID:',axis=0)

#Totals the number of X atom over the entire sample (# Cs, # Ns, etc...)
Total = dfcounts.sum(axis=1)
print("\nTotal:\n",Total)
#new dataframe using the Total, with appropriate column header of "Totals"
Total = pd.DataFrame({'Totals': Total})

#Joins/appends the Total column above with the dfcounts frame
incltotal = dfcounts.join(Total)

#Transposes dataframe, b/c having the atoms on column header made more sense
incT = incltotal.T

#NOTE: totalling caused even the ID number column to be 'totalled' (just 
#concatenated), but has no effect on data itself.

#Outputs to desired CSV file
incT.to_csv('countatoms_T.csv')
#incltotal.to_csv('countatoms_all.csv')