# -*- coding: utf-8 -*-
"""
Created on Tue Jan 16 14:07:18 2018

@author: lsprague

Trigger with:
    >python.exe SampleSelect.py <csv filename> <r or d>
NOTE: sys.argv[1] = <csv filename>
"""
import sys
#import random
#import matplotlib.pyplot as plt
import pandas as pd

#Some print triggering values
testing = False
details = True

if testing == True:
    print("Testing mode enabled. Verbose printing will occur. ")
    
#if/elif statements for checking for CSV filename. If none given, requests
#filename. If given as proper argument, will use it.
if len(sys.argv) == 1:
    print("Error: too few arguments, only ", len(sys.argv), "given. ")
    csvfile = input("Please give a CSV filename: ")
elif len(sys.argv) == 2 or len(sys.argv) == 3:
    csvfile = str(sys.argv[1])

#Reading csv file into data frame variable
csvcomps = pd.read_csv(csvfile)
if testing == True:
    print("Header of chosen file: ")
    print("\n",csvcomps.head())

#acquiring number of compounds/rows in csv file
ncs = csvcomps.shape[0]   #number of compounds/rows in csv

#lowest and highest masses in CSV file
lowm =  csvcomps['molmass'].min()   #lowest mass in molmass
highm =  csvcomps['molmass'].max()   #highest mass in molmass
print("\nLowest mass in CSV file is: ",lowm)
print("Highest mass in CSV file is: ",highm)

#optional command of 'r' for request all options/details, or 'd' (or nothing)
#for default settings to be run.
if len(sys.argv) == 3:
    if sys.argv[2] == 'r': #r for request options
        print("Requesting options... \n")
        print("NOTE: all lists should be values separated by commas (no spaces)")
        sampsize = float(input("How many largely separated samples? "))  #sample size desired from csv
        lrangelow =  float(input("\nLowest Mass? "))     #lowest mass in range
        lrangehigh = float(input("\nHighest Mass? ") )  #largest mass in range
        lspace = float(input("\nSeparation? ") )  #min diff in mw for large spacing
        #use a list for the desired small range start and end values
        ssampsize = float(input("\nHow many tightly separated samples (per range, one value)?"))
        srangelow_in = [input("\nList all starting MWs for tight cluster samples. ")]      #range(s) desired for tight spaced samples
        srangelow = list(map(float, srangelow_in[0].split(',')))  #splits up response with commas, into integer list
        srangehi_in = [input("\nList all ending MWs for tight cluster samples. ")]   #max range(s)
        srangehi = list(map(float, srangehi_in[0].split(','))) #splits with commas to integer vals
        sspace = float(input("\nDesired tight cluster spacing? ")) #min diff in molmass for tight spacing
        
    elif sys.argv[2] == 'd':
        print("Using default values... ")
        #All default spacing/sample size, min/max, etc... values are figured
        #out (via reading csv) or denoted here.
        sampsize = 600  #sample size desired from csv
        lrangelow =  150     #lowest mass in range
        lrangehigh = 750  #largest mass in range
        lspace = 1  #min diff in mw for large spacing
        #use a list for the desired small range start and end values
        ssampsize = 100
        srangelow =  [250,350,450]      #range(s) desired for tight spaced samples
        srangehi = [255,355,455]   #max range(s)
        sspace = 0.001 #min diff in molmass for tight spacing
        
else:
    print("\nNo alternate settings given. Default values will be used. \n")
    sampsize = 600  #sample size desired from csv
    lrangelow =  150     #lowest mass in range
    lrangehigh = 750  #largest mass in range
    lspace = 1  #min diff in mw for large spacing
    #use a list for the desired small range start and end values
    ssampsize = 100
    srangelow =  [250,350,450]      #range(s) desired for tight spaced samples
    srangehi = [255,355,455]   #max range(s)
    sspace = 0.001 #min diff in molmass for tight spacing

if details == True:   
    parameters = [ncs,sampsize,lrangelow,lrangehigh,lspace,ssampsize,srangelow,srangehi,sspace]
    print("Parameters are: \n")
    print("# Comps | L.Samp.Size | L.LowMW | L.HighMW | L.Spacing | S.Samp.Size | S.LowMW | S.HighMW | S.Spacing ")
    print(parameters)
    print("\n")

#Empty dataframe for storage of sample compounds
dfstorage = pd.DataFrame()
#print(storage.head())

##random starting position from first 100 compounds, for variety of spacing
#randstart = random.randint(0,100)

def getindex(bound,varchange):
    """Obtains the index in given CSV file for the molecular mass that is 
    greater than or equal to the given bound (i.e., low or high range), so 
    that the program has an index to begin/end its search for samples.
    """
    for index, row in csvcomps.iterrows():
        checkval = (row['molmass'])
        if checkval >= bound:
#            print("Checkval for mass range: ",checkval)
#            print("Index: ", index)
            varchange = (index,checkval) #the index of checkval
            break
        else:
            #print("nope")
            continue
#    print("Varchange: ", varchange)
    if varchange == 0:
        print("No checkval generated. Bound likely out of range. ")
#        print("Returning tuple of (0,0), and allowing IF statements to fix.")
        varchange = (0,0)
    elif len(varchange) == 2:
        print("Index found. No problems identified.")
#    print("Varchange updated (if req'd): ",varchange, "\n")
    return varchange


"""
Obtaining index of lowest mass in range from CSV, for large (spaced) sample
"""
#default row start if either lrangelow is out of range of dataset, or equiv.
rowiters = 0

lowestm = getindex(lrangelow,rowiters)
rowiters = lowestm[0]

if testing == True:
    print("Lowest Mass getindex return vals (index,checkval): ",lowestm)

if rowiters == 0 and lowestm[1] < lrangelow:
    print("\n!!!!\n")
    print("Error with given lrangelow. Setting rowiters start val to 0. \n")
print("Rowiters value: ",rowiters, "\n")

""" Done """

"""
Obtaining index of highest mass in range from CSV, for large (spaced) sample
"""
lrangehi_index = 0

highestm = getindex(lrangehigh, lrangehi_index)
#NOTE: the "-1" ensures that the index is within the bound, not one row over it
lrangehi_index = highestm[0] - 1

if testing == True:
    print("Highest Mass getindex return vals (index,checkval): ",highestm)

if lrangehi_index == 0 and highestm[1] < lrangehigh:
    print("\n!!!!\n")
    print("Error with given lrangehigh. Value likely larger than mass range of dataset. ")
    lrangehi_index = csvcomps.shape[0]
    print("Continuing with highest index of dataset... --> ",lrangehi_index)
print("Lrangehi_index value: ",lrangehi_index)


""" Done """

#Include the first row in the dataframe
print("\n\nIncluding starting row in sampling... \n")
dfstorage = dfstorage.append(csvcomps.loc[rowiters,:])
print("Done. \n")

""" 
Beginning the loop for sampling the large (spaced) dataset 
"""

#storing the differences calculated
diffsave = []

#iteration counter, gives number of runs in loop/number of samples gathered
#...starts at 1 because I've already included the first sample based on the 
#starting index.
loopiters = 1

#counter to stay at the row after compared row (rowiters)
nxtrow = rowiters + 1    

print("Beginning sampling from dataset... ")

"""
NOTE: Below while loop could DEFINITELY be a defined function... exactly the
same below in tightly clustered section with exchanced variables.
"""
#loops until the sample size is incremented upwards to desired number, and 
#checks the number of rows left so that it doesn't look for a higher indexed
#row than actually exists. Also cuts off at highest desired mass range.
while loopiters <= (sampsize-1) and nxtrow <= ncs-1 and rowiters <= lrangehi_index:
    m1 = csvcomps.loc[rowiters,'molmass']
    m2 = csvcomps.loc[nxtrow,'molmass']
    diff = m2 - m1
    if testing == True:
        print("Loop # = ",loopiters)
        print("Mass 1: ",m1)
        print("Mass 2: ",m2)
        print("Diff = ",diff)
        print("Next Row = ",nxtrow)
        
    if diff >= lspace:
        #append to list, set rowiters to nxtrow, and reset nxtrow
        #i.e., start again at new value, to find next val that fulfills spacing
        if testing == True:
            print("Accepted\n")
        diffsave.append(diff)
        dfstorage = dfstorage.append(csvcomps.loc[nxtrow,:])
        loopiters = loopiters + 1 #says, "i got another sample"
        rowiters = nxtrow
        nxtrow = rowiters + 1
    elif diff < lspace:
        #increment nxtrow by +1
        if testing == True:
            print("Declined\n")
        nxtrow = nxtrow + 1

#NOTE: the sampling still rolls over by 1, but i'm going to leave it...
if loopiters == sampsize:
    print("\nSuccessfully gathered all samples from dataset. ")
    print("Number of samples = ", loopiters, "\n")
elif loopiters != sampsize:
    print("\n!!!!\n")
    print("Complete sample size not reached. Check desired range and spacing.")
    print("Number of samples = ", loopiters, "\n")

""" 
End of large dataset sampling
"""

print("Some details on the large space sampling... \n")
print("Max difference calculated: ", max(diffsave))
print("Max difference in samples: ", dfstorage['molmass'].max() - dfstorage['molmass'].min())
print("Min difference calculated AND in samples: ", min(diffsave))
#print("\nDone. Check file " + ofile, "\n\n \t~~~~~~~~~~~~~~~~~~\n")
print("\nBeginning small data set inclusion with tighter spacing... \n")

"""
Beginning small (spaced) sampling process below....

~~~~~~~~~~~~~~~~~~~~~~~

Identify indices for all starting masses in srangelow
"""
startlist = []
for initials in srangelow:
    tstart = 0
    begin = getindex(initials,tstart)
    startlist.append(begin[0])
#    print(startlist, "\n")
    
print("These are the start indices for the chosen starting masses: ",startlist)
print("\n")


"""
Identify indices for all ending masses in srangehi
"""
endlist = []
for ends in srangehi:
    endpt = 0
    ending = getindex(ends,endpt)
    #"-1" is again to ensure ending row is within the mass boundary
    endlist.append(ending[0] - 1)
#    print(endlist, "\n")
    
print("These are the end indices for the chosen masses: ",endlist,"\n")


#zips start and endpoint together. Each tuple is a (start,end) pair.
rangelist = list(zip(startlist,endlist))
#print((rangelist))
##calling the first tuple's (first range's) starting index in the CSV file.
#print(rangelist[0][0]) #gives 0th tuple's 0th element in zipped list


"""
Loop through each defined range, check if diff is successful, check if 
ALREADY in dfstorage from large selection (to avoid overlap) and do not 
include in new dftight if so.
"""
diff2save = []
dftightstore = pd.DataFrame()
#for loop goes over all generated tuples, i.e. the start/end pairs
for tuples in rangelist:
    loopiters = 1
    nxtrow = loopiters + 1
    rowiters = tuples[0]
    dftightstore = dftightstore.append(csvcomps.loc[rowiters,:])
    srangehi_index = tuples[1]
    #rowiters will now use the starting points under the FOR loop, as tuples[0] 
    while loopiters <= ssampsize and nxtrow <= ncs-1 and rowiters <= srangehi_index:
        m1 = csvcomps.loc[rowiters,'molmass']
        m2 = csvcomps.loc[nxtrow,'molmass']
        diff = m2 - m1
        if testing == True:
            print("Loop # = ",loopiters)
            print("Mass 1: ",m1)
            print("Mass 2: ",m2)
            print("Diff = ",diff)
            print("Next Row = ",nxtrow)
        
        if diff >= sspace:
            #append to list, set rowiters to nxtrow, and reset nxtrow
            #i.e., start again at new value, to find next val that fulfills spacing
            if testing == True:
                print("Accepted\n")
            diff2save.append(diff)
            dftightstore = dftightstore.append(csvcomps.loc[nxtrow,:])
            loopiters = loopiters + 1 #says, "i got another sample"
            rowiters = nxtrow
            nxtrow = rowiters + 1
        elif diff < sspace:
            #increment nxtrow by +1
            if testing == True:
                print("Declined\n")
            nxtrow = nxtrow + 1

print("Tighter sampling finished.")
print("Total new samples (possibly including duplicates): ",len(diff2save), "\n")

print("Some details on the tighter sampling... \n")
print("Max difference calculated: ", max(diff2save))
print("Max difference in samples: ", dftightstore['molmass'].max() - dftightstore['molmass'].min())
print("Min difference calculated AND in samples: ", min(diff2save))


"""
Print data for newly created dftight, AND concatenated final dataframe????
"""
print("\nComparing large and tightly spaced dataframes now...")
print("Removing duplicates... ")

frames = [dfstorage,dftightstore]
cumeframe = pd.concat(frames)
finale = cumeframe.drop_duplicates(subset='molmass',keep='first')

print("Re-sorting the dataframe by molecular mass...")

finalframe = finale.sort_values('molmass')
print("Size of final dataframe: ",finalframe.shape)

splitfile = csvfile.split(sep='.') #cuts off extension
ofile = str(splitfile[0]) + '_samples.csv' #replaces with .csv
finalframe.to_csv(ofile)

print("\nDone. Check file " + ofile, "\n\n \t~~~~~~~~~~~~~~~~~~\n")
""" 
End of small dataset sampling
"""




