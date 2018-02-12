# -*- coding: utf-8 -*-
"""
Created on Thu Jan 18 15:23:19 2018

@author: lsprague
"""
import pandas as pd

print("Reading in CSV files within script... ")
df1 = pd.read_csv("NC_337553_sorted.csv")
df2 = pd.read_csv("IC_non-excl_89191_sorted.csv")
df3 = pd.read_csv("DC03_219356_sorted.csv")
df4 = pd.read_csv("DC02_350000_sorted.csv")
df5 = pd.read_csv("DC01_350000_sorted.csv")
df6 = pd.read_csv("DC_saltdata_not-available_122033_sorted.csv")
df7 = pd.read_csv("CDI_BBs_81676_sorted.csv")

print("\nDone, compiling into frames list... ")
frames = [df1,df2,df3,df4,df5,df6,df7]

#frames = []
cont = True
#print(cont)

if len(frames) == 0:
	while cont == True:
		csvfile = input("Got a csv file?")
		if csvfile in ('yes','y','Yes','Y'):
			df1 = pd.read_csv(input("Filename?"))
			frames.append(df1)
			print(cont)
		else:
			print("Okay, all set! Here we go! ")
			cont = False
			print(cont)
			break
else:
	print("\nUsing included frames array in script... ")

print("\nDone, now stacking... ")
stacked = pd.concat(frames, axis=0)
print("\nNow sorting... ")
sortstacked = stacked.sort_values('molmass')
print("\nConverting to CSV file... ")
filename = 'stackedtogether.csv'
sortstacked.to_csv(filename, index=False)
print("\nConverted. See file " + filename)
print("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
