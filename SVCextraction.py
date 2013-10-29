#!/usr/bin/env python
#########################################################################################
# SVCextraction.py extract data from .sig file (from SVC field-portable spectroradiometer)
# Note! There should be a directory name SVCreflectances within the main directory (SVC_spectra_date).
# Inputs:
# spectralist.txt: contains path to each .sig file as rows
# Outputs:
# Ex: L820130919R219T224spectra.txt: file containing only the columns (Only numbers)
# Created by: Javier A. Concha
# Rochester Institute of Technology
# 10/25/13
#########################################################################################
import os

date = raw_input("Enter the date (e.g. 130919):")

dirname = ''.join(('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/',date,'/SVC_spectra_',date,'/'))

# Create spectralist.txt
os.chdir(dirname)
outFile = open('./spectraList.txt', "w")
for files in os.listdir("."):
    if files.endswith(".sig"):
        outFile.writelines(files + '\n')
outFile.close()

# Create directory SVCreflectances
mypath = dirname + 'SVCreflectances/'
if not os.path.isdir(mypath):
   os.makedirs(mypath)

# Extract data
SpectraList = open(''.join((dirname,'spectralist.txt')))

for line in SpectraList:
    filename = line.replace('_','')
    filename = filename.replace(".sig","").strip()

    filepath = ''.join((dirname,line))
    inFile = open(filepath.strip())

    outFileName = ''.join((mypath,filename,'spectra.txt')).strip()
    outFile = open(outFileName, "w")

    key = "data="
    key_found = False

    for line2 in inFile:
        if key_found:

            # line2_stuff = line2.split() # for picking just some lines
            # outFile.write(' '.join((line2_stuff[0],' ', line2_stuff[3],'\n')))

            outFile.write(line2)

        elif line2.startswith(key):
            key_found = True

# Create SVCSpectraList.txt to be used in Matlab
os.chdir(mypath)
outFile = open('./SVCSpectraList.txt', "w")
for files in os.listdir("."):
    if files.endswith(".txt") and files.startswith("L"):
        outFile.writelines(files + '\n')
outFile.close()

# Close files
SpectraList.close()
inFile.close()