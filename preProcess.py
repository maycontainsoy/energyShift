# Code to clean up PHF and FCI output for use in Fortran for fit
import numpy             as np
import fileinput         as fi
import linecache         as lc
import matplotlib.pyplot as plt
import glob
import os

print('--------------------------------------------------------------')
print('                    Program: preProcess.py                    ')
print('                 Created by: Stephanie Lauber                 ')
print('                    Created August 2020                       ')
print('')
print(' preProcess.py takes .res files from PHF (LAMP) and FCI       ')
print(' (Bigstick) plus additional FCI files for increased M values  ')
print(' and matches states of angular momentum to create an output   ')
print(' file for use in rmsFortran.                                  ')
print('')
print(' Input for .res files does not require the .res to be added   ')
print(' by the user. This code assumes the additional FCI files for  ')
print(' higher M values have the form SsXX_FCI_*.res where Ss is the ')
print(' abbreviated name of the nucleus and XX is the total A.       ')
print('')
print(' *Currently recommended for use by CWJ nuclear research group ')
print(' only, not current set to work with other PHF or FCI codes.   ')
print('--------------------------------------------------------------')
print('')

# Gather input information from user
nucName      = input('Nucleus name (abbreviation):         ')
totalA       = int(input('Total A:                             '))
lampFileName = input('Enter PHF .res file:                 ') + '.res'
fciFileName  = input('Enter FCI .res file:                 ') + '.res'

multiFCI =     input('Are there multiple FCI files? (y/n): ')

if multiFCI == 'y':
    checkForFCI  = nucName + string(totalA) + 'FCI_*.res'
else:
    print('Not checking for additional FCI files.')

saveName     = input('Enter file output name:              ')

# Get information from PHF .res file output from LAMP
lampFileTemp = fi.input( lampFileName )

# Finds start and stop of energies in LAMP .res
# Assumes start of the energies output by LAMP are 2 lines below
# the labels for State, Energy, etc.
# Assumes the end of the energies output by LAMP are 1 line above
# the information on fractional occupation of states
for line in lampFileTemp:
    if line.find('  State') != -1:
        startLineLAMP = fi.lineno() + 2
    if line.find ('  Fraction') != -1:
        endLineLAMP = fi.lineno() - 1
        break

lampFileTemp.close()

# Calculate total number of LAMP states
numLampStates = int(endLineLAMP - startLineLAMP)

# Create array for LAMP energies
lampInfo = np.zeros(( numLampStates, 4 ))

# Populate lampInfo array with information from LAMP .res file
for i in range(numLampStates):
    lampInfo[i,0] = float(lc.getline( lampFileName, startLineLAMP + i ).split()[0]) # State number
    lampInfo[i,1] = float(lc.getline( lampFileName, startLineLAMP + i ).split()[1]) # Absolute energy
    lampInfo[i,2] = float(lc.getline( lampFileName, startLineLAMP + i ).split()[2]) # Excitation energy
    lampInfo[i,3] = float(lc.getline( lampFileName, startLineLAMP + i ).split()[3]) # Angular momentum

# Get information from FCI .res output from Bigstick
fciFileTemp = fi.input( fciFileName )

# Find start and stop of energies in Bigstick .res
# Assumes start of the energies output by Bigstick are 1 line below
# the labels for State, Energy, etc.
# Assumes the end of the energies output by Bigstick are 1 line above
# the information on total time to run
for line in fciFileTemp:
    if line.find('  State') != -1:
        startLineFCI = fi.lineno() + 1
    if line.find('  Total') != -1:
        endLineFCI = fi.lineno() - 1

fciFileTemp.close()

# Calculate total number of energies in Bigstick .res
numFCIStates = endLineFCI - startLineFCI

# Create array for FCI energies
fciInfo = np.zeros(( numFCIStates, 4 ))

# Populate fciInfo array with information from Bigstick .res file
for i in range(numFCIStates):
    fciInfo[i,0] = float(lc.getline( fciFileName, startLineFCI + i ).split()[0]) # State number
    fciInfo[i,1] = float(lc.getline( fciFileName, startLineFCI + i ).split()[1]) # Absolute energy
    fciInfo[i,2] = float(lc.getline( fciFileName, startLineFCI + i ).split()[2]) # Excitation energy
    fciInfo[i,3] = float(lc.getline( fciFileName, startLineFCI + i ).split()[3]) # Angular momentum

# Read in multiple .res files from Bigstick for increasing M
fileNameList = []

# Move to current working directory, can change if files are located elsewhere
# Addes all files that satisfy SsXX_FCI_*.res to list
if multiFCI == 'y':
    os.chdir('./')
    for file in glob.glob(checkForFCI):
        fileNameList.append(file)

    # For each additional FCI file in list
    if len(fileNameList) > 0:
        for i in range(len(fileNameList)):
            fciFileName = fileNameList[i]
            fciFileTemp = fi.input( fciFileName )

            # Finds start and stop of energies in Bigstick .res
            for line in fciFileTemp:
                if line.find('  State ') !=-1:
                    startLineFCI = fi.lineno() + 1
                if line.find('  Total') !=-1:
                    endLineFCI   = fi.lineno() - 2
                    break

            fciFileTemp.close()

            # Total number of FCI states
            numFCIStates = int((   lc.getline( fciFileName, endLineFCI )).split()[0])

            fciInfoAdd = np.zeros(( numFCIStates, 4 ))

            for i in range(numFCIStates):
                fciInfoAdd[i,0] = float(lc.getline( fciFileName, startLineFCI + i ).split()[0])
                fciInfoAdd[i,1] = float(lc.getline( fciFileName, startLineFCI + i ).split()[1])
                fciInfoAdd[i,2] = float(lc.getline( fciFileName, startLineFCI + i ).split()[2])
                fciInfoAdd[i,3] = float(lc.getline( fciFileName, startLineFCI + i ).split()[3])

            # Addes new FCI information to original fciInfo array
            fciInfo = np.concatenate((fciInfo,fciInfoAdd))

        i = 0
        j = 0

        # Checks for rows with identical absolute energy values
        while i in range(len(fciInfo)):
            j = i + 1
            while j in range(len(fciInfo)):
                if fciInfo[j,1] == fciInfo[i,1]:
                    fciInfo = np.delete(fciInfo,j,axis=0)
                j = j + 1
            i = i + 1

        # Recalcualte excitation energies based on absolute ground state
        # energy from lowest M value run (currently assumes that the first
        # file read from the FCI set is the lowest)
        for i in range(len(fciInfo)):
            fciInfo[i,2] = np.absolute(fciInfo[i,1]-fciInfo[0,1])

else:
    for i in range(len(fciInfo)):
        fciInfo[i,2] = np.absolute(fciInfo[i,1]-fciInfo[0,1])

# Sort arrays by angular momentum and absolute energy
lampSorted = lampInfo[np.lexsort((lampInfo[:,1], lampInfo[:,-1]))]
fciSorted  = fciInfo[np.lexsort(( fciInfo[:,1],  fciInfo[:,-1]))]

# Find common angular momentum values
commonAngMom = np.intersect1d( lampSorted[:,-1], fciSorted[:,-1])
diffAngMom   = np.setxor1d(    lampSorted[:,-1], fciSorted[:,-1])

# Mask values of differing angular momentum in both arrays
lampKeep = lampSorted[~np.in1d(lampSorted[:,-1], diffAngMom)]
fciKeep  = fciSorted[ ~np.in1d(fciSorted[:,-1],  diffAngMom)]

# Count instances of angular momenta for each array
uniqueLAMP, countsLAMP = np.unique(lampKeep[:,-1], return_counts=True)
uniqueFCI,  countsFCI  = np.unique(fciKeep[:,-1],  return_counts=True)

# Find the minimum number of instances for each shared value of angular
# momentum between FCI and PHF data sets
angMomKeep = np.minimum(countsLAMP,countsFCI)

angMomInst = np.zeros(len(angMomKeep))

for i in range(len(angMomKeep)):
    angMomInst[i] = np.minimum(countsLAMP[i],countsFCI[i])

idxArray = np.zeros((len(angMomInst),2))

for i in range(len(idxArray)):
    for j in range(len(lampKeep)):
        element = uniqueLAMP[i]
        if lampKeep[j,-1] == element:
            idxArray[i,0] = j
            break

for i in range(len(idxArray)):
    for j in range(len(fciKeep)):
        element = uniqueFCI[i]
        if fciKeep[j,-1] == element:
            idxArray[i,1] = j
            break

# Adds 10 additional rows to array as buffer (from Ti error)
# Not the most elegant fix but works
dimTemp = int(np.sum(angMomInst))+30
keepValues = np.zeros((dimTemp,5))

# Create array with pairs of energies from PHF and FCI
for i in range(len(idxArray)):
    keepValues[int(idxArray[i,0]):int(idxArray[i,0]+angMomInst[i]),0:2] = \
        lampKeep[int(idxArray[i,0]):int(idxArray[i,0]+angMomInst[i]),1:3]
    keepValues[int(idxArray[i,0]):int(idxArray[i,0]+angMomInst[i]),2:4] = \
        fciKeep[int(idxArray[i,1]):int(idxArray[i,1]+angMomInst[i]),1:3]
    keepValues[int(idxArray[i,0]):int(idxArray[i,0]+angMomInst[i]),-1]  = \
        fciKeep[int(idxArray[i,1]):int(idxArray[i,1]+angMomInst[i]),-1]

# Remove any used rows from buffer
keepValues = keepValues[~np.all(keepValues == 0, axis=1)]

# Save results to file for use in rmsFortran
np.savetxt(saveName,keepValues)
