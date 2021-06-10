#!/usr/bin/env python

import os, sys, math, glob, mmap
from collections import namedtuple

AtomNumberToSymbol = {1:'H', 6:'C', 7:'N', 8:'O', 9:'F', 16:'S', 17:'Cl', 35:'Br'}

#define a Conformer tuple [name, [freq], energy, Type]
Conformer = namedtuple('Conformer', ['name','freq','energy', 'Type'])

def getXYZMatrix(F):
#F is a file handle.  Rewind for safety
  F.seek(0, 0)
  StartKey = "Standard orientation"
  EndKey = 'Rotational constants'
  Ouput = []
#find last occurence of StartKey
  m = mmap.mmap(F.fileno(), 0, prot=mmap.PROT_READ)  
                          # prot argument is *nix only

  i = m.rfind(StartKey)   # search for last occurrence of 'word'
  m.seek(i)             # seek to the location
  line = m.readline()   # read to the end of the line
#skip 4 lines
  for i in range(4):
    line = m.readline
  for line in m:
     if EndKey in line: break
     words = line.strip().split()
     Elem = AtomNumberToSymbol[int(words[1])]
     X = words[3]
     Y = words[4]
     Z = words[5]
     Output.append([Elem, X, Y, Z])
  return Output
#----------------------------------------------
def getFrequencies(F, Type):
#F is a file handle.  Rewind for safety
  F.seek(0, 0)
  StartKey = " Frequencies -- "
  IRKey = " IR Inten    -- "
  EndKey = '- Thermochemistry -'
  Output = [0.0, 0.0]
  getIR = False
  freq = []
  IR = []
  hits = []
  Large = 10.0**16
  for line in F:
     if EndKey in line: break
     if StartKey in line[:len(StartKey)]: #search from start
# get the frequency data and the IR intensity data
       for w in line[len(StartKey):].split(): #split after key
         freq.append(float(w))
       getIR = True
     if IRKey in line[:len(IRKey)] and getIR:
# get the frequency data and the IR intensity data
       for w in line[len(IRKey):].split():
         IR.append(float(w))
       getIR = False
#Now we have arrays of the freqs and intensities
  if Type == 's':
   for I,fr in enumerate(freq):
#    if fr <-5.0:
#      return [Large, 0.0]
    if fr > 1600.0 and fr < 1740.0:
      if IR[I] > 100.0 and IR[I] > Output[1]:
        Output = [fr, IR[I]]
        hits.append(I)
  elif Type == 'u':
   for I,fr in enumerate(freq):
#    if fr <-5.0:
#      return [Large, 0.0]
    if fr>1700.0 and fr < 1950.0:
      if IR[I] > 200.0  and IR[I] > Output[1]:
        Output = [fr, IR[I]]
        hits.append(I)
  return Output
#  if len(hits) == 1:
#    return Output
#  else:
#    return [0-Large, 0-Large]
#----------------------------------------------
def isSameCompound(name1, name2):
  '''name will have a format like Sulf_234_cis.  Only the first two fields matter'''
  c1 = name1.split('_')
  c2 = name2.split('_')
  for i in [0,1]:
    if c1[i] != c2[i]: return False
  return True
    
#----------------------------------------------
def getEnergy(F):
#F is a file handle.  Rewind for safety
  F.seek(0, 0)
  StartKey = "Sum of electronic and thermal Free Energies= "
  for line in F:
    if StartKey in line:
       words = line[len(StartKey):].strip().split()
       return words[0]
  return '0.0'
#----------------------------------------------
def Boltz(ConfList):
  ''' take a list of conformers, then do Boltzmann weighting, returning the single
freq. value for the list'''
  Emin = 9999.0
  for i,C in enumerate(ConfList):
    if C.energy < Emin:
      Emin = C.energy
      IndMin = i
  weight = []
  Sum = 0.0
  for C in ConfList:
    Ekcal_norm = (C.energy-Emin)/0.0015936010974 #norm by Emin, convert from AU to kcal/mol
    contrib = math.exp(-Ekcal_norm/(0.001986*298.15)) # fix by RT
    weight.append(contrib)
    Sum = Sum + contrib
  NewFreq = [0.0, 0.0]
  for i,C in enumerate(ConfList):
    for j in [0,1]:
      NewFreq[j] = NewFreq[j] + C.freq[j]*weight[i]/Sum
  Parts = ConfList[0].name.split('_')
  NewName = ('_').join(Parts[0:2])
  return Parts[1],NewName, NewFreq, ConfList[0].Type
#----------------------------------------------
def ConfInit(File):
  FileIn = open(File, 'r')
  E = getEnergy(FileIn)
  Name = os.path.splitext(os.path.basename(FileIn.name))[0] # cut off the directory and extension
  Type = ''
  if 'Sulf' in Name:
    Type = 's'
  elif 'Urea' in Name:
    Type = 'u'
  else:
    return 'NULL'
  Freq = getFrequencies(FileIn, Type)
  FileIn.close()
  return Conformer(name=Name, freq=Freq, energy=float(E), Type=Type)
#----------------------------------------------
def Pred_In_kobs(f, A):
  if A == 's':
    return 0.0595*f-112.8
  elif A == 'u':
    return 0.0498*f-102.8

#----------------------------------------------
def main():
  FilePath = 'toste_py.txt'
  Type = 's'
  fw = open(FilePath, 'w')
  folder = '/da/CADD/lewisri2/toste'
  FileList = glob.glob('%s/Collect*/*.log'%folder)
  ConfList = []
  for File in FileList:
    # set up new conformer and append it
    ConfList.append(ConfInit(File))
# now group into families
  Family = []
  nFamily = -1
  for i,C in enumerate(ConfList):
    if C == 'NULL': continue #already captured
#    print(C.name)
    Family.append([C])
    nFamily = nFamily + 1
    for j in range(i+1,len(ConfList)):
      if ConfList[j] == 'NULL': continue #already captured
      if isSameCompound(C.name, ConfList[j].name):
        Family[nFamily].append(ConfList[j])
        ConfList[j] = 'NULL'
#  print(Family)
  Results = dict()
  for L in Family:
    key, name, freq, T = Boltz(L)
    if key not in Results:
      Results[key] = dict()
    Results[key][T] = [name, freq[0], Pred_In_kobs(freq[0], T)]
  fw.write('ID\tAv Kobs\tT\tFreq\tpred Kobs\tT\tFreq\tpred Kobs\n')
  for key in Results:
    av = 0.0
    S = ''
    for i,T in enumerate(['s', 'u']):
      if T not in Results[key]:
        S = S+ "\t \t \t "
      else:
        av = av + Results[key][T][2]
        S = S + "\t%s\t%.2f\t%.2f"%(Results[key][T][0], Results[key][T][1], Results[key][T][2])
    av = av / (i+1.0)
    fw.write("%s\t%.2f%s\n"%(key,av,S))

if __name__ == "__main__":
  main()
