import os
import sys
import string
import ROOT

expectedfiles = []
ExpNumbers = []
RealNumbers = []

directory = sys.argv[1]
MCorDATA = sys.argv[2]

file = open("samples.txt","r")
for line in file:
    cleaned = line.strip()
    cleaned = cleaned.replace(":","")
    cleaned = cleaned.split()
    expectedfiles.append(cleaned[0])
    ExpNumbers.append(int(cleaned[2]))

filenames= [os.path.normcase(f) for f in os.listdir(directory)]

for expfile in expectedfiles:
    RealNumbers.append(int(sum(expfile in x for x in filenames)))
 
#Check if all files are there
missing = 0

for i in range(0, len(expectedfiles)):
    if RealNumbers[i] == ExpNumbers[i]:
        print expectedfiles[i], ": OK"
    else:
        print expectedfiles[i], ":", ExpNumbers[i]-RealNumbers[i], "files are missing"
        missing = 1


#Show missing files
if missing:
    showmissingfiles = raw_input("Do you want to see which files are missing? (y/n) ")
    if showmissingfiles == "y":
        for name in expectedfiles:
            for i in range(0, ExpNumbers[expectedfiles.index(name)]):
                if str("uhh2.AnalysisModuleRunner." + MCorDATA + "." + name + "_" + str(i) + ".root") not in filenames:
                    print 'missing file:' , str(name + "_" + str(i) + ".root") 

#Rerun missing files locally
if missing:
    run = raw_input("Do you want run the missing files locally? (y/n) ")
    if run == "y":
        for name in expectedfiles:
            for i in range(0, ExpNumbers[expectedfiles.index(name)]):
                if str("uhh2.AnalysisModuleRunner." + MCorDATA + "." + name + "_" + str(i) + ".root") not in filenames:
                    print 'missing file:' , str(name + "_" + str(i) + ".root") , '   Rerun script' , i+1
                    os.system("sframe_main workdir/" + str(name) + "_" + str(i+1) + ".xml")

#Check if all files contain an input tree
checkinputtree = raw_input("Do you want check if all files contain an InputTree? (y/n) ")
if checkinputtree == "y":
    for i in range(0, ExpNumbers[expectedfiles.index(name)]):
        f = ROOT.TFile(directory + str(name + "_" + str(i) + ".root"))
        try:
            t = f.Get("AnalysisTree")
        
        except:
            print 'file ' + str(name + "_" + str(i) + ".root") + ' has no AnalysisTree' 
            rerun = raw_input("Do you want rerun the broken file locally? (y/n) ")
            if rerun == "y":
                os.system("sframe_main workdir/" + str(name) + "_" + str(i+1) + ".xml")
                f1 = ROOT.TFile(directory + str(name + "_" + str(i) + ".root"))
                try:
                    t = f1.Get("AnalysisTree")
                except:
                    print 'file ' + str(name + "_" + str(i) + ".root") + ' has still no AnalysisTree' 
                f1. Close()
        f.Close()



