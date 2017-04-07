import os
import sys
import string
import ROOT

directory = sys.argv[1]
Number = int(sys.argv[2])
MCorDATA = sys.argv[3]
Samples = sys.argv[4]

filenames_tmp = [os.path.normcase(f) for f in os.listdir(directory)]

filenames = [] 
expectedfiles = []
numberoffiles = []

for name in filenames_tmp:
    if Samples in name[0:name.rfind("_")]:
        filenames.append(name)

for name in filenames:
    if name[0:name.rfind("_")] not in expectedfiles and Samples in name[0:name.rfind("_")]:
        expectedfiles.append(name[0:name.rfind("_")])

for name in expectedfiles:
            numberoffiles.append(input("Number of files for " + name + ": "))

#print filenames
#print len(filenames)

if len(filenames) == Number:
    print 'All files are there.'
    for name in expectedfiles:
        for i in range(0, numberoffiles[expectedfiles.index(name)]):
            f = ROOT.TFile(directory + str(name + "_" + str(i) + ".root"))
            try:
                t = f.Get("AnalysisTree")
            except:
                os.system("sframe_main workdir/" + str(name.replace("uhh2.AnalysisModuleRunner."+ sys.argv[3] + ".","") + "_" + str(i+1) + ".xml"))
                f1 = ROOT.TFile(directory + str(name + "_" + str(i) + ".root"))
                try:
                    t = f1.Get("AnalysisTree")
                except:
                    print 'file ' + str(name + "_" + str(i) + ".root") + ' has no AnalysisTree' 
    
                    f. Close()
                    f1. Close()
                

else:
    if Number-len(filenames) < 0:
        print 'Number of expected files is too small!!!'
    if Number-len(filenames) == 1: 
        print Number-len(filenames) , ' file is missing.'
    if Number-len(filenames) > 1: 
        print Number-len(filenames) , ' files are missing.'


        for name in expectedfiles:
            for i in range(0, numberoffiles[expectedfiles.index(name)]):
                if str(name + "_" + str(i) + ".root") not in filenames:
                    print 'missing file: ' , str(name + "_" + str(i) + ".root") , '   run script ' , i+1
                    os.system("sframe_main workdir/" + str(name.replace("uhh2.AnalysisModuleRunner."+ sys.argv[3] + ".","") + "_" + str(i+1) + ".xml"))
                    f = ROOT.TFile(directory + str(name + "_" + str(i) + ".root"))
                    try:
                        t = f.Get("AnalysisTree")
                    except:
                        print 'file ' + str(name + "_" + str(i) + ".root") + ' has no AnalysisTree' 
                    f. Close()
                else:
                    f = ROOT.TFile(directory + str(name + "_" + str(i) + ".root"))
                    try:
                        t = f.Get("AnalysisTree")
                    except:
                        os.system("sframe_main workdir/" + str(name.replace("uhh2.AnalysisModuleRunner."+ sys.argv[3] + ".","") + "_" + str(i+1) + ".xml"))
                        f1 = ROOT.TFile(directory + str(name + "_" + str(i) + ".root"))
                        try:
                            t = f1.Get("AnalysisTree")
                        except:
                            print 'file ' + str(name + "_" + str(i) + ".root") + ' has no AnalysisTree' 
    
                    f. Close()
                    f1. Close()
                    

