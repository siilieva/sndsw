import ROOT
from argparse import ArgumentParser
import os

import SndlhcMuonReco

parser = ArgumentParser()
parser.add_argument("-f", "--inputFile", dest="inputFile", help="single input file", required=True)
parser.add_argument("-o", "--withOutput", dest="withOutput", help="persistent output", action='store_true',default=False)
parser.add_argument("-g", "--geoFile", dest="geoFile", help="geofile", required=True)
parser.add_argument("-n", "--nEvents", dest="nEvents",  type=int, help="number of events to process", default=100000)
parser.add_argument("-i", "--firstEvent",dest="firstEvent",  help="First event of input file to use", required=False,  default=0, type=int)
parser.add_argument("-t", "--tolerance", dest="tolerance",  type=float, help="How far away from Hough line hits assigned to the muon can be. In cm.", default=0.)

parser.add_argument("--hits_to_fit", dest = "hits_to_fit", type=str, help="Which detectors to use in the fit, in the format: vesfusds, where [ve] is veto, [sf] is Scifi, [us] is Upstream muon filter, and [ds] is downstream muon filter. Default: sfusds", default = "sfusds")
parser.add_argument("--hits_for_triplet", dest = "hits_for_triplet", type=str, help="Which detectors to use for the triplet condition. In the same format as --hits_to_fit. Default: ds", default = "ds")
parser.add_argument("-wX", "--withXmeas", dest="withXmeas", help="using X measurement by horizontal MuFi bars", action='store_true',default=False)

options = parser.parse_args()

import SndlhcGeo
geo = SndlhcGeo.GeoInterface(options.geoFile)

lsOfGlobals = ROOT.gROOT.GetListOfGlobals()
lsOfGlobals.Add(geo.modules['Scifi'])
lsOfGlobals.Add(geo.modules['MuFilter'])

x = options.inputFile
filename = x[x.rfind('/')+1:]
outFileName = filename.replace('.root','_muonReco.root')

fullPath = options.inputFile
if options.inputFile.find('/eos')==0:
     fullPath = os.environ['EOSSHIP']+options.inputFile
F = ROOT.TFile.Open(fullPath)

if options.withOutput:
  outFile = ROOT.TFile(outFileName, 'RECREATE')
else:
  outFile = ROOT.TMemFile(outFileName,'CREATE')

run = ROOT.FairRunAna()
print("Initialized FairRunAna")

source = ROOT.FairFileSource(F)
run.SetSource(source)

sink = ROOT.FairRootFileSink(outFile)
run.SetSink(sink)

# Hough transform
HoughTask = SndlhcMuonReco.MuonReco()    
#run.AddTask(HoughTask)
#run.Init()

# simple tracking
import SndlhcTracking
trackTask = SndlhcTracking.Tracking() 
trackTask.SetName('simpleTracking')
run.AddTask(trackTask)
run.Init()
#print('before trackTask.Init()')
#trackTask.Init()

# Run check of Nhits in each event
# Get input file
tchain = ROOT.TChain("cbmsim")
tchain.Add(options.inputFile)
print("Number of events", tchain.GetEntries())
Hits = {}
for i_event, event in enumerate(tchain) :
  #if i_event >11000 : break
  #if i_event%1000 ==0 :print(i_event)
  Hits[i_event]=[]
  Nve = 0
  Nus = 0
  Nsf = 0
  Nsfv = 0
  Nsfh = 0
  Ndsh = 0
  Ndsv = 0
  Nds = 0
  if False:
   for aHit in event.Digi_MuFilterHits:
     if aHit.GetSystem() == 1: Nve+=1
     if aHit.GetSystem() == 2: Nus+=1
     if aHit.GetSystem() == 3: 
       if aHit.isVertical() : Ndsv+=1
       else : Ndsh+=1
       Nds = Ndsv + Ndsh
   Nsf = len(event.Digi_ScifiHits)
   for aHit in event.Digi_ScifiHits:
     if aHit.isVertical() : Nsfv+=1
     else: Nsfh+=1
   Hits[i_event]=[Nve,Nsf,Nus,Nds, Nve+Nsf+Nus+Ndsh+Ndsv, Nsfv, Nsfh, Ndsh,Ndsv]

  if 1: #Hits[i_event][2] < 20 and Hits[i_event][1] < 20 and (Hits[i_event][3] < 10  or Hits[i_event][3] >= 10 ):
    #pass
    #print('now simple tracking')
    source.ReadEvent(i_event)
    trackTask.ExecuteTask()
    trackTask.FinishEvent() # fill
  else:
    # The following lines must be *after* run.Init()
    print('now Hough')
    HoughTask.SetTolerance(options.tolerance)
    HoughTask.SetHitsToFit("vesfusds")
    if Hits[i_event][3] >= 10 and Hits[i_event][2] >=20 :
      HoughTask.SetHitsForTriplet("sf")
    elif Hits[i_event][3] >= 10 and Hits[i_event][1] >= 20 :
      HoughTask.SetHitsForTriplet("us")
    else :  HoughTask.SetHitsForTriplet("ds")
    HoughTask.UseXmeas(options.withXmeas)
    HoughTask.Exec(options)  
#run.Run(options.firstEvent, options.firstEvent + options.nEvents)  
trackTask.FinishTask()
#HoughTask.FinishTask()
print("Done running")

