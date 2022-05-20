#!/usr/bin/env python
import ROOT, array
ROOT.gROOT.SetBatch(True)
import os,sys,subprocess
import ctypes
from array import array
import rootUtils as ut
from scipy import stats
import math
import numpy as np

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("-f", "--inputFile", dest="inputFile", help="single input file", required=False, default="sndLHC.Ntuple-TGeant4_dig.root")
parser.add_argument("-g", "--geoFile", dest="geoFile", help="geofile", required=True)
parser.add_argument("-d", "--inputDigiFile", dest="inputDigiFile", help="single  digi input file", required=False, default="sndLHC.Ntuple-TGeant4_dig.root")
parser.add_argument("-n", "--nEvents", dest="nEvents", help="number of events to process", default=100000)
parser.add_argument("-mc", "--mc", dest="MC", help="flag MC or data", default=1)
#parser.add_argument("-s", "--start", dest="start", type=int,help="file name with $*$", required=False,default=False)

options = parser.parse_args()
h={} # histogram storage

import SndlhcGeo
geo = SndlhcGeo.GeoInterface(options.geoFile)

detSize = {}
si = geo.snd_geo.Scifi
detSize[0] =[si.channel_width, si.channel_width, si.scifimat_z ]
mi = geo.snd_geo.MuFilter
detSize[1] =[mi.VetoBarX/2,                   mi.VetoBarY/2,            mi.VetoBarZ/2]
detSize[2] =[mi.UpstreamBarX/2,           mi.UpstreamBarY/2,    mi.UpstreamBarZ/2]
detSize[3] =[mi.DownstreamBarX_ver/2,mi.DownstreamBarY/2,mi.DownstreamBarZ/2]

minSipmMult = 1

nav = ROOT.gGeoManager.GetCurrentNavigator()

ut.bookHist(h,'Nsfv','Number of hits in Scifi vertical planes; N_sf V',1005,-5,1000)
ut.bookHist(h,'Nsfh','Number of hits in Scifi horizontal planes; N_sf H',1005,-5,1000)
ut.bookHist(h,'Nsf','Number of hits in Scifi; Nsf',1005,-5,1000)
ut.bookHist(h,'Nds','Number of hits in DS; Nds',1005,-5,1000)
ut.bookHist(h,'Ndsh','Number of hits in DS horizontal bars; Nds H',1005,-5,1000)
ut.bookHist(h,'Ndsv','Number of hits in DS vertical bars; Nds V',1005,-5,1000)
ut.bookHist(h,'Nve','Number of hits in Veto; Nve',55,-5,50)
ut.bookHist(h,'Nus','Number of hits in US; Nus',1005,-5,1000)

if options.MC:
     ut.bookHist(h,'Emuons','Energy of MC muons with points in the detector; E [GeV/c]',500, 0, 5000)
     ut.bookHist(h,'genEmuons','Energy of all MC muons at scoring plane; E [GeV/c]',500, 0, 5000)

ut.bookHist(h,'TallHitTo1st','Time btw 1st hit in event and rest of hits in event; T_hit -T_1Hit', 1000, 0, 1000)
 
# Get input file
tchain = ROOT.TChain("cbmsim")
tchain.Add(options.inputFile)

if options.inputDigiFile :
     tree_digi = ROOT.TChain("cbmsim")
     tree_digi.Add(options.inputDigiFile)
     tchain.AddFriend(tree_digi)

Hits={}
trks={}
eventList ={}
Emuon = {}

slopeArray={}
MCslopeArray={}
p={}
pe={}

for i_event, event in enumerate(tchain) :
     if i_event > int(options.nEvents) : break
     Hits[i_event]=[]
     Nve = 0
     Nus = 0
     Nsf = 0
     Nsfv = 0
     Nsfh = 0
     Ndsh = 0
     Ndsv = 0
     for aHit in event.Digi_MuFilterHits:
         if aHit.GetSystem() == 1: Nve+=1
         if aHit.GetSystem() == 2: Nus+=1
         if aHit.GetSystem() == 3: 
             if aHit.isVertical() : Ndsv+=1
             else : Ndsh+=1
     Nsf = len(event.Digi_ScifiHits)
     for aHit in event.Digi_ScifiHits:
         if aHit.isVertical() : Nsfv+=1
         else: Nsfh+=1
     Hits[i_event]=[Nve,Nsf,Nus,Ndsh,Ndsv, Nve+Nsf+Nus+Ndsh+Ndsv, Nsfv, Nsfh]   
     h['Nsf'].Fill(Hits[i_event][1])
     h['Nsfv'].Fill(Hits[i_event][6])
     h['Nsfh'].Fill(Hits[i_event][7])
     h['Ndsh'].Fill(Hits[i_event][3])
     h['Ndsv'].Fill(Hits[i_event][4])
     h['Nve'].Fill(Hits[i_event][0])
     h['Nds'].Fill(Hits[i_event][3]+Hits[i_event][4])
     h['Nus'].Fill(Hits[i_event][2])

     # SELECT INTERESTING EVENT TO INSPECT - BASED ON NHITS FOR NOW
     if (Hits[i_event][1] < 4  and Hits[i_event][1] !=0) : 
         if (Hits[i_event][3]+Hits[i_event][4]) < 3 and (Hits[i_event][3]+Hits[i_event][4])!=0 :
             if Hits[i_event][2] < 2 and Hits[i_event][2]!=0 : 
                 if Hits[i_event][5] < 4 and Hits[i_event][5] !=0 :
                     eventList[i_event] = [] 
                     eventList[i_event] = Hits[i_event]

     if i_event in eventList:
         ut.bookHist(h,'mc_ax_zx'+str(i_event),'MC points z-x view '+str(i_event)+'; z [cm]; x [cm]',100,250,600,100,-90,10)
         ut.bookHist(h,'mc_ax_zy'+str(i_event),'MC points z-y view '+str(i_event)+'; z [cm]; y [cm]',100,250,600,100,0,80)
         
     t_muon = 0. # reference T0
     if options.MC:
         # Loop over MCTracks
         for mctrack in event.MCTrack :      
             #primary muon
             if mctrack.GetMotherId()==-1:
                 h['genEmuons'].Fill(mctrack.GetEnergy())
                 t_muon = mctrack.GetStartT()
                 h['Emuons'].Fill(mctrack.GetEnergy())
                 if not i_event in Emuon: Emuon[i_event] = 0
                 Emuon[i_event]=mctrack.GetEnergy()                 
   
     trks[i_event]=[]
     trks[i_event]=[len(event.Reco_MuonTracks)]
     if trks[i_event][0] > 1 : print('Ntrks ', trks[i_event][0], ' in event ', i_event)
     
     # FOR EVENTS IN THE EVENTLIST, PLOT HITS
     if i_event not in eventList : continue
     p[i_event] = {}
     pe[i_event] = {}
     p[i_event]['xzX'] = array('d')
     p[i_event]['yzY'] = array('d')
     p[i_event]['xzZ'] = array('d')
     p[i_event]['yzZ'] = array('d')
     pe[i_event]['xzX'] = array('d')
     pe[i_event]['yzY'] = array('d')
     pe[i_event]['xzZ'] = array('d')
     pe[i_event]['yzZ'] = array('d')
     digis = []
     if event.FindBranch("Digi_ScifiHits"): digis.append(event.Digi_ScifiHits)
     if event.FindBranch("Digi_MuFilterHits"): digis.append(event.Digi_MuFilterHits)
     if event.FindBranch("Digi_MuFilterHit"): digis.append(event.Digi_MuFilterHit)
     empty = True
     for x in digis:
          if x.GetEntries()>0: empty = False
          if empty: continue
     systems = {1:'Veto',2:'US',3:'DS',0:'Scifi'}
     print(i_event)
     
     A,B = ROOT.TVector3(),ROOT.TVector3()
     trans2local = False
     t0 = 9e10
     tHit=[]
     for D in digis:
         for digi in D:      
             detID = digi.GetDetectorID()
             sipmMult = 1
             if digi.GetName()  == 'MuFilterHit':
                 #print('mu', digi.GetTime())
                 #tHit.append(digi.GetImpactT())
                 #if digi.GetImpactT() < t0 : t0 = digi.GetImpactT()
                 system = digi.GetSystem()
                 geo.modules['MuFilter'].GetPosition(detID,A,B)
                 sipmMult = len(digi.GetAllSignals())
                 if sipmMult<minSipmMult and (system==1 or system==2): continue
                 if trans2local:
                     curPath = nav.GetPath()
                     tmp = curPath.rfind('/')
                     nav.cd(curPath[:tmp])
             else:
                 #print('sf ', digi.GetTime())
                 tHit.append(digi.GetTime()-t_muon)# in ns
                 if (digi.GetTime(0) - t_muon) < t0 : t0 = digi.GetTime(0) - t_muon
                 geo.modules['Scifi'].GetSiPMPosition(detID,A,B)
                 if trans2local:
                     curPath = nav.GetPath()
                     tmp = curPath.rfind('/')
                     nav.cd(curPath[:tmp])
                 system = 0
             globA,locA = array('d',[A[0],A[1],A[2]]),array('d',[A[0],A[1],A[2]])
             if trans2local:   nav.MasterToLocal(globA,locA)
             Z = A[2]
             if digi.isVertical():
                 p[i_event]['xzZ'].append(Z)
                 p[i_event]['xzX'].append(locA[0])
                 pe[i_event]['xzX'].append(0)#detSize[system][0])
                 pe[i_event]['xzZ'].append(0)#detSize[system][2])               
             else:                         
                 p[i_event]['yzZ'].append(Z)
                 p[i_event]['yzY'].append(locA[1])
                 pe[i_event]['yzY'].append(0)#detSize[system][1])
                 pe[i_event]['yzZ'].append(0)#detSize[system][2])

     for T in tHit :
         #print('hit times', t0, T)
         h['TallHitTo1st'].Fill(T-t0)
         
file = ROOT.TFile('rec_muons.root', 'recreate')
if options.MC :
     h['genEmuons'].Write()
     h['Emuons'].Write()
h['Nsf'].Write()
h['Nsfh'].Write()
h['Nsfv'].Write()
h['Ndsh'].Write()
h['Ndsv'].Write()
h['Nve'].Write()
h['Nds'].Write()
h['Nus'].Write()
h['TallHitTo1st'].Write()

for i in eventList:
  ut.bookCanvas(h,'actionEvt '+str(i),' ',1024,768,2,2)
  h['actionEvt '+str(i)].cd(1)
  if len(p[i]['xzX']) > 0:
    grXZ = ROOT.TGraphErrors(len(p[i]['xzX']), p[i]['xzZ'], p[i]['xzX'], pe[i]['xzZ'], pe[i]['xzX'])
    grXZ.SetMarkerColor(4)
    grXZ.SetMarkerStyle(7)
    grXZ.SetTitle('Hits in x-z plane')
    grXZ.GetXaxis().SetRangeUser(250,600)
    grXZ.GetYaxis().SetRangeUser(-90,10)    
    grXZ.GetXaxis().SetTitle('z [cm]')
    grXZ.GetYaxis().SetTitle('x [cm]')
    grXZ.Draw('AP')    
  h['actionEvt '+str(i)].cd(2)
  if len(p[i]['yzY']) > 0:
   grYZ = ROOT.TGraphErrors(len(p[i]['yzY']), p[i]['yzZ'], p[i]['yzY'], pe[i]['yzZ'], pe[i]['yzY'])
   grYZ.SetMarkerColor(4)
   grYZ.SetMarkerStyle(7)
   grYZ.SetTitle('Hits in y-z plane')
   grYZ.GetXaxis().SetRangeUser(250,600)
   grYZ.GetYaxis().SetRangeUser(0,80)    
   grYZ.GetXaxis().SetTitle('z [cm]')
   grYZ.GetYaxis().SetTitle('y [cm]')
   grYZ.Draw('AP')   
  h['actionEvt '+str(i)].cd(3)  
  txt=ROOT.TLatex()
  txt.DrawLatexNDC(0.4,0.95,"Event %d"%i)
  txt.DrawLatexNDC(0.05,0.85,"E_{primary #mu } = %5.2f [GeV/c]"%Emuon[i])
  txt.DrawLatexNDC(0.05,0.45,"Digi_hits in event YZ: #color[2]{Veto %d } #color[4]{SciFi %d } #color[8]{US %d } #color[7]{DS %d }"%(Hits[i][0], Hits[i][7], Hits[i][2],Hits[i][3]))
  txt.DrawLatexNDC(0.05,0.65,"Digi_hits in event XZ: #color[4]{SciFi %d } #color[7]{DS %d }"%(Hits[i][6], Hits[i][4]))
  txt.Draw()
  h['actionEvt '+str(i)].Write() 
# Hits[i_event]=[Nve,Nsf,Nus,Ndsh,Ndsv, Nve+Nsf+Nus+Ndsh+Ndsv, Nsfv, Nsfh]   
print('Done running')