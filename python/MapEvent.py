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
#parser.add_argument("-d", "--inputDigiFile", dest="inputDigiFile", help="single  digi input file", required=False, default="sndLHC.Ntuple-TGeant4_dig.root")
#parser.add_argument("-n", "--nEvents", dest="nEvents", help="number of events to process", default=100000)
#parser.add_argument("-s", "--start", dest="start", type=int,help="file name with $*$", required=False,default=False)

# for fixing a root bug,  will be solved in the forthcoming 6.26 release.
ROOT.gInterpreter.Declare("""
#include "MuFilterHit.h"
#include "AbsMeasurement.h"
#include "TrackPoint.h"

void fixRoot(MuFilterHit& aHit,std::vector<int>& key,std::vector<float>& value, bool mask) {
   std::map<int,float> m = aHit.GetAllSignals(false);
   std::map<int, float>::iterator it = m.begin();
   while (it != m.end())
    {
        key.push_back(it->first);
        value.push_back(it->second);
        it++;
    }
}
void fixRootT(MuFilterHit& aHit,std::vector<int>& key,std::vector<float>& value, bool mask) {
   std::map<int,float> m = aHit.GetAllTimes(false);
   std::map<int, float>::iterator it = m.begin();
   while (it != m.end())
    {
        key.push_back(it->first);
        value.push_back(it->second);
        it++;
    }
}
void fixRoot(MuFilterHit& aHit, std::vector<TString>& key,std::vector<float>& value, bool mask) {
   std::map<TString, float> m = aHit.SumOfSignals();
   std::map<TString, float>::iterator it = m.begin();
   while (it != m.end())
    {
        key.push_back(it->first);
        value.push_back(it->second);
        it++;
    }
}
void fixRoot(std::vector<genfit::TrackPoint*>& points, std::vector<int>& d,std::vector<int>& k, bool mask) {
      for(std::size_t i = 0; i < points.size(); ++i) {
        genfit::AbsMeasurement*  m = points[i]->getRawMeasurement();
        d.push_back( m->getDetId() );
        k.push_back( int(m->getHitId()/1000) );
    }
}
""")

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

ut.bookHist(h,'kalman_angle_zx','Angle in ZX; [mrad]', 400,-200,200)
ut.bookHist(h,'kalman_angle_zy','Angle in ZY; [mrad]', 200,-200,200)
ut.bookHist(h,'kalman_res_zx','Slope resolution in ZX; MC true slopeXZ- rec slopeXZ [mrad]', 2000,-1000,1000)
ut.bookHist(h,'kalman_res_zy','Slope resolution in ZY; MC true slopeYZ- rec slopeYZ [mrad]', 2000,-1000,1000)
ut.bookHist(h,'slopes','Rec track slopes; slope X [mrad]; slope Y [mrad]',200,-100,100,200,-100,100)
ut.bookHist(h,'bs','beam spot; x[cm]; y[cm]',100,-90,10,100, -10, 90)
ut.bookHist(h,'ax_zx','Measured points on rec track z-x view; z [cm]; x [cm]',1000,250,600,1000,-90,10)
ut.bookHist(h,'ax_zy','Measured points on rec track z-y view; z [cm]; y [cm]',1000,250,600,1000,0,80)
ut.bookHist(h,'mc_ax_zx','MC points z-x view; z [cm]; x [cm]',1000,250,600,1000,-90,10)
ut.bookHist(h,'mc_ax_zy','MC points z-y view; z [cm]; y [cm]',1000,250,600,1000,0,80)
ut.bookHist(h,'Emuons','Energy of MC muons with points in the detector; E [GeV/c]',500, 0, 5000)
ut.bookHist(h,'genEmuons','Energy of all MC muons at scoring plane; E [GeV/c]',500, 0, 5000)
ut.bookHist(h,'Emuons_actionEvt','Action-on event: Energy of MC muons with points in the detector; E [GeV/c]',500, 0, 5000)
ut.bookHist(h,'genEmuons_actionEvt','Action-on event: Energy of MC muons at scoring plane; E [GeV/c]',500, 0, 5000)

ut.bookHist(h,'Taction','MC: Time difference btw primary muon start T and products startT; T_diff [ns]',2005,-5,2000)
ut.bookHist(h,'Taction_pid','MC: Time difference btw primary muon start T and products startT vs product PdgCode; daughter pid; T_diff [ns]',6000,-3000,3000,2005,-5,2000)
ut.bookHist(h,'Taction_procid','MC: Time difference btw primary muon start T and products startT vs TMCProcess id; procID; T_diff [ns]',50,-2,48, 2005,-5,2000)
ut.bookHist(h,'Taction_mufi','MC: Time difference btw MuFilter points; T_diff [ns]',2005,-5,2000)
ut.bookHist(h,'Taction_scifi','MC: Time difference btw SciFi points; T_diff [ns]',2005,-5,2000)
ut.bookHist(h,'allActionProcIds','TMCProcess id; TMCProcess id', 50, -2, 48)


Tkey  = ROOT.std.vector('TString')()
Ikey   = ROOT.std.vector('int')()
Value = ROOT.std.vector('float')()
def map2Dict(aHit,T='GetAllSignals',mask=True):
     if T=="SumOfSignals":
        key = Tkey
     elif T=="GetAllSignals" or T=="GetAllTimes":
        key = Ikey
     else: 
           print('use case not known',T)
           1/0
     key.clear()
     Value.clear()
     if T=="GetAllTimes": ROOT.fixRootT(aHit,key,Value,mask)
     else:                         ROOT.fixRoot(aHit,key,Value,mask)
     theDict = {}
     for k in range(key.size()):
         if T=="SumOfSignals": theDict[key[k].Data()] = Value[k]
         else: theDict[key[k]] = Value[k]
     return theDict
 
# Get input file
tchain = ROOT.TChain("cbmsim")
tchain.Add(options.inputFile)

#tree_digi = ROOT.TChain("cbmsim")
#tree_digi.Add(options.inputDigiFile)
#tchain.AddFriend(tree_digi)

MCpoints={}
Hits={}
trks={}
eventList ={}
pid = {}
procid = {}
counter = {}
Emuon = {}
Eothers = {}

slopeArray={}
MCslopeArray={}
why=0
yes=0
withPoints=0
withDSPoints=0
withDSHits=0
withVePoints=0
withVeHits=0
withVeDSPoints=0
withVeDSHits=0
withSfnoVePoints=0
withSfnoVeHits=0
withMuHits=0
th=0
th1=0
more=0

for i_event, event in enumerate(tchain) :
   ntrk=-1
   for aTrack in event.Reco_MuonTracks:
        ntrk+=1
        stateTrack = aTrack.getFittedState()
        posTrack   = stateTrack.getPos()
        rc = h['bs'].Fill(posTrack.x(),posTrack.y())
        mom = stateTrack.getMom()
        slopeX= mom.X()/mom.Z()
        slopeY= mom.Y()/mom.Z()
        h['slopes'].Fill(slopeX*1000,slopeY*1000)
        #h['theta'].Fill(asin((mom.)**2+mom.Y()**2)**0.5/mom.Mag()))
        h['kalman_angle_zx'].Fill(1000*math.atan(mom.X()/mom.Z()))
        h['kalman_angle_zy'].Fill(1000*math.atan(mom.Y()/mom.Z()))
        if i_event not in slopeArray : slopeArray[i_event] = []
        slopeArray[i_event]=[1000*math.atan(mom.X()/mom.Z()), 1000*math.atan(mom.Y()/mom.Z())]
        points = aTrack.getPoints() # 1st points or all
        points_z = []
        points_x = []
        points_y = []
        i = -1
        for Mpoint in aTrack.getPointsWithMeasurement():
            i+=1
            #RawHit = Mpoint.getRawMeasurement()
            #RawHitPos = RawHit.getRawHitCoords() 
            state = aTrack.getFittedState(i)
            pos    = state.getPos()
            points_z.append(pos[2])
            points_x.append(pos[0])
            points_y.append(pos[1])
            h['ax_zx'].Fill(pos[2], pos[0])
            h['ax_zy'].Fill(pos[2], pos[1])

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
   if Hits[i_event][1] > 40 : 
      eventList[i_event] = []
      eventList[i_event] = Hits[i_event]
   if (Hits[i_event][3]+Hits[i_event][4]) > 10 and not i_event in eventList :
      eventList[i_event] = [] 
      eventList[i_event] = Hits[i_event]
   if Hits[i_event][2] > 15 and not i_event in eventList : 
      eventList[i_event] = [] 
      eventList[i_event] = Hits[i_event]
   if Hits[i_event][0] > 10 and not i_event in eventList :
      eventList[i_event] = [] 
      eventList[i_event] = Hits[i_event]
   if i_event in eventList:
     ut.bookHist(h,'mc_ax_zx'+str(i_event),'MC points z-x view '+str(i_event)+'; z [cm]; x [cm]',100,250,600,100,-90,10)
     ut.bookHist(h,'mc_ax_zy'+str(i_event),'MC points z-y view '+str(i_event)+'; z [cm]; y [cm]',100,250,600,100,0,80)
     ut.bookHist(h,'mcmu_ax_zx'+str(i_event),'primary #mu: MC points z-x view '+str(i_event)+'; z [cm]; x [cm]',100,250,600,100,-90,10)
     ut.bookHist(h,'mcmu_ax_zy'+str(i_event),'primary #mu: MC points z-y view '+str(i_event)+'; z [cm]; y [cm]',100,250,600,100,0,80)
     #ut.bookHist(h,'pid '+str(i_event),'Pdg code; pdg code', 1000, 0, 1E6)
     #ut.bookHist(h,'procId '+str(i_event),'TMCProcess id; TMCProcess id', 50, -2, 48)
     #ut.bookHist(h,'procid_z '+str(i_event),'MC: TMCProcess vs zStart; zStart [cm]; procID',1000,250,600,51,-1,50)
     ut.bookHist(h,'mcnonmu_ax_zx'+str(i_event),'non-#mu particles in detector region: MCTrack start points z-x view '+str(i_event)+'; z [cm]; x [cm]',100,250,600,100,-90,10)
     ut.bookHist(h,'mcnonmu_ax_zy'+str(i_event),'non-#mu particles in detector region: MCTrack start points z-y view '+str(i_event)+'; z [cm]; y [cm]',100,250,600,100,0,80)

   MCpoints[i_event]=[]
   Nve = 0
   Nus = 0
   Nsf = 0
   Ndsh = 0
   Ndsv = 0
   hasPrMu = False
   prev = 9999999999.
   for mcpoint in event.MuFilterPoint:
     #print(mcpoint.GetDetectorID(), mcpoint.GetDetectorID()/10000)
     h['mc_ax_zx'].Fill(mcpoint.GetZ(), mcpoint.GetX())
     h['mc_ax_zy'].Fill(mcpoint.GetZ(), mcpoint.GetY())
     if int(mcpoint.GetDetectorID()/10000) == 1: Nve+=1
     if int(mcpoint.GetDetectorID()/10000) == 2: Nus+=1
     if int(mcpoint.GetDetectorID()/10000) == 3:
       #detID > 3xx60 are in vertical (YZ)
       if(mcpoint.GetDetectorID()%30000)%1000>120: print('exists')
       #print(mcpoint.GetDetectorID(), (mcpoint.GetDetectorID()%30000)%1000)
       if (mcpoint.GetDetectorID()%30000)%1000 > 60: Ndsv+=1
       else : Ndsh+=1
     if abs(mcpoint.PdgCode())==13 and mcpoint.GetTrackID()==0 : hasPrMu = True
     if i_event in eventList:
       #primary muon
       if abs(mcpoint.PdgCode())==13 and mcpoint.GetTrackID()==0 :
        h['mcmu_ax_zx'+str(i_event)].Fill(mcpoint.GetZ(), mcpoint.GetX())
        h['mcmu_ax_zy'+str(i_event)].Fill(mcpoint.GetZ(), mcpoint.GetY())
        if i_event ==2004: print("before" ,mcpoint.GetTime())
        if mcpoint.GetTime() < prev : 
           t_mufi=mcpoint.GetTime()
           prev = t_mufi
           if i_event==2004: print('Start', t_mufi)
       h['mc_ax_zx'+str(i_event)].Fill(mcpoint.GetZ(), mcpoint.GetX())
       h['mc_ax_zy'+str(i_event)].Fill(mcpoint.GetZ(), mcpoint.GetY())
       h['Taction_mufi'].Fill(t_mufi-mcpoint.GetTime())
       if i_event==2004: print(mcpoint.GetTime())         

   Nsf = len(event.ScifiPoint)
   prev = 9999999999.
   for mcpoint in event.ScifiPoint:
     h['mc_ax_zx'].Fill(mcpoint.GetZ(), mcpoint.GetX())
     h['mc_ax_zy'].Fill(mcpoint.GetZ(), mcpoint.GetY())
     if abs(mcpoint.PdgCode())==13 and mcpoint.GetTrackID()==0 : hasPrMu = True
     # even numbers are Y (horizontal plane), odd numbers X (vertical plane
     if (mcpoint.GetDetectorID()%1000+(mcpoint.GetDetectorID()%10000)//1000*128 )%2 == 0: Nsfh+=1
     else: Nsfv+=1
     if i_event in eventList:
        # primary muon
        if abs(mcpoint.PdgCode())==13 and mcpoint.GetTrackID()==0 :          
          h['mcmu_ax_zx'+str(i_event)].Fill(mcpoint.GetZ(), mcpoint.GetX())
          h['mcmu_ax_zy'+str(i_event)].Fill(mcpoint.GetZ(), mcpoint.GetY())
          if mcpoint.GetTime() < prev : 
            t_scifi=mcpoint.GetTime()
            prev = t_scifi
            #print(i_event, "once")
        h['mc_ax_zx'+str(i_event)].Fill(mcpoint.GetZ(), mcpoint.GetX())
        h['mc_ax_zy'+str(i_event)].Fill(mcpoint.GetZ(), mcpoint.GetY())
        h['Taction_scifi'].Fill(t_scifi-mcpoint.GetTime())
                          
   MCpoints[i_event]=[Nve,Nsf,Nus,Ndsh,Ndsv, Nve+Nsf+Nus+Ndsh+Ndsv, Nsfv, Nsfh]

   # Loop over MCTracks
   for mctrack in event.MCTrack :      
      #primary muon
      if mctrack.GetMotherId()==-1:
        h['genEmuons'].Fill(mctrack.GetEnergy())
        t_muon = mctrack.GetStartT()
        if hasPrMu : h['Emuons'].Fill(mctrack.GetEnergy())
        if i_event not in MCslopeArray : MCslopeArray[i_event] = []
        MCslopeArray[i_event]= [1000*math.atan(mctrack.GetPx()/mctrack.GetPz()), 1000*math.atan(mctrack.GetPy()/mctrack.GetPz())]        
        if i_event in eventList :
          if h['mcmu_ax_zx'+str(i_event)].GetEntries()+h['mcmu_ax_zy'+str(i_event)].GetEntries()>0 : h['Emuons_actionEvt'].Fill(mctrack.GetEnergy())
          else: print('action event with no mu in det', i_event)
          h['genEmuons_actionEvt'].Fill(mctrack.GetEnergy())
        if not i_event in Emuon: Emuon[i_event] = 0
        Emuon[i_event]=mctrack.GetEnergy()
        if not i_event in Eothers :         
          Eothers[i_event] = 0        
      else:
        h['Taction'].Fill(mctrack.GetStartT() - t_muon)   
        h['Taction_pid'].Fill(mctrack.GetPdgCode(), mctrack.GetStartT() - t_muon)   
        h['Taction_procid'].Fill(mctrack.GetProcID(), mctrack.GetStartT() - t_muon)      
        # take non-muons that are produced in detector region in Z
        if mctrack.GetMotherId()==0 : #and mctrack.GetStartZ()>200:
          Eothers[i_event] += mctrack.GetEnergy()
          if i_event==3319 and mctrack.GetEnergy()>100:print(ntrk, mctrack.GetMotherId(), mctrack.GetEnergy())          
        if i_event in eventList:
          #h['procid_z '+str(i_event)].Fill(mctrack.GetStartZ(), mctrack.GetProcID())
          h['mcnonmu_ax_zx'+str(i_event)].Fill(mctrack.GetStartZ(), mctrack.GetStartX())
          h['mcnonmu_ax_zy'+str(i_event)].Fill(mctrack.GetStartZ(), mctrack.GetStartY())          
        if not i_event in pid: 
          pid[i_event] = array('i')
          procid[i_event] = array('i')
        pid[i_event].append(mctrack.GetPdgCode())
        procid[i_event].append(mctrack.GetProcID())
        h['allActionProcIds'].Fill(mctrack.GetProcID())
        if i_event not in counter: counter[i_event] = {}          
        if pid[i_event][-1] not in counter[i_event]: 
          counter[i_event][pid[i_event][-1]]={}
        if procid[i_event][-1] not in counter[i_event][pid[i_event][-1]]:
          counter[i_event][pid[i_event][-1]][procid[i_event][-1]]= 0
        counter[i_event][pid[i_event][-1]][procid[i_event][-1]] += 1        
   if Emuon[i_event]<Eothers[i_event]: 
      more+=1
      print('cc', i_event)
   #if i_event==3319: print(Eothers[i_event] )
   
   trks[i_event]=[]
   trks[i_event]=[len(event.Reco_MuonTracks)]
print('Last eventN', i_event)   
for i in range(len(MCpoints)):
   if MCpoints[i][5]!=0:
     withPoints+=1
     if MCpoints[i][0]==0 :
       if MCpoints[i][1]!=0 : withSfnoVePoints+=1
   if MCpoints[i][0]!=0 :
     withVePoints+=1
     if MCpoints[i][3]>2 and MCpoints[i][4]>2:
       withVeDSPoints+=1
   if Hits[i][5]!=0: withMuHits+=1
   if Hits[i][0]!=0:
     withVeHits+=1
   if Hits[i][0]==0 :
     if Hits[i][1]!=0 : withSfnoVeHits+=1
   if MCpoints[i][3]>2 and MCpoints[i][4]>2:
     withDSPoints+=1
   if Hits[i][3]>2 and Hits[i][4]>2:
     withDSHits+=1     
     if Hits[i][0]!=0:
       withVeDSHits+=1
       withVeDSHits+=1
     if trks[i][0] < 1 : 
       why+=1#print(i, trks[i][0])
     else: yes+=1
#   if (MCpoints[i][3]>2 and MCpoints[i][4]>2) and not (Hits[i][3]>2 and Hits[i][4]>2):
#     print('why not digi', i)
#   if (Hits[i][3]>2 and Hits[i][4]>2) and (trks[i][0] < 1) :
#       print('why not rec', i)
   if trks[i][0] > 1 :
     print('more than one trk', trks[i][0], i)
   #if not len(event.Reco_MuonTracks)==0:
   #  print(i_event)
   #  print('MCPoints:', MCpoints[i_event])
   #  print('DigiHits', Hits[i_event])
     
   #TO BE DELETED   
   if (MCpoints[i][3]+MCpoints[i][4]+MCpoints[i][2])>40:
     th1+=1
   if (Hits[i][3]+Hits[i][4]+Hits[i][2])>40:     
     th+=1

p = {}
pe = {}
A,B = ROOT.TVector3(),ROOT.TVector3()
trans2local = False
for i, event in enumerate(tchain) :   
   for  aTrack in event.Reco_MuonTracks:
      S = aTrack.getFitStatus()
      if not S.isFitConverged(): print(i)
   if i not in eventList : continue
   print (i)
   p[i] = {}
   pe[i] = {}
   p[i]['xzX'] = array('d')
   p[i]['yzY'] = array('d')
   p[i]['xzZ'] = array('d')
   p[i]['yzZ'] = array('d')
   pe[i]['xzX'] = array('d')
   pe[i]['yzY'] = array('d')
   pe[i]['xzZ'] = array('d')
   pe[i]['yzZ'] = array('d')
   digis = []
   if event.FindBranch("Digi_ScifiHits"): digis.append(event.Digi_ScifiHits)
   if event.FindBranch("Digi_MuFilterHits"): digis.append(event.Digi_MuFilterHits)
   if event.FindBranch("Digi_MuFilterHit"): digis.append(event.Digi_MuFilterHit)
   empty = True
   for x in digis:
     if x.GetEntries()>0: empty = False
   if empty: continue
   systems = {1:'Veto',2:'US',3:'DS',0:'Scifi'}
   
   for D in digis:
      for digi in D:
         detID = digi.GetDetectorID()
         sipmMult = 1
         if digi.GetName()  == 'MuFilterHit':
            system = digi.GetSystem()
            geo.modules['MuFilter'].GetPosition(detID,A,B)
            sipmMult = len(digi.GetAllSignals())
            if sipmMult<minSipmMult and (system==1 or system==2): continue
            if trans2local:
                curPath = nav.GetPath()
                tmp = curPath.rfind('/')
                nav.cd(curPath[:tmp])
         else:
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
                   p[i]['xzZ'].append(Z)
                   p[i]['xzX'].append(locA[0])
                   pe[i]['xzX'].append(0)#detSize[system][0])
                   pe[i]['xzZ'].append(0)#detSize[system][2])
                   
         else:                         
                   p[i]['yzZ'].append(Z)
                   p[i]['yzY'].append(locA[1])
                   pe[i]['yzY'].append(0)#detSize[system][1])
                   pe[i]['yzZ'].append(0)#detSize[system][2])

for i in slopeArray :
  if i in MCslopeArray :
      h['kalman_res_zx'].Fill(MCslopeArray[i][0] - slopeArray[i][0])
      h['kalman_res_zy'].Fill(MCslopeArray[i][1] - slopeArray[i][1])
         
file = ROOT.TFile('rec_muons.root', 'recreate')
h['genEmuons'].Write()
h['Emuons'].Write()
h['genEmuons_actionEvt'].Write()
h['Emuons_actionEvt'].Write()
h['bs'].Write()
h['slopes'].Write()
#h['simslopes'].Write()
h['ax_zx'].Write()
h['ax_zy'].Write()
h['mc_ax_zx'].Write()
h['mc_ax_zy'].Write()
h['kalman_angle_zx'].Write()
h['kalman_angle_zy'].Write()
h['kalman_res_zx'].Write()
h['kalman_res_zy'].Write()
h['Nsf'].Write()
h['Nsfh'].Write()
h['Nsfv'].Write()
h['Ndsh'].Write()
h['Ndsv'].Write()
h['Nve'].Write()
h['Nds'].Write()
h['Nus'].Write()
h['Taction'].Write()
h['Taction_pid'].Write()
h['Taction_procid'].Write()
h['Taction_mufi'].Write()
h['Taction_scifi'].Write()
h['allActionProcIds'].Write()

for i in eventList:#range(500):
  #if i not in eventList: continue
  ut.bookCanvas(h,'actionEvt '+str(i),' ',1024,768,2,4)
  h['actionEvt '+str(i)].cd(1)
  #h['mc_ax_zx'+str(i)].Draw('colz')
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
  #h['mc_ax_zy'+str(i)].Draw('colz')
  if len(p[i]['yzY']) > 0:
   grYZ = ROOT.TGraphErrors(len(p[i]['yzY']), p[i]['yzZ'], p[i]['yzY'], pe[i]['yzZ'], pe[i]['yzY'])
   grYZ.SetMarkerColor(4)
   grYZ.SetMarkerStyle(7)
   grYZ.SetTitle('Hits in x-z plane')
   grYZ.GetXaxis().SetRangeUser(250,600)
   grYZ.GetYaxis().SetRangeUser(0,80)    
   grYZ.GetXaxis().SetTitle('z [cm]')
   grYZ.GetYaxis().SetTitle('Y [cm]')
   grYZ.Draw('AP')   
  h['actionEvt '+str(i)].cd(3)
  h['mcmu_ax_zx'+str(i)].SetMarkerStyle(7)
  h['mcmu_ax_zx'+str(i)].Draw('P')
  h['actionEvt '+str(i)].cd(4)
  h['mcmu_ax_zy'+str(i)].SetMarkerStyle(7)
  h['mcmu_ax_zy'+str(i)].Draw('P')
  h['actionEvt '+str(i)].cd(5)
  h['mcnonmu_ax_zx'+str(i)].Draw('colz')
  h['actionEvt '+str(i)].cd(6)
  h['mcnonmu_ax_zy'+str(i)].Draw('colz')
  h['actionEvt '+str(i)].cd(7)
  gr = ROOT.TGraph(len(pid[i]), pid[i], procid[i])
  gr.SetMarkerColor(ROOT.kRed)
  gr.SetMarkerStyle(20)
  gr.SetTitle('TMCProcess vs PID')
  gr.GetXaxis().SetTitle('Pdg Code')
  gr.GetYaxis().SetTitle('TMCProcess id')
  gr.Draw('AP')  
  h['actionEvt '+str(i)].cd(8)  
  txt=ROOT.TLatex()
  txt.DrawLatexNDC(0.4,0.95,"Event %d"%i)
  txt.DrawLatexNDC(0.05,0.85,"E_{primary #mu } = %5.2f [GeV/c]"%Emuon[i])
  if Emuon[i]<Eothers[i]:
    print('check that event',i)
    txt.DrawLatexNDC(0.5,0.85,"#Sigma(E_{primary #mu daughters }) = #color[2]{%5.2f} [GeV/c]"%Eothers[i])
  else: 
    txt.DrawLatexNDC(0.5,0.85,"#Sigma(E_{primary #mu daughters}) = %5.2f [GeV/c]"%Eothers[i])
  txt.DrawLatexNDC(0.05,0.65,"Digi_hits in event: #color[2]{Veto %d } #color[4]{SciFi %d } #color[8]{US %d } #color[7]{DS %d }"%(Hits[i][0], Hits[i][1], Hits[i][2], Hits[i][3]+Hits[i][4]))
  k=-1
  for ipid in counter[i] :
      for iproc in counter[i][ipid] :
        k+=1
        txt.DrawLatexNDC(0.5,0.7-k*0.1,"pid %d ProcID %d Number %d"%(ipid, iproc, counter[i][ipid][iproc]))
  txt.Draw()
  h['actionEvt '+str(i)].Write() 
    
print('N mctrk with points', withPoints, 'N mctrk with Veto points', withVePoints,'N mctrk with DS points', withDSPoints, 'N mctrk with Veto and DS points', withVeDSPoints, 'with Ve hits', withVeHits, 'with DS hits', withDSHits, 'with Veto and DS hits', withVeDSHits, 'not rec', why, 'rec', yes)
print('N mctrk with SciFi, but no Ve points', withSfnoVePoints, 'N with SciFi, but no Ve hits', withSfnoVeHits)
print('N events with muon hits in det', withMuHits)
print(th, th1)
print('N events with Eothers>Emu ', more)
print('N action events', len(eventList), 'of which primary muon left points in det in ', h['Emuons_actionEvt'].GetEntries())
print('Done')