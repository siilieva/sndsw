#!/usr/bin/env python
import ROOT, array
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

ut.bookHist(h,'Nsfv','Number of hits in Scifi vertical planes; N_sf V',1005,-5,1000)
ut.bookHist(h,'Nsfh','Number of hits in Scifi horizontal planes; N_sf H',1005,-5,1000)
ut.bookHist(h,'Nsf','Number of hits in Scifi; Nsf',1005,-5,1000)
ut.bookHist(h,'Nds','Number of hits in DS; Nds',1005,-5,1000)
ut.bookHist(h,'Ndsh','Number of hits in DS horizontal bars; Nds H',1005,-5,1000)
ut.bookHist(h,'Ndsv','Number of hits in DS vertical bars; Nds V',1005,-5,1000)
ut.bookHist(h,'Nve','Number of hits in Veto; Nve',55,-5,50)
ut.bookHist(h,'Nus','Number of hits in US; Nus',1005,-5,1000)

ut.bookHist(h,'kalman_angle_zx','Angle in ZX; [mrad]',400,-200,200)
ut.bookHist(h,'kalman_angle_zy','Angle in ZY; [mrad]',400,-200,200)
ut.bookHist(h,'slopes','Rec track slopes; slope X [mrad]; slope Y [mrad]',200,-100,100,200,-100,100)
ut.bookHist(h,'bs','beam spot; x[cm]; y[cm]',100,-80.,0.,100,0.,80.)
ut.bookHist(h,'ax_zx','z-x view; z [cm]; x [cm]',100,250,600,100,-80,10)
ut.bookHist(h,'ax_zy','z-y view; z [cm]; y [cm]',100,250,600,100,0,80)
ut.bookHist(h,'mc_ax_zx','MC z-x view; z [cm]; x [cm]',100,0,600,100,-80,10)
ut.bookHist(h,'mc_ax_zy','MC z-y view; z [cm]; y [cm]',100,250,600,100,0,80)

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
th=0
th1=0

for i_event, event in enumerate(tchain) :
   for mcpoint in event.MuFilterPoint:
         #if abs(mcpoint.PdgCode()) != 13 : continue
         #trackID = mcpoint.GetTrackID()
         #detID = mcpoint.GetDetectorID()
         h['mc_ax_zx'].Fill(mcpoint.GetZ(), mcpoint.GetX())
         h['mc_ax_zy'].Fill(mcpoint.GetZ(), mcpoint.GetY())

   for aTrack in event.Reco_MuonTracks:
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
   if (Hits[i_event][3]+Hits[i_event][4]) > 20 and not i_event in eventList :
      eventList[i_event] = [] 
      eventList[i_event] = Hits[i_event]
   if Hits[i_event][2] > 15 and not i_event in eventList : 
      eventList[i_event] = [] 
      eventList[i_event] = Hits[i_event]
   if Hits[i_event][0] > 10 and not i_event in eventList :
      eventList[i_event] = [] 
      eventList[i_event] = Hits[i_event]
   if i_event in eventList:
     ut.bookHist(h,'mc_ax_zx'+str(i_event),'MC z-x view'+str(i_event)+'; z [cm]; x [cm]',100,250,600,100,-80,10)
     ut.bookHist(h,'mc_ax_zy'+str(i_event),'MC z-y view'+str(i_event)+'; z [cm]; y [cm]',100,250,600,100,0,80)
     ut.bookHist(h,'mcmu_ax_zx'+str(i_event),'MC muons z-x view'+str(i_event)+'; z [cm]; x [cm]',100,250,600,100,-80,10)
     ut.bookHist(h,'mcmu_ax_zy'+str(i_event),'MC muons z-y view'+str(i_event)+'; z [cm]; y [cm]',100,250,600,100,0,80)
     ut.bookHist(h,'pid'+str(i_event),'Pdg code; pdg code', 1000, 0, 1E6)
     ut.bookHist(h,'procId'+str(i_event),'TMCProcess id; TMCProcess id', 50, -2, 48)
   
   MCpoints[i_event]=[]
   Nve = 0
   Nus = 0
   Nsf = 0
   Ndsh = 0
   Ndsv = 0
   for mcpoint in event.MuFilterPoint:
     #print(mcpoint.GetDetectorID(), mcpoint.GetDetectorID()/10000)
     if int(mcpoint.GetDetectorID()/10000) == 1: Nve+=1
     if int(mcpoint.GetDetectorID()/10000) == 2: Nus+=1
     if int(mcpoint.GetDetectorID()/10000) == 3:
       #detID > 3xx60 are in vertical (YZ)
       if(mcpoint.GetDetectorID()%30000)%1000>120: print('exists')
       #print(mcpoint.GetDetectorID(), (mcpoint.GetDetectorID()%30000)%1000)
       if (mcpoint.GetDetectorID()%30000)%1000 > 60: Ndsv+=1
       else : Ndsh+=1
       if i_event in eventList:
        h['mc_ax_zx'+str(i_event)].Fill(mcpoint.GetZ(), mcpoint.GetX())
        h['mc_ax_zy'+str(i_event)].Fill(mcpoint.GetZ(), mcpoint.GetY())
        if mcpoint.PdgCode()==13:
          h['mcmu_ax_zx'+str(i_event)].Fill(mcpoint.GetZ(), mcpoint.GetX())
          h['mcmu_ax_zy'+str(i_event)].Fill(mcpoint.GetZ(), mcpoint.GetY())
        #h['pid'+str(i_event)].Fill(mcpoint.PdgCode())
        if not i_event in pid: 
          pid[i_event] = array('i')
          procid[i_event] = array('i')
        pid[i_event].append(mcpoint.PdgCode())
        procid[i_event].append(event.MCTrack[mcpoint.GetTrackID()].GetProcID())
        #h['procId'+str(i_event)].Fill(event.MCTrack[mcpoint.GetTrackID()].GetProcID()) 
   Nsf = len(event.ScifiPoint)
   for mcpoint in event.ScifiPoint:
     # even numbers are Y (horizontal plane), odd numbers X (vertical plane
     if (mcpoint.GetDetectorID()%1000+(mcpoint.GetDetectorID()%10000)//1000*128 )%2 == 0: Nsfh+=1
     else: Nsfv+=1     
     if i_event in eventList:
        h['mc_ax_zx'+str(i_event)].Fill(mcpoint.GetZ(), mcpoint.GetX())
        h['mc_ax_zy'+str(i_event)].Fill(mcpoint.GetZ(), mcpoint.GetY())
        if mcpoint.PdgCode()==13:
          h['mcmu_ax_zx'+str(i_event)].Fill(mcpoint.GetZ(), mcpoint.GetX())
          h['mcmu_ax_zy'+str(i_event)].Fill(mcpoint.GetZ(), mcpoint.GetY())
        #h['pid'+str(i_event)].Fill(mcpoint.PdgCode())
        if not i_event in pid: 
          pid[i_event] = array('i')
          procid[i_event] = array('i')
        pid[i_event].append(mcpoint.PdgCode())
        procid[i_event].append(event.MCTrack[mcpoint.GetTrackID()].GetProcID())
        #counter[pid]
        #print(event.MCTrack[mcpoint.GetTrackID()].GetProcID())
        #h['procId'+str(i_event)].Fill(event.MCTrack[mcpoint.GetTrackID()].GetProcID())     
   MCpoints[i_event]=[Nve,Nsf,Nus,Ndsh,Ndsv, Nve+Nsf+Nus+Ndsh+Ndsv, Nsfv, Nsfh]
   
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
     #print(i)
     #print(Hits[i][1],i)
     th+=1
if False:     
 for i, event in enumerate(tchain) :
   if i!=9909: continue
   for aHit in event.Digi_MuFilterHits:
            # only DS hits
            if aHit.GetSystem() != 3: continue
            detID = aHit.GetDetectorID()            
            #if detID>31999: continue
            sumSignal = map2Dict(aHit,'SumOfSignals')
            print('sumSignals', sumSignal["Sum"], 'mcpoints:')
            for mc_point_i, _ in event.Digi_MuFilterHits2MCPoints[0].wList(detID) :
                print(aHit.isVertical(), detID, event.MuFilterPoint[mc_point_i].GetX(), event.MuFilterPoint[mc_point_i].GetY(),event.MuFilterPoint[mc_point_i].GetZ(),  event.MuFilterPoint[mc_point_i].PdgCode())
 for i, event in enumerate(tchain) :
   #if i!=9909: continue 
   for  aTrack in event.Reco_MuonTracks:
      S = aTrack.getFitStatus()
      if not S.isFitConverged(): print(i)   

file = ROOT.TFile('rec_muons.root', 'recreate')
#h['simE'].Draw()
#h['simE'].Write()
h['bs'].Write()
#h['us'].Write()
#h['ds'].Write()
h['slopes'].Write()
#h['simslopes'].Write()
#h['resslXZ'].Write()
#h['resslYZ'].Write()
h['ax_zx'].Write()
h['ax_zy'].Write()
h['kalman_angle_zx'].Write()
h['kalman_angle_zy'].Write()
h['Nsf'].Write()
h['Nsfh'].Write()
h['Nsfv'].Write()
h['Ndsh'].Write()
h['Ndsv'].Write()
h['Nve'].Write()
h['Nds'].Write()
h['Nus'].Write()
print(len(eventList))
for i in range(200):
  if i not in eventList: continue
  ut.bookCanvas(h,'strangeEvt'+str(i),' ',1024,1024,2,3)
  h['strangeEvt'+str(i)].cd(1)
  h['mc_ax_zx'+str(i)].Draw('colz')
  h['strangeEvt'+str(i)].cd(2)
  h['mc_ax_zy'+str(i)].Draw('colz')
  h['strangeEvt'+str(i)].cd(3)
  h['mcmu_ax_zx'+str(i)].Draw('*')
  h['strangeEvt'+str(i)].cd(4)
  h['mcmu_ax_zy'+str(i)].Draw('*')
  h['strangeEvt'+str(i)].cd(5)
  gr = ROOT.TGraph(len(pid[i]), pid[i], procid[i])
  gr.SetMarkerColor(4)
  gr.SetMarkerStyle(21)
  gr.SetTitle('TMCProcess vs PID')
  gr.GetXaxis().SetTitle('Pdg Code')
  gr.GetYaxis().SetTitle('TMCProcess id')
  gr.Draw('AP')
  #h['pid'+str(i)].Draw('l')
  #h['strangeEvt'+str(i)].cd(6)
  #h['procId'+str(i)].Draw('l') 
  h['strangeEvt'+str(i)].Write() 
     
print('N mctrk with points', withPoints, 'N mctrk with Veto points', withVePoints,'N mctrk with DS points', withDSPoints, 'N mctrk with Veto and DS points', withVeDSPoints, 'with Ve hits', withVeHits, 'with DS hits', withDSHits, 'with Veto and DS hits', withVeDSHits, 'not rec', why, 'rec', yes)
print('N mctrk with SciFi, but no Ve points', withSfnoVePoints, 'N with SciFi, but no Ve hits', withSfnoVeHits)
print(th, th1)
print('Done')
