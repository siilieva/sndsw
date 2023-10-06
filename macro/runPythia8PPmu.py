import ROOT
import rootUtils as ut
from array import array
import math
theSeed      = 0
ROOT.gRandom.SetSeed(theSeed)
ROOT.gSystem.Load("libpythia8")

mudict = {13: "mu", -13: "amu"}

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("-f", "--mutype",  action='append', type=int, dest='mu_type', help="muon and/or antimuon", default=None)
parser.add_argument("-b", "--heartbeat",  dest="heartbeat", type=int,  help="progress report", default=10000)
parser.add_argument("-n", "--pot",  dest="Np", type=int,  help="proton collisions", default=1000000)
parser.add_argument("-z", "--z",  dest="zStop", type=float,  help="position of the scoring plane [m]", default=19.)
parser.add_argument("-momMag", "--momMag",  dest="momMag", type=float,  help="momentum magnitude [GeV]", default=6800.)
parser.add_argument("-xsingTheta", "--xsingTheta",  dest="xsingTheta", type=float,  help="beam crossing half-angle [urad]", default=-160.)
parser.add_argument('-C', '--charm', action='store_true', dest='charm',help="ccbar production", default=False)
parser.add_argument('-B', '--beauty', action='store_true', dest='beauty',help="bbbar production", default=False)
parser.add_argument('-H', '--hard', action='store_true', dest='hard',help="all hard processes", default=False)
parser.add_argument('-X', '--PDFpSet',dest="PDFpSet",  type=str,  help="PDF pSet to use", default="13")

# for lhapdf, -X LHAPDF6:MMHT2014lo68cl (popular with LHC experiments, features LHC data till 2014)
# one PDF set, which is popular with IceCube, is HERAPDF15LO_EIG
# the default PDFpSet '13' is NNPDF2.3 QCD+QED LO

options = parser.parse_args()
X=ROOT.Pythia8Generator()

munames=''
L=len(options.mu_type)
for i in range(0, L):
  if i==L-1: munames +=mudict[options.mu_type[i]]
  else: munames +=mudict[options.mu_type[i]]+'_'

print(options.mu_type)

# Make pp events
generator = ROOT.Pythia8.Pythia()
generator.settings.mode("Next:numberCount",options.heartbeat)
generator.settings.mode("Beams:idA",  2212)
generator.settings.mode("Beams:idB",  2212)
#the beams are not back-to-back, and therefore the three-momentum of each incoming particle needs to be specified
generator.settings.mode("Beams:frameType", 3)
# Momentum of beam A - it is set using the momMag and beams xsing half-angle
generator.settings.parm("Beams:pxA", 0.)
generator.settings.parm("Beams:pyA", options.momMag*math.sin(1e-6*options.xsingTheta))
generator.settings.parm("Beams:pzA", options.momMag*math.cos(1e-6*options.xsingTheta))
# Momentum of beam B
generator.settings.parm("Beams:pxB", 0.)
generator.settings.parm("Beams:pyB", options.momMag*math.sin(1e-6*options.xsingTheta))
generator.settings.parm("Beams:pzB", -options.momMag*math.cos(1e-6*options.xsingTheta))

# The Monash 2013 tune (#14) is set by default for Pythia above v8.200. 
# This tune provides slightly higher Ds and Bs fractions, in better agreement with the data.
# Tune setting comes before PDF setting!
#generator.readString("Tune:pp = 14")
generator.readString("PDF:pSet = "+options.PDFpSet)
tag = 'nobias'
if options.charm:
     generator.readString("HardQCD:hardccbar = on")
     tag = 'ccbar'
elif options.beauty:
     generator.readString("HardQCD:hardbbbar = on")
     tag = 'bbbar'
elif options.hard:
     generator.readString("HardQCD:all = on")
     tag = 'hard'
else:
     generator.readString("SoftQCD:inelastic = on")
     generator.readString("SoftQCD:singleDiffractive = on")
     generator.readString("SoftQCD:doubleDiffractive = on")
     generator.readString("SoftQCD:centralDiffractive = on")
     generator.readString("PhotonCollision:gmgm2mumu = on")   
     
generator.init()

rc = generator.next()
processes = generator.info.codesHard()
hname = 'pythia8_'+tag+'_PDFpset'+options.PDFpSet+'_'+munames
hname = hname.replace('*','star')
hname = hname.replace('->','to')
hname = hname.replace('/','')

fout = ROOT.TFile(hname+".root","RECREATE")
dTree = ROOT.TTree('MuonTree', munames)
dAnc = ROOT.TClonesArray("TParticle")
# AncstrBranch will hold the neutrino at 0th TParticle entry.
# Neutrino ancestry is followed backwards in evolution up to the colliding TeV proton
# and saved as 1st, 2nd,... TParticle entries of AncstrBranch branch
dAncstrBranch = dTree.Branch("Ancstr",dAnc,32000,-1)
# EvtId will hold event id
evtId = array('i', [0])
dEvtId = dTree.Branch("EvtId", evtId, "evtId/I")

timer = ROOT.TStopwatch()
timer.Start()

nMade = 0
py = generator
for n in range(int(options.Np)):
  rc = py.next()
  for ii in range(1,py.event.size()):    
     evt = py.event[ii]
     # stop the simulation at the scoring plane
     if abs(evt.zProd()) > int(options.zStop)*1e+3: #// pythia units are mm, zStop is in [m]
        dAnc.Clear()
        # take muons only at the scoring plane
        if py.event[ii].id() in options.mu_type:
          print('got a MUON')
          neut = ROOT.TParticle(evt.id(), evt.status(),
                                evt.mother1(),evt.mother2(),0,0,
                                evt.px(),evt.py(),evt.pz(),evt.e(),
                                evt.xProd(),evt.yProd(),evt.zProd(),evt.tProd())
          dAnc[0] = neut
          evtId[0] = n
          gm = py.event[ii].mother1()
          # Chain all mothers (gm)
          while gm:
             evtM = py.event[gm]
             anc = ROOT.TParticle(evtM.id(),evtM.status(),
                                  evtM.mother1(),evtM.mother2(),evtM.daughter1(),evtM.daughter2(),
                                  evtM.px(),evtM.py(),evtM.pz(),evtM.e(),
                                  evtM.xProd(),evtM.yProd(),evtM.zProd(),evtM.tProd())
             nAnc = dAnc.GetEntries()
             if dAnc.GetSize() == nAnc: dAnc.Expand(nAnc+10)
             dAnc[nAnc] = anc
             gm = py.event[gm].mother1()
          dTree.Fill()
  nMade+=1
fout.cd() 
dTree.Write()
         

generator.stat()

timer.Stop()
rtime = timer.RealTime()
ctime = timer.CpuTime()
totalXsec = 0   # unit = mb,1E12 fb
processes = generator.info.codesHard()
for p in processes:
   totalXsec+=generator.info.sigmaGen(p)   
# nobias: 78.4mb, ccbar=4.47mb, bbbar=0.35mb

IntLumi = options.Np / totalXsec * 1E-12

print("Saving to output %s neutrino flavour(s) having PDG ID(s)"%(munames), options.mu_type)
print("simulated events = %i, equivalent to integrated luminosity of %5.2G fb-1. Real time %6.1Fs, CPU time %6.1Fs"%(options.Np,IntLumi,rtime,ctime))
# neutrino CC cross section about 0.7 E-38 cm2 GeV-1 nucleon-1, SND@LHC: 59mm tungsten 
# sigma_CC(100 GeV) = 4.8E-12  
print("corresponding to effective luminosity (folded with neutrino CC cross section at 100GeV) of %5.2G fb-1."%(IntLumi/4.8E-12))


def debugging(g):
   generator.settings.listAll()
