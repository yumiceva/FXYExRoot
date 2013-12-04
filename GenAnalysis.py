#! /usr/bin/env python

"""
   top.py

   Script to create histograms and cut flow table for lepton+jets.
   
   Francisco Yumiceva (yumiceva@gmail.com)
   Florida Institute of Technology, 2013
   
"""

class cTerm:
    PURPLE = '\033[95m'
    CYAN = '\033[96m'
    DARKCYAN = '\033[36m'
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    END = '\033[0m'                              

# External packages
import sys
import os
import math
import re
from optparse import OptionParser

# Check if pyROOT is available
try:
    from ROOT import *
    print cTerm.BOLD+"ROOT loaded"+cTerm.END
    
except:
    print cTerm.RED+"\nError: Cannot load PYROOT, make sure you have setup ROOT in the path"
    print "and pyroot library is also defined in the variable PYTHONPATH, try:\n"+cTerm.END
    if (os.getenv("PYTHONPATH")):
        print " setenv PYTHONPATH ${PYTHONPATH}:$ROOTSYS/lib\n"
    else:
        print " setenv PYTHONPATH $ROOTSYS/lib\n"
    print "Exit now\n"
    sys.exit()

def main():
    
    # OPTIONS
    usage = "usage: %prog [options] arg"
    parser = OptionParser(usage)
    parser.add_option("-b", "--batch",
                      action="store_true", dest="batch")
    parser.add_option("-s", "--sample",
                      default="madgraph",
                      help="input samples options are: madgraph, whizard [default: %default]")

    (options, args) = parser.parse_args()
    #print options
    #print args
    return options

if __name__ == '__main__':

    # parser
    options = main()
    
    if options.batch:
        ROOT.gROOT.SetBatch()
        print "running ROOT in batch mode."
    
    # Load ROOT libraries
    gSystem.Load('/uscms/home/yumiceva/work/sframe/CMSSW_5_3_3/src/ExRootAnalysis/lib/libExRootAnalysis.so')

    # Create output root file
    outname = "results_gen.root"
    if options.sample == "madgraph":
        outname = "results_gen_madgraph.root"
    if options.sample == "whizard":
        outname= "results_gen_whizard.root"
    
    outFile = TFile(outname,"RECREATE")
    
    # Create chain of root trees
    chain = TChain("LHEF")
    maxEntries = -1
    # MG files
    if options.sample = "madgraph":
        chain.Add("/uscms_data/d2/maravin/TTG_MG5/Two2Seven/ROOT/ttgamma_27_part1.root")
        chain.Add("/uscms_data/d2/maravin/TTG_MG5/Two2Seven/ROOT/ttgamma_27_part2.root")
    ##chain.Add("/uscms_data/d2/maravin/TTG_MG5/Two2Seven/ROOT/ttgamma_27_part3.root")
    ##chain.Add("/uscms_data/d2/maravin/TTG_MG5/Two2Seven/ROOT/ttgamma_27_part4.root")

    # Whizard files
    if options.sample = "whizard":
        chain.Add("/uscmst1b_scratch/lpc1/cmsroc/yumiceva/TTGamma/LHE/whizard/TTGamma_Whizard_2to7/ttgamma.root")
        maxEntries = 200000
        
    # setup ntuple object
    treeReader = ExRootTreeReader(chain)
    # number of entries
    numberOfEntries = treeReader.GetEntries()
    print "Total number of entries to be processed: " + str(numberOfEntries)

    # Get pointers to branches used in this analysis
    Particles = treeReader.UseBranch("Particle")

    # Book histograms
    h_nocut = {}
    h_cut = {}
    h_nocut['pt'] = TH1F("pt", "p_{T} [GeV]", 100, 0.0, 100.0)
    h_nocut['PID'] = TH1F("PID","Particle ID",50,-25,25)
    h_nocut['ttbar_prod'] = TH1F("ttbar_prod","ttbar production",4,0,4)
    h_nocut['ttbar_BR'] = TH1F("ttbar_BR","BR",6,0,6)
    h_nocut['photon_mom'] = TH1F("photon_mom","photon mother",50,-25,25)
    h_nocut['photon_pt'] = TH1F("photon_pt","Photon p_{T} [GeV]",100,0,100)
    h_nocut['photon_eta'] = TH1F("photon_eta","Photon #eta",100,-5,5)
    h_nocut['photon_phi'] = TH1F("photon_phi","Photon #phi",80,-3.2,3.2)
    h_nocut['photon_deltaRLep'] = TH1F('photon_deltaRLep',"#Delta R(#gamma,lepton)",100,0,6)
    h_nocut['N_photons'] = TH1F("N_photons","Number of Photons",5,0,5)
    h_nocut['top_pt'] = TH1F("top_pt","top p_{T} [GeV]",100,0,200)
    h_nocut['top_eta'] =TH1F("top_eta","top #eta",100,-5,5)
    h_nocut['top_phi'] =TH1F("top_phi","top #phi",80,-3.2,3.2)
    h_nocut['top_m'] = TH1F("top_m","top mass [GeV]",50,150,200)
    h_nocut['W_m'] = TH1F("W_m","W mass [GeV]",50,50,100)
    h_nocut['b_m'] = TH1F("b_m","b mass [GeV]",50,0,10)
    h_nocut['b_pt'] = TH1F("b_pt","b p_{T} [GeV]",100,0,200)
    h_nocut['b_eta'] = TH1F("b_eta","b #eta",100,-5,5)
    h_nocut['b_phi'] =TH1F("b_phi","b #phi",80,-3.2,3.2)
    h_nocut['q_pt'] = TH1F("q_pt","b p_{T} [GeV]",100,0,200)
    h_nocut['q_eta'] = TH1F("q_eta","b #eta",100,-5,5)
    h_nocut['q_phi'] =TH1F("q_phi","b #phi",80,-3.2,3.2)
    h_nocut['nu_pt'] = TH1F("nu_pt","Neutrino p_{T} [GeV]",100,0,200)
    h_nocut['nu_eta'] = TH1F("nu_eta","Neutrino #eta",100,-5,5)
    h_nocut['nu_phi'] =TH1F("nu_phi","Neutrino #phi",80,-3.2,3.2)
    h_nocut['lep_pt'] = TH1F("lep_pt","Lepton p_{T} [GeV]",100,0,200)
    h_nocut['lep_eta'] = TH1F("lep_eta","Lepton #eta",100,-5,5)
    h_nocut['lep_phi'] =TH1F("lep_phi","Lepton #phi",80,-3.2,3.2)
            
    for key in h_nocut.keys():
        h_cut[key] = h_nocut[key].Clone(h_nocut[key].GetName())
        h_cut[key].Sumw2()
        h_nocut[key].Sumw2()
        h_cut[key].SetTitle( h_cut[key].GetTitle() )
        h_nocut[key].SetTitle( h_nocut[key].GetTitle() )
        
    # Loop over all events
    for entry in xrange(0, numberOfEntries):

        # Load selected branches with data from specified event
        treeReader.ReadEntry(entry)

        if entry%1000 == 0:
            print "entry=",entry
            #### 
            if maxEntries > entry: break
        # Check ttbar production mechanism
        ttbarMech = 0
        if Particles[0].PID == 21 and Particles[1].PID == 21: ttbarMech = 1
        elif Particles[0].PID == 21 and Particles[1].PID != 21: ttbarMech = 2
        elif Particles[0].PID != 21 and Particles[1].PID == 21: ttbarMech = 2
        else: ttbarMech = 3
        h_nocut['ttbar_prod'].Fill( ttbarMech )
        
        # TRootLHEFParticle
        index = 0
        p4Electron = TLorentzVector()
        p4Muon = TLorentzVector()
        p4Photon = TLorentzVector()
        p4Lepton = TLorentzVector()
        
        numElectrons = 0
        numMuons = 0
        numTaus = 0
        numPhotons =0
        
        for p in Particles:
            index += 1
            
            h_nocut['PID'].Fill( p.PID )
            
            # MG Status code: -1 initial, 2 intermediate, 1 final
            
            # Electrons
            if math.fabs(p.PID) == 11 and p.Status == 1:
                if numElectrons == 0:
                    p4Electron.SetPtEtaPhiE( p.PT, p.Eta, p.Phi, p.E)
                    p4Lepton.SetPtEtaPhiE( p.PT, p.Eta, p.Phi, p.E)
                numElectrons += 1
            # Muons
            if math.fabs(p.PID) == 13 and p.Status == 1:
                if numMuons == 0:
                    p4Muon.SetPtEtaPhiE( p.PT, p.Eta, p.Phi, p.E)
                    p4Lepton.SetPtEtaPhiE( p.PT, p.Eta, p.Phi, p.E)
                numMuons += 1
            # Taus
            if math.fabs(p.PID) == 15 and p.Status == 1:
                #p4Muon.SetPtEtaPhiE( p.PT, p.Eta, p.Phi, p.E)
                numTaus += 1
            # Photons
            if p.PID == 22 and p.Status == 1:
                if numPhotons==0:
                    p4Photon.SetPtEtaPhiE( p.PT, p.Eta, p.Phi, p.E)
                numPhotons += 1
                motherIndex = p.Mother1 - 1
                pMom = Particles[motherIndex]
                h_nocut['photon_mom'].Fill( pMom.PID )
            # top
            if math.fabs(p.PID) == 6:
                h_nocut['top_m'].Fill( p.M )
                h_nocut['top_pt'].Fill( p.PT )
                h_nocut['top_eta'].Fill( p.Eta )
                h_nocut['top_phi'].Fill( p.Phi )
            # W
            if math.fabs(p.PID) == 24:
                h_nocut['W_m'].Fill( p.M )
            # b quark
            if math.fabs(p.PID) == 5:
                h_nocut['b_m'].Fill( p.M )
                h_nocut['b_pt'].Fill( p.PT )
                h_nocut['b_eta'].Fill( p.Eta )
                h_nocut['b_phi'].Fill( p.Phi )
            # light quarks
            if math.fabs(p.PID) > 0 and math.fabs(p.PID) < 5:
                h_nocut['q_pt'].Fill( p.PT )
                h_nocut['q_eta'].Fill( p.Eta )
                h_nocut['q_phi'].Fill( p.Phi )
            # neutrino
            if (math.fabs(p.PID) == 12 or math.fabs(p.PID) == 14) and p.Status == 1:
                p4Nu.SetPtEtaPhiE( p.PT, p.Eta, p.Phi, p.E )
                                
            # test histogram
            h_nocut['pt'].Fill( p.PT )
            #print p.PT
        # end loop over Particles

        # Check BR
        BR_code = 0
        if (numElectrons > 1 and numMuons ==0 and numTaus ==0) or \
               (numElectrons ==0 and numMuons > 1 and numTaus ==0) or \
               (numElectrons == 1 and numMuons == 1 and numTaus==0) or \
               (numElectrons ==0 and numMuons ==0 and numTaus > 1) or \
               (numElectrons ==1 and numMuons ==0 and numTaus==1) or \
               (numElectrons ==0 and numMuons==1 and numTaus ==1):
            # Dileptons
            BR_code = 1
        elif numElectrons == 1 and numMuons == 0 and numTaus ==0:
            # e+jets
            BR_code = 2
        elif numElectrons == 0 and numMuons == 1 and numTaus ==0:
            # mu+jets
            BR_code = 3
        elif numElectrons == 0 and numMuons == 0 and numTaus ==0:
            # All jets
            BR_code = 5
        else:
            # tau+jets
            BR_code = 4

        h_nocut['ttbar_BR'].Fill( BR_code )

        h_nocut['N_photons'].Fill(numPhotons)
        if numPhotons > 0 and (BR_code == 2 or BR_code == 3):
            h_nocut['photon_pt'].Fill(p4Photon.Pt())
            h_nocut['photon_eta'].Fill(p4Photon.Eta())
            h_nocut['photon_phi'].Fill(p4Photon.Phi())
            # DeltaR
            deltaRwLep = p4Photon.DeltaR( p4Lepton )
            h_nocut['photon_deltaRLep'].Fill( deltaRwLep )

            h_nocut['lep_pt'].Fill( p4Lepton.Pt() )
            h_nocut['lep_eta'].Fill( p4Lepton.Eta() )
            h_nocut['lep_phi'].Fill( p4Lepton.Phi() )

            h_nocut['nu_pt'].Fill( p.PT )
            h_nocut['nu_eta'].Fill( p.Eta )
            h_nocut['nu_phi'].Fill( p.Phi )
            
    # END loop over entries
    h_nocut['pt'].Draw()

    outFile.cd()
    for key in h_nocut.keys():
        
        if h_nocut[key].GetEntries() > 0:
            h_nocut[key].Write()
        if h_cut[key].GetEntries() > 0:
            h_cut[key].Write()
        
    
    # Wait
    rep = ''
    while not rep in [ 'q', 'Q', '.q', 'qq' 'p']:
        rep = raw_input( '\nenter: ["q",".q" to quit] ["p" or "print" to print all canvas]: ' )
        if 0<len(rep):
            if rep=='quit': rep = 'q'
            #if rep=='p' or rep=='print':

    #del(treeReader)
    #del(chain)
    del treeReader
    print "done."
