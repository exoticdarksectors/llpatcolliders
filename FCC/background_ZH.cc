// main381.cc is a part of the PYTHIA event generator.
// Copyright (C) 2025 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Keywords: Higgs; electron-positron;

// Authors: Torbjörn Sjostrand <torbjorn.sjostrand@fysik.lu.se>

// Simple example of Higgs pruduction at future e+e- colliders.

#include "Pythia8/Pythia.h"

//include ROOT functions to generate histograms
#include "TH1D.h"
#include "TFile.h"
#include "TApplication.h"
#include "TVirtualPad.h"
#include "TMath.h"
#include "TTree.h"

//include FastJet3 functions to perform jet clustering
#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"



using namespace Pythia8;

//==========================================================================



int main() {
    
    // Number of events.
    int nEvent = 50000;
    
    // Generator. Incoming beams. (Switch off iniial-state photon radiation.)
    Pythia pythia;
    pythia.readString("Beams:idA = -11");
    pythia.readString("Beams:idB = 11");
    pythia.readString("Beams:eCM = 240.");
    //pythia.readString("PDF:lepton = off");

    // All Higgs production channels.
    pythia.readString("HiggsSM:all = on");
    pythia.readString("25:onMode = 0");
    //pythia.readString("23:onMode = 1");
    pythia.readString("25:onIfAny = 5");
    pythia.readString("23:onMode = off");
    pythia.readString("23:onIfAny = 1 2 3 4 5");
    //pythia.readString("PartonLevel:ISR = off");
    //pythia.readString("PartonLevel:FSR = off");
    //pythia.readString("TimeShower:QEDshowerByL=off");
    //pythia.readString("PDF:lepton = off");

    // If Pythia fails to initialize, exit with error.
    if (!pythia.init()) return 1;
    
    // Create file on which histogram(s) can be saved.
    TFile* outFile = new TFile("bg_ZH_2b.root", "RECREATE");

    TTree *t1 = new TTree("t1","t1");
    TTree *t_jet = new TTree("t_jet","t_jet");
    
    //Set up the TTree and define all the variables needed
    
    int pid;
    int MC_event;
    double x,y,z,t,energy,phi,theta,px,py,pz;
    std::vector<int> MotherList,DaughterListRec,DaughterList,SisterList;
    int list,status;
    bool isFinal,isCharged,isCarbonCopy;
    int daughter1,daughter2;
    int jetAssignment;
    int final_carbon_copy;
    
    
    t1->Branch("energy",&energy,"energy/D");
    t1->Branch("x",&x,"x/D");
    t1->Branch("y",&y,"y/D");
    t1->Branch("z",&z,"z/D");
    t1->Branch("t",&t,"t/D");
    t1->Branch("pid",&pid,"pid/I");
    t1->Branch("phi",&phi,"phi/D");
    t1->Branch("theta",&theta,"theta/D");
    t1->Branch("px",&px,"px/D");
    t1->Branch("py",&py,"py/D");
    t1->Branch("pz",&pz,"pz/D");
    t1->Branch("MC_event",&MC_event,"MC_event/I");
    t1->Branch("MotherList",&MotherList);
    t1->Branch("list",&list);
    t1->Branch("status",&status);
    t1->Branch("isFinal",&isFinal);
    t1->Branch("isCharged",&isCharged);
    t1->Branch("DaughterListRec",&DaughterListRec);
    t1->Branch("DaughterList",&DaughterList);
    t1->Branch("daughter1",&daughter1);
    t1->Branch("daughter2",&daughter2);
    t1->Branch("jetAssignment",&jetAssignment);
    t1->Branch("isCarbonCopy",&isCarbonCopy);
    t1->Branch("final_carbon_copy",&final_carbon_copy);
    t1->Branch("SisterList",&SisterList);
    
    int num_event;
    int num_jet;
    double jet_e,jet_px,jet_py,jet_pz;
    std::vector<int> constituents;
    
    t_jet->Branch("num_event",&num_event);
    t_jet->Branch("num_jet",&num_jet);
    t_jet->Branch("jet_e",&jet_e);
    t_jet->Branch("jet_px",&jet_px);
    t_jet->Branch("jet_py",&jet_py);
    t_jet->Branch("jet_pz",&jet_pz);
    t_jet->Branch("constituents",&constituents);
    
    
    //set up for fast jet
    double d_cut = 10.0;
    double theta_cut = 0.154;
    int nSel =2;
    
    fastjet::JetDefinition          jetDef = fastjet::JetDefinition(fastjet::ee_kt_algorithm);
    
    //fastjet input
    std::vector <fastjet::PseudoJet> fjInputs;
    
    bool firstEvent=true;
    int previous_size=0;
    
    
    // Begin event loop. Generate event. Skip if error.
    for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
        if (!pythia.next()) continue;
        
        
        for (int ijet = 0; ijet < pythia.event.size();++ijet){
            if(!pythia.event[ijet].isFinal()) continue;
            
            if ( pythia.event[ijet].idAbs() == 12 || pythia.event[ijet].idAbs() == 14
                || pythia.event[ijet].idAbs() == 16 ) continue;
            
            if ((pythia.event[ijet].theta() < theta_cut) || ((TMath::Pi()-pythia.event[ijet].theta() < theta_cut))) continue;
            
            fastjet::PseudoJet fastjet_input_particle;
            fastjet_input_particle = fastjet::PseudoJet( pythia.event[ijet].px(),
                                                        pythia.event[ijet].py(), pythia.event[ijet].pz(), pythia.event[ijet].e() );
            fastjet_input_particle.set_user_index(ijet);
            fjInputs.push_back(fastjet_input_particle);
        }
        
        //run fast jet algorithm
        vector <fastjet::PseudoJet> exclusiveJets, sortedJets;
        fastjet::ClusterSequence clustSeq(fjInputs, jetDef);
        
        if (firstEvent) {
          cout << "Ran " << jetDef.description() << endl;
          cout << "Strategy adopted by FastJet was "
               << clustSeq.strategy_string() << endl << endl;
          firstEvent = false;
        }
        
        // Extract inclusive jets sorted by E(d_cut=10.0 GeV mode).
        exclusiveJets = clustSeq.exclusive_jets(pow(d_cut,2));
        sortedJets    = sorted_by_E(exclusiveJets);
        fjInputs.resize(0);
        
        
        for (int jet_tree=0; jet_tree<sortedJets.size(); ++jet_tree ){
            num_event = previous_size;
            jet_e = sortedJets[jet_tree].e();
            jet_px = sortedJets[jet_tree].px();
            jet_py = sortedJets[jet_tree].py();
            jet_pz = sortedJets[jet_tree].pz();
            vector<fastjet::PseudoJet> jet_constituents = sortedJets[jet_tree].constituents();
            num_jet = jet_constituents.size();
            for (int j=0; j<jet_constituents.size();++j){
                constituents.push_back(jet_constituents[j].user_index());
            }
            t_jet->Fill();
            constituents.resize(0);
            
        }
        
        previous_size+=pythia.event.size();
 
        
        for (int ipt = 0; ipt < pythia.event.size();++ipt){ //store MC information
            pid = pythia.event[ipt].id();
            x = pythia.event[ipt].xProd();
            y = pythia.event[ipt].yProd();
            z = pythia.event[ipt].zProd();
            t = pythia.event[ipt].tProd();
            phi = pythia.event[ipt].phi();
            theta = pythia.event[ipt].theta();
            energy = pythia.event[ipt].e();
            px = pythia.event[ipt].px();
            py = pythia.event[ipt].py();
            pz = pythia.event[ipt].pz();
            list = ipt;
            status = pythia.event[ipt].status();
            MC_event = iEvent;
            MotherList = pythia.event[ipt].motherList();
            DaughterListRec = pythia.event[ipt].daughterListRecursive();
            isFinal = pythia.event[ipt].isFinal();
            isCharged = pythia.event[ipt].isCharged();
            SisterList = pythia.event[ipt].sisterList();
            
            if(pythia.event[ipt].daughter1()==pythia.event[ipt].daughter2() && pythia.event[ipt].daughter1()>0){
                int q;
                q = pythia.event[ipt].daughter1();
                while(pythia.event[q].daughter1()==pythia.event[q].daughter2() && pythia.event[q].daughter1()>0){
                    q = pythia.event[q].daughter1();
                }
                daughter1 = pythia.event[q].daughter1();
                daughter2 = pythia.event[q].daughter2();
                DaughterList = pythia.event[q].daughterList();
                isCarbonCopy = 1;
                final_carbon_copy = q;
            }
            
            else{
                daughter1 = pythia.event[ipt].daughter1();
                daughter2 = pythia.event[ipt].daughter2();
                DaughterList = pythia.event[ipt].daughterList();
                isCarbonCopy = 0;
                final_carbon_copy = ipt;
            }
            
            t1->Fill();
        }
        
        
        
        
        
        
    }
    // Statistics on event generation.
    pythia.stat();

    // Write everything into a root file
    t1->Write();
    t_jet->Write();
    outFile->Close();
    
    delete outFile;
    return 0;
    
}
