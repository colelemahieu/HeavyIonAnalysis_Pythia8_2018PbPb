//miniAOD_mcskim.C

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include "TString.h"
#include "TTree.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TBranch.h"
#include "TFile.h"
#include <TTree.h>
#include <TChain.h>
#include "TCanvas.h"
#include "THStack.h"
#include "TLine.h"
#include "TMath.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TLorentzVector.h"
#include <cmath>
using namespace std;

void miniAOD_mcSkim(string inputfile, string outfile)
{    
  TChain* calochain = new TChain("ggHiNtuplizer/EventTree");
  //TChain* jetchain = new TChain("ak2PFJetAnalyzer/t");
  TChain* jetchain = new TChain("ak4PFJetAnalyzer/t");
  //TChain* jetchain = new TChain("ak6PFJetAnalyzer/t");
  //TChain* jetchain = new TChain("ak8PFJetAnalyzer/t");
  TChain* hiEvtchain = new TChain("hiEvtAnalyzer/HiTree");
  TChain* pfchain = new TChain("particleFlowAnalyser/pftree");

  //calochain->Add(inputfile.c_str());
  jetchain->Add(inputfile.c_str());
  hiEvtchain->Add(inputfile.c_str());
  pfchain->Add(inputfile.c_str());

  //int caloentries = calochain->GetEntries();
  int jetentries = jetchain->GetEntries();
  int hiEvtEntries = hiEvtchain->GetEntries();
  //int trackEntries = trackchain->GetEntries();
  //cout << "Entries in calo chain = " << caloentries << endl;
  cout << "Entries in jet chain = " << jetentries << endl;
  cout << "Entries in hiEvt chain = " << hiEvtEntries << endl;
  //cout << "Entries in track chain = " << trackEntries << endl;

  // Tree Variables
  int nRef=0, nGen=0, nTrk=0, nTower=0, n=0, nPFpart=0, vtxFilter=0, clusterFilter=0;
  unsigned int run=0, lumis=0;
  ULong64_t event=0;
  float vtx_z=0;
  float trkPt[10000], trkEta[10000];
  float jetPt[200]={0}, jetEta[200]={0}, jetPhi[200]={0}, jetM[200]={0};
  float genPt[200]={0}, genEta[200]={0}, genPhi[200]={0}, genM[200]={0};

  vector<float> *jetPt_vec=0, *jetEta_vec=0, *jetPhi_vec=0;
  vector<float> *caloEta=0, *caloPhi=0, *caloEn=0, *calo_hadE=0, *calo_emE=0;
  vector<float> *trkPt_vec=0, *trkEta_vec=0;
  vector<float> *pfEnergy=0, *pfEta=0, *pfPhi=0, *pfPt=0, *pfEt=0, *pfE=0, *pfTrkEta=0, *pfTrkPt=0;
  vector<int> *pfId=0;
  //vector<bool> *highPurity=0;

  // miniAOD read Variables
  float hiHE_pfle=0, hiHEMinus_pfle=0, hiHEPlus_pfle=0;
  float hiEE_pfle=0, hiEEMinus_pfle=0, hiEEPlus_pfle=0; 

  hiEvtchain->SetBranchAddress("vz",&vtx_z);

  jetchain->SetBranchAddress("nref", &nRef);
  jetchain->SetBranchStatus("*",0);
  jetchain->SetBranchStatus("jtpt",1);
  jetchain->SetBranchStatus("jteta",1);
  jetchain->SetBranchStatus("jtphi",1);
  jetchain->SetBranchStatus("jtm",1);
  jetchain->SetBranchStatus("ngen",1);
  jetchain->SetBranchStatus("genpt",1);
  jetchain->SetBranchStatus("geneta",1);
  jetchain->SetBranchStatus("genphi",1);
  jetchain->SetBranchStatus("genm",1);
  jetchain->SetBranchAddress("jtpt", &jetPt);
  jetchain->SetBranchAddress("jteta", &jetEta);
  jetchain->SetBranchAddress("jtphi", &jetPhi);
  jetchain->SetBranchAddress("jtm", &jetM);
  jetchain->SetBranchAddress("ngen", &nGen);
  jetchain->SetBranchAddress("genpt", &genPt);
  jetchain->SetBranchAddress("geneta", &genEta);
  jetchain->SetBranchAddress("genphi", &genPhi);
  jetchain->SetBranchAddress("genm", &genM);

  /*calochain->SetBranchStatus("*",0);
  calochain->SetBranchStatus("run",1);
  calochain->SetBranchStatus("event",1);
  calochain->SetBranchStatus("lumis",1);
  calochain->SetBranchStatus("CaloTower_e", 1);*/
  /*calochain->SetBranchStatus("CaloTower_eta",1);
  calochain->SetBranchStatus("CaloTower_phi",1);
  calochain->SetBranchStatus("CaloTower_hadE",1);
  calochain->SetBranchStatus("CaloTower_emE",1);
  calochain->SetBranchAddress("run",&run);
  calochain->SetBranchAddress("event",&event);
  calochain->SetBranchAddress("lumis",&lumis);
  calochain->SetBranchAddress("CaloTower_e", &caloEn);
  calochain->SetBranchAddress("CaloTower_eta", &caloEta);
  calochain->SetBranchAddress("CaloTower_phi", &caloPhi);
  calochain->SetBranchAddress("CaloTower_hadE", &calo_hadE);
  calochain->SetBranchAddress("CaloTower_emE", &calo_emE);*/

  /*trackchain->SetBranchStatus("*",0);
  trackchain->SetBranchStatus("nTrk",1);
  trackchain->SetBranchStatus("trkPt",1);
  trackchain->SetBranchStatus("trkEta",1);
  trackchain->SetBranchAddress("nTrk", &nTrk);
  trackchain->SetBranchAddress("trkPt", &trkPt_vec);
  trackchain->SetBranchAddress("trkEta", &trkEta_vec);*/
  
  pfchain->SetBranchStatus("*",0);
  pfchain->SetBranchStatus("nPF",1);
  pfchain->SetBranchStatus("pfPt",1);
  pfchain->SetBranchStatus("pfEt",1);
  pfchain->SetBranchStatus("pfE",1);
  pfchain->SetBranchStatus("pfEta",1);
  pfchain->SetBranchStatus("pfPhi",1);
  pfchain->SetBranchStatus("pfId",1);
  pfchain->SetBranchStatus("pfTrkEta",1);
  pfchain->SetBranchStatus("pfTrkPt",1);
  //pfchain->SetBranchStatus("highPurity",1);
  pfchain->SetBranchAddress("nPF", &nPFpart);
  pfchain->SetBranchAddress("pfPt", &pfPt);
  pfchain->SetBranchAddress("pfEt", &pfEt);
  pfchain->SetBranchAddress("pfE", &pfE);
  pfchain->SetBranchAddress("pfEta", &pfEta);
  pfchain->SetBranchAddress("pfPhi", &pfPhi);
  pfchain->SetBranchAddress("pfId", &pfId);
  pfchain->SetBranchAddress("pfTrkEta", &pfTrkEta);
  pfchain->SetBranchAddress("pfTrkPt", &pfTrkPt);
  //pfchain->SetBranchAddress("highPurity", &highPurity);


  // Clone Trees
  TFile newfile(outfile.c_str(),"recreate");
  TTree *pfClone = pfchain->CloneTree(0);
  //TTree *trkClone = trackchain->CloneTree(0);
  TTree *jetClone = jetchain->CloneTree(0);
  //TTree *caloClone = calochain->CloneTree(0);
  

  // New ROOT File
  TTree *hiEvent = new TTree("hiEvent","hiEvent");
  hiEvent->Branch("Vertex_Z", &vtx_z);

  int goodEvt=0;
  int pfGoodTrk, pfHadronYes;

  // Event Loop for All events
  //for (int i=0; i<jetentries; i++)
  int halfEntries = ceil(jetentries/2);
  int thirdEntries = ceil(jetentries/3);
  //for (int i=0; i<jetentries; i++)
  for (int i=0; i<thirdEntries; i++)
    {
      jetchain->GetEntry(i);
      hiEvtchain->GetEntry(i);
      pfchain->GetEntry(i);
     
      pfGoodTrk=0, pfHadronYes=0;

      // Cuts
      if (fabs(vtx_z) > 20) continue;

      
      for (int pf=0; pf<nPFpart; pf++)
      {
        if (pfId->at(pf)==1) pfHadronYes = pfHadronYes + 1;
      }     

      if (pfHadronYes==0) continue; 

      
      // RECO
      if (nRef==0) continue;
      if (nRef==1) continue;
      if ((fabs(jetEta[0]) > 2.4) || (fabs(jetEta[1]) > 2.4)) continue;
      if (jetPt[0]<30 || jetPt[1]<20) continue;

      // find which jet really is leading

      // find distance between gen0 and reco0, reco1
      /*float phiDiff[2]={0}, etaDiff[2]={0};
      float jetDist[2]={0};
      int indexMin=0;

      for (int nJet=0; nJet<2; nJet++)
	{
	  phiDiff[nJet] = genPhi[0] - jetPhi[nJet];
	  etaDiff[nJet] = genEta[0] - jetEta[nJet];
	  jetDist[nJet] = sqrt(phiDiff[nJet]*phiDiff[nJet] + etaDiff[nJet]*etaDiff[nJet]);
	  if (jetDist[nJet] < jetDist[0]) indexMin=nJet;
	}

      // indexMin=0
      if (indexMin==0)
	{
          if (jetPt[0]<30 || jetPt[1]<20) continue;
        }
      // indexMin=1
      if (indexMin==1)
	{
          if (jetPt[1]<30 || jetPt[0]<20) continue;
        }*/
     
      // Fill Branches
      pfClone->Fill();     
      //trkClone->Fill(); 
      jetClone->Fill();
      //caloClone->Fill();
      hiEvent->Fill();
   
      goodEvt=  goodEvt+ 1;
    }
 
  cout << "goodEvt=" << goodEvt << endl;
  newfile.Write();
  newfile.Close();
}
