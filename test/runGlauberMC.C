#include <TRandom3.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TTree.h>
#include "NA6PNuclearShape.h"
#include "NA6PGlauber.h"

void runGlauberMC(int tottry = 20000,
                  float bmin = -1.,
                  float bmax = -1.)
{

  NA6PNuclearShape nucA, nucB;
  nucA.setNucleusPb();
  nucB.setNucleusPb();
  float sigmb = 30.;
  nucA.setMinimumDistanceBetweenNucleons(0.4);
  nucB.setMinimumDistanceBetweenNucleons(0.4);

  NA6PGlauber gl(nucA, nucB, sigmb);

  TH1F* hNColl = new TH1F("hNColl", ";N_{coll};entries", 100, 0, 2000);
  TH1F* hNPart = new TH1F("hNPart", ";N_{part};entries", 100, 0, 500);
  TGraph* gNColl = new TGraph(0);
  TGraph* gNPart = new TGraph(0);
  TTree* tree = new TTree("treeEvents", "Nuclear collision events");
  float b, phib, ecc;
  int nPart, nSpectA, nSpectB, nColl;
  tree->Branch("b", &b, "b/F");
  tree->Branch("phib", &phib, "phib/F");
  tree->Branch("nPart", &nPart, "nPart/I");
  tree->Branch("nSpectA", &nSpectA, "nSpectA/I");
  tree->Branch("nSpectB", &nSpectB, "nSpectB/I");
  tree->Branch("nColl", &nColl, "nColl/I");
  tree->Branch("eccentr", &ecc, "ecc/F");

  for (int itry = 0; itry < tottry; itry++) {
    NA6PGlauber::GlauberEvent ev = gl.generateEvent(bmin, bmax);
    nColl = ev.nCollisions;
    nPart = ev.nParticipants;
    nSpectA = ev.nSpectatorsA;
    nSpectB = ev.nSpectatorsB;
    b = ev.impactParameter;
    phib = ev.phiImpPar;
    ecc = ev.eccentricity;
    if ((itry % 2000) == 0)
      printf("Event %d\n", itry);
    hNColl->Fill(nColl);
    hNPart->Fill(nPart);
    gNColl->SetPoint(itry, b, nColl);
    gNPart->SetPoint(itry, b, nPart);
    tree->Fill();
  }

  TFile* filout = new TFile("glauberMC.root", "recreate");
  filout->cd();
  tree->Write();
  filout->Close();

  TCanvas* c1 = new TCanvas("c1", "", 1200, 600);
  c1->Divide(2, 1);
  c1->cd(1);
  gPad->SetLogy();
  hNPart->Draw();
  c1->cd(2);
  gPad->SetLogy();
  hNColl->Draw();

  gNColl->SetMarkerStyle(7);
  gNPart->SetMarkerStyle(7);

  TCanvas* c2 = new TCanvas("c2", "", 1400, 700);
  c2->Divide(2, 1);
  c2->cd(1);
  gNColl->Draw("AP");
  gNColl->GetXaxis()->SetTitle("Impact parameter (fm)");
  gNColl->GetYaxis()->SetTitle("N_{coll}");
  c2->cd(2);
  gNPart->Draw("AP");
  gNPart->GetXaxis()->SetTitle("Impact parameter (fm)");
  gNPart->GetYaxis()->SetTitle("N_{part}");
}
