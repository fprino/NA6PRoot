#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TParticle.h>
#include "Propagator.h"
#include "NA6PVertex.h"
#include "NA6PVertexerTracks.h"
#include "NA6PMCEventHeader.h"
#endif

void runVertexerTracks(const char* dirSimu = ".")
{
  if (!Propagator::loadField() || !Propagator::loadGeometry(Form("%s/geometry.root", dirSimu))) {
    return;
  }

  TFile* ft = new TFile(Form("%s/TracksVerTel.root", dirSimu));
  if (!ft)
    return;
  TTree* trTree = (TTree*)ft->Get("tracksVerTel");
  std::vector<NA6PTrack>* trArr = nullptr;
  trTree->SetBranchAddress("VerTel", &trArr);

  TFile* fk = new TFile(Form("%s/MCKine.root", dirSimu));
  TTree* mcTree = (TTree*)fk->Get("mckine");
  std::vector<TParticle>* mcArr = nullptr;
  mcTree->SetBranchAddress("tracks", &mcArr);
  NA6PMCEventHeader* mcHead = nullptr;
  mcTree->SetBranchAddress("header", &mcHead);

  TH1F* hncontr = new TH1F("hncontr", ";N_{contributors}", 100, -0.5, 999.5);
  TH1F* hnvert = new TH1F("hnvert", ";Number of reco vertices", 11, -0.5, 10.5);
  TH2F* hxrecgen = new TH2F("hxrecgen", ";x_{gen} (cm);x_{rec} (cm)", 1000, -0.5, 0.5, 1000, -0.5, 0.5);
  TH2F* hyrecgen = new TH2F("hyrecgen", ";y_{gen} (cm);y_{rec} (cm)", 1000, -0.5, 0.5, 1000, -0.5, 0.5);
  TH2F* hzrecgen = new TH2F("hzrecgen", ";z_{gen} (cm);z_{rec} (cm)", 100, -6., 1., 100, -6., 1.);
  TH1F* hxrec = new TH1F("hxrec", ";x_{rec,main} (cm)", 100, -0.1, 0.1);
  TH1F* hyrec = new TH1F("hyrec", ";y_{rec,main} (cm)", 100, -0.1, 0.1);
  TH1F* hzrec = new TH1F("hzrec", ";z_{rec,main} (cm)", 200, -6., 1.);
  TH1F* hdx = new TH1F("hdx", ";x_{rec} - x_{gen} (#mum)", 100, -40., 40.);
  TH1F* hdy = new TH1F("hdy", ";y_{rec} - y_{gen} (#mum)", 100, -40., 40.);
  TH1F* hdz = new TH1F("hdz", ";z_{rec} - z_{gen} (#mum)", 100, -200., 200.);
  TH1F* hsigx = new TH1F("hsigx", ";#sigma_{x} (#mum)", 100, 0., 50.);
  TH1F* hsigy = new TH1F("hsigy", ";#sigma_{y} (#mum)", 100, 0., 50.);
  TH1F* hsigz = new TH1F("hsigz", ";#sigma_{z} (#mum)", 100, 0., 250.);
  TH1F* hpullx = new TH1F("hpullx", ";(x_{rec} - x_{gen}) / #sigma_{x}", 100, -5., 5.);
  TH1F* hpully = new TH1F("hpully", ";(y_{rec} - y_{gen}) / #sigma_{y}", 100, -5., 5.);
  TH1F* hpullz = new TH1F("hpullz", ";(z_{rec} - z_{gen}) / #sigma_{z}", 100, -5., 5.);
  TH1F* hzRecDiff = new TH1F("hzRecDiff", ";z_{rec}^{pil} - z_{rec}^{main} (cm)", 100, -5, 5);
  TH1F* hcontRatio = new TH1F("hcontRatio", ";n_{Contrib}^{pil} / n_{Contrib}^{main}", 101, 0., 1.01);
  TH2F* hzdrc = new TH2F("hzdrc", ";n_{Contrib}^{A} / n_{Contrib}^{B};|z_{rec}^{A} - z_{rec}^{B}| (cm)", 101, 0., 1.01, 100, 0, 5);

  NA6PVertexerTracks* vertxr = new NA6PVertexerTracks();
  vertxr->setVerbosity(true);

  int nEv = trTree->GetEntries();
  printf("Number of events = %d\n", nEv);
  for (int jEv = 0; jEv < nEv; jEv++) {
    mcTree->GetEvent(jEv);
    trTree->GetEvent(jEv);
    int nPart = mcArr->size();
    int nTracks = trArr->size();
    double xVertGen = mcHead->getVX();
    double yVertGen = mcHead->getVY();
    double zVertGen = mcHead->getVZ();

    printf("Event %d particles = %d tracks = %d true z vertex position = %f\n", jEv, nPart, nTracks, zVertGen);

    vertxr->setBeamX(xVertGen);
    vertxr->setBeamY(yVertGen);
    vertxr->createTracksPool(*trArr);
    std::vector<NA6PVertex> vertices;
    vertxr->findVertices(vertices);
    int nVertices = vertices.size();
    hnvert->Fill(nVertices);
    int nv = vertices.size();
    std::vector<int> multSort(nv); // sort time indices in multiplicity
    std::iota(multSort.begin(), multSort.end(), 0);
    std::sort(multSort.begin(), multSort.end(), [vertices](int i, int j) {
      return vertices[i].getNContributors() > vertices[j].getNContributors();
    });
    auto& vt0 = vertices[0];
    for (int im = 0; im < nv; im++) { // loop from highest multiplicity to lowest one
      int it1 = multSort[im];
      auto& vtI = vertices[it1];
      double xRec = vtI.getX();
      double yRec = vtI.getY();
      double zRec = vtI.getZ();
      double sigx = vtI.getSigmaX();
      double sigy = vtI.getSigmaY();
      double sigz = vtI.getSigmaZ();
      if (im == 0) {
        hxrec->Fill(xRec);
        hyrec->Fill(yRec);
        hzrec->Fill(zRec);
        hxrecgen->Fill(xVertGen, xRec);
        hyrecgen->Fill(yVertGen, yRec);
        hzrecgen->Fill(zVertGen, zRec);
        hdx->Fill((xRec - xVertGen) * 1e4);
        hdy->Fill((yRec - yVertGen) * 1e4);
        hdz->Fill((zRec - zVertGen) * 1e4);
        hsigx->Fill(sigx * 1e4);
        hsigy->Fill(sigy * 1e4);
        hsigz->Fill(sigz * 1e4);
        hpullx->Fill((xRec - xVertGen) / sigx);
        hpully->Fill((yRec - yVertGen) / sigy);
        hpullz->Fill((zRec - zVertGen) / sigz);
        hncontr->Fill(vtI.getNContributors());
      }
      if (im > 0) {
        float zDiff0 = vtI.getZ() - vt0.getZ();
        float rCont0 = float(vtI.getNContributors()) / float(vt0.getNContributors());
        hzRecDiff->Fill(zDiff0);
        hcontRatio->Fill(rCont0);
      }
      for (int in = im + 1; in < nv; in++) { // loop from highest multiplicity to lowest one
        int it2 = multSort[in];
        auto& vtJ = vertices[it2];
        float zDiff = vtI.getZ() - vtJ.getZ();
        float rCont = float(vtJ.getNContributors()) / float(vtI.getNContributors());
        hzdrc->Fill(rCont, std::abs(zDiff));
      }
    }
  }

  TCanvas* cv = new TCanvas("cv", "", 1400, 800);
  cv->Divide(3, 2);
  cv->cd(1);
  hnvert->SetLineWidth(2);
  hnvert->Draw();
  cv->cd(2);
  hncontr->SetLineWidth(2);
  hncontr->Draw();
  cv->cd(4);
  hxrec->SetFillColor(kBlue - 9);
  hxrec->SetFillStyle(1001);
  hxrec->Draw();
  cv->cd(5);
  hyrec->SetFillColor(kBlue - 9);
  hyrec->SetFillStyle(1001);
  hyrec->Draw();
  cv->cd(6);
  hzrec->SetFillColor(kBlue - 9);
  hzrec->SetFillStyle(1001);
  hzrec->Draw();
  cv->SaveAs(Form("%s/vertexerSummary.png", dirSimu));

  TCanvas* cres = new TCanvas("cres", "", 1400, 800);
  cres->Divide(3, 2);
  cres->cd(1);
  hxrecgen->Draw("colz");
  cres->cd(2);
  hyrecgen->Draw("colz");
  cres->cd(3);
  hzrecgen->Draw("colz");
  cres->cd(4);
  hdx->SetFillColor(kOrange + 2);
  hdx->SetFillStyle(1001);
  hdx->Draw();
  cres->cd(5);
  hdy->SetFillColor(kOrange + 2);
  hdy->SetFillStyle(1001);
  hdy->Draw();
  cres->cd(6);
  hdz->SetFillColor(kOrange + 2);
  hdz->SetFillStyle(1001);
  hdz->Draw();
  cres->SaveAs(Form("%s/vertexerResiduals.png", dirSimu));

  TCanvas* cpull = new TCanvas("cpull", "", 1400, 800);
  cpull->Divide(3, 2);
  cpull->cd(1);
  hsigx->SetLineWidth(2);
  hsigx->Draw();
  cpull->cd(2);
  hsigy->SetLineWidth(2);
  hsigy->Draw();
  cpull->cd(3);
  hsigz->SetLineWidth(2);
  hsigz->Draw();
  cpull->cd(4);
  hpullx->SetFillColor(kTeal - 7);
  hpullx->SetFillStyle(1001);
  hpullx->Draw();
  cpull->cd(5);
  hpully->SetFillColor(kTeal - 7);
  hpully->SetFillStyle(1001);
  hpully->Draw();
  cpull->cd(6);
  hpullz->SetFillColor(kTeal - 7);
  hpullz->SetFillStyle(1001);
  hpullz->Draw();
  cpull->SaveAs(Form("%s/vertexerPulls.png", dirSimu));

  TCanvas* cpil = new TCanvas("cz", "", 1200, 800);
  cpil->Divide(2, 2);
  cpil->cd(1);
  hzRecDiff->Draw();
  cpil->cd(2);
  hcontRatio->Draw();
  cpil->cd(3);
  hzdrc->Draw("colz");
  cpil->SaveAs(Form("%s/vertexerPileup.png", dirSimu));
}
