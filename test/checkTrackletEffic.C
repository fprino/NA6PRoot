#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>
#include <TVector3.h>
#include <TParticle.h>
#include <TCanvas.h>
#include <TPaveStats.h>
#include <TLegend.h>
#include "NA6PVerTelCluster.h"
#include "NA6PVertex.h"
#include "NA6PTrackerCA.h"
#include "MagneticField.h"
#include "NA6PVerTelHit.h"
#include "NA6PMCTruthContainer.h"
#include "NA6PReconstruction.h"
#endif

void checkTrackletEffic(int firstLay = 0,
                        int secondLay = 1,
                        int firstEv = 0,
                        int lastEv = 9999,
                        const char* dirSimu = ".")
{

  auto magField = new MagneticField();
  if (!Propagator::loadField() || !Propagator::loadGeometry(Form("%s/geometry.root", dirSimu))) {
    return;
  }

  int nMomBins = 40;
  TH1F* hMomGen = new TH1F("hMomGen", ";p (GeV/c);counts", nMomBins, 0., 10.);
  TH1F* hEtaGen = new TH1F("hEtaGen", ";#eta;counts", 20, 1., 5.);
  TH1F* hMomReco = new TH1F("hMomReco", ";p (GeV/c);counts", nMomBins, 0., 10.);
  TH1F* hEtaReco = new TH1F("hEtaReco", ";#eta;counts", 20, 1., 5.);
  TH1F* hMomGoodReco = new TH1F("hMomGoodReco", ";p (GeV/c);counts", nMomBins, 0., 10.);
  TH1F* hEtaGoodReco = new TH1F("hEtaGoodReco", ";#eta;counts", 20, 1., 5.);
  TH1F* hInnerLay = new TH1F("hInnerLay", "inner;layer;counts", 15, -0.5, 14.5);
  TH1F* hOuterLay = new TH1F("hOuterLay", "outer;layer;counts", 15, -0.5, 14.5);

  std::unique_ptr<NA6PTrackerCA> tracker = std::make_unique<NA6PTrackerCA>();
  tracker->configureFromRecoParamVT();
  tracker->setNumberOfIterations(1);
  tracker->setIterationParams(0, 0.05, 0.1, 4., 0.4, 0.02, 0.002, 5000, 5000, 5000, 5);
  tracker->setIterationParams(1, 0.10, 0.2, 9., 0.6, 0.05, 0.0075, 1000, 1000, 1000, 3);
  tracker->printConfiguration();

  tracker->setVerbosity(true);
  TFile* fk = new TFile(Form("%s/MCKine.root", dirSimu));
  TTree* mcTree = (TTree*)fk->Get("mckine");
  std::vector<TParticle>* mcArr = nullptr;
  mcTree->SetBranchAddress("tracks", &mcArr);

  TFile* fh = new TFile(Form("%s/HitsVerTel.root", dirSimu));
  TTree* th = (TTree*)fh->Get("hitsVerTel");
  std::vector<NA6PVerTelHit> vtHits, *vtHitsPtr = &vtHits;
  th->SetBranchAddress("VerTel", &vtHitsPtr);

  TFile* fc = new TFile(Form("%s/ClustersVerTel.root", dirSimu));
  printf("Open cluster file: %s\n", fc->GetName());
  TTree* tc = (TTree*)fc->Get("clustersVerTel");
  std::vector<NA6PVerTelCluster> vtClus, *vtClusPtr = &vtClus;
  NA6PMCTruthContainer vtCluMCLabels, *vtCluMCLabelsPtr = &vtCluMCLabels;
  tc->SetBranchAddress("VerTel", &vtClusPtr);
  tc->SetBranchAddress("VerTelMCTruth", &vtCluMCLabelsPtr);
  int nEv = tc->GetEntries();
  if (lastEv > nEv || lastEv < 0)
    lastEv = nEv;
  if (firstEv < 0)
    firstEv = 0;

  NA6PVertex primVert;
  NA6PReconstruction rec("dummy");
  std::vector<NA6PMCComposedLabel> trkMCLabs;

  int nIterationsCA = tracker->getNIterations();

  std::vector<bool> isTrackable;
  std::vector<int> nFoundTracklets;
  for (int jEv = firstEv; jEv < lastEv; jEv++) {
    mcTree->GetEvent(jEv);
    th->GetEvent(jEv);
    int nPart = mcArr->size();
    isTrackable.resize(nPart);
    nFoundTracklets.resize(nPart);
    double zvert = 0;
    // get primary vertex position from the Kine Tree
    for (int jp = 0; jp < nPart; jp++) {
      auto curPart = mcArr->at(jp);
      if (curPart.IsPrimary()) {
        zvert = curPart.Vz();
        break;
      }
    }
    primVert.setXYZ(0., 0., zvert);
    uint nHits = vtHits.size();
    for (int jp = 0; jp < nPart; jp++) {
      nFoundTracklets[jp] = 0;
      isTrackable[jp] = false;
      auto curPart = mcArr->at(jp);
      int maskHits = 0;
      for (size_t jHit = 0; jHit < nHits; ++jHit) {
        const auto& hit = vtHits.at(jHit);
        int idPart = hit.getTrackID();
        if (idPart == jp) {
          int nLay = hit.getDetectorID() / 4;
          maskHits |= (1 << nLay);
        }
      }
      if ((maskHits & (1 << firstLay)) && (maskHits & (1 << secondLay))) {
        isTrackable[jp] = true;
        double pxPart = curPart.Px();
        double pyPart = curPart.Py();
        double pzPart = curPart.Pz();
        double momPart = curPart.P();
        double phiPart = curPart.Phi();
        double thetaPart = std::acos(pzPart / momPart);
        double etaPart = -std::log(std::tan(thetaPart / 2.));
        hMomGen->Fill(momPart);
        hEtaGen->Fill(etaPart);
      }
    }
    tc->GetEvent(jEv);
    for (size_t i = 0; i < vtClus.size(); ++i) {
      vtClus[i].setClusterIndex(static_cast<int>(i));
    }
    tracker->setClusterMCTruth(&vtCluMCLabels);
    std::vector<std::pair<NA6PVerTelCluster, NA6PVerTelCluster>> trklts = tracker->findTracklets(firstLay, secondLay, vtClus, &primVert);
    int nTrklets = trklts.size();
    printf("Event %d nHits = %d nClusters = %d nMCTruth entries = %d nTracklets = %d\n", jEv, nHits, (int)vtClus.size(), (int)vtCluMCLabels.getNElements(), nTrklets);
    for (int jT = 0; jT < nTrklets; jT++) {
      NA6PVerTelCluster clu1 = trklts[jT].first;
      NA6PVerTelCluster clu2 = trklts[jT].second;
      hInnerLay->Fill(clu1.getLayer());
      hOuterLay->Fill(clu2.getLayer());
      if (clu1.getLayer() != firstLay || clu2.getLayer() != secondLay) {
        printf("ERROR in tracklet layers: %d %d\n", clu1.getLayer(), clu2.getLayer());
        continue;
      }
      int cluID[2] = {clu1.getClusterIndex(), clu2.getClusterIndex()};
      std::vector<std::pair<NA6PMCComposedLabel, int>> countLabs;
      for (int jClu = 0; jClu < 2; ++jClu) {
        std::span labels = vtCluMCLabels.getLabels(cluID[jClu]);
        int nLabels = labels.size();
        for (int jLab = 0; jLab < nLabels; jLab++) {
          NA6PMCComposedLabel lbl = labels[jLab];
          bool found = false;
          for (auto& p : countLabs) {
            if (p.first == lbl) {
              ++p.second;
              found = true;
              break;
            }
          }
          if (!found) {
            countLabs.push_back({lbl, 1});
          }
        }
      }
      NA6PMCComposedLabel lblTrack;
      lblTrack.unset();
      int maxCountLabs = 0;
      for (const auto& p : countLabs) {
        if (p.second > maxCountLabs) {
          maxCountLabs = p.second;
          lblTrack = p.first;
        }
      }
      if (lblTrack.isSet() && maxCountLabs == 2) {
        int idPartTrack = lblTrack.getTrackID();
        if (!isTrackable[idPartTrack]) {
          printf("Mismatch!!!\n");
        }
        nFoundTracklets[idPartTrack]++;
        if (nFoundTracklets[idPartTrack] > 1) {
          printf("Particle %d had already another tracklet, tot tracklets = %d\n", idPartTrack, nFoundTracklets[idPartTrack]);
          continue;
        }
        auto curPart = mcArr->at(idPartTrack);
        double pxPart = curPart.Px();
        double pyPart = curPart.Py();
        double pzPart = curPart.Pz();
        double momPart = curPart.P();
        double phiPart = curPart.Phi();
        double thetaPart = std::acos(pzPart / momPart);
        double etaPart = -std::log(std::tan(thetaPart / 2.));
        hMomReco->Fill(momPart);
        hEtaReco->Fill(etaPart);
      }
    }
  }

  TCanvas* clay = new TCanvas("clay", "", 1200, 400);
  clay->Divide(2, 1);
  clay->cd(1);
  hInnerLay->Draw();
  clay->cd(2);
  hOuterLay->Draw();

  TCanvas* cef = new TCanvas("cef", "", 1400, 800);
  cef->Divide(2, 2);
  cef->cd(1);
  hMomGen->SetLineColor(kGray + 1);
  hMomGen->SetLineWidth(3);
  hMomGen->Draw();
  gPad->Update();
  TPaveStats* st1 = (TPaveStats*)hMomGen->GetListOfFunctions()->FindObject("stats");
  if (st1) {
    st1->SetY1NDC(0.72);
    st1->SetY2NDC(0.92);
    st1->SetTextColor(hMomGen->GetLineColor());
  }
  hMomReco->SetLineWidth(2);
  hMomReco->Draw("sames");
  gPad->Update();
  TPaveStats* st2 = (TPaveStats*)hMomReco->GetListOfFunctions()->FindObject("stats");
  if (st2) {
    st2->SetY1NDC(0.51);
    st2->SetY2NDC(0.71);
    st2->SetTextColor(hMomReco->GetLineColor());
  }
  gPad->Modified();
  TLegend* leg = new TLegend(0.25, 0.6, 0.89, 0.8);
  leg->AddEntry(hMomGen, Form("Generated particles with hits in layers %d and %d", firstLay, secondLay));
  leg->AddEntry(hMomReco, Form("Tracklets in layers %d and %d", firstLay, secondLay));
  leg->Draw();
  cef->cd(2);
  TH1F* hEffMom = (TH1F*)hMomReco->Clone("hEffMom");
  hEffMom->Divide(hMomReco, hMomGen, 1., 1., "B");
  hEffMom->GetYaxis()->SetTitle("Efficiency");
  hEffMom->SetStats(0);
  hEffMom->Draw();
  cef->cd(3);
  hEtaGen->SetLineColor(kGray + 1);
  hEtaGen->SetLineWidth(3);
  hEtaGen->Draw();
  gPad->Update();
  TPaveStats* st3 = (TPaveStats*)hEtaGen->GetListOfFunctions()->FindObject("stats");
  if (st3) {
    st3->SetY1NDC(0.72);
    st3->SetY2NDC(0.92);
    st3->SetTextColor(hEtaGen->GetLineColor());
  }
  hEtaReco->SetLineWidth(2);
  hEtaReco->Draw("sames");
  gPad->Update();
  TPaveStats* st4 = (TPaveStats*)hEtaReco->GetListOfFunctions()->FindObject("stats");
  if (st4) {
    st4->SetY1NDC(0.51);
    st4->SetY2NDC(0.71);
    st4->SetTextColor(hEtaReco->GetLineColor());
  }
  gPad->Modified();
  cef->cd(4);
  TH1F* hEffEta = (TH1F*)hEtaReco->Clone("hEffEta");
  hEffEta->Divide(hEtaReco, hEtaGen, 1., 1., "B");
  hEffEta->GetYaxis()->SetTitle("Efficiency");
  hEffEta->SetStats(0);
  hEffEta->Draw();
  cef->SaveAs(Form("TrackletEffic-VT%d%d.png", firstLay, secondLay));
}
