#include "NA6PMCTruthContainer.h"
#include "NA6PVerTelHit.h"
#include "NA6PVerTelDigit.h"
#include "NA6PVerTelClusterizer.h"
#include <TFile.h>
#include <TTree.h>
#include <TParticle.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TStyle.h>

void plotVerTelClusters(const char* dirSimu = ".")
{

  TFile* fk = new TFile(Form("%s/MCKine.root", dirSimu));
  TTree* mcTree = (TTree*)fk->Get("mckine");
  std::vector<TParticle>* mcArr = nullptr;
  mcTree->SetBranchAddress("tracks", &mcArr);

  TFile* fh = new TFile(Form("%s/HitsVerTel.root", dirSimu));
  TTree* th = (TTree*)fh->Get("hitsVerTel");
  std::vector<NA6PVerTelHit> vtHits, *vtHitsPtr = &vtHits;
  th->SetBranchAddress("VerTel", &vtHitsPtr);

  TFile* fd = TFile::Open(Form("%s/DigitsVerTel.root", dirSimu));
  TTree* td = (TTree*)fd->Get("digitsVerTel");
  std::vector<NA6PVerTelDigit> vtDigits, *vtDigitsPtr = &vtDigits;
  NA6PMCTruthContainer vtDigMCLabels, *vtDigMCLabelsPtr = &vtDigMCLabels;
  td->SetBranchAddress("VerTel", &vtDigitsPtr);
  td->SetBranchAddress("VerTelMCTruth", &vtDigMCLabelsPtr);

  TFile* fc = TFile::Open(Form("%s/ClustersVerTel.root", dirSimu));
  TTree* tc = (TTree*)fc->Get("clustersVerTel");
  std::vector<NA6PVerTelCluster> vtClusters, *vtClustersPtr = &vtClusters;
  NA6PMCTruthContainer vtCluMCLabels, *vtCluMCLabelsPtr = &vtCluMCLabels;
  tc->SetBranchAddress("VerTel", &vtClustersPtr);
  tc->SetBranchAddress("VerTelMCTruth", &vtCluMCLabelsPtr);

  int nEv = tc->GetEntries();

  TH1F* hMCLabelsPerCluster = new TH1F("hMCLabelsPerCluster", ";n MC labels per cluster; counts", 11, -0.5, 10.5);
  TH1F* hCluSiz = new TH1F("hCluSiz", ";cluster size; counts", 11, -0.5, 10.5);
  TH1F* hNCluPerHit = new TH1F("hNCluPerHit", ";n_{clusters}/n_{hits}", 100, 0.5, 1.);
  TH1F* hNCluPerDigit = new TH1F("hNCluPerDigit", ";n_{clusters}/n_{digits}", 100, 0.5, 1.);
  TH1F** hNCluPerEvent = new TH1F*[5];
  TH2F** hXYClu = new TH2F*[5];
  TH1F** hZClu = new TH1F*[5];
  TH2F** hXYHit = new TH2F*[5];
  TH1F** hZHit = new TH1F*[5];
  for (int jLay = 0; jLay < 5; jLay++) {
    hNCluPerEvent[jLay] = new TH1F(Form("hNCluPerEventLay%d", jLay), Form(";n_{Clusters, lay %d}/event;counts", jLay), 100, 0., 1000.);
    hXYClu[jLay] = new TH2F(Form("hXYCluLay%d", jLay), Form("clusters layer %d", jLay), 100, -15., 15., 100, -15., 15.);
    hZClu[jLay] = new TH1F(Form("hZCluLay%d", jLay), Form("clusters layer %d", jLay), 100, 0., 40.);
    hXYHit[jLay] = new TH2F(Form("hXYHitLay%d", jLay), Form("hits layer %d", jLay), 100, -15., 15., 100, -15., 15.);
    hZHit[jLay] = new TH1F(Form("hZHitLay%d", jLay), Form("hits layer %d", jLay), 100, 0., 40.);
  }
  TH1F* hDeltaXAll = new TH1F("hDeltaXAll", "All clusters;x_{glo}^{clu}-x_{glo}^{hit} (#mum);counts", 100, -300., 300.);
  TH1F* hDeltaYAll = new TH1F("hDeltaYAll", "All clusters;y_{glo}^{clu}-y_{glo}^{hit} (#mum);counts", 100, -300., 300.);
  TH1F* hDeltaXPrim = new TH1F("hDeltaXPrim", "Primary-particle clusters;x_{glo}^{clu}-x_{glo}^{hit} (#mum);counts", 100, -300., 300.);
  TH1F* hDeltaYPrim = new TH1F("hDeltaYPrim", "Primary-particle clusters;y_{glo}^{clu}-y_{glo}^{hit} (#mum);counts", 100, -300., 300.);
  TH1F* hDeltaXSingle = new TH1F("hDeltaXSingle", "Single-hit clusters;x_{glo}^{clu}-x_{glo}^{hit} (#mum);counts", 100, -300., 300.);
  TH1F* hDeltaYSingle = new TH1F("hDeltaYSingle", "Single-hit clusters;y_{glo}^{clu}-y_{glo}^{hit} (#mum);counts", 100, -300., 300.);
  TH2F** hDeltaXVsX = new TH2F*[20];
  TH2F** hDeltaXVsY = new TH2F*[20];
  TH2F** hDeltaYVsX = new TH2F*[20];
  TH2F** hDeltaYVsY = new TH2F*[20];
  for (int jMod = 0; jMod < 20; jMod++) {
    hDeltaXVsX[jMod] = new TH2F(Form("hDeltaXVsXMod%d", jMod), ";x (cm); x_{glo}^{clu}-x_{glo}^{hit} (#mum);counts", 30, -15., 15., 100, -300000, 300000);
    hDeltaXVsY[jMod] = new TH2F(Form("hDeltaXVsYMod%d", jMod), ";y (cm); x_{glo}^{clu}-x_{glo}^{hit} (#mum);counts", 30, -15., 15., 100, -300000, 300000);
    hDeltaYVsX[jMod] = new TH2F(Form("hDeltaYVsXMod%d", jMod), ";x (cm); y_{glo}^{clu}-y_{glo}^{hit} (#mum);counts", 30, -15., 15., 100, -300000, 300000);
    hDeltaYVsY[jMod] = new TH2F(Form("hDeltaYVsYMod%d", jMod), ";y (cm); y_{glo}^{clu}-x_{glo}^{hit} (#mum);counts", 30, -15., 15., 100, -300000, 300000);
  }
  for (int jEv = 0; jEv < nEv; jEv++) {
    mcTree->GetEvent(jEv);
    td->GetEvent(jEv);
    th->GetEvent(jEv);
    tc->GetEvent(jEv);
    int nPart = mcArr->size();
    int nHits = vtHits.size();
    int nDigits = vtDigits.size();
    int nDigMClabels = vtDigMCLabels.getNElements();
    int nClusters = vtClusters.size();
    int nCluMClabels = vtCluMCLabels.getNElements();
    printf("Event %d particles = %d hits = %d digits = %d digMClabels = %d clusters = %d cluMClabels = %d\n", jEv, nPart, nHits, nDigits, nDigMClabels, nClusters, nCluMClabels);
    hNCluPerHit->Fill((float)nClusters / (float)nHits);
    hNCluPerDigit->Fill((float)nClusters / (float)nDigits);
    for (const auto& hit : vtHits) {
      int nDet = hit.getDetectorID();
      int lay = nDet / 4;
      if (lay < 0 || lay >= 5) {
        printf("ERROR wrong layer %d\n", lay);
        continue;
      }
      hXYHit[lay]->Fill(hit.getX(), hit.getY());
      hZHit[lay]->Fill(hit.getZ());
    }
    int nCluPerLay[5] = {0, 0, 0, 0, 0};
    for (int jClu = 0; jClu < nClusters; ++jClu) {
      const auto& clu = vtClusters.at(jClu);
      int lay = clu.getLayer();
      int detID = clu.getDetectorID();
      if (lay < 0 || lay >= 5) {
        printf("ERROR wrong layer %d\n", lay);
        continue;
      }
      nCluPerLay[lay]++;
      hXYClu[lay]->Fill(clu.getX(), clu.getY());
      hZClu[lay]->Fill(clu.getZ());
      hCluSiz->Fill(clu.getClusterSize());
      std::span labels = vtCluMCLabels.getLabels(jClu);
      int nLabels = labels.size();
      hMCLabelsPerCluster->Fill(nLabels);
      if (nLabels == 1 && clu.getClusterSize() == 1) {
        NA6PMCComposedLabel lbl = labels[0];
        int trackID = lbl.getTrackID();
        auto curPart = mcArr->at(trackID);
        int nFoundHits = 0;
        for (const auto& hit : vtHits) {
          int idPart = hit.getTrackID();
          int nDet = hit.getDetectorID();
          if (idPart == trackID && detID == nDet) {
            float dx = (clu.getX() - hit.getX()) * 1e4;
            float dy = (clu.getY() - hit.getY()) * 1e4;
            hDeltaXAll->Fill(dx);
            hDeltaYAll->Fill(dy);
            if (curPart.IsPrimary()) {
              hDeltaXPrim->Fill(dx);
              hDeltaYPrim->Fill(dy);
            }
            nFoundHits++;
          }
        }
        if (nFoundHits == 1) {
          for (const auto& hit : vtHits) {
            int idPart = hit.getTrackID();
            int nDet = hit.getDetectorID();
            if (idPart == trackID && detID == nDet) {
              float dx = (clu.getX() - hit.getX()) * 1e4;
              float dy = (clu.getY() - hit.getY()) * 1e4;
              hDeltaXSingle->Fill(dx);
              hDeltaYSingle->Fill(dy);
              hDeltaXVsX[nDet]->Fill(clu.getX(), dx);
              hDeltaXVsY[nDet]->Fill(clu.getY(), dx);
              hDeltaYVsX[nDet]->Fill(clu.getX(), dy);
              hDeltaYVsY[nDet]->Fill(clu.getY(), dy);
            }
          }
        }
      }
    }
    for (int jLay = 0; jLay < 5; jLay++) {
      hNCluPerEvent[jLay]->Fill(nCluPerLay[jLay]);
    }
  }

  TCanvas* cclu = new TCanvas("cclu", "Clusters", 1400, 800);
  cclu->Divide(3, 2);
  cclu->cd(1);
  hNCluPerEvent[0]->Draw();
  cclu->cd(2);
  hNCluPerEvent[4]->Draw();
  cclu->cd(3);
  hCluSiz->Draw();
  cclu->cd(4);
  hNCluPerHit->Draw();
  cclu->cd(5);
  hNCluPerDigit->Draw();
  cclu->cd(6);
  gPad->SetLogy();
  hMCLabelsPerCluster->Draw();

  TCanvas* cxyh = new TCanvas("cxyh", "Hits XY", 1400, 800);
  cxyh->Divide(3, 2);
  for (int jLay = 0; jLay < 5; jLay++) {
    cxyh->cd(jLay + 1);
    hXYHit[jLay]->Draw("colz");
  }
  TCanvas* cxyc = new TCanvas("cxyc", "Clusters XY", 1400, 800);
  cxyc->Divide(3, 2);
  for (int jLay = 0; jLay < 5; jLay++) {
    cxyc->cd(jLay + 1);
    hXYClu[jLay]->Draw("colz");
  }

  TCanvas* czh = new TCanvas("czh", "Hits Z", 1400, 800);
  czh->Divide(3, 2);
  for (int jLay = 0; jLay < 5; jLay++) {
    czh->cd(jLay + 1);
    hZHit[jLay]->Draw();
  }
  TCanvas* czc = new TCanvas("czc", "Clusters Z", 1400, 800);
  czc->Divide(3, 2);
  for (int jLay = 0; jLay < 5; jLay++) {
    czc->cd(jLay + 1);
    hZClu[jLay]->Draw("colz");
  }
  gStyle->SetOptStat(111111);
  TCanvas* cdel = new TCanvas("cdel", "Residuals", 1400, 800);
  cdel->Divide(3, 2);
  cdel->cd(1);
  gPad->SetLogy();
  hDeltaXAll->Draw();
  cdel->cd(2);
  gPad->SetLogy();
  hDeltaXPrim->Draw();
  cdel->cd(3);
  gPad->SetLogy();
  hDeltaXSingle->Draw();
  cdel->cd(4);
  gPad->SetLogy();
  hDeltaYAll->Draw();
  cdel->cd(5);
  gPad->SetLogy();
  hDeltaYPrim->Draw();
  cdel->cd(6);
  gPad->SetLogy();
  hDeltaYSingle->Draw();
}
