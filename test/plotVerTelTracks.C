#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TParticle.h>
#include <TCanvas.h>
#include <TPaveStats.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include "NA6PTrack.h"
#include "MagneticField.h"
#include "Propagator.h"
#include "NA6PVerTelHit.h"
#include "NA6PMCComposedLabel.h"
#endif

void fillMeanAndRms(TH2F* h2d, TH1F* hMean, TH1F* hRms, TH1F* hSig)
{
  if (!h2d || !hMean || !hRms || !hSig) {
    printf("One of pointers is not set: h2d=%p hMean=%p hRms=%p hSig=%p\n", h2d, hMean, hRms, hSig);
    return;
  }
  int nMomBins = h2d->GetXaxis()->GetNbins();
  hMean->SetStats(0);
  hRms->SetStats(0);
  hSig->SetStats(0);
  for (int jp = 1; jp <= nMomBins; jp++) {
    TH1D* htmp1 = h2d->ProjectionY("htmp1", jp, jp);
    double mean = htmp1->GetMean();
    double emean = htmp1->GetMeanError();
    double rms = htmp1->GetRMS();
    double erms = htmp1->GetRMSError();
    hMean->SetBinContent(jp, mean);
    hMean->SetBinError(jp, emean);
    hRms->SetBinContent(jp, rms);
    hRms->SetBinError(jp, erms);
    double gw = rms;
    double egw = erms;
    if (htmp1->GetEntries() > 20. && htmp1->Fit("gaus", "Q", "", -3. * rms, 3. * rms)) {
      TF1* fg = (TF1*)htmp1->GetListOfFunctions()->FindObject("gaus");
      if (!fg)
        continue;
      gw = fg->GetParameter(2);
      egw = fg->GetParError(2);
    }
    hSig->SetBinContent(jp, gw);
    hSig->SetBinError(jp, egw);
    delete htmp1;
  }
}

void superposHistos(TH1* h1, TH1* h2, TH1* h3, TH1* h4 = nullptr, int col1 = 1, int col2 = kGreen + 1, int col3 = 2, int col4 = 4)
{
  h1->SetLineColor(col1);
  h1->SetLineWidth(2);
  if (h2->GetMaximum() > h1->GetMaximum())
    h1->SetMaximum(1.02 * h2->GetMaximum());
  if (h3->GetMaximum() > h1->GetMaximum())
    h1->SetMaximum(1.02 * h3->GetMaximum());
  h1->SetMinimum(0.5);
  h1->Draw();
  h2->SetLineColor(col2);
  h2->SetLineWidth(2);
  h2->Draw("sames");
  h3->SetLineColor(col3);
  h3->SetLineWidth(2);
  h3->Draw("sames");
  if (h4) {
    h4->SetLineColor(col4);
    h4->SetLineWidth(2);
    h4->Draw("sames");
  }
  gPad->Update();
  TPaveStats* st1 = (TPaveStats*)h1->GetListOfFunctions()->FindObject("stats");
  if (st1) {
    st1->SetY1NDC(0.72);
    st1->SetY2NDC(0.92);
    st1->SetTextColor(h1->GetLineColor());
  }
  TPaveStats* st2 = (TPaveStats*)h2->GetListOfFunctions()->FindObject("stats");
  if (st2) {
    st2->SetY1NDC(0.51);
    st2->SetY2NDC(0.71);
    st2->SetTextColor(h2->GetLineColor());
  }
  TPaveStats* st3 = (TPaveStats*)h3->GetListOfFunctions()->FindObject("stats");
  if (st3) {
    st3->SetY1NDC(0.30);
    st3->SetY2NDC(0.50);
    st3->SetTextColor(h3->GetLineColor());
  }
  if (h4) {
    TPaveStats* st4 = (TPaveStats*)h4->GetListOfFunctions()->FindObject("stats");
    if (st4) {
      st4->SetY1NDC(0.09);
      st4->SetY2NDC(0.29);
      st4->SetTextColor(h4->GetLineColor());
    }
  }
  gPad->Modified();
}

void plotVerTelTracks(const char* dirSimu = ".")
{
  auto magField = new MagneticField();
  magField->loadField();
  magField->setAsGlobalField();

  TFile* ft = new TFile(Form("%s/TracksVerTel.root", dirSimu));
  if (!ft)
    return;
  TTree* trTree = (TTree*)ft->Get("tracksVerTel");
  std::vector<NA6PTrack>* trArr = nullptr;
  std::vector<NA6PMCComposedLabel>* labArr = nullptr;
  trTree->SetBranchAddress("VerTel", &trArr);
  trTree->SetBranchAddress("VerTelMCTruth", &labArr);

  TFile* fk = new TFile(Form("%s/MCKine.root", dirSimu));
  TTree* mcTree = (TTree*)fk->Get("mckine");
  std::vector<TParticle>* mcArr = nullptr;
  mcTree->SetBranchAddress("tracks", &mcArr);

  TFile* fh = new TFile(Form("%s/HitsVerTel.root", dirSimu));
  TTree* th = (TTree*)fh->Get("hitsVerTel");
  std::vector<NA6PVerTelHit> vtHits, *vtHitsPtr = &vtHits;
  th->SetBranchAddress("VerTel", &vtHitsPtr);

  int nMomBins = 20;
  double maxP = 10.;
  TH1F* hMomGen = new TH1F("hMomGen", ";p_{gen} (GeV/c);counts", nMomBins, 0., maxP);
  TH1F* hEtaGen = new TH1F("hEtaGen", ";#eta_{gen};counts", 20, 1., 5.);
  TH1F* hMomTrackable = new TH1F("hMomTrackable", ";p_{gen} (GeV/c);counts", nMomBins, 0., maxP);
  TH1F* hEtaTrackable = new TH1F("hEtaTrackable", ";#eta_{gen};counts", 20, 1., 5.);
  TH1F* hMomTracked5clu = new TH1F("hMomTracked5clu", ";p_{gen} (GeV/c);counts", nMomBins, 0., maxP);
  TH1F* hEtaTracked5clu = new TH1F("hEtaTracked5clu", ";#eta_{gen};counts", 20, 1., 5.);
  TH1F* hMomAllReco5clu = new TH1F("hMomAllReco5clu", ";p (GeV/c); counts", nMomBins, 0., maxP);
  TH1F* hMomGoodReco5clu = new TH1F("hMomGoodReco5clu", ";p (GeV/c); counts", nMomBins, 0., maxP);
  TH1F* hMomFakeReco5clu = new TH1F("hMomFakeReco5clu", ";p (GeV/c); counts", nMomBins, 0., maxP);
  TH1F* hEtaAllReco5clu = new TH1F("hEtaAllReco5clu", ";#eta; counts", 20, 1., 5.);
  TH1F* hEtaGoodReco5clu = new TH1F("hEtaGoodReco5clu", ";#eta; counts", 20, 1., 5.);
  TH1F* hEtaFakeReco5clu = new TH1F("hEtaFakeReco5clu", ";#eta; counts", 20, 1., 5.);
  TH1F* hMomTrackedge4clu = new TH1F("hMomTrackedge4clu", ";p_{gen} (GeV/c);counts", nMomBins, 0., maxP);
  TH1F* hEtaTrackedge4clu = new TH1F("hEtaTrackedge4clu", ";#eta_{gen};counts", 20, 1., 5.);
  TH1F* hMomAllRecoge4clu = new TH1F("hMomAllRecoge4clu", ";p (GeV/c); counts", nMomBins, 0., maxP);
  TH1F* hMomGoodRecoge4clu = new TH1F("hMomGoodRecoge4clu", ";p (GeV/c); counts", nMomBins, 0., maxP);
  TH1F* hMomFakeRecoge4clu = new TH1F("hMomFakeRecoge4clu", ";p (GeV/c); counts", nMomBins, 0., maxP);
  TH1F* hEtaAllRecoge4clu = new TH1F("hEtaAllRecoge4clu", ";#eta; counts", 20, 1., 5.);
  TH1F* hEtaGoodRecoge4clu = new TH1F("hEtaGoodRecoge4clu", ";#eta; counts", 20, 1., 5.);
  TH1F* hEtaFakeRecoge4clu = new TH1F("hEtaFakeRecoge4clu", ";#eta; counts", 20, 1., 5.);
  TH1F* hNclu = new TH1F("hNclu", ";n_{ITSclus};counts", 7, -0.5, 6.5);
  double impMax = 500.;  // in microns
  double deltaMax = 0.5; // in GeV/c
  TH2F* hImpParXVsP = new TH2F("hImpParXVsP", "5-cluster tracks;p (GeV/c);Track Imp. Par. X (#mum)};counts", nMomBins, 0., maxP, 100, -impMax, impMax);
  TH2F* hImpParYVsP = new TH2F("hImpParYVsP", "5-cluster tracks;p (GeV/c);Track Imp. Par. Y (#mum)};counts", nMomBins, 0., maxP, 100, -impMax, impMax);
  TH2F* hDeltaPxVsP = new TH2F("hDeltaPxVsP", "5-cluster tracks;p (GeV/c);p_{x}^{rec}-p_{x}^{gen} (GeV/c);counts", nMomBins, 0., maxP, 500, -deltaMax, deltaMax);
  TH2F* hDeltaPyVsP = new TH2F("hDeltaPyVsP", "5-cluster tracks;p (GeV/c);p_{y}^{rec}-p_{y}^{gen} (GeV/c);counts", nMomBins, 0., maxP, 500, -deltaMax, deltaMax);
  TH2F* hDeltaPzVsP = new TH2F("hDeltaPzVsP", "5-cluster tracks;p (GeV/c);p_{z}^{rec}-p_{z}^{gen} (GeV/c);counts", nMomBins, 0., maxP, 500, -deltaMax, deltaMax);
  TH2F* hRelDeltaPxVsP = new TH2F("hRelDeltaPxVsP", "5-cluster tracks;p (GeV/c);(p_{x}^{rec}-p_{x}^{gen})/p_{x}^{gen};counts", nMomBins, 0., maxP, 100, -deltaMax, deltaMax);
  TH2F* hRelDeltaPyVsP = new TH2F("hRelDeltaPyVsP", "5-cluster tracks;p (GeV/c);(p_{y}^{rec}-p_{y}^{gen})/p_{y}^{gen};counts", nMomBins, 0., maxP, 100, -deltaMax, deltaMax);
  TH2F* hRelDeltaPzVsP = new TH2F("hRelDeltaPzVsP", "5-cluster tracks;p (GeV/c);(p_{z}^{rec}-p_{z}^{gen})/p_{z}^{gen};counts", nMomBins, 0., maxP, 100, -deltaMax, deltaMax);

  int nEv = trTree->GetEntries();
  printf("Number of events = %d\n", nEv);
  for (int jEv = 0; jEv < nEv; jEv++) {
    mcTree->GetEvent(jEv);
    trTree->GetEvent(jEv);
    th->GetEvent(jEv);
    int nPart = mcArr->size();
    int nTracks = trArr->size();
    int nHits = vtHits.size();
    printf("Event %d particles = %d tracks = %d\n", jEv, nPart, nTracks);
    double xvert = 0;
    double yvert = 0;
    double zvert = 0;

    for (int jp = 0; jp < nPart; jp++) {
      auto curPart = mcArr->at(jp);
      if (curPart.IsPrimary()) {
        xvert = curPart.Vx();
        yvert = curPart.Vy();
        zvert = curPart.Vz();
      }
      int jDau = curPart.GetFirstDaughter();
      double zDecay = 999.;
      double zOrig = curPart.Vz();
      if (jDau >= 0) {
        TParticle dauPart = mcArr->at(jDau);
        zDecay = dauPart.Vz();
      }
      double pxPart = curPart.Px();
      double pyPart = curPart.Py();
      double pzPart = curPart.Pz();
      double momPart = curPart.P();
      double phiPart = curPart.Phi();
      double thetaPart = std::acos(pzPart / momPart);
      double etaPart = -std::log(std::tan(thetaPart / 2.));
      if (zOrig < 7 && zDecay > 40 && etaPart > 1 && etaPart < 5) {
        hMomGen->Fill(momPart);
        hEtaGen->Fill(etaPart);
      }
      int maskHits = 0;
      for (int jHit = 0; jHit < nHits; ++jHit) {
        const auto& hit = vtHits.at(jHit);
        int idPart = hit.getTrackID();
        if (idPart == jp) {
          int nLay = hit.getDetectorID() / 4;
          maskHits |= (1 << nLay);
        }
      }
      if (maskHits == 31 || maskHits == 30 || maskHits == 15) {
        hMomTrackable->Fill(momPart);
        hEtaTrackable->Fill(etaPart);
      }
    }
    for (int jTr = 0; jTr < nTracks; ++jTr) {
      NA6PTrack tr = trArr->at(jTr);
      Propagator::PropOpt propOpt;
      propOpt.matCorr = Propagator::MatCorrType::USEMatCorrNONE;
      Propagator::Instance()->propagateToZ(tr, zvert, propOpt);
      NA6PMCComposedLabel mcCompLabel = labArr->at(jTr);
      int nClusters = tr.getNHits();
      hNclu->Fill(nClusters);
      if (nClusters < 4)
        continue;
      auto pxyzReco = tr.getPXYZ<double>();
      double pxReco = pxyzReco[0];
      double pyReco = pxyzReco[1];
      double pzReco = pxyzReco[2];
      double ptReco = std::sqrt(pxyzReco[0] * pxyzReco[0] + pxyzReco[1] * pxyzReco[1]);
      double momReco = tr.getP();
      double thetaReco = std::acos(pxyzReco[2] / momReco);
      double etaReco = -std::log(std::tan(thetaReco / 2.));
      double impparX = tr.getX() - xvert;
      double impparY = tr.getY() - yvert;
      hMomAllRecoge4clu->Fill(momReco);
      hEtaAllRecoge4clu->Fill(etaReco);
      if (nClusters == 5) {
        hImpParXVsP->Fill(momReco, impparX * 1e4);
        hImpParYVsP->Fill(momReco, impparY * 1e4);
        hMomAllReco5clu->Fill(momReco);
        hEtaAllReco5clu->Fill(etaReco);
      }
      int mcLabel = mcCompLabel.getTrackID();
      if (!mcCompLabel.isFake()) {
        hMomGoodRecoge4clu->Fill(momReco);
        hEtaGoodRecoge4clu->Fill(etaReco);
        if (nClusters == 5) {
          hMomGoodReco5clu->Fill(momReco);
          hEtaGoodReco5clu->Fill(etaReco);
        }
      } else {
        hMomFakeRecoge4clu->Fill(momReco);
        hEtaFakeRecoge4clu->Fill(etaReco);
        if (nClusters == 5) {
          hMomFakeReco5clu->Fill(momReco);
          hEtaFakeReco5clu->Fill(etaReco);
        }
      }

      TParticle part = mcArr->at(mcLabel);
      double pxPart = part.Px();
      double pyPart = part.Py();
      double pzPart = part.Pz();
      double momPart = part.P();
      double phiPart = part.Phi();
      double thetaPart = std::acos(pzPart / momPart);
      double etaPart = -std::log(std::tan(thetaPart / 2.));
      hMomTrackedge4clu->Fill(momPart);
      hEtaTrackedge4clu->Fill(etaPart);
      if (nClusters == 5) {
        hMomTracked5clu->Fill(momPart);
        hEtaTracked5clu->Fill(etaPart);
        hDeltaPxVsP->Fill(momReco, pxReco - pxPart);
        hDeltaPyVsP->Fill(momReco, pyReco - pyPart);
        hDeltaPzVsP->Fill(momReco, pzReco - pzPart);
        hRelDeltaPxVsP->Fill(momReco, (pxReco - pxPart) / pxPart);
        hRelDeltaPyVsP->Fill(momReco, (pyReco - pyPart) / pyPart);
        hRelDeltaPzVsP->Fill(momReco, (pzReco - pzPart) / pzPart);
      }
    }
  }

  TH1F* hPurityMom5clu = (TH1F*)hMomGoodReco5clu->Clone("hPurityMom5clu");
  hPurityMom5clu->GetYaxis()->SetTitle("purity");
  TH1F* hPurityMomge4clu = (TH1F*)hMomGoodRecoge4clu->Clone("hPurityMomge4clu");
  hPurityMomge4clu->GetYaxis()->SetTitle("purity");
  for (int iBin = 1; iBin <= hMomGoodReco5clu->GetNbinsX(); iBin++) {
    double cg5 = hMomGoodReco5clu->GetBinContent(iBin);
    double ct5 = hMomAllReco5clu->GetBinContent(iBin);
    if (ct5 == 0) {
      hPurityMom5clu->SetBinContent(iBin, 0.);
      hPurityMom5clu->SetBinError(iBin, 0.);
    } else {
      double p = cg5 / ct5;
      double ep = std::sqrt(p * (1 - p) / ct5);
      hPurityMom5clu->SetBinContent(iBin, p);
      hPurityMom5clu->SetBinError(iBin, ep);
    }
    double cg4 = hMomGoodRecoge4clu->GetBinContent(iBin);
    double ct4 = hMomAllRecoge4clu->GetBinContent(iBin);
    if (ct4 == 0) {
      hPurityMomge4clu->SetBinContent(iBin, 0.);
      hPurityMomge4clu->SetBinError(iBin, 0.);
    } else {
      double p = cg4 / ct4;
      double ep = std::sqrt(p * (1 - p) / ct4);
      hPurityMomge4clu->SetBinContent(iBin, p);
      hPurityMomge4clu->SetBinError(iBin, ep);
    }
  }
  TH1F* hPurityEta5clu = (TH1F*)hEtaGoodReco5clu->Clone("hPurityEta5clu");
  hPurityEta5clu->GetYaxis()->SetTitle("purity");
  TH1F* hPurityEtage4clu = (TH1F*)hEtaGoodRecoge4clu->Clone("hPurityEtage4clu");
  hPurityEtage4clu->GetYaxis()->SetTitle("purity");
  for (int iBin = 1; iBin <= hEtaGoodReco5clu->GetNbinsX(); iBin++) {
    double cg5 = hEtaGoodReco5clu->GetBinContent(iBin);
    double ct5 = hEtaAllReco5clu->GetBinContent(iBin);
    if (ct5 == 0) {
      hPurityEta5clu->SetBinContent(iBin, 0.);
      hPurityEta5clu->SetBinError(iBin, 0.);
    } else {
      double p = cg5 / ct5;
      double ep = std::sqrt(p * (1 - p) / ct5);
      hPurityEta5clu->SetBinContent(iBin, p);
      hPurityEta5clu->SetBinError(iBin, ep);
    }
    double cg4 = hEtaGoodRecoge4clu->GetBinContent(iBin);
    double ct4 = hEtaAllRecoge4clu->GetBinContent(iBin);
    if (ct4 == 0) {
      hPurityEtage4clu->SetBinContent(iBin, 0.);
      hPurityEtage4clu->SetBinError(iBin, 0.);
    } else {
      double p = cg4 / ct4;
      double ep = std::sqrt(p * (1 - p) / ct4);
      hPurityEtage4clu->SetBinContent(iBin, p);
      hPurityEtage4clu->SetBinError(iBin, ep);
    }
  }

  TH1F* hImpParXMean = new TH1F("hImpParXMean", ";p (GeV/c);<Imp Par X> (#mum)", nMomBins, 0., maxP);
  TH1F* hImpParYMean = new TH1F("hImpParYMean", ";p (GeV/c);<Imp Par Y> (#mum)", nMomBins, 0., maxP);
  TH1F* hImpParXRms = new TH1F("hImpParXRms", ";p (GeV/c);rms (Imp Par X) (#mum)", nMomBins, 0., maxP);
  TH1F* hImpParYRms = new TH1F("hImpParYRms", ";p (GeV/c);rms (Imp Par Y) (#mum)", nMomBins, 0., maxP);
  TH1F* hImpParXSig = new TH1F("hImpParXSig", ";p (GeV/c);#sigma(Imp Par X) (#mum)", nMomBins, 0., maxP);
  TH1F* hImpParYSig = new TH1F("hImpParYSig", ";p (GeV/c);#sigma(Imp Par Y) (#mum)", nMomBins, 0., maxP);
  fillMeanAndRms(hImpParXVsP, hImpParXMean, hImpParXRms, hImpParXSig);
  fillMeanAndRms(hImpParYVsP, hImpParYMean, hImpParYRms, hImpParYSig);

  TH1F* hDeltaPxMean = new TH1F("hDeltaPxMean", ";p (GeV/c);<p_{x}^{rec}-p_{x}^{gen}> (GeV/c)", nMomBins, 0., maxP);
  TH1F* hDeltaPxRms = new TH1F("hDeltaPxRms", ";p (GeV/c);rms (p_{x}^{rec}-p_{x}^{gen}) (GeV/c)", nMomBins, 0., maxP);
  TH1F* hDeltaPxSig = new TH1F("hDeltaPxSig", ";p (GeV/c);#sigma(p_{x}^{rec}-p_{x}^{gen}) (GeV/c)", nMomBins, 0., maxP);
  TH1F* hRelDeltaPxMean = new TH1F("hRelDeltaPxMean", ";p (GeV/c);<(p_{x}^{rec}-p_{x}^{gen})/p_{x}^{gen}>", nMomBins, 0., maxP);
  TH1F* hRelDeltaPxRms = new TH1F("hRelDeltaPxRms", ";p (GeV/c);rms (p_{x}^{rec}-p_{x}^{gen})/p_{x}^{gen}", nMomBins, 0., maxP);
  TH1F* hRelDeltaPxSig = new TH1F("hRelDeltaPxSig", ";p (GeV/c);#sigma(p_{x}^{rec}-p_{x}^{gen})/p_{x}^{gen} (GeV/c)", nMomBins, 0., maxP);
  fillMeanAndRms(hDeltaPxVsP, hDeltaPxMean, hDeltaPxRms, hDeltaPxSig);
  fillMeanAndRms(hRelDeltaPxVsP, hRelDeltaPxMean, hRelDeltaPxRms, hRelDeltaPxSig);

  TH1F* hDeltaPyMean = new TH1F("hDeltaPyMean", ";p (GeV/c);<p_{y}^{rec}-p_{y}^{gen}> (GeV/c)", nMomBins, 0., maxP);
  TH1F* hDeltaPyRms = new TH1F("hDeltaPyRms", ";p (GeV/c);rms (p_{y}^{rec}-p_{y}^{gen}) (GeV/c)", nMomBins, 0., maxP);
  TH1F* hDeltaPySig = new TH1F("hDeltaPySig", ";p (GeV/c);#sigma(p_{y}^{rec}-p_{y}^{gen}) (GeV/c)", nMomBins, 0., maxP);
  TH1F* hRelDeltaPyMean = new TH1F("hRelDeltaPyMean", ";p (GeV/c);<(p_{y}^{rec}-p_{y}^{gen})/p_{y}^{gen}>", nMomBins, 0., maxP);
  TH1F* hRelDeltaPyRms = new TH1F("hRelDeltaPyRms", ";p (GeV/c);rms (p_{y}^{rec}-p_{y}^{gen})/p_{y}^{gen}", nMomBins, 0., maxP);
  TH1F* hRelDeltaPySig = new TH1F("hRelDeltaPySig", ";p (GeV/c);#sigma(p_{y}^{rec}-p_{y}^{gen})/p_{y}^{gen} (GeV/c)", nMomBins, 0., maxP);
  fillMeanAndRms(hDeltaPyVsP, hDeltaPyMean, hDeltaPyRms, hDeltaPySig);
  fillMeanAndRms(hRelDeltaPyVsP, hRelDeltaPyMean, hRelDeltaPyRms, hRelDeltaPySig);

  TH1F* hDeltaPzMean = new TH1F("hDeltaPzMean", ";p (GeV/c);<p_{z}^{rec}-p_{z}^{gen}> (GeV/c)", nMomBins, 0., maxP);
  TH1F* hDeltaPzRms = new TH1F("hDeltaPzRms", ";p (GeV/c);rms (p_{z}^{rec}-p_{z}^{gen}) (GeV/c)", nMomBins, 0., maxP);
  TH1F* hDeltaPzSig = new TH1F("hDeltaPzSig", ";p (GeV/c);#sigma(p_{z}^{rec}-p_{z}^{gen}) (GeV/c)", nMomBins, 0., maxP);
  TH1F* hRelDeltaPzMean = new TH1F("hRelDeltaPzMean", ";p (GeV/c);<(p_{z}^{rec}-p_{z}^{gen})/p_{z}^{gen}>", nMomBins, 0., maxP);
  TH1F* hRelDeltaPzRms = new TH1F("hRelDeltaPzRms", ";p (GeV/c);rms (p_{z}^{rec}-p_{z}^{gen})/p_{z}^{gen}", nMomBins, 0., maxP);
  TH1F* hRelDeltaPzSig = new TH1F("hRelDeltaPzSig", ";p (GeV/c);#sigma(p_{z}^{rec}-p_{z}^{gen})/p_{z}^{gen} (GeV/c)", nMomBins, 0., maxP);
  fillMeanAndRms(hDeltaPzVsP, hDeltaPzMean, hDeltaPzRms, hDeltaPzSig);
  fillMeanAndRms(hRelDeltaPzVsP, hRelDeltaPzMean, hRelDeltaPzRms, hRelDeltaPzSig);

  TCanvas* cef = new TCanvas("cef", "Efficiency", 1400, 800);
  cef->Divide(2, 2);
  cef->cd(1);
  superposHistos(hMomTrackable, hMomTracked5clu, hMomTrackedge4clu, 0x0, kMagenta + 1, 1, kBlue + 1);
  cef->cd(2);
  TH1F* hEffMom5clu = (TH1F*)hMomTracked5clu->Clone("hEffMom5clu");
  hEffMom5clu->Divide(hMomTracked5clu, hMomTrackable, 1., 1., "B");
  hEffMom5clu->GetYaxis()->SetTitle("Efficiency");
  hEffMom5clu->SetStats(0);
  hEffMom5clu->SetMaximum(1.05);
  hEffMom5clu->SetMinimum(0.);
  hEffMom5clu->Draw();
  TH1F* hEffMomge4clu = (TH1F*)hMomTrackedge4clu->Clone("hEffMomge4clu");
  hEffMomge4clu->Divide(hMomTrackedge4clu, hMomTrackable, 1., 1., "B");
  hEffMomge4clu->GetYaxis()->SetTitle("Efficiency");
  hEffMomge4clu->SetStats(0);
  hEffMomge4clu->Draw("same");
  TLegend* lege = new TLegend(0.3, 0.15, 0.85, 0.3);
  lege->SetMargin(0.1);
  lege->AddEntry(hEffMom5clu, "5 cluster tracks", "L")->SetTextColor(hEffMom5clu->GetLineColor());
  lege->AddEntry(hEffMomge4clu, ">= 4 cluster tracks", "L")->SetTextColor(hEffMomge4clu->GetLineColor());
  lege->Draw();
  cef->cd(3);
  superposHistos(hEtaTrackable, hEtaTracked5clu, hEtaTrackedge4clu, 0x0, kMagenta + 1, 1, kBlue + 1);
  cef->cd(4);
  TH1F* hEffEta5clu = (TH1F*)hEtaTracked5clu->Clone("hEffEta5clu");
  hEffEta5clu->Divide(hEtaTracked5clu, hEtaTrackable, 1., 1., "B");
  hEffEta5clu->GetYaxis()->SetTitle("Efficiency");
  hEffEta5clu->SetStats(0);
  hEffEta5clu->SetMaximum(1.05);
  hEffEta5clu->SetMinimum(0.);
  hEffEta5clu->Draw();
  TH1F* hEffEtage4clu = (TH1F*)hEtaTrackedge4clu->Clone("hEffEtage4clu");
  hEffEtage4clu->Divide(hEtaTrackedge4clu, hEtaTrackable, 1., 1., "B");
  hEffEtage4clu->GetYaxis()->SetTitle("Efficiency");
  hEffEtage4clu->SetStats(0);
  hEffEtage4clu->Draw("same");
  cef->SaveAs("TrackingEfficVT-4and5clu.png");

  TCanvas* cpu = new TCanvas("cpu", "Purity", 1400, 800);
  cpu->Divide(3, 2);
  cpu->cd(1);
  gPad->SetLogy();
  superposHistos(hMomAllReco5clu, hMomGoodReco5clu, hMomFakeReco5clu);
  cpu->cd(2);
  gPad->SetLogy();
  superposHistos(hMomAllRecoge4clu, hMomGoodRecoge4clu, hMomFakeRecoge4clu);
  cpu->cd(3);
  hPurityMom5clu->GetYaxis()->SetTitle("Purity");
  hPurityMom5clu->SetMinimum(0.8);
  hPurityMom5clu->SetStats(0);
  hPurityMom5clu->SetLineWidth(2);
  hPurityMom5clu->Draw();
  hPurityMomge4clu->GetYaxis()->SetTitle("Purity");
  hPurityMomge4clu->SetMinimum(0.8);
  hPurityMomge4clu->SetStats(0);
  hPurityMomge4clu->SetLineWidth(2);
  hPurityMomge4clu->SetLineColor(kBlue + 1);
  hPurityMomge4clu->Draw("same");
  lege->Draw();
  cpu->cd(4);
  superposHistos(hEtaAllReco5clu, hEtaGoodReco5clu, hEtaFakeReco5clu);
  cpu->cd(5);
  superposHistos(hEtaAllRecoge4clu, hEtaGoodRecoge4clu, hEtaFakeRecoge4clu);
  cpu->cd(6);
  hPurityEta5clu->GetYaxis()->SetTitle("Purity");
  hPurityEta5clu->SetMinimum(0.8);
  hPurityEta5clu->SetStats(0);
  hPurityEta5clu->SetLineWidth(2);
  hPurityEta5clu->Draw();
  hPurityEtage4clu->GetYaxis()->SetTitle("Purity");
  hPurityEtage4clu->SetMinimum(0.8);
  hPurityEtage4clu->SetStats(0);
  hPurityEtage4clu->SetLineWidth(2);
  hPurityEtage4clu->SetLineColor(kBlue + 1);
  hPurityEtage4clu->Draw("same");
  cpu->SaveAs("TrackingPurityVT-4and5clu.png");

  TCanvas* cip = new TCanvas("cip", "Impact Parameter", 1200, 800);
  cip->Divide(2, 2);
  cip->cd(1);
  gPad->SetTickx();
  gPad->SetTicky();
  hImpParXMean->SetMarkerStyle(25);
  hImpParXMean->SetMarkerColor(1);
  hImpParXMean->SetLineColor(1);
  hImpParXMean->SetMinimum(-20);
  hImpParXMean->SetMaximum(20);
  hImpParXMean->Draw("P");
  cip->cd(2);
  gPad->SetTickx();
  gPad->SetTicky();
  hImpParXSig->SetMarkerStyle(25);
  hImpParXSig->SetMarkerColor(1);
  hImpParXSig->SetLineColor(1);
  hImpParXSig->SetMinimum(0);
  hImpParXSig->SetMaximum(120);
  hImpParXSig->Draw("P");
  cip->cd(3);
  gPad->SetTickx();
  gPad->SetTicky();
  hImpParYMean->SetMarkerStyle(25);
  hImpParYMean->SetMarkerColor(1);
  hImpParYMean->SetLineColor(1);
  hImpParYMean->SetMinimum(-20);
  hImpParYMean->SetMaximum(20);
  hImpParYMean->Draw("P");
  cip->cd(4);
  gPad->SetTickx();
  gPad->SetTicky();
  hImpParYSig->SetMarkerStyle(25);
  hImpParYSig->SetMarkerColor(1);
  hImpParYSig->SetLineColor(1);
  hImpParYSig->SetMinimum(0);
  hImpParYSig->SetMaximum(120);
  hImpParYSig->Draw("P");

  TCanvas* cmom = new TCanvas("cmom", "Momentum Resolution", 1400, 900);
  cmom->Divide(3, 3);
  cmom->cd(1);
  gPad->SetTickx();
  gPad->SetTicky();
  hDeltaPxMean->SetMarkerStyle(25);
  hDeltaPxMean->SetMarkerColor(1);
  hDeltaPxMean->SetLineColor(1);
  hDeltaPxMean->SetMinimum(-0.05);
  hDeltaPxMean->SetMaximum(0.05);
  hDeltaPxMean->Draw("P");
  cmom->cd(2);
  gPad->SetTickx();
  gPad->SetTicky();
  hDeltaPxSig->SetMarkerStyle(25);
  hDeltaPxSig->SetMarkerColor(1);
  hDeltaPxSig->SetLineColor(1);
  hDeltaPxSig->SetMinimum(0);
  hDeltaPxSig->SetMaximum(0.05);
  hDeltaPxSig->Draw("P");
  cmom->cd(3);
  gPad->SetTickx();
  gPad->SetTicky();
  hRelDeltaPxSig->SetMarkerStyle(25);
  hRelDeltaPxSig->SetMarkerColor(1);
  hRelDeltaPxSig->SetLineColor(1);
  hRelDeltaPxSig->SetMinimum(0);
  hRelDeltaPxSig->SetMaximum(0.1);
  hRelDeltaPxSig->Draw("P");
  cmom->cd(4);
  gPad->SetTickx();
  gPad->SetTicky();
  hDeltaPyMean->SetMarkerStyle(25);
  hDeltaPyMean->SetMarkerColor(1);
  hDeltaPyMean->SetLineColor(1);
  hDeltaPyMean->SetMinimum(-0.05);
  hDeltaPyMean->SetMaximum(0.05);
  hDeltaPyMean->Draw("P");
  cmom->cd(5);
  gPad->SetTickx();
  gPad->SetTicky();
  hDeltaPySig->SetMarkerStyle(25);
  hDeltaPySig->SetMarkerColor(1);
  hDeltaPySig->SetLineColor(1);
  hDeltaPySig->SetMinimum(0);
  hDeltaPySig->SetMaximum(0.05);
  hDeltaPySig->Draw("P");
  cmom->cd(6);
  gPad->SetTickx();
  gPad->SetTicky();
  hRelDeltaPySig->SetMarkerStyle(25);
  hRelDeltaPySig->SetMarkerColor(1);
  hRelDeltaPySig->SetLineColor(1);
  hRelDeltaPySig->SetMinimum(0);
  hRelDeltaPySig->SetMaximum(0.1);
  hRelDeltaPySig->Draw("P");
  cmom->cd(7);
  gPad->SetTickx();
  gPad->SetTicky();
  hDeltaPzMean->SetMarkerStyle(25);
  hDeltaPzMean->SetMarkerColor(1);
  hDeltaPzMean->SetLineColor(1);
  hDeltaPzMean->SetMinimum(-0.05);
  hDeltaPzMean->SetMaximum(0.05);
  hDeltaPzMean->Draw("P");
  cmom->cd(8);
  gPad->SetTickx();
  gPad->SetTicky();
  hDeltaPzSig->SetMarkerStyle(25);
  hDeltaPzSig->SetMarkerColor(1);
  hDeltaPzSig->SetLineColor(1);
  hDeltaPzSig->SetMinimum(0);
  hDeltaPzSig->SetMaximum(0.2);
  hDeltaPzSig->Draw("P");
  cmom->cd(9);
  hRelDeltaPzSig->SetMarkerStyle(25);
  hRelDeltaPzSig->SetMarkerColor(1);
  hRelDeltaPzSig->SetLineColor(1);
  hRelDeltaPzSig->SetMinimum(0);
  hRelDeltaPzSig->SetMaximum(0.05);
  hRelDeltaPzSig->Draw("P");

  TCanvas* c2d = new TCanvas("c2d", "", 1500, 500);
  c2d->Divide(3, 1);
  c2d->cd(1);
  hDeltaPxVsP->Draw("colz");
  c2d->cd(2);
  hDeltaPyVsP->Draw("colz");
  c2d->cd(3);
  hDeltaPzVsP->Draw("colz");

  TFile* outRoot = new TFile("TrackingPerformance.root", "recreate");
  hImpParXVsP->Write();
  hImpParYVsP->Write();
  hDeltaPxVsP->Write();
  hDeltaPyVsP->Write();
  hDeltaPzVsP->Write();
  hRelDeltaPxVsP->Write();
  hRelDeltaPyVsP->Write();
  hRelDeltaPzVsP->Write();
  hMomGen->Write();
  hEtaGen->Write();
  hMomTrackable->Write();
  hEtaTrackable->Write();
  hMomTracked5clu->Write();
  hEtaTracked5clu->Write();
  hMomAllReco5clu->Write();
  hMomGoodReco5clu->Write();
  hMomFakeReco5clu->Write();
  hEtaAllReco5clu->Write();
  hEtaGoodReco5clu->Write();
  hEtaFakeReco5clu->Write();
  hMomTrackedge4clu->Write();
  hEtaTrackedge4clu->Write();
  hMomAllRecoge4clu->Write();
  hMomGoodRecoge4clu->Write();
  hMomFakeRecoge4clu->Write();
  hEtaAllRecoge4clu->Write();
  hEtaGoodRecoge4clu->Write();
  hEtaFakeRecoge4clu->Write();
  outRoot->Close();
}
