#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TPaveStats.h>
#include <TParticle.h>
#include "NA6PTrack.h"
#include "NA6PMCEventHeader.h"
#include "NA6PDCAFitterN.h"
#include "MagneticField.h"
#include "ConfigurableParam.h"
#endif

void FillMeanAndRms(TH2F* h2D, TH1F* hMean, TH1F* hRms, TH1F* hSig)
{
  if (!h2D || !hMean || !hRms || !hSig) {
    printf("One of pointers is not set: h2D=%p hMean=%p hRms=%p hSig=%p\n", h2D, hMean, hRms, hSig);
    return;
  }
  int nBins = h2D->GetXaxis()->GetNbins();
  hMean->SetStats(0);
  hRms->SetStats(0);
  hSig->SetStats(0);
  for (int jp = 1; jp <= nBins; jp++) {
    TH1D* htmp1 = h2D->ProjectionY("htmp1", jp, jp);
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
    if (htmp1->GetEntries() > 20.) {
      int fitStatus = htmp1->Fit("gaus", "N0Q", "", mean - 3. * rms, mean + 3. * rms);
      if (fitStatus == 0) {
        TF1* fg = (TF1*)htmp1->GetListOfFunctions()->FindObject("gaus");
        if (fg) {
          gw = fg->GetParameter(2);
          egw = fg->GetParError(2);
        }
      }
    }
    hSig->SetBinContent(jp, gw);
    hSig->SetBinError(jp, egw);
    delete htmp1;
  }
}

void printDecay(TParticle& mothPart,
                std::vector<TParticle>* mcArr)
{

  int nDau = mothPart.GetNDaughters();
  int pdg = mothPart.GetPdgCode();
  printf("%d (zprod = %f) ->", pdg, mothPart.Vz());
  for (int jd = 0; jd < nDau; ++jd) {
    int indexDau = mothPart.GetFirstDaughter() + jd;
    auto dauPart = mcArr->at(indexDau);
    int pdgDau = dauPart.GetPdgCode();
    printf(" %d (zprod = %f)", pdgDau, dauPart.Vz());
    int nGranDau = dauPart.GetNDaughters();
    if (nGranDau > 0) {
      printf(" [ ");
      printDecay(dauPart, mcArr);
      printf(" ] ");
    }
  }
}

bool CheckCharmDecay(TParticle& mothPart,
                     std::vector<TParticle>* mcArr,
                     int pdgProngs[5],
                     int idDau[5],
                     float decvert[3])
{

  bool okDau[5] = {false, false, false, false, false};
  int pdg = mothPart.GetPdgCode();
  int nDau = mothPart.GetNDaughters();
  int finalProngs = nDau;
  int nExpectedProngs = 0;
  for (int i = 0; i < 5; ++i) {
    if (pdgProngs[i] != 0)
      nExpectedProngs++;
  }
  bool isDecVertSet = false;
  for (int jd = 0; jd < nDau; ++jd) {
    bool isExpected = false;
    int indexDau = mothPart.GetFirstDaughter() + jd;
    auto& dauPart = mcArr->at(indexDau);
    int pdgDau = dauPart.GetPdgCode();
    int absPdgDau = std::abs(pdgDau);
    for (int ld = 0; ld < 5; ++ld) {
      if (pdgDau == pdgProngs[ld]) {
        isExpected = true;
        break;
      }
    }
    if (absPdgDau == 313 || absPdgDau == 333)
      isExpected = true;
    if (isExpected && !isDecVertSet) {
      decvert[0] = dauPart.Vx();
      decvert[1] = dauPart.Vy();
      decvert[2] = dauPart.Vz();
      isDecVertSet = true;
    }
  }
  if (!isDecVertSet)
    return false;
  for (int jd = 0; jd < nDau; ++jd) {
    int indexDau = mothPart.GetFirstDaughter() + jd;
    auto& dauPart = mcArr->at(indexDau);
    int pdgDau = dauPart.GetPdgCode();
    int absPdgDau = std::abs(pdgDau);
    float dist2 = (dauPart.Vx() - decvert[0]) * (dauPart.Vx() - decvert[0]) +
                  (dauPart.Vy() - decvert[1]) * (dauPart.Vy() - decvert[1]) +
                  (dauPart.Vz() - decvert[2]) * (dauPart.Vz() - decvert[2]);

    if (pdgDau == 11 && dist2 > 1e-8) {
      // delta electron, should not count as daughter
      finalProngs--;
      continue;
    }
    bool found = false;
    for (int ld = 0; ld < 5; ++ld) {
      if (pdgDau == pdgProngs[ld] && okDau[ld] == false) {
        okDau[ld] = true;
        idDau[ld] = indexDau;
        found = true;
        break;
      }
    }
    if (!found) {
      if (absPdgDau == 313 || absPdgDau == 333) {
        int nDauRes = dauPart.GetNDaughters();
        finalProngs = finalProngs - 1 + nDauRes;
        for (int kd = 0; kd < nDauRes; ++kd) {
          int indexDau2 = dauPart.GetFirstDaughter() + kd;
          auto dauPart2 = mcArr->at(indexDau2);
          int pdgDau2 = dauPart2.GetPdgCode();
          for (int ld = 0; ld < 5; ++ld) {
            if (pdgDau2 == pdgProngs[ld] && okDau[ld] == false) {
              okDau[ld] = true;
              idDau[ld] = indexDau2;
              break;
            }
          }
        }
      }
    }
  }

  if (finalProngs != nExpectedProngs)
    return false;

  for (int ld = 0; ld < 5; ++ld) {
    if (pdgProngs[ld] != 0 && !okDau[ld])
      return false;
  }
  return true;
}

void CheckDecayVertices(const char* part = "D0",
                        const char* dirSimu = ".")
{
  na6p::conf::ConfigurableParam::updateFromFile(Form("%s/na6pLayout.ini",dirSimu), "", true);
  if (TGeoGlobalMagField::Instance()->GetField() == nullptr) {
    auto magField = new MagneticField();
    magField->loadField();
    magField->setAsGlobalField();
  }

  int pdgPart = 421;
  int nProngs = 2;
  int pdgProngs[5] = {0, 0, 0, 0, 0};
  float maxDL = 5000.;
  std::string decay;
  if (strcmp(part, "D0") == 0) {
    pdgPart = 421;
    nProngs = 2;
    pdgProngs[0] = 211;
    pdgProngs[1] = -321;
    maxDL = 1000.;
    decay = "D^{0} #rightarrow K^{-} #pi^{+}";
  } else if (strcmp(part, "D0bar") == 0) {
    pdgPart = -421;
    nProngs = 2;
    pdgProngs[0] = -211;
    pdgProngs[1] = 321;
    maxDL = 1000.;
    decay = "#bar{D^{0}} #rightarrow K^{+} #pi^{-}";
  } else if (strcmp(part, "Dplus") == 0) {
    pdgPart = 411;
    nProngs = 3;
    pdgProngs[0] = 211;
    pdgProngs[1] = 211;
    pdgProngs[2] = -321;
    maxDL = 3000.;
    decay = "D^{+} #rightarrow K^{-} #pi^{+} #pi^{+}";
  } else if (strcmp(part, "Dminus") == 0) {
    pdgPart = -411;
    nProngs = 3;
    pdgProngs[0] = -211;
    pdgProngs[1] = -211;
    pdgProngs[2] = 321;
    maxDL = 3000.;
    decay = "D^{-} #rightarrow K^{+} #pi^{-} #pi^{-}";
  } else if (strcmp(part, "Ds") == 0) {
    pdgPart = 431;
    nProngs = 3;
    pdgProngs[0] = 211;
    pdgProngs[1] = 321;
    pdgProngs[2] = -321;
    maxDL = 1000.;
    decay = "D_{s}^{+} #rightarrow K^{-} K^{+} #pi^{+}";
  } else if (strcmp(part, "DsMinus") == 0) {
    pdgPart = -431;
    nProngs = 3;
    pdgProngs[0] = -211;
    pdgProngs[1] = 321;
    pdgProngs[2] = -321;
    maxDL = 1000.;
    decay = "D_{s}^{-} #rightarrow K^{-} K^{+} #pi^{-}";
  } else if (strcmp(part, "LambdaC") == 0 || strcmp(part, "LcPlus") == 0) {
    pdgPart = 4122;
    nProngs = 3;
    pdgProngs[0] = 2212;
    pdgProngs[1] = 211;
    pdgProngs[2] = -321;
    maxDL = 1000.;
    decay = "#Lambda_{c}^{+} #rightarrow p K^{-} #pi^{+}";
  } else if (strcmp(part, "antiLambdaC") == 0 || strcmp(part, "LcMinus") == 0) {
    pdgPart = -4122;
    nProngs = 3;
    pdgProngs[0] = -2212;
    pdgProngs[1] = -211;
    pdgProngs[2] = 321;
    maxDL = 1000.;
    decay = "#bar{#Lambda}_{c}^{-} #rightarrow #bar{p} K^{+} #pi^{-}";
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

  NA6PDCAFitterN<2> df2;
  df2.setPropagateToPCA(true);
  df2.setUseAbsDCA(true);
  NA6PDCAFitterN<3> df3;
  df3.setPropagateToPCA(true);
  df3.setUseAbsDCA(true);

  TH1F* hGenVsRap = new TH1F("hGenVsRap", ";y;entries", 20, 1., 5.);
  TH1F* hGenVsPt = new TH1F("hGenVsPt", ";p_{T} (GeV/c);entries", 20, 0., 5.);
  TH1F* hRecoVsRap = new TH1F("hRecoVsRap", ";y;entries", 20, 1., 5.);
  TH1F* hRecoVsPt = new TH1F("hRecoVsPt", ";p_{T} (GeV/c);entries", 20, 0., 5.);
  TH1F* hDecLen = new TH1F("hDecLen", "; decay length L_{xyz} (#mum); entries", 100, 0., 5 * maxDL);
  TH1F* hPseudoPropDecLen = new TH1F("hPseudoPropDecLen", "; L_{xyz} M / p (#mum); entries", 100, 0., maxDL);
  TH1F* hInvMass = new TH1F("hInvMass", ";M (GeV/c^{2});entries", 200, 1.7, 2.5);
  TH2F* hInvMassVsRap = new TH2F("hInvMassVsRap", ";y;M (GeV/c^{2});entries", 20, 1., 5., 200, 1.7, 2.5);
  TH1F* hSecVx = new TH1F("hSecVx", ";x_{sec vert} (cm);entries", 100, -1., 1.);
  TH1F* hSecVy = new TH1F("hSecVy", ";y_{sec vert} (cm);entries", 100, -1., 1.);
  TH1F* hSecVz = new TH1F("hSecVz", ";z_{sec vert} (cm);entries", 100, -5., 5.);
  TH1F* hDeltaVx = new TH1F("hDeltaVx", ";x_{rec} - x_{true} (#mum);entries", 100, -500., 500.);
  TH1F* hDeltaVy = new TH1F("hDeltaVy", ";y_{rec} - y_{true} (#mum);entries", 100, -500., 500.);
  TH1F* hDeltaVz = new TH1F("hDeltaVz", ";z_{rec} - z_{true} (#mum);entries", 100, -500., 500.);
  TH2F* hDeltaVxVsRap = new TH2F("hDeltaVxVsRap", ";y;x_{rec} - x_{true} (#mum);entries", 20, 1., 5., 100, -500., 500.);
  TH2F* hDeltaVyVsRap = new TH2F("hDeltaVyVsRap", ";y;y_{rec} - y_{true} (#mum);entries", 20, 1., 5., 100, -500., 500.);
  TH2F* hDeltaVzVsRap = new TH2F("hDeltaVzVsRap", ";y;z_{rec} - z_{true} (#mum);entries", 20, 1., 5., 100, -500., 500.);

  int nGenerated = 0, nGoodNdau = 0, nGoodDecay = 0, nReconstructed = 0;

  int nEv = trTree->GetEntries();
  printf("Number of events = %d\n", nEv);
  for (int jEv = 0; jEv < nEv; jEv++) {
    mcTree->GetEvent(jEv);
    trTree->GetEvent(jEv);
    int nPart = mcArr->size();
    int nTracks = trArr->size();
    double xPrimVertGen = mcHead->getVX();
    double yPrimVertGen = mcHead->getVY();
    double zPrimVertGen = mcHead->getVZ();
    double bxyz[3];
    double xyz[3] = {xPrimVertGen, yPrimVertGen, zPrimVertGen};
    TGeoGlobalMagField::Instance()->Field(xyz, bxyz);
    df2.setBy(bxyz[1]);
    df3.setBy(bxyz[1]);
    if ((jEv % 1000) == 0)
      printf("--- Event %d particles %d ---\n", jEv, nPart);
    for (int jp = 0; jp < nPart; jp++) {
      auto curPart = mcArr->at(jp);
      int pdg = curPart.GetPdgCode();
      if (pdg != pdgPart)
        continue;
      float prodvert[3] = {(float)curPart.Vx(), (float)curPart.Vy(), (float)curPart.Vz()};
      float genmom = curPart.P();
      float genpt = curPart.Pt();
      float geny = curPart.Y();
      float mass = curPart.GetMass();
      ++nGenerated;
      int nDau = curPart.GetNDaughters();
      if (nDau == nProngs) {
        ++nGoodNdau;
      }
      int idDau[5] = {-1, -1, -1, -1, -1};
      float decvert[3] = {0., 0., -9999.};
      bool decOk = CheckCharmDecay(curPart, mcArr, pdgProngs, idDau, decvert);
      if (!decOk) {
        // printf(" Decay not ok\n");
        // printDecay(curPart, mcArr);
        // printf("\n");
        continue;
      }
      ++nGoodDecay;
      hGenVsRap->Fill(geny);
      hGenVsPt->Fill(genpt);
      float decLen = std::sqrt((decvert[0] - prodvert[0]) * (decvert[0] - prodvert[0]) +
                               (decvert[1] - prodvert[1]) * (decvert[1] - prodvert[1]) +
                               (decvert[2] - prodvert[2]) * (decvert[2] - prodvert[2]));

      hDecLen->Fill(decLen * 1e4);
      float ppdl = decLen * 1e4 * mass / genmom;
      hPseudoPropDecLen->Fill(ppdl);

      std::vector<NA6PTrack> prongs;
      int nRecoDau = 0;
      for (int kd = 0; kd < nProngs; ++kd) {
        for (int jTr = 0; jTr < nTracks; ++jTr) {
          NA6PTrack tr = trArr->at(jTr);
          int nClusters = tr.getNVTHits();
          if (nClusters < 5)
            continue;
          int mcLabel = tr.getParticleID();
          if (mcLabel == idDau[kd]) {
            prongs.push_back(tr);
            ++nRecoDau;
          }
        }
      }
      if (nRecoDau != nProngs)
        continue;
      ++nReconstructed;
      hRecoVsRap->Fill(geny);
      hRecoVsPt->Fill(genpt);
      int n = 0;
      if (nProngs == 2) {
        n = df2.process(prongs[0], prongs[1]);
      } else if (nProngs == 3) {
        n = df3.process(prongs[0], prongs[1], prongs[2]);
      }
      if (n == 0)
        continue;
      const auto& secondaryVertex = (nProngs == 2) ? df2.getPCACandidate() : df3.getPCACandidate();
      float dvx = (secondaryVertex[0] - decvert[0]) * 1e4;
      float dvy = (secondaryVertex[1] - decvert[1]) * 1e4;
      float dvz = (secondaryVertex[2] - decvert[2]) * 1e4;

      float sumPx = 0., sumPy = 0., sumPz = 0., sumE = 0.;
      for (int kd = 0; kd < nProngs; ++kd) {
        double pxyz[3];
        prongs[kd].getPXYZ(pxyz);
        sumPx += pxyz[0];
        sumPy += pxyz[1];
        sumPz += pxyz[2];
        float mom = prongs[kd].getP();
        auto mcPart = mcArr->at(prongs[kd].getParticleID());
        float mass = mcPart.GetMass();
        float energy = std::sqrt(mass * mass + mom * mom);
        sumE += energy;
      }
      float mom2 = sumPx * sumPx + sumPy * sumPy + sumPz * sumPz;
      float invMass = std::sqrt(sumE * sumE - mom2);

      hInvMass->Fill(invMass);
      hInvMassVsRap->Fill(geny, invMass);
      hSecVx->Fill(secondaryVertex[0]);
      hSecVy->Fill(secondaryVertex[1]);
      hSecVz->Fill(secondaryVertex[2]);
      hDeltaVx->Fill(dvx);
      hDeltaVy->Fill(dvy);
      hDeltaVz->Fill(dvz);
      hDeltaVxVsRap->Fill(geny, dvx);
      hDeltaVyVsRap->Fill(geny, dvy);
      hDeltaVzVsRap->Fill(geny, dvz);
    }
  }
  printf("Number of generated %s = %d\n", part, nGenerated);
  printf("Number of %s with %d daughters = %d\n", part, nProngs, nGoodNdau);
  printf("Number of %s with selected decay = %d\n", part, nGoodDecay);
  printf("Number of %s with reconstructed daughters = %d\n", part, nReconstructed);

  TH1F* hEffVsRap = (TH1F*)hRecoVsRap->Clone("hEffVsRap");
  hEffVsRap->Divide(hRecoVsRap, hGenVsRap, 1., 1., "B");
  hEffVsRap->GetYaxis()->SetTitle("Acceptance x Efficiency");
  hEffVsRap->SetStats(0);
  TH1F* hEffVsPt = (TH1F*)hRecoVsPt->Clone("hEffVsPt");
  hEffVsPt->Divide(hRecoVsPt, hGenVsPt, 1., 1., "B");
  hEffVsPt->GetYaxis()->SetTitle("Acceptance x Efficiency");
  hEffVsPt->SetStats(0);

  TH1F* hInvMassMean = new TH1F("hInvMassMean", ";y;<M> (GeV/c^{2})", hInvMassVsRap->GetNbinsX(), hInvMassVsRap->GetXaxis()->GetXmin(), hInvMassVsRap->GetXaxis()->GetXmax());
  TH1F* hInvMassRms = new TH1F("hInvMassRms", ";y;rms (M) (GeV/c^{2})", hInvMassVsRap->GetNbinsX(), hInvMassVsRap->GetXaxis()->GetXmin(), hInvMassVsRap->GetXaxis()->GetXmax());
  TH1F* hInvMassSig = new TH1F("hInvMassSig", ";y;#sigma(M) (GeV/c^{2})", hInvMassVsRap->GetNbinsX(), hInvMassVsRap->GetXaxis()->GetXmin(), hInvMassVsRap->GetXaxis()->GetXmax());
  FillMeanAndRms(hInvMassVsRap, hInvMassMean, hInvMassRms, hInvMassSig);

  TH1F* hDeltaVxMean = new TH1F("hDeltaVxMean", ";y;<x_{reco} - x_{gen}> (#mum)", hDeltaVxVsRap->GetNbinsX(), hDeltaVxVsRap->GetXaxis()->GetXmin(), hDeltaVxVsRap->GetXaxis()->GetXmax());
  TH1F* hDeltaVxRms = new TH1F("hDeltaVxRms", ";y;rms (x_{reco} - x_{gen}) (#mum)", hDeltaVxVsRap->GetNbinsX(), hDeltaVxVsRap->GetXaxis()->GetXmin(), hDeltaVxVsRap->GetXaxis()->GetXmax());
  TH1F* hDeltaVxSig = new TH1F("hDeltaVxSig", ";y;#sigma(x_{reco} - x_{gen}) (#mum)", hDeltaVxVsRap->GetNbinsX(), hDeltaVxVsRap->GetXaxis()->GetXmin(), hDeltaVxVsRap->GetXaxis()->GetXmax());
  TH1F* hDeltaVyMean = new TH1F("hDeltaVyMean", ";y;<y_{reco} - y_{gen}> (#mum)", hDeltaVyVsRap->GetNbinsX(), hDeltaVyVsRap->GetXaxis()->GetXmin(), hDeltaVyVsRap->GetXaxis()->GetXmax());
  TH1F* hDeltaVyRms = new TH1F("hDeltaVyRms", ";y;rms (y_{reco} - y_{gen}) (#mum)", hDeltaVyVsRap->GetNbinsX(), hDeltaVyVsRap->GetXaxis()->GetXmin(), hDeltaVyVsRap->GetXaxis()->GetXmax());
  TH1F* hDeltaVySig = new TH1F("hDeltaVySig", ";y;#sigma(y_{reco} - y_{gen}) (#mum)", hDeltaVyVsRap->GetNbinsX(), hDeltaVyVsRap->GetXaxis()->GetXmin(), hDeltaVyVsRap->GetXaxis()->GetXmax());
  TH1F* hDeltaVzMean = new TH1F("hDeltaVzMean", ";y;<z_{reco} - z_{gen}> (#mum)", hDeltaVzVsRap->GetNbinsX(), hDeltaVzVsRap->GetXaxis()->GetXmin(), hDeltaVzVsRap->GetXaxis()->GetXmax());
  TH1F* hDeltaVzRms = new TH1F("hDeltaVzRms", ";y;rms (z_{reco} - z_{gen}) (#mum)", hDeltaVzVsRap->GetNbinsX(), hDeltaVzVsRap->GetXaxis()->GetXmin(), hDeltaVzVsRap->GetXaxis()->GetXmax());
  TH1F* hDeltaVzSig = new TH1F("hDeltaVzSig", ";y;#sigma(z_{reco} - z_{gen}) (#mum)", hDeltaVzVsRap->GetNbinsX(), hDeltaVzVsRap->GetXaxis()->GetXmin(), hDeltaVzVsRap->GetXaxis()->GetXmax());
  FillMeanAndRms(hDeltaVxVsRap, hDeltaVxMean, hDeltaVxRms, hDeltaVxSig);
  FillMeanAndRms(hDeltaVyVsRap, hDeltaVyMean, hDeltaVyRms, hDeltaVySig);
  FillMeanAndRms(hDeltaVzVsRap, hDeltaVzMean, hDeltaVzRms, hDeltaVzSig);

  TLatex* tdec = new TLatex(0.45, 0.83, decay.data());
  tdec->SetNDC();

  TCanvas* cdl = new TCanvas("cdl", "Decay Length (gen)", 1400, 600);
  cdl->Divide(2, 1);
  cdl->cd(1);
  gPad->SetLogy();
  hDecLen->Draw();
  cdl->cd(2);
  gPad->SetLogy();
  hPseudoPropDecLen->Draw();
  tdec->Draw();

  TCanvas* ceff = new TCanvas("ceff", "Acc x Eff", 800, 800);
  ceff->Divide(2, 2);
  ceff->cd(1);
  gPad->SetTickx();
  gPad->SetTicky();
  hGenVsRap->SetMarkerStyle(21);
  hGenVsRap->SetMarkerColor(1);
  hGenVsRap->SetLineColor(1);
  hRecoVsRap->SetMarkerStyle(20);
  hRecoVsRap->SetMarkerColor(kGreen + 2);
  hRecoVsRap->SetLineColor(kGreen + 2);
  hGenVsRap->SetMinimum(0.);
  hGenVsRap->Draw("E");
  gPad->Update();
  TPaveStats* stg = (TPaveStats*)hGenVsRap->GetListOfFunctions()->FindObject("stats");
  if (stg) {
    stg->SetY1NDC(0.72);
    stg->SetY2NDC(0.92);
    stg->SetTextColor(hGenVsRap->GetMarkerColor());
  }
  hRecoVsRap->Draw("Esames");
  gPad->Update();
  TPaveStats* str = (TPaveStats*)hRecoVsRap->GetListOfFunctions()->FindObject("stats");
  if (str) {
    str->SetY1NDC(0.51);
    str->SetY2NDC(0.71);
    str->SetTextColor(hRecoVsRap->GetMarkerColor());
  }
  hRecoVsRap->Draw("Esames");
  gPad->Modified();
  ceff->cd(2);
  gPad->SetTickx();
  gPad->SetTicky();
  gPad->SetLogy();
  hGenVsPt->SetMarkerStyle(21);
  hGenVsPt->SetMarkerColor(1);
  hGenVsPt->SetLineColor(1);
  hRecoVsPt->SetMarkerStyle(20);
  hRecoVsPt->SetMarkerColor(kGreen + 2);
  hRecoVsPt->SetLineColor(kGreen + 2);
  hGenVsPt->Draw("E");
  gPad->Update();
  TPaveStats* stg2 = (TPaveStats*)hGenVsPt->GetListOfFunctions()->FindObject("stats");
  if (stg2) {
    stg2->SetY1NDC(0.72);
    stg2->SetY2NDC(0.92);
    stg2->SetTextColor(hGenVsPt->GetMarkerColor());
  }
  hRecoVsPt->Draw("Esames");
  gPad->Update();
  TPaveStats* str2 = (TPaveStats*)hRecoVsPt->GetListOfFunctions()->FindObject("stats");
  if (str2) {
    str2->SetY1NDC(0.51);
    str2->SetY2NDC(0.71);
    str2->SetTextColor(hRecoVsPt->GetMarkerColor());
  }
  gPad->Modified();
  tdec->Draw();
  ceff->cd(3);
  gPad->SetTickx();
  gPad->SetTicky();
  hEffVsRap->SetMarkerStyle(20);
  hEffVsRap->SetMarkerColor(kGreen + 2);
  hEffVsRap->SetLineColor(kGreen + 2);
  hEffVsRap->SetMinimum(0.);
  hEffVsRap->SetMaximum(1.);
  hEffVsRap->Draw();
  ceff->cd(4);
  gPad->SetTickx();
  gPad->SetTicky();
  hEffVsPt->SetMarkerStyle(20);
  hEffVsPt->SetMarkerColor(kGreen + 2);
  hEffVsPt->SetLineColor(kGreen + 2);
  hEffVsPt->SetMinimum(0.);
  hEffVsPt->SetMaximum(1.);
  hEffVsPt->Draw();
  ceff->SaveAs(Form("%s-AccEff.png", part));

  TCanvas* cm = new TCanvas("cm", "Inv Mass", 1400, 500);
  cm->Divide(3, 1);
  cm->cd(1);
  gPad->SetLeftMargin(0.13);
  gPad->SetRightMargin(0.07);
  gPad->SetTickx();
  gPad->SetTicky();
  hInvMass->Draw();
  cm->cd(2);
  gPad->SetLeftMargin(0.13);
  gPad->SetRightMargin(0.07);
  gPad->SetTickx();
  gPad->SetTicky();
  hInvMassMean->SetStats(0);
  hInvMassMean->SetMarkerStyle(20);
  hInvMassMean->SetMarkerColor(1);
  hInvMassMean->SetLineColor(1);
  hInvMassMean->SetMinimum(hInvMass->GetMean() - 5. * hInvMass->GetRMS());
  hInvMassMean->SetMaximum(hInvMass->GetMean() + 5. * hInvMass->GetRMS());
  hInvMassMean->GetYaxis()->SetTitleOffset(1.7);
  hInvMassMean->Draw();
  tdec->Draw();
  cm->cd(3);
  gPad->SetLeftMargin(0.13);
  gPad->SetRightMargin(0.07);
  gPad->SetTickx();
  gPad->SetTicky();
  hInvMassSig->SetStats(0);
  hInvMassSig->SetMarkerStyle(20);
  hInvMassSig->SetMarkerColor(1);
  hInvMassSig->SetLineColor(1);
  hInvMassSig->GetYaxis()->SetTitleOffset(1.7);
  hInvMassSig->SetMinimum(0.);
  hInvMassSig->SetMaximum(0.04);
  hInvMassSig->Draw();
  cm->SaveAs(Form("%s-InvMass.png", part));

  TCanvas* cps = new TCanvas("cps", "Dec Vert Pos", 1400, 800);
  cps->Divide(3, 2);
  cps->cd(1);
  hSecVx->Draw();
  cps->cd(2);
  hSecVy->Draw();
  cps->cd(3);
  hSecVz->Draw();
  cps->cd(4);
  hDeltaVx->Draw();
  cps->cd(5);
  hDeltaVy->Draw();
  cps->cd(6);
  hDeltaVz->Draw();

  TCanvas* csig = new TCanvas("csig", "Dec Vert Resid", 1400, 800);
  csig->Divide(3, 2);
  csig->cd(1);
  gPad->SetTickx();
  gPad->SetTicky();
  hDeltaVxMean->SetStats(0);
  hDeltaVxMean->SetMarkerStyle(20);
  hDeltaVxMean->SetMarkerColor(1);
  hDeltaVxMean->SetLineColor(1);
  hDeltaVxMean->SetMinimum(-5.);
  hDeltaVxMean->SetMaximum(5.);
  hDeltaVxMean->GetYaxis()->SetTitleOffset(1.2);
  hDeltaVxMean->Draw();
  csig->cd(2);
  gPad->SetTickx();
  gPad->SetTicky();
  hDeltaVyMean->SetStats(0);
  hDeltaVyMean->SetMarkerStyle(20);
  hDeltaVyMean->SetMarkerColor(1);
  hDeltaVyMean->SetLineColor(1);
  hDeltaVyMean->SetMinimum(-5.);
  hDeltaVyMean->SetMaximum(5.);
  hDeltaVyMean->GetYaxis()->SetTitleOffset(1.2);
  hDeltaVyMean->Draw();
  tdec->Draw();
  csig->cd(3);
  gPad->SetTickx();
  gPad->SetTicky();
  hDeltaVzMean->SetStats(0);
  hDeltaVzMean->SetMarkerStyle(20);
  hDeltaVzMean->SetMarkerColor(1);
  hDeltaVzMean->SetLineColor(1);
  hDeltaVzMean->SetMinimum(-20.);
  hDeltaVzMean->SetMaximum(20.);
  hDeltaVzMean->GetYaxis()->SetTitleOffset(1.2);
  hDeltaVzMean->Draw();
  csig->cd(4);
  gPad->SetTickx();
  gPad->SetTicky();
  hDeltaVxSig->SetStats(0);
  hDeltaVxSig->SetMarkerStyle(20);
  hDeltaVxSig->SetMarkerColor(1);
  hDeltaVxSig->SetLineColor(1);
  hDeltaVxSig->SetMinimum(0.);
  hDeltaVxSig->SetMaximum(40.);
  hDeltaVxSig->GetYaxis()->SetTitleOffset(1.2);
  hDeltaVxSig->Draw();
  csig->cd(5);
  gPad->SetTickx();
  gPad->SetTicky();
  hDeltaVySig->SetStats(0);
  hDeltaVySig->SetMarkerStyle(20);
  hDeltaVySig->SetMarkerColor(1);
  hDeltaVySig->SetLineColor(1);
  hDeltaVySig->SetMinimum(0.);
  hDeltaVySig->SetMaximum(40.);
  hDeltaVySig->GetYaxis()->SetTitleOffset(1.2);
  hDeltaVySig->Draw();
  csig->cd(6);
  gPad->SetTickx();
  gPad->SetTicky();
  hDeltaVzSig->SetStats(0);
  hDeltaVzSig->SetMarkerStyle(20);
  hDeltaVzSig->SetMarkerColor(1);
  hDeltaVzSig->SetLineColor(1);
  hDeltaVzSig->SetMinimum(0.);
  hDeltaVzSig->SetMaximum(200.);
  hDeltaVzSig->GetYaxis()->SetTitleOffset(1.2);
  hDeltaVzSig->Draw();
  csig->SaveAs(Form("%s-DecayVertPos.png", part));

  TFile* fout = new TFile(Form("%s-DecayVert.root", part), "recreate");
  hDecLen->Write();
  hPseudoPropDecLen->Write();
  hInvMassVsRap->Write();
  hInvMassMean->Write();
  hInvMassRms->Write();
  hInvMassSig->Write();
  hSecVx->Write();
  hSecVy->Write();
  hSecVz->Write();
  hDeltaVxVsRap->Write();
  hDeltaVyVsRap->Write();
  hDeltaVzVsRap->Write();
  hDeltaVxMean->Write();
  hDeltaVxRms->Write();
  hDeltaVxSig->Write();
  hDeltaVyMean->Write();
  hDeltaVyRms->Write();
  hDeltaVySig->Write();
  hDeltaVzMean->Write();
  hDeltaVzRms->Write();
  hDeltaVzSig->Write();
  hGenVsRap->Write();
  hRecoVsRap->Write();
  hEffVsRap->Write();
  hGenVsPt->Write();
  hRecoVsPt->Write();
  hEffVsPt->Write();
  fout->Close();
}
