#include <Riostream.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TAxis.h>
#include <TArc.h>
#include <TNtuple.h>
#include <TH2.h>
#include "NA6PNuclearShape.h"
#include "NA6PGlauber.h"

void runOpticalGlauber(float bmin = 0,
                       float bmax = 18,
                       int nstep = 72,
                       bool doProfiles = false)
{

  NA6PNuclearShape nucA, nucB;
  nucA.setNucleusPb();
  nucB.setNucleusPb();
  float sigmb = 60.;

  NA6PGlauber gl(nucA, nucB, sigmb);
  TNtuple* nt = new TNtuple("nt", "nt", "b:Ta:Tb:TaTb:sigmain:xmed:ymed:xrms:yrms:ecc:part:coll:pdensmed:pdensmax:L");
  float array[20];
  FILE* fil;
  fil = fopen("glauber_optical.txt", "w");
  fprintf(fil,
          " %7s %9s %9s %8s %8s %8s %8s %9s %9s %9s %9s %10s %11s %10s %11s %9s\n",
          "b", "Tab", "ProbInt", "<x>", "<y>", "<x2>", "<y2>",
          "eccentr", "<Npart>", "<Npart>_i", "<Ncoll>",
          "<PartDens>", "PartDensMax", "<CollDens>",
          "CollDensMax", "<L>");

  float APro = nucA.getA();
  float ATgt = nucB.getA();

  float db = (bmax - bmin) / nstep;
  TGraph* gpar = new TGraph();
  TGraph* gcol = new TGraph();

  std::vector<float> bForProfiles;

  for (int ib = 0; ib < nstep; ++ib) {
    float b = bmin + ib * db;
    if (std::abs(std::floor(b) - b) < 0.001) {
      bForProfiles.push_back(b);
    }
    gl.doOpticalIntegrals(b);
    float tta = gl.getTa();
    float ttb = gl.getTb();
    float ttab = gl.getTAB();
    float sigmaIn = gl.getInteractionProbability();
    float nPartInt = gl.getAverageNPartInInteractions();
    float nPart = gl.getAverageNPart();
    float partDensAve = gl.getAveragePartDens();
    float partDensMax = gl.getMaximumPartDens();
    float collDensAve = gl.getAverageCollDens();
    float collDensMax = gl.getMaximumCollDens();
    float x = gl.getMeanX();
    float y = gl.getMeanY();
    float xx = gl.getVarX();
    float yy = gl.getVarY();
    float ecc = gl.getEccentricity();
    float nColl = gl.getAverageNColl();
    float L = gl.getRhoL() / 0.17;
    if (b * 2 == int(b * 2)) {
      cout << "b= " << b << " Npart=" << nPart << " Ncoll=" << nColl << " L=" << L << endl;
    }
    gpar->SetPoint(ib, b, nPart);
    gcol->SetPoint(ib, b, nColl);

    fprintf(fil,
            " %7.3f %9.4g %9.4f %8.2f %8.2f %8.2f %8.2f %9.3f %9.3f %9.3f %9.3f %10.3f %11.3f %10.3f %11.3f %9.3f\n",
            b, ttab, sigmaIn, x, y, xx, yy, ecc, nPart, nPartInt, nColl, partDensAve, partDensMax, collDensAve, collDensMax, L);
    array[0] = b;
    array[1] = tta;
    array[2] = ttb;
    array[3] = ttab;
    array[4] = sigmaIn;
    array[5] = x;
    array[6] = y;
    array[7] = xx;
    array[8] = yy;
    array[9] = ecc;
    array[10] = nPart;
    array[11] = nColl;
    array[12] = partDensAve;
    array[13] = partDensMax;
    array[14] = L;
    nt->Fill(array);
  }
  fclose(fil);

  TCanvas* c1 = new TCanvas("c1", "", 700, 900);
  c1->Divide(1, 2);
  c1->cd(1);
  gpar->Draw("AL");
  gpar->GetXaxis()->SetTitle("b (fm)");
  gpar->GetYaxis()->SetTitle("Number of Participants");
  c1->cd(2);
  gcol->Draw("AL");
  gcol->GetXaxis()->SetTitle("b (fm)");
  gcol->GetYaxis()->SetTitle("Number of Collisions");
  c1->Update();

  // Draw density profiles
  int nProfs = bForProfiles.size();
  TH2F** hProfC = new TH2F*[nProfs];
  TH2F** hProfP = new TH2F*[nProfs];

  if (doProfiles) {
    double maxrad = 15.;
    double dx = 0.1, dy = 0.1;
    double nptx = (2. * maxrad) / dx + 1.00001;
    double npty = (2. * maxrad) / dy + 1.00001;
    double xmin = -maxrad;
    double ymin = -maxrad;
    double xmax = maxrad;
    double ymax = maxrad;

    for (int j = 0; j < nProfs; ++j) {
      float b = bForProfiles.at(j);
      hProfC[j] = new TH2F(Form("hCollProf%d", (int)b), Form("Collision density b = %.0f fm; x (fm); y (fm)", b), nptx, xmin, xmax, npty, ymin, ymax);
      hProfC[j]->SetStats(0);
      hProfP[j] = new TH2F(Form("hPartProf%d", (int)b), Form("Paticipant density b = %.0f fm; x (fm); y (fm)", b), nptx, xmin, xmax, npty, ymin, ymax);
      hProfP[j]->SetStats(0);
      for (int ix = 0; ix < nptx; ++ix) {
        double x = xmin + ix * dx;
        for (int iy = 0; iy < npty; ++iy) {
          double y = ymin + iy * dy;
          double colld = gl.calcCollDensFromCoord(x, y, b);
          hProfC[j]->SetBinContent(iy + 1, ix + 1, colld);
          double partd = gl.calcPartDensFromCoord(x, y, b);
          hProfP[j]->SetBinContent(iy + 1, ix + 1, partd);
        }
      }
    }

    if (nProfs > 0) {
      TCanvas* c2 = new TCanvas("c2", "", 1200, 600);
      if (nProfs > 1) {
        c2->Divide(2, 1);
        c2->cd(1);
      }
      hProfP[0]->Draw("colz");
      float b = bForProfiles.at(0);
      double rA = nucA.getRadius();
      double rB = nucB.getRadius();
      TArc* nA = new TArc(-b / 2., 0., rA);
      TArc* nB = new TArc(b / 2., 0., rB);
      nA->SetLineColor(6);
      nB->SetLineColor(6);
      nA->SetFillStyle(0);
      nB->SetFillStyle(0);
      nA->SetLineWidth(2);
      nB->SetLineWidth(2);
      nA->Draw();
      nB->Draw();
      if (nProfs > 1) {
        c2->cd(2);
        hProfP[nProfs - 1]->Draw("colz");
        b = bForProfiles.at(nProfs - 1);
        TArc* nA2 = new TArc(-b / 2., 0., rA);
        TArc* nB2 = new TArc(b / 2., 0., rB);
        nA2->SetLineColor(6);
        nB2->SetLineColor(6);
        nA2->SetFillStyle(0);
        nB2->SetFillStyle(0);
        nA2->SetLineWidth(2);
        nB2->SetLineWidth(2);
        nA2->Draw();
        nB2->Draw();
      }
    }
  }

  TFile* fnt = new TFile("glauber_optical.root", "recreate");
  fnt->cd();
  nt->Write();
  if (doProfiles) {
    for (int j = 0; j < nProfs; ++j) {
      hProfC[j]->Write();
      hProfP[j]->Write();
    }
  }
}
