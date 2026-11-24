// NA6PCCopyright

#include "NA6PGenParam.h"
#include "NA6PBeamParam.h"
#include "NA6PGenCutParam.h"
#include "NA6PMCStack.h"
#include <TDatabasePDG.h>
#include <TRandom.h>
#include <TMath.h>
#include <cmath>

NA6PGenParam::NA6PGenParam(const std::string& name, int pdg, float mult, const std::string& parTrans, const std::string& parLong,
                           float ptmin, float ptmax, float lmin, float lmax, bool istransPt, bool islongY, bool isPoisson) : NA6PGenerator(name), mPDGCode(pdg), mMult(mult), mPoisson(isPoisson), mTransIsPt(istransPt), mLongIsY(islongY)
{
  setFunTrans(parTrans, ptmin, ptmax);
  setFunLong(parLong, lmin, lmax);
}

void NA6PGenParam::setFunTrans(const std::string& fnameTrans, float ptmin, float ptmax)
{
  // requires name of interpreted function
  if (!mTransIsPt) { // convert to Mt
    if (!TDatabasePDG::Instance()->GetParticle(mPDGCode)) {
      LOGP(fatal, "Particle PDG={} is not defined", mPDGCode);
    }
    auto m = TDatabasePDG::Instance()->GetParticle(mPDGCode)->Mass();
    ptmin = std::sqrt(m * m + ptmin * ptmin);
    ptmax = std::sqrt(m * m + ptmax * ptmax);
  }
  mFunTrans = std::make_unique<TF1>(fnameTrans.c_str(), fnameTrans.c_str(), ptmin, ptmax);
}

void NA6PGenParam::setFunLong(const std::string& fnameLong, float lmin, float lmax)
{
  // requires name of interpreted function
  mFunLong = std::make_unique<TF1>(fnameLong.c_str(), fnameLong.c_str(), lmin, lmax);
}

void NA6PGenParam::init()
{
  if (isInitDone()) {
    return;
  }
  initMultiplicity();

  NA6PGenerator::init();
}

void NA6PGenParam::initMultiplicity()
{
  if (!TDatabasePDG::Instance()->GetParticle(mPDGCode)) {
    LOGP(fatal, "Particle PDG={} is not defined", mPDGCode);
  }
  auto getLbl = [this](int i) { return i == 0 ? (mTransIsPt ? "Pt" : "Mt") : (mLongIsY ? "Y" : "Eta"); };
  TF1* ff[2] = {mFunTrans.get(), mFunLong.get()};
  for (int i = 0; i < 2; i++) {
    auto f = ff[i];
    if (!f) {
      LOGP(fatal, "{}-parametrization for PDG={} is not set", getLbl(i), mPDGCode);
    }
    if (f->GetNpar() > 0 && !f->TestBit(ParamsSetBit)) {
      LOGP(fatal, "{}-parametrization {} for PDG={} needs {} parameters which are not set", getLbl(i), f->GetName(), mPDGCode, f->GetNpar());
    }
    if (!f->IsValid()) {
      LOGP(fatal, "{}-parametrization {} for PDG={} is not valid", getLbl(i), f->GetName(), mPDGCode);
    }
  }
  if (mdNdLong > 0) { // dneta or dndy at was provided, calculate multiplicity to generate
    const auto& beamParam = NA6PBeamParam::Instance();
    auto ycm = beamParam.getYCM();
    // 1st check which fraction of pTs we generate
    double scalePt = 1., scaleLong = 1.;
    {
      auto tmin = mFunTrans->GetXmin(), tmax = mFunTrans->GetXmax();
      double tmin0 = 0, tmax0 = 20.;
      if (!mTransIsPt) { // convert to Mt
        auto m = TDatabasePDG::Instance()->GetParticle(mPDGCode)->Mass();
        tmin0 = std::sqrt(m * m + tmin0 * tmin0);
        tmax0 = std::sqrt(m * m + tmax0 * tmax0);
      }
      mFunTrans->SetRange(tmin0, tmax0);
      mFunTrans->SetNpx(std::max(100, int(100 * (tmax0 - tmin0))));
      auto integTransTot = mFunTrans->Integral(tmin0, tmax0);
      auto integTransRange = mFunTrans->Integral(tmin, tmax);
      mFunTrans->SetRange(tmin, tmax);
      scalePt = integTransRange / integTransTot;
    }
    // check which fraction of mid-rapidity we generate
    {
      auto tmin = mFunLong->GetXmin(), tmax = mFunLong->GetXmax();
      double tmin0 = ycm - 0.5, tmax0 = ycm + 0.5;
      mFunLong->SetRange(tmin0, tmax0);
      mFunLong->SetNpx(100);
      auto integLongMid = mFunLong->Integral(tmin0, tmax0);
      auto integLongRange = mFunLong->Integral(tmin, tmax);
      mFunLong->SetRange(tmin, tmax);
      scaleLong = integLongRange / integLongMid;
    }
    mMult = mdNdLong * scalePt * scaleLong;
    LOGP(info, "Renormalizing dN/d{} = {:.2f} to mean multiplicity {} in requested phase space (factors {:.2f} and {:.2f} from transverse and longitudinal distributions)",
         mLongIsY ? "Y" : "Eta", mdNdLong, mMult, scalePt, scaleLong);
  }
  if (mMult <= 0) {
    LOGP(fatal, "Multiplicity for particle PDG={} is {}", mPDGCode, mMult);
  }
  mFunTrans->SetNpx(std::max(100, int(100 * (mFunTrans->GetXmax() - mFunTrans->GetXmin()))));
  mFunLong->SetNpx(std::max(100, int(100 * (mFunLong->GetXmax() - mFunLong->GetXmin()))));
  LOGP(info, "Initialized generator {} for PDG={}, {} ({}) particles in the requested phase-space", getName(), mPDGCode, mMult, isPoisson() ? "Poisson" : "fixed");
  LOGP(info, "Transverse distribution in {:.2f}<{}<{:.2f}", mFunTrans->GetXmin(), isTransPt() ? "Pt" : "Mt", mFunTrans->GetXmax());
  mFunTrans->Print();
  LOGP(info, "Longitudinal distribution in {:.2f}<{}<{:.2f}", mFunLong->GetXmin(), isLongY() ? "Y" : "Eta", mFunLong->GetXmax());
  mFunLong->Print();
}

void NA6PGenParam::setParametersTrans(const std::vector<float>& v)
{
  int npf = mFunTrans->GetNpar();
  if (npf != int(v.size())) {
    LOGP(fatal, "{}-parametrization {} needs {} parameters, {} provided", mTransIsPt ? "Pt" : "Mt", mFunTrans->GetName(), npf, v.size());
  }
  for (int i = 0; i < npf; i++) {
    mFunTrans->SetParameter(i, v[i]);
  }
  mFunTrans->SetBit(ParamsSetBit);
}

void NA6PGenParam::setParametersLong(const std::vector<float>& v)
{
  int npf = mFunLong->GetNpar();
  if (npf != int(v.size())) {
    LOGP(fatal, "{}-parametrization {} needs {} parameters, {} provided", mLongIsY ? "Y" : "Eta", mFunTrans->GetName(), npf, v.size());
  }
  for (int i = 0; i < npf; i++) {
    mFunLong->SetParameter(i, v[i]);
  }
  mFunLong->SetBit(ParamsSetBit);
}

void NA6PGenParam::generate()
{
  generatePrimaryVertex(); // will generate if not generated yet, otherwise, will set the origin from the MCHeader of the stack.
  const auto& param = NA6PGenCutParam::Instance();
  auto m = TDatabasePDG::Instance()->GetParticle(mPDGCode)->Mass();
  int nPart = mPoisson ? gRandom->Poisson(mMult) : TMath::Nint(mMult);
  for (int i = 0; i < nPart; i++) {
    int ntrial = 0;
    do {
      double vTrans, vLong, phi, pt2, pt, mt2, mt, pZ, en;
      try {
        vTrans = mFunTrans->GetRandom();
        vLong = mFunLong->GetRandom();
        phi = gRandom->Uniform(0, TMath::TwoPi());
        pt2 = vTrans * vTrans;
        mt2 = pt2;
        pt = vTrans;
        mt = pt;
        if (!mTransIsPt) { // mT
          pt2 -= m * m;
          pt = std::sqrt(pt2);
        } else {
          mt2 += m * m;
          mt = std::sqrt(mt2);
        }
        if (mLongIsY) {
          pZ = mt * std::sinh(vLong);
        } else { // eta
          auto tht = 2. * std::atan(std::exp(-vLong));
          pZ = pt / std::tan(tht);
        }
        en = std::sqrt(mt2 + pZ * pZ);
      } catch (std::exception const& e) {
        LOGP(info, "Exception at trial {} of generating particle PDG={}, VTrans={}, vLong={}, error: {}", ntrial, mPDGCode, vTrans, vLong, e.what());
        if (ntrial++ > param.maxTrialsPerParticle) {
          LOGP(fatal, "Failed to generate track PDG:{} after max allowed {} trials, check generator", mPDGCode, ntrial);
        }
        continue;
      }
      ntrial = 0;
      int dummy = 0;
      float sn = std::sin(phi), cs = std::cos(phi);
      getStack()->PushTrack(true, -1, mPDGCode, pt * cs, pt * sn, pZ, en, getOriginX(), getOriginY(), getOriginZ(),
                            0., 0., 0., 0., TMCProcess::kPPrimary, dummy, 1., 0);
    } while (ntrial != 0);
  }
  // register generatot header in the MCHeader
  auto mcHead = getStack()->getEventHeader();
  static std::string info = fmt::format("{}_pdg{}x{}", getName(), mPDGCode, getMultiplicity());
  mcHead->getGenHeaders().emplace_back(nPart, 0, mcHead->getNPrimaries(), 0, info);
  mcHead->incNPrimaries(nPart);
}

void NA6PGenParam::setdNdY(float v)
{
  if (v <= 0.) {
    LOGP(fatal, "dNdY cannot be negative, {} requested", v);
  }
  mdNdLong = v;
}

void NA6PGenParam::setdNdEta(float v)
{
  if (v <= 0.) {
    LOGP(fatal, "dNdeta cannot be negative, {} requested", v);
  }
  mdNdLong = v;
}

void NA6PGenParam::setLongIsY(bool v)
{
  if (mFunLong) {
    LOGP(fatal, "setLongIsY(..) must be called before setting the longitudinal parametrization");
  }
  mLongIsY = v;
}

void NA6PGenParam::setTransIsPt(bool v)
{
  if (mFunTrans) {
    LOGP(fatal, "setTransIsPt(..) must be called before setting the transverse parametrization");
  }
  mTransIsPt = v;
}
