// mcp_dy.cc -- Drell-Yan production of a heavy fractionally charged fermion
//
//   q qbar -> gamma* -> chi chibar        (PURE photon exchange, no Z)
//
// Implemented as a Pythia semi-internal 2 -> 2 process so that
//   (a) the mediator is strictly the photon (no Z pole, no gamma*/Z
//       interference), matching a photon-only reference cross section for a
//       fractionally charged particle, whose Z coupling is model dependent;
//   (b) the final-state mass enters BOTH the matrix element and the phase
//       space (Pythia's built-in ffbar2ffbar(s:gmZ) treats the final state as
//       massless in the phase space, so it cannot be used for a heavy chi).
//
// With MCP:includeZ = on the Z is included as well, with the couplings an MCP
// actually has when the kinetic mixing is with hypercharge (as it must be above
// the EW scale): after EWSB the chi behaves as a fermion with T3 = 0 and
// electric charge Q_chi, i.e. a PURE VECTOR Z coupling
//     v_chi = -4 sin^2(thetaW) Q_chi,   a_chi = 0,
// equivalently g_Z = Q_chi e tan(thetaW).  Because a_chi = 0 the final-state
// Lorentz structure is unchanged (and there is no forward-backward asymmetry),
// so the Z enters only by reweighting the sHat spectrum.
//
// The chi is given electric charge 1 in the matrix element; the cross section
// for charge eps is obtained by scaling with eps^2.  The eta spectrum is
// independent of eps.
//
// Usage: ./mcp_dy <mass/GeV> <nEvent> <seed> <eCM/GeV> <outfile>

#include "Pythia8/Pythia.h"
#include <cstdio>

using namespace Pythia8;

const int ID_MCP = 6000015;

//==========================================================================

class Sigma2qqbar2mcp : public Sigma2Process {

public:

  Sigma2qqbar2mcp(int idMCPin, double qMCPin)
    : idMCP(idMCPin), qMCP(qMCPin), useZ(false), m2Res(0.), GamMRat(0.),
      thetaWRat(0.), vChi(0.), gamProp(0.), intProp(0.), resProp(0.),
      angFac(0.) {}

  virtual void initProc() {
    useZ      = flag("MCP:includeZ");
    double mZ = particleDataPtr->m0(23);
    m2Res     = mZ * mZ;
    GamMRat   = particleDataPtr->mWidth(23) / mZ;
    thetaWRat = 1. / (16. * coupSMPtr->sin2thetaW()
              * (1. - coupSMPtr->sin2thetaW()));
    // T3 = 0, charge Q_chi  ->  a_chi = 0, v_chi = -4 sin^2(thetaW) Q_chi.
    vChi      = -4. * coupSMPtr->sin2thetaWbar() * qMCP;
  }

  // sHat-dependent (flavour-independent) part of d(sigmaHat)/d(tHat).
  virtual void sigmaKin() {
    // Velocity of chi in the qqbar rest frame, and the scattering angle.
    double beta2 = 1. - 4. * s3 / sH;
    double beta  = sqrtpos(beta2);
    double cThe  = (beta > 0.) ? (tH - uH) / (sH * beta) : 0.;
    double c2    = cThe * cThe;
    // The chi coupling is pure vector for both mediators, so the transverse and
    // longitudinal pieces share one angular factor:
    //   (1 + cos^2) + (1 - beta^2)(1 - cos^2),
    // and there is no cos(theta) (forward-backward) term.
    angFac  = (1. + c2) + (1. - beta2) * (1. - c2);
    gamProp = M_PI * pow2(alpEM) / sH2;
    if (useZ) {
      double den = pow2(sH - m2Res) + pow2(sH * GamMRat);
      intProp = gamProp * 2. * thetaWRat * sH * (sH - m2Res) / den;
      resProp = gamProp * pow2(thetaWRat * sH) / den;
    } else {
      intProp = 0.;
      resProp = 0.;
    }
  }

  // Incoming-flavour-dependent part: the quark photon and Z couplings.
  // The 1/3 is the colour average for a colour-singlet final state.
  virtual double sigmaHat() {
    int    idA = abs(id1);
    double ei  = coupSMPtr->ef(idA);
    double vi  = coupSMPtr->vf(idA);
    double ai  = coupSMPtr->af(idA);
    double coef = ei * ei * gamProp * (qMCP * qMCP)
                + ei * vi * intProp * (qMCP * vChi)
                + (vi * vi + ai * ai) * resProp * (vChi * vChi);
    return coef * angFac / 3.;
  }

  // Outgoing flavours and colour flow (q qbar annihilation to singlet).
  virtual void setIdColAcol() {
    int id3New = (id1 > 0) ? idMCP : -idMCP;
    setId( id1, id2, id3New, -id3New);
    setColAcol( 1, 0, 0, 1, 0, 0, 0, 0);
    if (id1 < 0) swapColAcol();
  }

  virtual string name()     const {return "q qbar -> gamma* -> chi chibar";}
  virtual int    code()     const {return 9901;}
  virtual string inFlux()   const {return "qqbarSame";}
  // Make the phase space use the chi mass rather than treating it massless.
  virtual int    id3Mass()  const {return idMCP;}
  virtual int    id4Mass()  const {return idMCP;}
  // Declare the s-channel Z so that PhaseSpace adds a Breit-Wigner sampling
  // term at sHat = mZ^2.  Without this the tau sampling is only 1/tau + 1/tau^2
  // and the pole is badly under-sampled whenever it sits well above the
  // 2 m_chi threshold (i.e. at low mass).  Pythia drops the term by itself once
  // mHatMin rises above the pole.
  virtual int    resonanceA() const {return useZ ? 23 : 0;}

private:
  int    idMCP;
  double qMCP;
  bool   useZ;
  double m2Res, GamMRat, thetaWRat, vChi, gamProp, intProp, resProp, angFac;
};

//==========================================================================

int main(int argc, char* argv[]) {

  if (argc < 6) {
    cout << "usage: mcp_dy <mass> <nEvent> <seed> <eCM> <outfile> [extra.cmnd]" << endl;
    return 1;
  }
  double mChi   = atof(argv[1]);
  int    nEvent = atoi(argv[2]);
  int    seed   = atoi(argv[3]);
  double eCM    = atof(argv[4]);
  string outName = argv[5];

  Pythia pythia;

  // Photon only by default; set "MCP:includeZ = on" to add the Z with the
  // hypercharge-kinetic-mixing couplings described above.
  pythia.settings.addFlag("MCP:includeZ", false);

  SigmaProcessPtr sigmaMCP = make_shared<Sigma2qqbar2mcp>(ID_MCP, 1.0);
  pythia.addSigmaPtr(sigmaMCP);

  // Define chi: spin 1/2, colour singlet, stable.  Electric charge is set to
  // zero in the particle data (the millicharge lives in the matrix element)
  // so that Pythia does not shower QED radiation off it; that radiation is
  // suppressed by eps^2 and is negligible for a millicharged particle.
  char buf[512];
  snprintf(buf, sizeof(buf),
    "%d:new = chi chibar 2 0 0 %.6f 0.0 0.0 0.0 0.0", ID_MCP, mChi);
  pythia.readString(buf);
  snprintf(buf, sizeof(buf), "%d:mayDecay = off", ID_MCP);
  pythia.readString(buf);
  snprintf(buf, sizeof(buf), "%d:isResonance = off", ID_MCP);
  pythia.readString(buf);

  pythia.readString("Beams:idA = 2212");
  pythia.readString("Beams:idB = 2212");
  snprintf(buf, sizeof(buf), "Beams:eCM = %.1f", eCM);
  pythia.readString(buf);

  // Drell-Yan-like scale choice: mu_R = mu_F = mHat.
  pythia.readString("SigmaProcess:renormScale2 = 4");
  pythia.readString("SigmaProcess:factorScale2 = 4");

  // Showers on (they set the chi pT and slightly broaden eta); MPI and
  // hadronisation are irrelevant for the chi kinematics and cost time.
  pythia.readString("PartonLevel:ISR = on");
  pythia.readString("PartonLevel:FSR = on");
  pythia.readString("PartonLevel:MPI = off");
  pythia.readString("HadronLevel:all = off");

  pythia.readString("Random:setSeed = on");
  snprintf(buf, sizeof(buf), "Random:seed = %d", seed);
  pythia.readString(buf);

  pythia.readString("Print:quiet = on");
  pythia.readString("Init:showChangedSettings = on");
  pythia.readString("Next:numberCount = 0");

  // Optional extra settings file (used by the validation runs).
  if (argc > 6) pythia.readFile(argv[6]);

  if (!pythia.init()) return 1;

  // Full-eta histogram of the single-particle spectrum, for the
  // "normalise over all eta" convention.
  const int    NFULL = 240;
  const double ETAFULL = 6.0;
  vector<long> fullHist(NFULL, 0);

  // Validation histograms: cos(theta*) of chi w.r.t. the incoming quark, in
  // the hard-process CM, split into a near-threshold and a relativistic
  // sHat window; plus the mHat spectrum.
  const int NCOS = 20;
  vector<long> cosNear(NCOS, 0), cosFast(NCOS, 0);
  const int NMHAT = 100;
  const double MHATMAX = 8. * mChi;
  vector<long> mhatHist(NMHAT, 0);

  long nChiAll = 0;
  vector<float> etaOut, bgOut, ptOut;
  // Per-event chi pair, for per-event (rather than per-particle) acceptance.
  vector<float> pairEta1, pairEta2;

  int nAcc = 0;
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (!pythia.next()) continue;
    ++nAcc;

    // Hard-process validation quantities (before showering).
    {
      Event& pr = pythia.process;
      int iq = (pr[3].id() > 0) ? 3 : 4;         // the incoming quark
      int ic = (pr[5].id() > 0) ? 5 : 6;         // the outgoing chi
      Vec4 pPair = pr[5].p() + pr[6].p();
      double mHatNow = pPair.mCalc();
      int ibm = int(mHatNow / MHATMAX * NMHAT);
      if (ibm >= 0 && ibm < NMHAT) ++mhatHist[ibm];
      Vec4 pChi = pr[ic].p(), pQ = pr[iq].p();
      pChi.bstback(pPair);
      pQ.bstback(pPair);
      double cosStar = costheta(pChi, pQ);
      int ibc = int((cosStar + 1.) / 2. * NCOS);
      if (ibc >= 0 && ibc < NCOS) {
        double beta = sqrtpos(1. - 4. * mChi * mChi / (mHatNow * mHatNow));
        if (beta < 0.5)      ++cosNear[ibc];
        else if (beta > 0.9) ++cosFast[ibc];
      }
    }
    vector<double> etaPair;
    for (int i = 0; i < pythia.event.size(); ++i) {
      if (!pythia.event[i].isFinal()) continue;
      if (pythia.event[i].idAbs() != ID_MCP) continue;
      Particle& p = pythia.event[i];
      etaPair.push_back(p.eta());
      double eta = p.eta();
      double bg  = p.pAbs() / mChi;
      ++nChiAll;
      int ib = int((eta + ETAFULL) / (2. * ETAFULL) * NFULL);
      if (ib >= 0 && ib < NFULL) ++fullHist[ib];
      if (fabs(eta) < 1.5) {
        etaOut.push_back(eta);
        bgOut.push_back(bg);
        ptOut.push_back(p.pT());
      }
    }
    // Keep the pair only if at least one chi can have non-zero A(eta).
    if (etaPair.size() == 2
        && (fabs(etaPair[0]) < 1.2 || fabs(etaPair[1]) < 1.2)) {
      pairEta1.push_back(etaPair[0]);
      pairEta2.push_back(etaPair[1]);
    }
  }

  pythia.stat();

  FILE* f = fopen(outName.c_str(), "w");
  fprintf(f, "# mass_GeV %.6f\n", mChi);
  fprintf(f, "# eCM_GeV %.1f\n", eCM);
  fprintf(f, "# nEventRequested %d\n", nEvent);
  fprintf(f, "# nEventAccepted %d\n", nAcc);
  fprintf(f, "# sigma_mb %.8e\n", pythia.info.sigmaGen());
  fprintf(f, "# sigmaErr_mb %.8e\n", pythia.info.sigmaErr());
  fprintf(f, "# nChiAll %ld\n", nChiAll);
  fprintf(f, "# fullHist_etaMax %.1f nbins %d\n", ETAFULL, NFULL);
  fprintf(f, "FULLHIST");
  for (int i = 0; i < NFULL; ++i) fprintf(f, " %ld", fullHist[i]);
  fprintf(f, "\n");
  fprintf(f, "# cos/mhat validation histograms\n");
  fprintf(f, "COSNEAR");
  for (int i = 0; i < NCOS; ++i) fprintf(f, " %ld", cosNear[i]);
  fprintf(f, "\nCOSFAST");
  for (int i = 0; i < NCOS; ++i) fprintf(f, " %ld", cosFast[i]);
  fprintf(f, "\n# mhatHist_max %.3f nbins %d\n", MHATMAX, NMHAT);
  fprintf(f, "MHATHIST");
  for (int i = 0; i < NMHAT; ++i) fprintf(f, " %ld", mhatHist[i]);
  fprintf(f, "\n");
  fprintf(f, "# nPairsKept %zu  (events with >=1 chi at |eta|<1.2)\n",
          pairEta1.size());
  for (size_t i = 0; i < pairEta1.size(); ++i)
    fprintf(f, "PAIR %.6g %.6g\n", pairEta1[i], pairEta2[i]);
  fprintf(f, "# columns: eta betagamma pt (|eta|<1.5)\n");
  for (size_t i = 0; i < etaOut.size(); ++i)
    fprintf(f, "%.6g %.6g %.6g\n", etaOut[i], bgOut[i], ptOut[i]);
  fclose(f);

  cout << "wrote " << outName << " with " << etaOut.size()
       << " chi in |eta|<1.5 out of " << nChiAll << endl;
  return 0;
}
