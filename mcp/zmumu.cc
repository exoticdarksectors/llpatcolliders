// Independent cross-check: muon pseudorapidity from standard Z production.
#include "Pythia8/Pythia.h"
using namespace Pythia8;
int main(int argc, char* argv[]) {
  Pythia pythia;
  pythia.readString("Beams:idA = 2212");
  pythia.readString("Beams:idB = 2212");
  pythia.readString("Beams:eCM = 13600.");
  // argv[2] = "22" selects the full 2 -> 2 matrix element, which carries the
  // correct decay angular distribution; default is the 2 -> 1 resonance.
  bool use22 = (argc > 2 && string(argv[2]) == "22");
  pythia.readString(use22 ? "WeakSingleBoson:ffbar2ffbar(s:gmZ) = on"
                          : "WeakSingleBoson:ffbar2gmZ = on");
  pythia.readString("23:onMode = off");
  pythia.readString("23:onIfAny = 13");
  pythia.readString("PhaseSpace:mHatMin = 60.");
  pythia.readString("PhaseSpace:mHatMax = 120.");
  pythia.readString("SigmaProcess:renormScale1 = 1");
  pythia.readString("PartonLevel:MPI = off");
  pythia.readString("HadronLevel:all = off");
  pythia.readString("Print:quiet = on");
  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = 4242");
  if (!pythia.init()) return 1;
  int n = (argc > 1) ? atoi(argv[1]) : 200000;
  long nMu = 0, nMuEta1 = 0, nEv = 0, nAny = 0, nBoth = 0, nZy1 = 0;
  for (int i = 0; i < n; ++i) {
    if (!pythia.next()) continue;
    ++nEv;
    vector<double> eta; Vec4 pZ;
    for (int j = 0; j < pythia.event.size(); ++j)
      if (pythia.event[j].isFinal() && pythia.event[j].idAbs() == 13) {
        eta.push_back(pythia.event[j].eta());
        pZ += pythia.event[j].p();
      }
    if (eta.size() != 2) continue;
    if (fabs(pZ.rap()) < 1.) ++nZy1;
    for (double e : eta) { ++nMu; if (fabs(e) < 1.) ++nMuEta1; }
    bool a = fabs(eta[0]) < 1., b = fabs(eta[1]) < 1.;
    if (a || b) ++nAny;
    if (a && b) ++nBoth;
  }
  printf("Z -> mu mu (60 < mHat < 120 GeV), %ld events\n", nEv);
  printf("  per-muon  f(|eta|<1) = %.4f\n", double(nMuEta1)/nMu);
  printf("  P(>=1 mu at |eta|<1) = %.4f\n", double(nAny)/nEv);
  printf("  P(both  at |eta|<1) = %.4f\n", double(nBoth)/nEv);
  printf("  f(|y_Z| < 1)         = %.4f   [the boson itself]\n",
         double(nZy1)/nEv);
  return 0;
}
