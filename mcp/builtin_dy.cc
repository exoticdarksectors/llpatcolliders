// Reference: Pythia's built-in pure-gamma* Drell-Yan to muons, for
// cross-section validation of the custom mcp process.
#include "Pythia8/Pythia.h"
using namespace Pythia8;
int main(int argc, char* argv[]) {
  Pythia pythia;
  pythia.readString("Beams:idA = 2212");
  pythia.readString("Beams:idB = 2212");
  pythia.readString("Beams:eCM = 14000.");
  pythia.readString("WeakSingleBoson:ffbar2ffbar(s:gmZ) = on");
  pythia.readString("WeakZ0:gmZmode = 1");          // pure photon
  pythia.readString("23:onMode = off");
  pythia.readString("23:onIfAny = 13");             // gamma* -> mu mu only
  pythia.readString("SigmaProcess:renormScale2 = 4");
  pythia.readString("SigmaProcess:factorScale2 = 4");
  pythia.readString("PartonLevel:all = off");
  pythia.readString("Print:quiet = on");
  if (argc > 1) pythia.readString(string("PhaseSpace:mHatMin = ") + argv[1]);
  if (argc > 3) pythia.readFile(argv[3]);
  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = 991");
  if (!pythia.init()) return 1;
  int n = (argc > 2) ? atoi(argv[2]) : 20000;
  for (int i = 0; i < n; ++i) pythia.next();
  printf("SIGMA_mb %.6e +- %.3e   (mHatMin=%s)\n",
         pythia.info.sigmaGen(), pythia.info.sigmaErr(), argc>1?argv[1]:"def");
  return 0;
}
