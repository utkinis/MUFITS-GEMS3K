#include "usereos.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

// Various constants

// Module name appears in MUFITS logs
const char *moduleName = "GEMS Test module";
const int filenameLength = 256;
const int titleLength = 80;

// We model solid and fluid phases as slightly compressible
const double fluidCompressibility = 1e-7;
const double solidCompressibility = 1e-8;

// Used only if there is no solid phase present for numerical reasons
const double solidDensity = 3000;

// Solid viscosity is set to a very high value to model precipitation
// Actual values don't matter because in this example pressure distribution becomes stationary after
// first time step
const double fluidViscosity = 1e-3;
const double solidViscosity = 1e7;

// Pressure is fixed at 1 bar
const double systemPressure = 1e5;

// Instance of GemsModule
static GemsModule mod;

void ReadConfigurationFile(char *filename, int *ierr, char *title) {
  std::string cFilename(filename, filename + filenameLength);
  size_t spacePos = cFilename.find(' ');
  if (spacePos != std::string::npos) {
    cFilename[cFilename.find(' ')] = '\0';
  }
  memset(title, ' ', size_t(titleLength));
  memcpy(title, moduleName, strlen(moduleName));

  *ierr = mod.init(cFilename);
}

void GetDimensions(int *nComponents, int *nPhaseMax, int *nAux) {
  *nComponents = mod.nComponents();
  *nPhaseMax = mod.maxPhases();
  *nAux = mod.nEndMembers() + mod.nComponents();
}

void GetGlobalParameters(char *cmpNames, double *molWeights, char *phNames, char *auxNames,
                         char *auxUnits, int8_t *opt) {
  mod.getParams(cmpNames, molWeights, phNames, auxNames);
  char unit[] = "NODIM   ";
  for (int idx = 0; idx < mod.nEndMembers() + mod.nComponents(); ++idx) {
    memcpy(auxUnits + size_t(idx) * 8, unit, 8);
  }
  opt[0] = 0;
  opt[1] = 0;
}

void PhaseEquilibrium(double *pres, double *temp, double *z, int *nPhase, double *props,
                      int8_t *phaseId, double *auxArray, int8_t *mode) {
  mod.calculateEquilibrium(*pres, *temp, z, *nPhase, props, phaseId, auxArray, *mode);
}

int GemsModule::init(const std::string &filename) {
  // Create GEMS3K structures
  node = new TNode();
  dch = node->pCSD();
  dbr = node->pCNode();

  // Initialize GEMS3K
  auto status = node->GEM_init(filename.c_str());
  if (status != 0) {
    return status;
  }

  compNames.resize(size_t(nComponents()));
  molarWeights.resize(size_t(nComponents()));
  phaseNames.resize(size_t(maxPhases()));
  auxiliaryNames.resize(size_t(nEndMembers() + nComponents()));

  for (size_t idx = 0; idx < size_t(nComponents()); ++idx) {
    memset(&compNames[idx], ' ', 3);
    memcpy(&compNames[idx], dch->ICNL[dch->xic[idx]],
           std::min(strlen(dch->ICNL[dch->xic[idx]]), size_t(3)));
  }

  for (size_t idx = 0; idx < size_t(dch->nPSb); ++idx) {
    memset(&phaseNames[idx], ' ', 3);
    memcpy(&phaseNames[idx], dch->PHNL[dch->xph[idx]],
           std::min(strlen(dch->PHNL[dch->xph[idx]]), size_t(3)));
  }
  memcpy(&phaseNames[size_t(dch->nPSb)], "SOL", 3);

  for (size_t idx = 0; idx < size_t(nEndMembers()); ++idx) {
    std::string aux("#");
    aux.append(dch->DCNL[dch->xdc[idx]]);
    memset(&auxiliaryNames[idx], ' ', 8);
    memcpy(&auxiliaryNames[idx], aux.c_str(), std::min(aux.size(), size_t(8)));
  }

  for (size_t idx = 0; idx < size_t(nComponents()); ++idx) {
    std::string aux("#");
    aux.append(dch->ICNL[dch->xic[idx]]);
    memset(&auxiliaryNames[idx + nEndMembers()], ' ', 8);
    memcpy(&auxiliaryNames[idx + nEndMembers()], aux.c_str(), std::min(aux.size(), size_t(8)));
    molarWeights[idx] = dch->ICmm[dch->xic[idx]];
  }

  return 0;
}

void GemsModule::getParams(char cmpNames[], double molWeights[], char phNames[], char auxNames[]) {
  memcpy(cmpNames, compNames.data(), 3 * size_t(nComponents()));
  memcpy(phNames, phaseNames.data(), 3 * size_t(maxPhases()));
  if (auxNames != nullptr) {
    memcpy(auxNames, auxiliaryNames.data(), 8 * size_t(nEndMembers() + nComponents()));
  }
  std::copy(molarWeights.begin(), molarWeights.end(), molWeights);
}

bool GemsModule::calculateEquilibrium(const double P, const double T, const double z[], int &nPhase,
                                      double props[], int8_t phaseId[], double auxArray[], int8_t) {
  dbr->NodeStatusCH = NEED_GEM_AIA;
  dbr->IterDone = 0;
  dbr->P = systemPressure;
  dbr->TK = T;
  dbr->Vs = 0;
  dbr->Mi = 0;
  dbr->Vi = 0;

  // Specify system composition
  for (size_t idx = 0; idx < size_t(nComponents()); ++idx) {
    dbr->bIC[idx] = z[idx] / molarWeights[idx];
  }

  nPhase = maxPhases();

  // Run GEM algorithm
  auto status = node->GEM_run(false);
  if (status != OK_GEM_AIA) {
    logFailure(P, T, z, status);
  }

  size_t offset = size_t(6 + nComponents());
  std::fill(props, props + maxPhases() * offset, 0.0);

  double fluidVol = 0.0;
  double fluidMass = 0.0;
  double fluidEnthalpy = 0.0;
  for (int iPh = 0; iPh < dch->nPSb; ++iPh) {
    fluidVol += dbr->vPS[iPh];
    fluidMass += dbr->mPS[iPh];
    fluidEnthalpy += node->Ph_Enthalpy(iPh);
    for (size_t iC = 0; iC < size_t(nComponents()); ++iC) {
      props[6 + iC] += dbr->bPS[iPh * dch->nICb + iC] * molarWeights[iC];
    }
  }

  if (fluidMass > 0.0) {
    for (int iC = 0; iC < nComponents(); ++iC) {
      props[6 + iC] /= fluidMass;
    }
  }

  double fluidDensity = 0.0;
  if (fluidVol > 0.0) {
    fluidDensity = fluidMass / fluidVol;
  }

  if (fluidMass > 0.0) {
    fluidEnthalpy /= fluidMass;
  }

  double totVol = fluidVol;

  phaseId[0] = 1;

  if (mod.dch->nPHb > mod.dch->nPSb) {
    double solVol = 0.0;
    double solMass = 0.0;
    double solEnthalpy = 0.0;
    for (int iPh = dch->nPSb; iPh < mod.dch->nPHb; ++iPh) {
      solVol += mod.node->Ph_Volume(iPh);
      solMass += mod.node->Ph_Mass(iPh);
      solEnthalpy += mod.node->Ph_Enthalpy(iPh);
    }

    double solDensity = solidDensity;
    if (solMass > 0.0) {
      if (solVol > 0.0) {
        solDensity = solMass / solVol;
        totVol += solVol;
      }
      props[nFluidPhases() * offset + 0] = solDensity + solidCompressibility * (P - systemPressure);
      props[nFluidPhases() * offset + 1] = solEnthalpy / solMass;
      props[nFluidPhases() * offset + 2] = solidViscosity;
      props[nFluidPhases() * offset + 3] = solVol / totVol;
      for (size_t iC = 0; iC < size_t(mod.nComponents()); ++iC) {
        props[nFluidPhases() * offset + 6 + iC] = mod.dbr->bSP[iC] * mod.molarWeights[iC] / solMass;
      }
    } else {
      props[nFluidPhases() * offset + 0] = solDensity + solidCompressibility * (P - systemPressure);
      props[nFluidPhases() * offset + 2] = solidViscosity;
    }
    phaseId[nFluidPhases()] = int8_t(nFluidPhases() + 1);
  }

  props[0] = fluidDensity + fluidCompressibility * (P - systemPressure);
  props[1] = fluidEnthalpy;
  props[2] = fluidViscosity;
  props[3] = fluidVol / totVol;

  for (int idx = 0; idx < nEndMembers(); ++idx) {
    auxArray[idx] = mod.dbr->xDC[idx] / fluidMass;
  }

  for (int idx = 0; idx < nComponents(); ++idx) {
    auxArray[idx + nEndMembers()] = dbr->bPS[idx] / fluidMass;
  }

  return true;
}

void GemsModule::logFailure(const double P, const double T, const double z[], long status) {
  std::cerr << std::right << std::setfill('-') << std::setw(50);
  std::cerr << " !!! GEMS ERROR !!! " << std::setfill('-') << std::setw(30) << '-' << std::endl;
  std::cerr << std::setfill(' ') << std::left;
  std::cerr << setw(79);
  if (status == BAD_GEM_AIA) {
    std::cerr << "| ERROR CODE: BAD_GEM_AIA";
  } else if (status == ERR_GEM_AIA) {
    std::cerr << "| ERROR CODE: ERR_GEM_AIA";
  }
  std::cerr << "|\n";

  node->GEM_print_ipm("GEMS_fail.txt");
  std::cerr << std::setw(79) << "| GEMS REPORT SAVED TO FILE 'GEMS_fail.txt'" << '|' << std::endl;
  std::cerr << "| PRES: " << std::setprecision(8) << std::scientific << P << " PA";
  std::cerr << std::setw(55) << std::right << '|' << std::endl;
  std::cerr << "| TEMP: " << std::setprecision(8) << std::scientific << T << " K ";
  std::cerr << std::setw(55) << std::right << '|' << std::endl;
  std::cerr << "| MASS: " << std::setprecision(8) << std::scientific << dbr->Ms << " KG";
  std::cerr << std::setw(55) << std::right << '|' << std::endl;

  std::cerr << std::setw(79) << std::left << "| BULK COMPOSITION:" << '|' << std::endl;
  std::cerr << "|          ";
  for (size_t idx = 0; idx < size_t(nComponents()); ++idx) {
    char name[4] = {};
    memcpy(name, compNames[idx].data(), 3);
    std::cerr << std::setw(16) << std::left << name;
  }
  std::cerr << std::setw(5) << std::right << '|' << std::endl;

  std::cerr << "| MASS:  ";
  for (int idx = 0; idx < nComponents(); ++idx) {
    std::cerr << setw(16) << std::setprecision(8) << std::scientific << z[idx];
  }
  std::cerr << std::setw(7) << std::right << '|' << std::endl;

  std::cerr << "| MOLAR: ";
  for (int idx = 0; idx < nComponents(); ++idx) {
    std::cerr << setw(16) << std::setprecision(8) << std::scientific << dbr->bIC[idx];
  }
  std::cerr << std::setw(7) << std::right << '|' << std::endl;
  std::cerr << std::right << setw(80) << std::setfill('-') << '-' << std::endl << std::left;
}
