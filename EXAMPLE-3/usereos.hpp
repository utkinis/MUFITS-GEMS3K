#pragma once

#define IPMGEMPLUGIN
#define NOPARTICLEARRAY
#define NODEARRAYLEVEL
#include <GEMS3K/node.h>

#include <cstdint>

#define MF_API extern "C" __declspec(dllexport)

// Exported functions for USEREOS

MF_API void ReadConfigurationFile(char *filename, int *ierr, char *title);

MF_API void GetDimensions(int *nComponents, int *nPhaseMax, int *nAux);

MF_API void GetGlobalParameters(char *cmpNames, double *molWeights, char *phNames, char *auxNames,
                                char *auxUnits, int8_t *opt);

MF_API void PhaseEquilibrium(double *pres, double *temp, double *z, int *nPhase, double *props,
                             int8_t *phaseId, double *auxArray, int8_t *mode);

// Class GemsModule provides functionality related to coupling GEMS3K and MUFITS
class GemsModule {
public:
  // Initialize GEMS3K
  int init(const std::string &filename);

  // Pass fluid params to USEREOS arrays
  void getParams(char cmpNames[], double molWeights[], char phNames[], char auxNames[]);

  // Calculate chemical equilibrium and store results in arrays passed from USEREOS
  bool calculateEquilibrium(const double P, const double T, const double z[], int &nPhase,
                            double props[], int8_t phaseId[], double auxArray[], int8_t mode);

  // Number of independent components (ICs)
  int nComponents() const {
    return dch->nICb - 1; // Minus charge
  }

  // Number of dependent components (DCs)
  int nEndMembers() const { return dch->nDCb; }

  // Number of multicomponent fluid phases
  int nFluidPhases() const { return dch->nPSb; }

  // Maximum number of phases including multicomponent fluid phases and solid phases
  int maxPhases() const {
    if (dch->nPSb == dch->nPHb) {
      return dch->nPSb;
    } else {
      return dch->nPSb + 1; // One for the solid phase
    }
  }

private:
  void logFailure(const double P, const double T, const double z[], long status);

  // GEMS3K node
  TNode *node;

  // GEMS3K structures
  DATABR *dbr;
  DATACH *dch;

  // Molar weights of ICs
  std::vector<double> molarWeights;

  // Names of components
  std::vector<std::array<char, 3>> compNames;

  // Names of phases
  std::vector<std::array<char, 3>> phaseNames;

  // Names of auxiliary arrays for USEREOS
  std::vector<std::array<char, 8>> auxiliaryNames;
};
