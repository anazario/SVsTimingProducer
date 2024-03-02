#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

enum LepType {kW, kZ, kTau, kConversion, kLight, kHeavy, kSusy, kUnmatched};

bool isSignalGenElectron(const reco::GenParticle &genElectron);
std::vector<int> MomIDs(const reco::GenParticle &genElectron);
LepType ClassifyGenElectron(const std::vector<int> &motherIDs);
LepType AssignLeptonMomType(const int motherID);
void PrintGenElectronInfo(const reco::GenParticle &genElectron);
void PrintMother(const LepType &momType);

// Implementations

bool isSignalGenElectron(const reco::GenParticle &genElectron) {

  LepType momType = ClassifyGenElectron(MomIDs(genElectron));
  return (momType == kZ || momType == kSusy);

}

std::vector<int> MomIDs(const reco::GenParticle &genElectron) {

  auto mother = genElectron.mother();

  std::vector<int> motherIDs;
  while(mother->pt() > 0) {
    const int motherID = mother->pdgId();

    mother = mother->mother();

    if(motherID == mother->pdgId())
      continue;

    motherIDs.push_back(motherID);
  }

  return motherIDs;
}

LepType ClassifyGenElectron(const std::vector<int> &motherIDs) {

  LepType momType = kUnmatched;
  for(auto const& id : motherIDs) {
    momType = AssignLeptonMomType(id);

    if(momType != kUnmatched)
      break;
  }
  return momType;
}

LepType AssignLeptonMomType(const int motherID) {

  LepType type = kUnmatched;

  if(abs(motherID) == 24)
    type = kW;
  else if(motherID == 23 || motherID == 36)
    type = kZ;
  else if(abs(motherID) == 15)
    type = kTau;
  else if((abs(motherID%1000) > 100 && abs(motherID%1000) < 400)
          || (abs(motherID%1000) > 1000 && abs(motherID%1000) < 4000)
          || (abs(motherID) > 0 && abs(motherID) < 4)
          || motherID == 21)
    type = kLight;
  else if((abs(motherID%1000) > 400 && abs(motherID%1000) < 600)
          || (abs(motherID%1000) > 4000 && abs(motherID%1000) < 6000)
          || (abs(motherID) > 3 && abs(motherID) < 7))
    type = kHeavy;
  else if(motherID == 22)
    type = kConversion;
  else if(motherID == 1000022)
    type = kSusy;

  return type;
}

void PrintGenElectronInfo(const reco::GenParticle &genElectron) {

  const int genCharge(genElectron.charge());
  const float gnX(genElectron.vx() );
  const float gnY(genElectron.vy() );
  const float gnZ(genElectron.vz() );
  const float distXY(sqrt(gnX*gnX + gnY*gnY));
  const float genPt(genElectron.pt());
  const float genEta(genElectron.eta());
  const float genPhi(genElectron.phi());
  const float genEnergy(genElectron.energy());

  std::cout << "  vertex position: (" << gnX << ", " << gnY <<  ", " << gnZ << ")" << std::endl;
  std::cout << "  Total transverse displacement: " << distXY << std::endl;
  std::cout << "  pT: " << genPt << ", eta: " << genEta << ", phi: " << genPhi << ", energy: " << genEnergy << std::endl;
  std::cout << "  gen charge: " << genCharge <<  std::endl;

  for(auto const& id : MomIDs(genElectron))
    std::cout << "  motherID: " << id << std::endl;
  PrintMother(ClassifyGenElectron(MomIDs(genElectron)));
}

void PrintMother(const LepType &momType) {

  std::cout << "mother: ";
  if (momType == kW) std::cout << "W boson" << std::endl;
  else if (momType == kZ) std::cout << "Z boson" << std::endl;
  else if (momType == kTau)  std::cout << "tau lepton" << std::endl;
  else if (momType == kConversion) std::cout << "photon (conversion)" << std::endl;
  else if (momType == kLight) std::cout << "light quark" << std::endl;
  else if (momType == kHeavy) std::cout << "heavy quark" << std::endl;
  else if (momType == kSusy) std::cout << "prompt from Susy particle" << std::endl;
  else if (momType == kUnmatched) std::cout << "unmatched (this shouldn't happen!)" << std::endl;

}

