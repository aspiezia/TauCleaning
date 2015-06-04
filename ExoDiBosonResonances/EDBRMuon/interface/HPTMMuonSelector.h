#ifndef HPTMMUONSELECTOR_63FF6402B78B442EB96CDD86A0595875_H
#define HPTMMUONSELECTOR_63FF6402B78B442EB96CDD86A0595875_H
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include <bitset>


/// High Pt Muon 
namespace hptm {

  /// Two tipes of muon ID: the "old" high pt muon id,
  /// and the custom one we developed.
  enum MuonIDType {HIGHPT_OLD, TRACKER};

  /// Auxiliary class to select muons in the context of high-pt muon
  /// analysis. Could be extended to other analyses, of course.
  class MuonSelector {

  private: 

  public:
    MuonSelector(){}
    ~MuonSelector(){} 
    
    /// Typedef for three-dimensional point
    typedef math::XYZPoint Point;
    
    /// Returns a bitset where each bit indicates if the muon passed
    /// that cut or not.  Could be extended to become a full-fledged
    /// self-descriptive structure. For the time being, we simply
    /// describe it here. A given bit is true if the condition is true.
    /// Notice that for some bits, the condition is true depending
    /// on the chosen muon ID.
    ///
    /// [0]: Global muon
    /// [1]: Tracker muon
    /// [2]: Number of muon chamber hits cut
    /// [3]: Number of matched stations cut
    /// [4]: dB cut
    /// [5]: dZ cut
    /// [6]: Number of pixel hits cut
    /// [7]: Number of tracker layers cut
    std::bitset<8> muonBits(const reco::Muon&, const Point&, MuonIDType);
  
    /// Checks if a given muon passes the muon ID or not.
    bool checkMuonID(const reco::Muon&, const Point&,MuonIDType);
  };
}
#endif
