#ifndef SimDataFormats_CaloParticle_h
#define SimDataFormats_CaloParticle_h

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include "SimDataFormats/CaloAnalysis/interface/SimClusterFwd.h"
#include "SimDataFormats/CaloHit/interface/PCaloHit.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include <vector>

class SimTrack;
class EncodedEventId;

class CaloParticle {
  friend std::ostream &operator<<(std::ostream &s, CaloParticle const &tp);

public:
  typedef int Charge;                                       ///< electric charge type
  typedef math::XYZTLorentzVectorD LorentzVector;           ///< Lorentz vector
  typedef math::PtEtaPhiMLorentzVector PolarLorentzVector;  ///< Lorentz vector
  typedef math::XYZPointD Point;                            ///< point in the space
  typedef math::XYZVectorD Vector;                          ///< point in the space

  /// reference to reco::GenParticle
  typedef reco::GenParticleRefVector::iterator genp_iterator;
  typedef std::vector<SimTrack>::const_iterator g4t_iterator;
  typedef SimClusterRefVector::iterator sc_iterator;

  CaloParticle();

  CaloParticle(const SimTrack &simtrk);
  CaloParticle(EncodedEventId eventID, uint32_t particleID);  // for PU

  // destructor
  ~CaloParticle();

  /** @brief PDG ID.
   *
   * Returns the PDG ID of the first associated gen particle. If there are no
   * gen particles associated then it returns type() from the first SimTrack. */
  int pdgId() const {
    if (genParticles_.empty())
      return g4Tracks_[0].type();
    else
      return (*genParticles_.begin())->pdgId();
  }

  /** @brief Signal source, crossing number.
   *
   * Note this is taken from the first SimTrack only, but there shouldn't be any
   * SimTracks from different crossings in the CaloParticle. */
  EncodedEventId eventId() const { return event_; }

  uint64_t particleId() const { return particleId_; }

  // Setters for G4 and reco::GenParticle
  void addGenParticle(const reco::GenParticleRef &ref) { genParticles_.push_back(ref); }
  void addSimCluster(const SimClusterRef &ref) { simClusters_.push_back(ref); }
  void addG4Track(const SimTrack &t) { g4Tracks_.push_back(t); }
  /// iterators
  genp_iterator genParticle_begin() const { return genParticles_.begin(); }
  genp_iterator genParticle_end() const { return genParticles_.end(); }
  g4t_iterator g4Track_begin() const { return g4Tracks_.begin(); }
  g4t_iterator g4Track_end() const { return g4Tracks_.end(); }
  sc_iterator simCluster_begin() const { return simClusters_.begin(); }
  sc_iterator simCluster_end() const { return simClusters_.end(); }

  // Getters for Embd and Sim Tracks
  const reco::GenParticleRefVector &genParticles() const { return genParticles_; }
  const SimClusterRefVector &simClusters() const { return simClusters_; }
  // Only for clusters from the signal vertex
  const std::vector<SimTrack> &g4Tracks() const { return g4Tracks_; }

  void clearSimClusters() { simClusters_.clear(); }

  /// @brief Electric charge. Note this is taken from the first SimTrack only.
  float charge() const { return g4Tracks_[0].charge(); }
  /// Gives charge in unit of quark charge (should be 3 times "charge()")
  int threeCharge() const { return lrintf(3.f * charge()); }

  /// @brief Four-momentum Lorentz vector. Note this is taken from the first
  /// SimTrack only.
  const math::XYZTLorentzVectorF &p4() const { return theMomentum_; }

  /// @brief spatial momentum vector
  math::XYZVectorF momentum() const { return p4().Vect(); }

  /// @brief Vector to boost to the particle centre of mass frame.
  math::XYZVectorF boostToCM() const { return p4().BoostToCM(); }

  /// @brief Magnitude of momentum vector. Note this is taken from the first
  /// SimTrack only.
  float p() const { return p4().P(); }

  /// @brief Energy. Note this is taken from the first SimTrack only.
  float energy() const { return p4().E(); }

  /// @brief Transverse energy. Note this is taken from the first SimTrack only.
  float et() const { return p4().Et(); }

  /// @brief Mass. Note this is taken from the first SimTrack only.
  float mass() const { return p4().M(); }

  /// @brief Mass squared. Note this is taken from the first SimTrack only.
  float massSqr() const { return pow(mass(), 2); }

  /// @brief Transverse mass. Note this is taken from the first SimTrack only.
  float mt() const { return p4().Mt(); }

  /// @brief Transverse mass squared. Note this is taken from the first SimTrack
  /// only.
  float mtSqr() const { return p4().Mt2(); }

  /// @brief x coordinate of momentum vector. Note this is taken from the first
  /// SimTrack only.
  float px() const { return p4().Px(); }

  /// @brief y coordinate of momentum vector. Note this is taken from the first
  /// SimTrack only.
  float py() const { return p4().Py(); }

  /// @brief z coordinate of momentum vector. Note this is taken from the first
  /// SimTrack only.
  float pz() const { return p4().Pz(); }

  /// @brief Transverse momentum. Note this is taken from the first SimTrack
  /// only.
  float pt() const { return p4().Pt(); }

  /// @brief Momentum azimuthal angle. Note this is taken from the first
  /// SimTrack only.
  float phi() const { return p4().Phi(); }

  /// @brief Momentum polar angle. Note this is taken from the first SimTrack
  /// only.
  float theta() const { return p4().Theta(); }

  /// @brief Momentum pseudorapidity. Note this is taken from the simtrack
  /// before the calorimeter
  float eta() const { return p4().Eta(); }

  /// @brief Rapidity. Note this is taken from the simtrack before the
  /// calorimeter
  float rapidity() const { return p4().Rapidity(); }

  /// @brief Same as rapidity().
  float y() const { return rapidity(); }

  /** @brief Status word.
   *
   * Returns status() from the first gen particle, or -99 if there are no gen
   * particles attached. */
  int status() const { return genParticles_.empty() ? -99 : (*genParticles_[0]).status(); }

  static const unsigned int longLivedTag;  ///< long lived flag

  /// is long lived?
  bool longLived() const { return status() & longLivedTag; }

  /** @brief Gives the total number of SimHits, in the cluster */
  int numberOfSimHits() const { return nsimhits_; }

  /** @brief Gives the total number of SimHits, in the cluster */
  int numberOfRecHits() const { return hits_.size(); }

  /** @brief returns the time in ns of the caloparticle */
  float simTime() const { return time_; }

  /** @brief add rechit with fraction */
  void addRecHitAndFraction(uint32_t hit, float fraction) {
    hits_.emplace_back(hit);
    fractions_.emplace_back(fraction);
  }

  /** @brief Returns list of rechit IDs and fractions for this CaloParticle */
  std::vector<std::pair<uint32_t, float>> hits_and_fractions() const {
    std::vector<std::pair<uint32_t, float>> result;
    for (size_t i = 0; i < hits_.size(); ++i) {
      result.emplace_back(hits_[i], fractions_[i]);
    }
    return result;
  }

  /** @brief returns the accumulated sim energy in the cluster */
  float simEnergy() const { return simhit_energy_; }

  /** @brief add simhit's energy to cluster */
  void addSimHit(const PCaloHit &hit) {
    simhit_energy_ += hit.energy();
    ++nsimhits_;
  }

  /** @brief add vertex time to the caloparticle */
  void setSimTime(const float time) { time_ = time; }

protected:
  uint64_t nsimhits_{0};
  EncodedEventId event_;

  uint32_t particleId_{0};
  float simhit_energy_{0.f};
  float time_{std::numeric_limits<float>::lowest()};
  std::vector<uint32_t> hits_;
  std::vector<float> fractions_;

  math::XYZTLorentzVectorF theMomentum_;

  /// references to G4 and reco::GenParticle tracks
  std::vector<SimTrack> g4Tracks_;
  reco::GenParticleRefVector genParticles_;

  SimClusterRefVector simClusters_;
};

#endif  // SimDataFormats_CaloParticle_H
