import FWCore.ParameterSet.Config as cms

hiEvtAnalyzer = cms.EDAnalyzer('HiEvtAnalyzer',
   CentralitySrc    = cms.InputTag("hiCentrality"),
   CentralityBinSrc = cms.InputTag("centralityBin","HFtowers"),
   pfCandidateSrc   = cms.InputTag('packedPFCandidates'),
   EvtPlane         = cms.InputTag("hiEvtPlane"),
   EvtPlaneFlat     = cms.InputTag("hiEvtPlaneFlat",""),
   HiMC             = cms.InputTag("heavyIon"),
   Vertex           = cms.InputTag("offlineSlimmedPrimaryVertices"),
   HFfilters = cms.InputTag("hiHFfilters","hiHFfilters"),
   doCentrality     = cms.bool(True),
   doEvtPlane       = cms.bool(True),
   doEvtPlaneFlat   = cms.bool(True),
   doVertex         = cms.bool(True),
   doMC             = cms.bool(True),
   doHiMC           = cms.bool(True),
   useHepMC         = cms.bool(False),
   doHFfilters      = cms.bool(True),
   evtPlaneLevel    = cms.int32(0),

   minHFEnergy     = cms.untracked.double(4.0),
   minAbsEtaHF     = cms.untracked.double(3.0),
   maxAbsEtaHF     = cms.untracked.double(5.2),
)

#from HeavyIonsAnalysis.EventAnalysis.hievtanalyzer_data_cfi import *

#hiEvtAnalyzer.HFfilters = cms.InputTag("hiHFfilters","hiHFfilters","DQM"),

#hiEvtAnalyzer.doMC   = True
#hiEvtAnalyzer.doHiMC = True
