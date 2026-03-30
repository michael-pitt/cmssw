import FWCore.ParameterSet.Config as cms

def customiseExclusive(process):
    """
    Customisation function to add MPI/Exclusive track footprint 
    tables to the standard NanoAOD sequence.
    """
    
    # Define your producer and attach it to the process
    process.objectTrackCount = cms.EDProducer("ObjectTrackCountProducer",
        vtxs = cms.InputTag("offlineSlimmedPrimaryVertices"),
        jets = cms.InputTag("slimmedJets"),
        muons = cms.InputTag("slimmedMuons"),
        electrons = cms.InputTag("slimmedElectrons"),
        pfcands = cms.InputTag("packedPFCandidates"),
        dzCut = cms.double(0.1),
        etaCut = cms.double(2.1),
        drJet = cms.double(0.4),
        drLep = cms.double(0.3)
    )

    # Add it to the standard NanoAOD Common Table Task
    if hasattr(process, 'nanoTableTaskCommon'):
        process.nanoTableTaskCommon.add(process.objectTrackCount)
    else:
        # Fallback if running a very specific custom sequence
        process.objectTrackCountTask = cms.Task(process.objectTrackCount)
        process.schedule.associate(process.objectTrackCountTask)

    return process