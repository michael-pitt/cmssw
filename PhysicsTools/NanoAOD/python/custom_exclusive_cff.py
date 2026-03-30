import FWCore.ParameterSet.Config as cms

def customiseExclusive(process):
    """
    Customisation function to add MPI/Exclusive track footprint 
    tables to the standard NanoAOD sequence.
    """
    
    # Define your producer and attach it to the process
    process.objectTrackCountTable = cms.EDProducer("ObjectTrackCountProducer",
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

    if hasattr(process, 'nanoAOD_step'):
        process.nanoAOD_step += process.objectTrackCountTable
    else:
        # Fallback
        process.objectTrackCountTask = cms.Task(process.objectTrackCountTable)
        process.schedule.associate(process.objectTrackCountTask)
    
    for output in ["NANOAODoutput", "NANOAODSIMoutput"]:
        if hasattr(process, output):
            out_module = getattr(process, output)
            out_module.outputCommands.append("keep *_objectTrackCountTable_*_*")

    return process
