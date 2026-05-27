import FWCore.ParameterSet.Config as cms

def customiseTightVertices(process):
    """
    Adds Tighter Vertex counters to the standard NanoAOD PV table.
    """
    process.tightVertexTable = cms.EDProducer("TightVertexTableProducer",
        pvSrc = cms.InputTag("offlineSlimmedPrimaryVertices"),
        goodPvCut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2")
    )
    
    if hasattr(process, 'nanoAOD_step'):
        process.nanoAOD_step += process.tightVertexTable
    else:
        process.tightVertexTask = cms.Task(process.tightVertexTable)
        process.schedule.associate(process.tightVertexTask)

    return process
    
def customiseExclusive(process):
    """
    Customisation function to add MPI/Exclusive track footprint 
    tables to the standard NanoAOD sequence.
    """
       
    # Define your producer and attach it to the process
    process.objectTrackCountTable = cms.EDProducer("ObjectTrackCountProducer",
        vtxs = cms.InputTag("offlineSlimmedPrimaryVertices"),
        jets = cms.InputTag("linkedObjects", "jets"),
        muons = cms.InputTag("linkedObjects", "muons"),
        electrons = cms.InputTag("linkedObjects", "electrons"),
        pfcands = cms.InputTag("packedPFCandidates"),
        dzOutCut = cms.double(1.0),
        dzCut = cms.double(0.1),
        etaCut = cms.double(2.1),
        drJet = cms.double(0.4),
        drEl = cms.double(0.2),
        drMu = cms.double(0.1),
    )

    if hasattr(process, 'nanoAOD_step'):
        process.nanoAOD_step += process.objectTrackCountTable
    else:
        # Fallback
        process.objectTrackCountTask = cms.Task(process.objectTrackCountTable)
        process.schedule.associate(process.objectTrackCountTask)
    
    return process

def customiseLowPUAnalysis(process):
    """
    The complete customisation chain for the Low-PU Diffractive Analysis.
    Sequentially applies the tight vertex counters and the footprint logic.
    """
    # Apply the vertex customisation
    process = customiseTightVertices(process)
    
    # Apply the exclusive footprint customisation
    process = customiseExclusive(process)
    
    return process

