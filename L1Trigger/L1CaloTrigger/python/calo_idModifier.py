import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Modifier_phase2_hgcalV11_cff import phase2_hgcalV11

option_egid_select = {'new_tight':1,'new_loose':2 ,'old_tight':3 ,'old_loose':4 }
option_pfid_select = {'':0, 'WithPUNewNoPionID':1,'WithPUNewPionOldID':2}

def ModifyEGID(process, option):
  if(option==1):
    phase2_hgcalV11.toModify(process.hgcalBackEndLayer2Producer.ProcessorParameters.C3d_parameters.histoMax_C3d_clustering_parameters.EGIdentification,
        Weights=cms.vstring(
        'L1Trigger/L1THGCal/data/egid_electron_200PU_vs_neutrino_200PU_full_Histomaxvardr_loweta.xml',
        'L1Trigger/L1THGCal/data/egid_electron_200PU_vs_neutrino_200PU_full_Histomaxvardr_higheta.xml'
        ),
        WorkingPoints=cms.vdouble( 0.8808619, 0.9958049) #tightID
    )
  elif(option==2):
    phase2_hgcalV11.toModify(process.hgcalBackEndLayer2Producer.ProcessorParameters.C3d_parameters.histoMax_C3d_clustering_parameters.EGIdentification,
        Weights=cms.vstring(
        'L1Trigger/L1THGCal/data/egid_electron_200PU_vs_neutrino_200PU_full_Histomaxvardr_loweta.xml',
        'L1Trigger/L1THGCal/data/egid_electron_200PU_vs_neutrino_200PU_full_Histomaxvardr_higheta.xml'
        ),
        WorkingPoints=cms.vdouble(-0.6905887, 0.9806464) #looseID
    )
  elif(option==3):  # see https://cmssdt.cern.ch/lxr/source/L1Trigger/L1THGCal/python/egammaIdentification.py
    phase2_hgcalV11.toModify(process.hgcalBackEndLayer2Producer.ProcessorParameters.C3d_parameters.histoMax_C3d_clustering_parameters.EGIdentification,
        Weights=cms.vstring(
        'L1Trigger/L1THGCal/data/egamma_id_histomax_3151_loweta_v0.xml',
        'L1Trigger/L1THGCal/data/egamma_id_histomax_3151_higheta_v0.xml'
        ),
        WorkingPoints=cms.vdouble( 0.8292287, 0.9932326) #tightID
    )
  elif(option==4):
    phase2_hgcalV11.toModify(process.hgcalBackEndLayer2Producer.ProcessorParameters.C3d_parameters.histoMax_C3d_clustering_parameters.EGIdentification,
        Weights=cms.vstring(
        'L1Trigger/L1THGCal/data/egamma_id_histomax_3151_loweta_v0.xml',
        'L1Trigger/L1THGCal/data/egamma_id_histomax_3151_higheta_v0.xml'
        ),
        WorkingPoints=cms.vdouble(-0.7099538, 0.9611762) #looseID
    )
  else:
    print('WARNING: Called ModifyEGID with wrong option, keep it default (old tight configuration)')

#from https://github.com/VourMa/FastPUPPI/blob/8e031c44f70558f96f7c1a8155eb0e7bb1a8fa51/NtupleProducer/python/runPerformanceNTuple.py#L661
def ModifyPFID(process, option):
  process.load("L1Trigger.Phase2L1ParticleFlow.l1ParticleFlow_cff")
  process.runPF_step = cms.Path(process.l1ParticleFlow)
  process.extraPFStuff = cms.Task()
  nameID = ''

  if(option):
    #variablesPU_New weightsPU_New WPPU_New from 
    #https://github.com/VourMa/FastPUPPI/blob/8e031c44f70558f96f7c1a8155eb0e7bb1a8fa51/NtupleProducer/python/runPerformanceNTuple.py#L769
    variablesPU = cms.VPSet(
      cms.PSet(name = cms.string("eMaxOverE"), value = cms.string("eMax()/energy()")),
      cms.PSet(name = cms.string("sigmaRRTot"), value = cms.string("sigmaRRMean()")),
      cms.PSet(name = cms.string("zBarycenter"), value = cms.string("zBarycenter()")),
      cms.PSet(name = cms.string("layer90percent"), value = cms.string("layer90percent()")),
      cms.PSet(name = cms.string("triggerCells67percent"), value = cms.string("triggerCells67percent()")),
    )
    weightsPU = "L1Trigger/L1THGCal/data/Train_DoublePhoton_FlatPt-1To100_PU200_SinglePion_PT0to200_PU200_vs_PU_110X_v1_5var/weights/TMVA_BDT.weights.xml"
    WPPU = "-0.01"
     
    #variablesPion_Old weightsPion_OldOld WPPion_OldOld from 
    #https://github.com/VourMa/FastPUPPI/blob/8e031c44f70558f96f7c1a8155eb0e7bb1a8fa51/NtupleProducer/python/runPerformanceNTuple.py#L810
    variablesPion = cms.VPSet(
        cms.PSet(name = cms.string("fabs(eta)"), value = cms.string("abs(eta())")),
        cms.PSet(name = cms.string("coreShowerLength"), value = cms.string("coreShowerLength()")),
        cms.PSet(name = cms.string("maxLayer"), value = cms.string("maxLayer()")),
        cms.PSet(name = cms.string("hOverE"), value = cms.string("hOverE()")),
        cms.PSet(name = cms.string("sigmaZZ"), value = cms.string("sigmaZZ()")),
    )
    weightsPion = "L1Trigger/Phase2L1ParticleFlow/data/hgcal_egID/Photon_vs_Pion_BDTweights.xml.gz"
    WPPion = "0.01"
    
    inv_map = {v: k for k, v in option_pfid_select.iteritems()}
    nameID = inv_map[option]
      
    #Clone the pfClustersFromHGC3DClusters and modify them to include whatever ID you want
    pfClustersFromHGC3DClustersID = process.pfClustersFromHGC3DClusters.clone(
      emVsPUID = cms.PSet(
        isPUFilter = cms.bool(True),
        preselection = cms.string(""),
        method = cms.string("BDT"),
        variables = variablesPU,
        weightsFile = cms.string(weightsPU),
        wp = cms.string(WPPU)
      ),
      emVsPionID = cms.PSet(
        isPUFilter = cms.bool(False),
        preselection = cms.string(""),
        method = cms.string("BDT"),
        variables = variablesPion,
        weightsFile = cms.string(weightsPion),
        wp = cms.string(WPPion)
      ),
    )
    setattr(process, 'pfClustersFromHGC3DClusters'+nameID, pfClustersFromHGC3DClustersID)
    
    #Clone the l1pfProducerHGCal and modify it to use the modified version of the pfClustersFromHGC3DClusters that you created
    l1pfProducerHGCalID = process.l1pfProducerHGCal.clone(hadClusters = [ cms.InputTag("pfClustersFromHGC3DClusters"+nameID) ])
    setattr(process, 'l1pfProducerHGCal'+nameID, l1pfProducerHGCalID)
    
    #Clone and create a modified l1pfProducerHGCalNoTKID from the modified l1pfProducerHGCal
    l1pfProducerHGCalNoTKID = l1pfProducerHGCalID.clone(
        regions = cms.VPSet(
            cms.PSet(
                etaBoundaries = cms.vdouble(-3,-2.5),
                phiSlices = cms.uint32(1),
                etaExtra = cms.double(0.3),
                phiExtra = cms.double(0.0)
            ),
            cms.PSet(
                etaBoundaries = cms.vdouble(2.5,3),
                phiSlices = cms.uint32(1),
                etaExtra = cms.double(0.3),
                phiExtra = cms.double(0.0)
            ),
        )
    )  
    setattr(process, 'l1pfProducerHGCalNoTK'+nameID, l1pfProducerHGCalNoTKID)
    
    #Clone and modify the l1pfCandidates to use the modified versions of the producers created
    l1pfCandidatesHGCalID = process.l1pfCandidates.clone(
        pfProducers = cms.VInputTag(
            cms.InputTag("l1pfProducerBarrel"),
            cms.InputTag("l1pfProducerHGCal"+nameID),
            cms.InputTag("l1pfProducerHGCalNoTK"+nameID),
            cms.InputTag("l1pfProducerHF")
        ),
    )
    setattr(process, 'l1pfCandidatesHGCal'+nameID, l1pfCandidatesHGCalID)
    
    process.extraPFStuff.add(getattr(process,'pfClustersFromHGC3DClusters'+nameID),getattr(process,'l1pfProducerHGCal'+nameID),getattr(process,'l1pfProducerHGCalNoTK'+nameID),getattr(process,'l1pfCandidatesHGCal'+nameID))
  
  process.runPF_step.associate(process.extraPFStuff)
  print(process.runPF_step)
  
  

  