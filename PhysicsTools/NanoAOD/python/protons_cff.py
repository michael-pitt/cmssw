import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *
from PhysicsTools.NanoAOD.nano_eras_cff import *
from RecoCTPPS.ProtonReconstruction.ppsFilteredProtonProducer_cfi import *

from Validation.CTPPS.simu_config.base_cff import *

# proton reco settings
#from Validation.CTPPS.simu_config.year_2017_preTS2_cff import ctppsDirectProtonSimulation as ppsDirectSim2017_preTS2
#from Validation.CTPPS.simu_config.year_2017_postTS2_cff import ctppsDirectProtonSimulation as ppsDirectSim2017_postTS2
#from Validation.CTPPS.simu_config.year_2018_cff import ctppsDirectProtonSimulation as ppsDirectSim2018
#from Validation.CTPPS.simu_config.year_2017_cff import reco_local as reco_local_2017
#from Validation.CTPPS.simu_config.year_2018_cff import reco_local as reco_local_2018

def SetAperture(ctppsDirectProtonSimulation, year):
  if year=='2016_preTS2':
    ctppsDirectProtonSimulation.empiricalAperture45_xi0_int = 0.111
    ctppsDirectProtonSimulation.empiricalAperture45_xi0_slp = 0.0
    ctppsDirectProtonSimulation.empiricalAperture45_a_int = 127.0
    ctppsDirectProtonSimulation.empiricalAperture45_a_slp = -0.000
    ctppsDirectProtonSimulation.empiricalAperture56_xi0_int = 0.138
    ctppsDirectProtonSimulation.empiricalAperture56_xi0_slp = 0.0
    ctppsDirectProtonSimulation.empiricalAperture56_a_int = 191.6
    ctppsDirectProtonSimulation.empiricalAperture56_a_slp = 0.0
  elif year=='2016_preTS2':
    ctppsDirectProtonSimulation.empiricalAperture45_xi0_int = 0.104
    ctppsDirectProtonSimulation.empiricalAperture45_xi0_slp = 0.0
    ctppsDirectProtonSimulation.empiricalAperture45_a_int = 116.4
    ctppsDirectProtonSimulation.empiricalAperture45_a_slp = 0.0
    ctppsDirectProtonSimulation.empiricalAperture56_xi0_int = 0.110
    ctppsDirectProtonSimulation.empiricalAperture56_xi0_slp = 0.0
    ctppsDirectProtonSimulation.empiricalAperture56_a_int = 150.0
    ctppsDirectProtonSimulation.empiricalAperture56_a_slp = 0.0
  elif year=='2017_preTS2':
    ctppsDirectProtonSimulation.empiricalAperture45_xi0_int = 0.066
    ctppsDirectProtonSimulation.empiricalAperture45_xi0_slp = 3.536E-04
    ctppsDirectProtonSimulation.empiricalAperture45_a_int = 47.7
    ctppsDirectProtonSimulation.empiricalAperture45_a_slp = 0.447
    ctppsDirectProtonSimulation.empiricalAperture56_xi0_int = 0.062
    ctppsDirectProtonSimulation.empiricalAperture56_xi0_slp = 5.956E-04
    ctppsDirectProtonSimulation.empiricalAperture56_a_int = -31.9
    ctppsDirectProtonSimulation.empiricalAperture56_a_slp = 1.323
  elif year=='2017_postTS2':
    ctppsDirectProtonSimulation.empiricalAperture45_xi0_int = 0.073
    ctppsDirectProtonSimulation.empiricalAperture45_xi0_slp = 4.107E-04
    ctppsDirectProtonSimulation.empiricalAperture45_a_int = 39.0
    ctppsDirectProtonSimulation.empiricalAperture45_a_slp = 0.768
    ctppsDirectProtonSimulation.empiricalAperture56_xi0_int = 0.067
    ctppsDirectProtonSimulation.empiricalAperture56_xi0_slp = 6.868E-04
    ctppsDirectProtonSimulation.empiricalAperture56_a_int = -50.2
    ctppsDirectProtonSimulation.empiricalAperture56_a_slp = 1.740
  elif year=='2018':
    ctppsDirectProtonSimulation.empiricalAperture45_xi0_int = 0.079
    ctppsDirectProtonSimulation.empiricalAperture45_xi0_slp = 4.211E-04
    ctppsDirectProtonSimulation.empiricalAperture45_a_int = 42.8
    ctppsDirectProtonSimulation.empiricalAperture45_a_slp = 0.669
    ctppsDirectProtonSimulation.empiricalAperture56_xi0_int = 0.074
    ctppsDirectProtonSimulation.empiricalAperture56_xi0_slp = 6.604E-04
    ctppsDirectProtonSimulation.empiricalAperture56_a_int = -22.7
    ctppsDirectProtonSimulation.empiricalAperture56_a_slp = 1.600
    
  return ctppsDirectProtonSimulation

#optics
from CalibPPS.ESProducers.ctppsOpticalFunctionsESSource_cfi import *
from CalibPPS.ESProducers.ctppsInterpolatedOpticalFunctionsESSource_cfi import *
from Validation.CTPPS.simu_config.year_2016_preTS2_cff import config_2016_preTS2
from Validation.CTPPS.simu_config.year_2016_postTS2_cff import config_2016_postTS2
from Validation.CTPPS.simu_config.year_2017_cff import config_2017
from Validation.CTPPS.simu_config.year_2018_cff import config_2018
  
#alignment
from CalibPPS.ESProducers.ctppsRPAlignmentCorrectionsDataESSourceXML_cfi import *
from Validation.CTPPS.simu_config.year_2016_preTS2_cff import alignmentFile as alignmentFile_2016_preTS2
from Validation.CTPPS.simu_config.year_2016_postTS2_cff import alignmentFile as alignmentFile_2016_postTS2
from Validation.CTPPS.simu_config.year_2017_preTS2_cff import alignmentFile as alignmentFile_2017_preTS2
from Validation.CTPPS.simu_config.year_2017_postTS2_cff import alignmentFile as alignmentFile_2017_postTS2
from Validation.CTPPS.simu_config.year_2018_cff import alignmentFile as alignmentFile_2018

from Geometry.VeryForwardGeometry.geometryRPFromDD_2017_cfi import XMLIdealGeometryESSource_CTPPS, ctppsGeometryESModule

singleRPProtons = True

filteredProtons = ppsFilteredProtonProducer.clone(
    protons_single_rp = cms.PSet(
        include = cms.bool(singleRPProtons)
    )
)

protonTable = cms.EDProducer("ProtonProducer",
                             tagRecoProtonsMulti  = cms.InputTag("filteredProtons", "multiRP"),
                             tagTrackLite         = cms.InputTag("ctppsLocalTrackLiteProducer"),
                             storeSingleRPProtons = cms.bool(singleRPProtons)
)
protonTable.tagRecoProtonsSingle = cms.InputTag("filteredProtons" if singleRPProtons else "ctppsProtons","singleRP")


multiRPTable = cms.EDProducer("SimpleProtonTrackFlatTableProducer",
    src = cms.InputTag("filteredProtons","multiRP"),
    cut = cms.string(""),
    name = cms.string("Proton_multiRP"),
    doc  = cms.string("bon"),
    singleton = cms.bool(False),
    extension = cms.bool(False),
    skipNonExistingSrc = cms.bool(True),
    variables = cms.PSet(
        xi = Var("xi",float,doc="xi or dp/p",precision=12),
        thetaX = Var("thetaX",float,doc="theta x",precision=13),
        thetaY = Var("thetaY",float,doc="theta y",precision=13),
        t = Var("t",float,doc="Mandelstam variable t",precision=13),
        time = Var("time()",float,doc="time",precision=16),
        timeUnc = Var("timeError",float,doc="time uncertainty",precision=13),
    ),
    externalVariables = cms.PSet(
        arm = ExtVar("protonTable:arm",int,doc="0 = sector45, 1 = sector56"),
    ),
)

singleRPTable = cms.EDProducer("SimpleProtonTrackFlatTableProducer",
    src = cms.InputTag("filteredProtons","singleRP"),
    cut = cms.string(""),
    name = cms.string("Proton_singleRP"),
    doc  = cms.string("bon"),
    singleton = cms.bool(False),
    extension = cms.bool(False),
    skipNonExistingSrc = cms.bool(True),
    variables = cms.PSet(
        xi = Var("xi",float,doc="xi or dp/p",precision=12),
        thetaY = Var("thetaY",float,doc="th y",precision=10),
    ),
    externalVariables = cms.PSet(
        decRPId = ExtVar("protonTable:protonRPId",int,doc="Detector ID",precision=8), 
    ),
)

protonTables = cms.Sequence(  
    filteredProtons
    +protonTable
    +multiRPTable
)

# setup proton simulation
beamDivergenceVtxGenerator.src = cms.InputTag("")
beamDivergenceVtxGenerator.srcGenParticle = cms.VInputTag(
    #cms.InputTag("genPUProtons",""), # should be used in private production
    #cms.InputTag("genPUProtons","genPUProtons"),  # works with step2_premix modifier
    cms.InputTag("prunedGenParticles"), # when ~premix_stage2 signal protons proporate to genPUProtons
)

ctppsLocalTrackLiteProducer.includeStrips = False
ctppsLocalTrackLiteProducer.includePixels = False 
reco_local = cms.Sequence(ctppsLocalTrackLiteProducer)

#setup postTS2 2017 configuration

#ctppsDirectProtonSimulation = ppsDirectSim2017_postTS2
#reco_local = reco_local_2017
ctppsDirectProtonSimulation.useEmpiricalApertures = True
ctppsDirectProtonSimulation=SetAperture(ctppsDirectProtonSimulation,'2017_postTS2')
ctppsLocalTrackLiteProducer.includeStrips = True
ctppsLocalTrackLiteProducer.includePixels = True 
reco_local.insert(0,totemRPUVPatternFinder*totemRPLocalTrackFitter*ctppsPixelLocalTracks)

ctppsOpticalFunctionsESSource.configuration.append(config_2017)
ctppsRPAlignmentCorrectionsDataESSourceXML.MisalignedFiles = [alignmentFile_2017_postTS2]
ctppsRPAlignmentCorrectionsDataESSourceXML.RealFiles = [alignmentFile_2017_postTS2]

# https://twiki.cern.ch/twiki/bin/viewauth/CMS/TaggedProtonsDirectSimulation#Configuring_the_beam_divergence
# do not apply vertex smearing again
ctppsBeamParametersESSource.vtxStddevX = 0
ctppsBeamParametersESSource.vtxStddevY = 0
ctppsBeamParametersESSource.vtxStddevZ = 0
  
#undo CMS vertex shift
ctppsBeamParametersESSource.vtxOffsetX45 = +0.2475 * 1E-1
ctppsBeamParametersESSource.vtxOffsetY45 = -0.6924 * 1E-1
ctppsBeamParametersESSource.vtxOffsetZ45 = -8.1100 * 1E-1

# crossing angle distribution
from CalibPPS.ESProducers.ctppsLHCInfoRandomXangleESSource_cfi import *
ctppsLHCInfoRandomXangleESSource.xangleHistogramObject = cms.string("hxang")
ctppsLHCInfoRandomXangleESSource.xangleHistogramFile = cms.string("/afs/cern.ch/work/m/mpitt/PROPOG/NANOAOD/DirectSimulation/CrossingAngles2017H.root")
ctppsLHCInfoRandomXangleESSource.beamEnergy = cms.double(6500.)
ctppsLHCInfoRandomXangleESSource.betaStar = cms.double(0.4)
esPreferLHCInfo = cms.ESPrefer("CTPPSLHCInfoRandomXangleESSource", "ctppsLHCInfoRandomXangleESSource")

# for multiRP fit, set if you want to use x* and y* as free parameters or set them to zero
ctppsProtons.fitVtxY = True
#if false then ndof=1 and chi2 values will be big (filteredProton container will be empty)        
      
protonTablesMC = cms.Sequence(  
    beamDivergenceVtxGenerator
    +ctppsDirectProtonSimulation
    +reco_local
    +ctppsProtons
    +filteredProtons
    +protonTable
    +multiRPTable
)

if singleRPProtons: protonTables.insert(protonTables.index(multiRPTable),singleRPTable)

(run2_nanoAOD_92X | run2_miniAOD_80XLegacy | run2_nanoAOD_94X2016 | run2_nanoAOD_94X2016 | \
    run2_nanoAOD_94XMiniAODv1 | run2_nanoAOD_94XMiniAODv2 | \
    run2_nanoAOD_102Xv1 | ( run2_nanoAOD_106Xv1 & ~run2_nanoAOD_devel) ).toReplaceWith(protonTables, cms.Sequence())


