import FWCore.ParameterSet.Config as cms
from SimPPS.DirectSimProducer.profile_base_cff import profile_base as _base

from CalibPPS.ESProducers.ctppsOpticalFunctions_non_DB_cff import optics_2025 as _optics

_base_2025 = _base.clone(
    ctppsOpticalFunctions = _base.ctppsOpticalFunctions.clone(
        opticalFunctions = _optics.opticalFunctions,
        scoringPlanes = _optics.scoringPlanes,
    ),
    ctppsDirectSimuData = _base.ctppsDirectSimuData.clone(
        empiricalAperture45 = "1.e3*([xi] - 0.20)",
        empiricalAperture56 = "1.e3*([xi] - 0.20)"
    )
)

profile_2025_default = _base_2025.clone(
    L_int = 1.,
    ctppsLHCInfo = _base_2025.ctppsLHCInfo.clone(
        # We can update this to 2025/2026 histograms when they are ready, using 2021 for now is fine
        xangleBetaStarHistogramObject = "2021/h2_betaStar_vs_xangle"
    ),
    ctppsRPAlignmentCorrectionsDataXML = _base_2025.ctppsRPAlignmentCorrectionsDataXML.clone(
        MisalignedFiles = ["Validation/CTPPS/alignment/alignment_2022.xml"],
        RealFiles = ["Validation/CTPPS/alignment/alignment_2022.xml"]
    ),
    ctppsDirectSimuData = _base_2025.ctppsDirectSimuData.clone(
        timeResolutionDiamonds45 = "0.200",
        timeResolutionDiamonds56 = "0.200"
    )
)

