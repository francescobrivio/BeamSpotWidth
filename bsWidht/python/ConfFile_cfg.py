import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = '80X_dataRun2_2016LegacyRepro_Candidate_v1'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) ) # -1 for all events

process.TFileService = cms.Service("TFileService",
								   fileName = cms.string("tracksFile.root"),
								   closeFileFast = cms.untracked.bool(False)
								   )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        #'file:myfile.root'
		'/store/data/Run2016G/ZeroBias/ALCARECO/TkAlMinBias-PromptReco-v1/000/278/820/00000/16619F74-7364-E611-84C6-FA163EF92E66.root',
		#'/store/data/Run2016G/ZeroBias/ALCARECO/TkAlMinBias-PromptReco-v1/000/278/820/00000/16A528EF-7164-E611-AF07-02163E0133EA.root',
		#'/store/data/Run2016G/ZeroBias/ALCARECO/TkAlMinBias-PromptReco-v1/000/278/820/00000/26CD48AE-6464-E611-9E44-02163E0144D7.root',
		#'/store/data/Run2016G/ZeroBias/ALCARECO/TkAlMinBias-PromptReco-v1/000/278/820/00000/44AA15D0-6664-E611-83D9-02163E012695.root',
		#'/store/data/Run2016G/ZeroBias/ALCARECO/TkAlMinBias-PromptReco-v1/000/278/820/00000/4A72AAEE-6C64-E611-9E56-FA163E9FAA3A.root',
		#'/store/data/Run2016G/ZeroBias/ALCARECO/TkAlMinBias-PromptReco-v1/000/278/820/00000/609614C7-6D64-E611-BCFD-02163E012991.root',
		#'/store/data/Run2016G/ZeroBias/ALCARECO/TkAlMinBias-PromptReco-v1/000/278/820/00000/8826C85C-6764-E611-ADEB-02163E01446C.root',
		#'/store/data/Run2016G/ZeroBias/ALCARECO/TkAlMinBias-PromptReco-v1/000/278/820/00000/90020AD8-7064-E611-BAAF-02163E014683.root',
		#'/store/data/Run2016G/ZeroBias/ALCARECO/TkAlMinBias-PromptReco-v1/000/278/820/00000/AC24F5AA-6F64-E611-AB62-FA163E205596.root',
		#'/store/data/Run2016G/ZeroBias/ALCARECO/TkAlMinBias-PromptReco-v1/000/278/820/00000/AEC41B44-6864-E611-9E86-FA163E8B3345.root',
		#'/store/data/Run2016G/ZeroBias/ALCARECO/TkAlMinBias-PromptReco-v1/000/278/820/00000/B4A4B2E5-C164-E611-A291-FA163E1851FB.root',
		#'/store/data/Run2016G/ZeroBias/ALCARECO/TkAlMinBias-PromptReco-v1/000/278/820/00000/BAD9998F-6964-E611-9BD6-02163E011BD0.root',
		#'/store/data/Run2016G/ZeroBias/ALCARECO/TkAlMinBias-PromptReco-v1/000/278/820/00000/BAF3292F-6B64-E611-A9B1-02163E0141CD.root',
		#'/store/data/Run2016G/ZeroBias/ALCARECO/TkAlMinBias-PromptReco-v1/000/278/820/00000/D89E6BA4-6D64-E611-B3FF-02163E01382C.root',
		#'/store/data/Run2016G/ZeroBias/ALCARECO/TkAlMinBias-PromptReco-v1/000/278/820/00000/DC1C5D1F-6164-E611-B7A5-FA163EA0486D.root',
		#'/store/data/Run2016G/ZeroBias/ALCARECO/TkAlMinBias-PromptReco-v1/000/278/820/00000/DC4CA2A2-7264-E611-A27C-02163E013927.root',
		#'/store/data/Run2016G/ZeroBias/ALCARECO/TkAlMinBias-PromptReco-v1/000/278/820/00000/E2D1751A-7964-E611-BF7F-02163E011979.root',
		#'/store/data/Run2016G/ZeroBias/ALCARECO/TkAlMinBias-PromptReco-v1/000/278/820/00000/E867F108-6464-E611-A6BA-FA163E3C23FF.root',
		#'/store/data/Run2016G/ZeroBias/ALCARECO/TkAlMinBias-PromptReco-v1/000/278/820/00000/F220C7B0-6E64-E611-BECA-02163E011C19.root',
		#'/store/data/Run2016G/ZeroBias/ALCARECO/TkAlMinBias-PromptReco-v1/000/278/820/00000/F8C8F41F-6A64-E611-9FA4-02163E01288C.root',
		#'/store/data/Run2016G/ZeroBias/ALCARECO/TkAlMinBias-PromptReco-v1/000/278/820/00000/FE716F77-6664-E611-AFBE-02163E01386E.root',
		#'/store/data/Run2016G/ZeroBias/ALCARECO/TkAlMinBias-PromptReco-v1/000/278/821/00000/0C24573A-6A64-E611-997D-02163E01374E.root',
		#'/store/data/Run2016G/ZeroBias/ALCARECO/TkAlMinBias-PromptReco-v1/000/278/822/00000/06012C5C-7964-E611-9B34-02163E0145A6.root',
		#'/store/data/Run2016G/ZeroBias/ALCARECO/TkAlMinBias-PromptReco-v1/000/278/822/00000/0846B306-6F64-E611-A280-02163E011C84.root',
		#'/store/data/Run2016G/ZeroBias/ALCARECO/TkAlMinBias-PromptReco-v1/000/278/822/00000/0EC56827-7364-E611-B00A-02163E014643.root',
		#'/store/data/Run2016G/ZeroBias/ALCARECO/TkAlMinBias-PromptReco-v1/000/278/822/00000/18B36F5F-7864-E611-8F89-FA163E8E8695.root',
		#'/store/data/Run2016G/ZeroBias/ALCARECO/TkAlMinBias-PromptReco-v1/000/278/822/00000/28753889-7264-E611-8420-02163E014599.root',
		#'/store/data/Run2016G/ZeroBias/ALCARECO/TkAlMinBias-PromptReco-v1/000/278/822/00000/5896ECCB-7064-E611-8D80-02163E0134DF.root',
		#'/store/data/Run2016G/ZeroBias/ALCARECO/TkAlMinBias-PromptReco-v1/000/278/822/00000/96FB494E-7564-E611-9C71-02163E0138BB.root',
		#'/store/data/Run2016G/ZeroBias/ALCARECO/TkAlMinBias-PromptReco-v1/000/278/822/00000/9E5F6081-7464-E611-BE31-02163E01351B.root',
		#'/store/data/Run2016G/ZeroBias/ALCARECO/TkAlMinBias-PromptReco-v1/000/278/822/00000/B0EF5C3F-8B64-E611-9BEF-02163E0134D3.root',
		#'/store/data/Run2016G/ZeroBias/ALCARECO/TkAlMinBias-PromptReco-v1/000/278/822/00000/B2683426-9864-E611-85BD-02163E013707.root',
		#'/store/data/Run2016G/ZeroBias/ALCARECO/TkAlMinBias-PromptReco-v1/000/278/822/00000/C0E8A2A5-7664-E611-BB9E-02163E0125FD.root',
		#'/store/data/Run2016G/ZeroBias/ALCARECO/TkAlMinBias-PromptReco-v1/000/278/822/00000/C881DD0C-7064-E611-AA6E-FA163E951746.root',
		#'/store/data/Run2016G/ZeroBias/ALCARECO/TkAlMinBias-PromptReco-v1/000/278/822/00000/DA0A2910-7C64-E611-9E40-FA163E7FD727.root',
		#'/store/data/Run2016G/ZeroBias/ALCARECO/TkAlMinBias-PromptReco-v1/000/278/822/00000/F6900F72-6E64-E611-A88C-02163E011882.root', 
		#'/store/data/Run2016G/ZeroBias/ALCARECO/TkAlMinBias-PromptReco-v1/000/278/822/00000/FC4B08C4-7164-E611-87E1-FA163E87BDC8.root',
    )
)

process.demo = cms.EDAnalyzer('bsWidht',
							  vtxCollection   = cms.InputTag("offlinePrimaryVertices"),
							  beamSpot        = cms.InputTag("offlineBeamSpot"),
							  trackCollection = cms.InputTag("ALCARECOTkAlMinBias")
							  )


process.p = cms.Path(process.demo)
