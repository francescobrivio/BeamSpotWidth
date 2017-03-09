
from CRABClient.UserUtilities import config
config = config()

#config.General.requestName  = '2016B_Fills'
config.General.requestName  = 'Fill_5199'
#config.General.requestName  = 'Run_274959_vtx_trks'

config.General.workArea     = 'beamSpot'

config.JobType.pluginName   = 'Analysis'
config.JobType.psetName     = '/afs/cern.ch/work/f/fbrivio/beamSpot/CMSSW_8_0_10_patch1/src/dxy_phi/dxy_phi/python/dxy_phi_cfg.py'
config.JobType.outputFiles  = ['tracksFile.root']

config.Data.inputDataset    = '/ZeroBias/Run2016G-TkAlMinBias-PromptReco-v1/ALCARECO'
config.Data.inputDBS        = 'global'
#config.Data.lumiMask        = 'Cert_246908-251883_13TeV_PromptReco_Collisions15_JSON_MuonPhys_v4.txt'
config.Data.lumiMask	    = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/DCSOnly/json_DCSONLY.txt'
config.Data.splitting       = 'LumiBased'
config.Data.unitsPerJob     = 30

#config.Data.runRange 	    = '274420-274422' #Fill 4988
#config.Data.runRange 	    = '274440-274443' #Fill 4990
#config.Data.runRange 	    = '274954-274959' #Fill 5005
#config.Data.runRange 	    = '274966-274971' #Fill 5013
#config.Data.runRange 	    = '274998-275001' #Fill 5017
#config.Data.runRange 	    = '275059-275074' #Fill 5020
#config.Data.runRange 	    = '275124-275125' #Fill 5021
#config.Data.runRange 	    = '275282-275293' #Fill 5024
#config.Data.runRange 	    = '274959' #Run 274959
#config.Data.runRange 	    = '275311' #Run 275311

config.Data.runRange		= '278820-278822' #Fill 5199



config.Data.outLFNDirBase = '/store/user/fbrivio/BeamSpot/'
config.Data.publication = False
config.Site.storageSite = 'T3_IT_MIB' #'T2_CH_CERN'
