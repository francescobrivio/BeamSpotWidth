
from CRABClient.UserUtilities import config
config = config()

config.General.requestName  = 'Fill_5199_bsWidth_good'

config.General.workArea     = 'beamSpot'

config.JobType.pluginName   = 'Analysis'
config.JobType.psetName     = '/afs/cern.ch/work/f/fbrivio/beamSpot/CMSSW_8_0_18/src/bsWidth/bsWidht/python/ConfFile_cfg.py'
config.JobType.outputFiles  = ['tracksFile.root']

config.Data.inputDataset    = '/ZeroBias/Run2016G-TkAlMinBias-PromptReco-v1/ALCARECO'
config.Data.inputDBS        = 'global'
#config.Data.lumiMask        = 'Cert_246908-251883_13TeV_PromptReco_Collisions15_JSON_MuonPhys_v4.txt'
config.Data.lumiMask	    = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/DCSOnly/json_DCSONLY.txt'
config.Data.splitting       = 'LumiBased'
config.Data.unitsPerJob     = 30


config.Data.runRange		= '278820-278822' #Fill 5199



config.Data.outLFNDirBase = '/store/user/fbrivio/BeamSpot/'
config.Data.publication = False
config.Site.storageSite = 'T3_IT_MIB' #'T2_CH_CERN'
