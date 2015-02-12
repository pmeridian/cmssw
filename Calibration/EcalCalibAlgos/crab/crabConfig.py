from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'DYJetsToLL_Phys14DR_zeeSkim'
config.General.transferLogs = False

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'skimZee.py'
#config.JobType.inputFiles = ['./ecalcovariances_popcon_data.db','./ecalcovariances_popcon_mc.db','./ecalnoisecorrelation_popcon_data.db','./ecalnoisecorrelation_popcon_mc.db','./ecaltemplates_popcon_data.db','./ecaltemplates_popcon_mc.db']

config.section_("Data")
config.Data.inputDataset = '/DYJetsToLL_M-50_13TeV-madgraph-pythia8/Phys14DR-PU20bx25_PHYS14_25_V1-v1/AODSIM'
config.Data.inputDBS = 'global'
#config.Data.splitting = 'LumiBased'
#config.Data.unitsPerJob = 25 # 200
#config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions12/8TeV/Prompt/Cert_190456-208686_8TeV_PromptReco_Collisions12_JSON.txt'
#config.Data.lumiMask = 'Cert_190456-208686_8TeV_PromptReco_Collisions12_JSON.txt' # if you downloaded the file in the working directory
#config.Data.runRange = '211760-211831' # '193093-194075'
config.Data.splitting = "FileBased"
config.Data.unitsPerJob = 1
config.Data.publication = False
#config.Data.publishDBS = 'phys03'
#config.Data.publishDataName = 'ECAL-multifit-50ns'
config.Data.outLFN = "/store/group/dpg_ecal/alca_ecalcalib/crovelli/zeeSkim"

config.section_("Site")
config.Site.storageSite = "T2_CH_CERN"
