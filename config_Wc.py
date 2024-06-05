import pocket_coffea
c = pocket_coffea.utils.configurator.Configurator
from pocket_coffea.utils.configurator import Configurator
from pocket_coffea.lib.cut_definition import Cut
from pocket_coffea.lib.cut_functions import get_nObj_min, get_HLTsel
from pocket_coffea.parameters.cuts import *
from pocket_coffea.parameters.histograms import *
import CoffeaBTVProcessor
from CoffeaBTVProcessor import CommBTVBaseProcessor

import CommonSelectors
from CommonSelectors import *

import cloudpickle
cloudpickle.register_pickle_by_value(CoffeaBTVProcessor)
cloudpickle.register_pickle_by_value(CommonSelectors)

import os
localdir = os.path.dirname(os.path.abspath(__file__))

# Loading default parameters
from pocket_coffea.parameters import defaults
default_parameters = defaults.get_default_parameters()
defaults.register_configuration_dir("config_dir", localdir+"/params")

parameters = defaults.merge_parameters_from_files(default_parameters,
                                                  f"{localdir}/params/object_preselection.yaml",
                                                  f"{localdir}/params/triggers.yaml",
                                                  update=True)



cfg = Configurator(
    parameters = parameters,
    datasets = {
        "jsons": [f"{localdir}/Run2UL2017_MC_VJets.json",
                  f"{localdir}/datasets/Run2UL2017_DATA.json"],
        
        
        "filter" : {
            "samples": ["WJetsToLNu_MLM","DATA_DoubleMuon",
                
            ],
            "samples_exclude" : [],
            "year": ["2017"]
            
        }
    },

    workflow = CommBTVBaseProcessor,

    #skim = [get_HLTsel(primaryDatasets=["SingleMuon","SingleEle"])],
    skim = [get_HLTsel(primaryDatasets=["SingleMuon"])],

    preselections = [ll_2j],
    categories = {
        "baseline_Wc": [passthrough],
        "Added cuts Wc": [req_mtw, ptratiobtv, req_dilepveto, mujetselbtv, reqmuon, jetselbtv, reqsoftmuon, req_dilepmass ],
    },

    weights = {
        "common": {
            "inclusive": ["signOf_genWeight","lumi","XS",
                          "pileup",
                          "sf_mu_id","sf_mu_iso",
                          ],
            "bycategory" : {
                
            }
        },
    },
    
    variations = {
        "weights": {
            "common": {
                "inclusive": [
                    "pileup",
                    "sf_mu_id", "sf_mu_iso",
                ]
            },
            "bysample": {
            }
        },
    },

    variables = {
        **lepton_hists(coll="LeptonGood", pos=0),
        **lepton_hists(coll="LeptonGood", pos=1),
        **lepton_hists(coll="MuonGood", pos=0),
        **lepton_hists(coll="SoftMuonGood", pos=0),
        **count_hist(name="nElectronGood", coll="ElectronGood",bins=5, start=0, stop=5),
        **count_hist(name="nMuonGood", coll="MuonGood",bins=5, start=0, stop=5),
        **count_hist(name="nJets", coll="JetGood",bins=8, start=0, stop=8),
        #**count_hist(name="nBJets", coll="BJetGood",bins=8, start=0, stop=8),
        **jet_hists(coll="JetGood", pos=0),
        **jet_hists(coll="JetGood", pos=1),
        
        "nJet": HistConf( [Axis(field="nJet", bins=15, start=0, stop=15, label=r"nJet direct from NanoAOD")] ),
        "nMuon": HistConf( [Axis(field="nMuonGood", bins=15, start=0, stop=15, label=r"nMuon direct from NanoAOD")] ),
        "nSoftMuon": HistConf( [Axis(field="nSoftMuonGood", bins=15, start=0, stop=15, label=r"nSoftMuon direct from NanoAOD")] ),
        "hl_ptratio" : HistConf( [Axis(field="hl_ptratio", bins=100, start=0, stop=8, label=r"$\frac{p_T(\mu)}{p_T(j)}$")] ),
        "soft_l_ptratio" : HistConf( [Axis(field="soft_l_ptratio", bins=100, start=0, stop=8, label=r"$\frac{p_T(\mu_{soft})}{p_T(j)}$")] ),
        "dr_soft_l_jet" : HistConf( [Axis(field="dr_soft_l_jet", bins=100, start=0, stop=8, label=r"$\Delta R_{\mu_{soft}, j}")] ),
        "dr_l_jet" : HistConf( [Axis(field="dr_l_jet", bins=100, start=0, stop=8, label=r"$\Delta R_{\mu, j}")] ),
        "dr_l_soft_l" : HistConf( [Axis(field="dr_l_soft_l", bins=100, start=0, stop=8, label=r"$\Delta R_{\mu_{soft}, \mu}")] ),
        "Z_pt" : HistConf( [Axis(field="Z_pt", bins=100, start=0, stop=200, label=r"$p_T{ll}")] ),
        "Z_mass" : HistConf( [Axis(field="Z_mass", bins=100, start=0, stop=150, label=r"$m_{ll}")] ),
        "Z_eta" : HistConf( [Axis(field="Z_eta", bins=100, start=0, stop=10, label=r"$\eta_{ll}")] ),
        "Z_phi" : HistConf( [Axis(field="Z_phi", bins=100, start=0, stop=200, label=r"$\phi_{ll}")] ),
        "W_pt" : HistConf( [Axis(field="W_pt", bins=100, start=0, stop=200, label=r"$p_T{l\nu}")] ),
        "W_mass" : HistConf( [Axis(field="W_mass", bins=100, start=0, stop=150, label=r"$m_{l\nu}")] ),
        "W_eta" : HistConf( [Axis(field="W_eta", bins=100, start=0, stop=10, label=r"$\eta_{l\nu}")] ),
        "W_phi" : HistConf( [Axis(field="W_phi", bins=100, start=0, stop=200, label=r"$\phi_{l\nu}")] ),
        "MET_pt" : HistConf( [Axis(field="MET_pt", bins=100, start=0, stop=200, label=r"$p_T{MET}")] ),
        "MET_phi" : HistConf( [Axis(field="MET_phi", bins=100, start=0, stop=200, label=r"$\phi_{MET}")] ),
    }
)



