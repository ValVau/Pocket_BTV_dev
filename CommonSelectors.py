import awkward as ak
import numpy as np
from pocket_coffea.lib.cut_definition import Cut

def diLepton(events, params, year, sample, **kwargs):

    # Masks for same-flavor (SF) and opposite-sign (OS)
    EE = ((events.nMuonGood == 0) & (events.nElectronGood == 2))
    MuMu = ((events.nMuonGood == 2) & (events.nElectronGood == 0))

    OS = events.ll.charge == 0

    mask = (
        (MuMu | EE) & OS
        & (ak.firsts(events.LeptonGood.pt) > params["pt_leading_lep"])
        & (events.ll.mass > params["mll"]["low"])
        & (events.ll.mass < params["mll"]["high"])
    )
    
    # Pad None values with False
    return ak.where(ak.is_none(mask), False, mask)

def TwoMuons(events, **kwargs):
    mask = (events.nMuonGood >= 2)
    return ak.where(ak.is_none(mask), False, mask)

def TwoElectrons(events, **kwargs):
    mask = (events.nElectronGood >= 2)
    return ak.where(ak.is_none(mask), False, mask)

def TwoJets(events, **kwargs):
    mask = (events.nJetGood >= 2)
    return ak.where(ak.is_none(mask), False, mask)

def TwoLepTwoJets(events, params, **kwargs):
    mask = ( (events.nJetGood >= 2)
             & ( ( (params["lep_flav"]=="mu") & (events.nMuonGood>=2) ) |
                 ( (params["lep_flav"]=="el") & (events.nElectronGood>=2) )  |
                 ( (params["lep_flav"]=="both") & (events.nLeptonGood>=2) )
                )
             & (events.ll.pt > params["pt_dilep"])
             & (events.ll.mass > params["mll"]["low"])
             & (events.ll.mass < params["mll"]["high"])
            )
    return ak.where(ak.is_none(mask), False, mask)

def OneLeptonPlusMet(events, params, **kwargs):
    #mask = (events.nLeptonGood == 1 )
    mask = ( (events.nLeptonGood == 1 )
             & (ak.firsts(events.LeptonGood.pt) > params["pt_lep"])
             & (events.MET.pt > params["pt_met"])
            )
    return ak.where(ak.is_none(mask), False, mask)

def LepMetTwoJets(events, params, **kwargs):
    mask = ( (events.nLeptonGood == 1 )
             & (ak.firsts(events.LeptonGood.pt) > params["pt_lep"])
             & (events.MET.pt > params["pt_met"])
             & (events.nJetGood >= 2)
            )
    return ak.where(ak.is_none(mask), False, mask)

def MetTwoJetsNoLep(events, params, **kwargs):    
    mask = ( (events.nLeptonGood == 0 )
             & (events.MET.pt > params["pt_met"])
             & (events.nJetGood >= 2)
             & (ak.firsts(events.JetGood.pt) > params["pt_jet1"])
             & (ak.pad_none(events.JetGood.pt, 2, axis=1)[:,1] > params["pt_jet2"])
            )
    return ak.where(ak.is_none(mask), False, mask)

def WLNuTwoJets(events, params, **kwargs):

    fields = {
        "pt": events.MET.pt,
	"eta": ak.zeros_like(events.MET.pt),
        "phi": events.MET.phi,
        "mass": ak.zeros_like(events.MET.pt),
        "charge": ak.zeros_like(events.MET.pt),
    }

    METs = ak.zip(fields, with_name="PtEtaPhiMCandidate")
    LepPlusMet = METs + ak.firsts(events.LeptonGood)
    mask = ( (events.nJetGood >= 2)
             & ( ( (params["lep_flav"]=="mu") & (events.nMuonGood==1) ) |
                 ( (params["lep_flav"]=="el") & (events.nElectronGood==1) )  |
                 (params["lep_flav"]=="both") 
                )
             & (LepPlusMet.pt > params["pt_w"])
        )
    return ak.where(ak.is_none(mask), False, mask)

def ctag(events, params, **kwargs):
    #print(events.JetsCvsL.btagDeepFlavCvL[:, 0]>0.2)
    mask = (events.JetsCvsL.btagDeepFlavCvL[:,0]>0.2)
    return ak.where(ak.is_none(mask), False, mask)

def CvsLsorted(jets, ctag):
    return jets[ak.argsort(jets[ctag["tagger"]], axis=1, ascending=False)]

def DiJetPtCut(events, params, **kwargs):
    mask = (  (events.nJetGood >= 2)
              & (events.dijet.pt > params["pt_dijet"])
              #& (events.dijet_csort.pt > params["pt_dijet"])
            )
    return ak.where(ak.is_none(mask), False, mask)

def DeltaPhiJetMetCut(events, params, **kwargs):
    mask = ( (events['deltaPhi_jet1_MET'] > params["jet_met_dphi_cut"])
             & (events['deltaPhi_jet2_MET'] > params["jet_met_dphi_cut"])
            )
    return ak.where(ak.is_none(mask), False, mask)

def MuonSelBTV(events, params, **kwargs):
    mask = ak.num(events.JetGood[( ((events.JetGood.muonIdx1 != -1) | (events.JetGood.muonIdx2 != -1))
                & ((events.JetGood.muEF + events.JetGood.neEmEF) < 0.7)
            )].pt, axis=1) >= 1
    return ak.where(ak.is_none(mask), False, mask)

def JetSelBTV(events, params, **kwargs):
    mask = (ak.num(events.JetGood.pt) >= 1) & (ak.num(events.JetGood.pt) <= 3)
    return ak.where(ak.is_none(mask), False, mask)

def reqmu(events, params, **kwargs):
    mask = ak.count(events.MuonGood.pt, axis=1) == 1
    return ak.where(ak.is_none(mask), False, mask)

def reqsoftmu(events, params, **kwargs):
    mask = ak.count(events.SoftMuonGood.pt, axis=1) >= 1
    return ak.where(ak.is_none(mask), False, mask)

def JetMuPtratio(events, params, **kwargs):
    first_muon_pt = ak.firsts(events.SoftMuonGood.pt)
    first_jet_pt = ak.firsts(events.JetGood.pt)
    mask = (first_muon_pt / first_jet_pt) < 0.4
    return ak.where(ak.is_none(mask), False, mask)

def QCDVetoratio(events, params, **kwargs):
    first_muon_pt = ak.firsts(events.MuonGood.pt)
    first_jet_pt = ak.firsts(events.JetGood.pt)
    
    mask = (  (first_muon_pt / first_jet_pt > 0.75))
    return ak.where(ak.is_none(mask), False, mask)

def QCDVeto(events, params, **kwargs):
    mask = (  (events.MuonGood.pfRelIso04_all < 0.05) 
           & (abs(events.MuonGood.dz) < 0.01) & (abs(events.MuonGood.dxy) < 0.002) & (events.MuonGood.sip3d < 2))
    print(mask)
    return ak.where(ak.is_none(mask), False, mask)

def dimuonBTV(events, params, **kwargs):
    mask = ( (events.MuonGood.pt > 12) 
            & (abs(events.MuonGood.eta) < 2.4) & (events.MuonGood.tightId > 0.5) & (events.MuonGood.pfRelIso04_all <= 0.15))
    return ak.where(ak.is_none(mask), False, mask)

def dieleBTV(events, params, **kwargs):
    mask = ( (events.ElectronGood.pt > 12) 
            & (abs(events.ElectronGood.eta) < 1.4442) | ((events.ElectronGood.eta < 2.5) & (events.ElectronGood.eta > 1.566)) & (events.ElectronGood.mvaFall17V2Iso_WP80 > 0.5))
    return ak.where(ak.is_none(mask), False, mask)

def dilepveto(events, params, **kwargs):
    mask = ak.count( events.Electron[(events.ElectronGood.pt > 12) 
            & (abs(events.ElectronGood.eta) < 1.4442) | ((events.ElectronGood.eta < 2.5) & (events.ElectronGood.eta > 1.566)) & 
            (events.ElectronGood.mvaFall17V2Iso_WP80 > 0.5)].pt, axis = 1) + ak.count( events.MuonGood[(events.MuonGood.pt > 12) 
            & (abs(events.MuonGood.eta) < 2.4) & (events.MuonGood.tightId > 0.5) & (events.MuonGood.pfRelIso04_all <= 0.15)].pt, axis = 1) != 2
    return ak.where(ak.is_none(mask), False, mask)

def dilepmass(events, params, **kwargs):
    mask = ((ak.firsts(events.MuonGood) + ak.firsts(events.SoftMuonGood)).mass > 12) & (((ak.firsts(events.MuonGood) + ak.firsts(events.SoftMuonGood)).mass < 80) | ((ak.firsts(events.MuonGood) + ak.firsts(events.SoftMuonGood)).mass > 100))
    return ak.where(ak.is_none(mask), False, mask)

def mtw(events, params, **kwargs):
    mask = (np.sqrt(2 * ak.firsts(events.MuonGood).pt * events.MET.pt * (1 - np.cos(ak.firsts(events.MuonGood).delta_phi(events.MET)))) > 55)
    return ak.where(ak.is_none(mask), False, mask)        


mujetselbtv = Cut(
    name = 'mujetselbtv',
    function=MuonSelBTV,
    params=None
)

jetselbtv = Cut(
    name = 'jetselbtv',
    function=JetSelBTV,
    params=None
)

reqmuon = Cut(
    name = 'reqmuon',
    function=reqmu,
    params=None
)

reqsoftmuon = Cut(
    name = 'reqsoftmuon',
    function=reqsoftmu,
    params=None
)

ptratiobtv = Cut(
    name = 'ptratiobtv',
    function=JetMuPtratio,
    params=None
)

qcdveto = Cut(
    name = 'qcdveto',
    function=QCDVeto,
    params=None
)
qcdvetoratio = Cut(
    name = 'qcdvetoratio',
    function=QCDVetoratio,
    params=None
)

dimubtv = Cut(
    name = 'dimubtv',
    function=dimuonBTV,
    params=None
)

dielebtv = Cut(
    name = 'dielebtv',
    function=dieleBTV,
    params=None
)

req_dilepveto = Cut(
    name = 'req_dilepveto',
    function=dilepveto,
    params=None
)

req_dilepmass = Cut(
    name = 'req_dilepmass',
    function=dilepmass,
    params=None
)

req_mtw = Cut(
    name = 'req_mtw',
    function=mtw,
    params=None
)


dijet_pt_cut = Cut(
    name="dijet_pt_cut",
    function=DiJetPtCut,
    params={
	"pt_dijet": 120,
    },
)

jet_met_dphi_cut = Cut(
    name='jet_met_dphi_cut',
    function=DeltaPhiJetMetCut,
    params={
	"jet_met_dphi_cut": 0.6,
    },
)
# Cuts for 0-Lep channel

met_2jets_0lep = Cut(
    name="met_2jets_0lep",
    function=MetTwoJetsNoLep,
    params={
        "pt_met": 170,
        "pt_jet1": 60,
        "pt_jet2": 35,
    },
)

# Cuts for 1-Lep channel

onelep_plus_met = Cut(
    name="onelep_plus_met",
    function=OneLeptonPlusMet,
    params={
        "pt_lep": 33,
        "pt_met": 10,
    },
)

lep_met_2jets = Cut(
    name="lep_met_2jets",
    function=LepMetTwoJets,
    params={
        "pt_lep": 33,
        "pt_met": 10,
    },
)

wlnu_plus_2j = Cut(
    name="w_plus_2j",
    function=WLNuTwoJets,
    params={
        "lep_flav": "both",
        "pt_w": 100
    }
)

wmunu_plus_2j = Cut(
    name="w_plus_2j",
    function=WLNuTwoJets,
    params={
        "lep_flav": "mu",
        "pt_w": 100
    }
)

welnu_plus_2j = Cut(
    name="w_plus_2j",
    function=WLNuTwoJets,
    params={
        "lep_flav": "el",
        "pt_w": 100
    }
)

# Cuts for 1-Lep channel

dilepton = Cut(
    name="dilepton",
    function=diLepton,
    params={
        "pt_leading_lep": 33,
        "mll": {'low': 60, 'high': 120},
    },
)

mumu_channel = Cut(
    name = 'mumu',
    function=TwoMuons,
    params=None
)
ee_channel = Cut(
    name = 'ee',
    function=TwoElectrons,
    params=None
)

ll_2j = Cut(
    name = 'll_2j',
    function=TwoLepTwoJets,
    params={"lep_flav": "both",
            "pt_dilep": 60,
            "mll": {'low': 70, 'high': 120}
            }
)
mumu_2j = Cut(
    name = 'mumu_2j',
    function=TwoLepTwoJets,
    params={"lep_flav": "mu",
            "pt_dilep": 60,
            "mll": {'low': 70, 'high': 120}
        }
)
ee_2j = Cut(
    name = 'ee_2j',
    function=TwoLepTwoJets,
    params={"lep_flav": "el",
            "pt_dilep": 60,
            "mll": {'low': 70, 'high': 120}
        }
)
