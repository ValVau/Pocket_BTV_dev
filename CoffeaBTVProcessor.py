import awkward as ak
import numpy as np
import warnings
import os
import uproot
import pandas as pd

from pocket_coffea.workflows.base import BaseProcessorABC
from pocket_coffea.utils.configurator import Configurator
from pocket_coffea.lib.hist_manager import Axis
from pocket_coffea.lib.objects import (
    jet_correction,
    lepton_selection,
    soft_lepton_selection,
    jet_selection,
    btagging,
    get_dilepton
)


class CommBTVBaseProcessor(BaseProcessorABC):
    def __init__(self, cfg: Configurator):
        super().__init__(cfg)

        # self.proc_type = self.params["proc_type"]
        # self.isArray=self.cfg["isArray"]
        #  self.campaign["run_options]/
        print(self.params,self.cfg)
    def apply_object_preselection(self, variation):
        '''
        
        '''
        # Include the supercluster pseudorapidity variable
        electron_etaSC = self.events.Electron.eta + self.events.Electron.deltaEtaSC
        self.events["Electron"] = ak.with_field(
            self.events.Electron, electron_etaSC, "etaSC"
        )
        # Build masks for selection of muons, electrons, jets, fatjets
        self.events["MuonGood"] = lepton_selection(
            self.events, "Muon", self.params
        )
        
        self.events["SoftMuonGood"] = soft_lepton_selection(
            self.events, "Muon", self.params
        )
        
        self.events["ElectronGood"] = lepton_selection(
            self.events, "Electron", self.params
        )
        leptons = ak.with_name(
           ak.concatenate((self.events.MuonGood, self.events.ElectronGood), axis=1),
            name='PtEtaPhiMCandidate',
        )
        self.events["LeptonGood"] = leptons[ak.argsort(leptons.pt, ascending=False)]

        self.events["ll"] = get_dilepton(
            self.events.ElectronGood, self.events.MuonGood
        )

        self.events["JetGood"], self.jetGoodMask = jet_selection(
            self.events, "Jet", self.params, "LeptonGood"
        )
        
        #print(self.params)
        # self.events["BJet_XXT"] = btagging(
            # self.events["JetGood"], self.params.bctagging[self.campaign].b.DeepFlav.XXT)
        # self.events["BJet_XT"] = btagging(
        #     self.events["JetGood"], self.params.btagging.T[self.campaign])
        # self.events["BJet_T"] = btagging(
        #     self.events["JetGood"], self.params.btagging.working_point[self.campaign])
        

    def count_objects(self, variation):
        self.events["nMuonGood"] = ak.num(self.events.MuonGood)
        self.events["nSoftMuonGood"] = ak.num(self.events.SoftMuonGood)
        self.events["nElectronGood"] = ak.num(self.events.ElectronGood)
        self.events["nLeptonGood"] = ak.num(self.events.LeptonGood)

        self.events["nJet"] = ak.num(self.events.Jet)
        self.events["nJetGood"] = ak.num(self.events.JetGood)
        # self.events["nBJetGood"] = ak.num(self.events.BJetGood)
        # self.events["nfatjet"]   = ak.num(self.events.FatJetGood)

    # Function that defines common variables employed in analyses and save them as attributes of `events`
    def define_common_variables_before_presel(self, variation):
        self.events["JetGood_pt"] = self.events.JetGood.pt
    def define_common_variables_after_presel(self, variation):
        self.events["MET"] = ak.zip(
                {
                    "pt": self.events.PuppiMET.pt,
                    "eta": ak.zeros_like(self.events.PuppiMET.pt),
                    "phi": self.events.PuppiMET.phi,
                    "mass": ak.zeros_like(self.events.PuppiMET.pt),
                },
                with_name="PtEtaPhiMLorentzVector",
            )
        
        
        self.events["Z"] = get_dilepton(
            self.events.SoftMuonGood , self.events.MuonGood
        )
        self.events["W"] = ak.zip(
                {
                    "pt": (self.events.MuonGood + self.events.MET).pt,
                    "eta": (self.events.MuonGood + self.events.MET).eta,
                    "phi": (self.events.MuonGood + self.events.MET).phi,
                    "mass": (self.events.MuonGood + self.events.MET).mass,
                },
                with_name="PtEtaPhiMLorentzVector",
            )
        
        self.events["hl_ptratio"] = ak.firsts(self.events.MuonGood.pt)/ak.firsts(self.events.JetGood.pt)
        self.events["soft_l_ptratio"] = ak.firsts(self.events.SoftMuonGood.pt)/ak.firsts(self.events.JetGood.pt)
        self.events["dr_soft_l_jet"] = ak.firsts(self.events.JetGood).delta_r(ak.firsts(self.events.SoftMuonGood))
        self.events["dr_l_jet"] = ak.firsts(self.events.JetGood).delta_r(ak.firsts(self.events.MuonGood))
        self.events["dr_l_soft_l"] = ak.firsts(self.events.MuonGood).delta_r(ak.firsts(self.events.SoftMuonGood))
        self.events["Z_pt"] = self.events.Z.pt
        self.events["Z_mass"] = self.events.Z.mass
        self.events["Z_eta"] = self.events.Z.eta
        self.events["Z_phi"] = self.events.Z.phi
        self.events["W_pt"] = self.events.W.pt
        self.events["W_mass"] = self.events.W.mass
        self.events["W_eta"] = self.events.W.eta
        self.events["W_phi"] = self.events.W.phi
        self.events["MET_pt"] = self.events.MET.pt
        self.events["MET_phi"] = self.events.MET.phi
        variables_to_save = ak.zip({
            "nJet": self.events.nJet,
            "nMuon": self.events.nMuonGood,
            "nSoftMuon": self.events.nSoftMuonGood,
            "hl_ptratio": self.events.hl_ptratio,
            "soft_l_ptratio": self.events.soft_l_ptratio,
            "dr_soft_l_jet": self.events.dr_soft_l_jet,
            "dr_l_jet": self.events.dr_l_jet,
            "dr_l_soft_l": self.events.dr_l_soft_l,
            "Z_pt": self.events.Z_pt,
            "Z_mass": self.events.Z_mass,
            "Z_eta": self.events.Z_eta,
            "Z_phi": self.events.Z_phi,
            "W_pt": self.events.W_pt,
            "W_mass": self.events.W_mass,
            "W_eta": self.events.W_eta,
            "W_phi": self.events.W_phi,
            "MET_pt": self.events.MET_pt,
            "MET_phi": self.events.MET_phi,
                })
        
        
        #osss = ak.values_astype(self.events.MuonGood * self.events.SoftMuonGood.charge * -1, int)
        
        self.events["JetGood_pt"] = self.events.JetGood.pt
        # self.events["dilep_deltaR"] = self.events.ll.deltaR
        with warnings.catch_warnings():
            # Suppress FutureWarning
            warnings.filterwarnings("ignore", category=FutureWarning)
            
            # Check if the directory exists
            if not os.path.exists(f"Saved_root_files/{self.events.metadata['dataset']}"):
                # If not, create it
                os.system(f"mkdir -p Saved_root_files/{self.events.metadata['dataset']}")
            
            df = ak.to_pandas(variables_to_save)    
            # Write the events to a ROOT file
            #with uproot.recreate(f"Saved_root_files/{self.events.metadata['dataset']}/{self.events.metadata['filename'].split('/')[-1].replace('.root','')}_{int(self.events.metadata['entrystart'])}_{int(self.events.metadata['entrystop'])}.root") as f: 
            #    for col in df.columns:
            #        f[col] = np.array(df[col])

            # Write the events to a Parquet file
            df.to_parquet(f"Saved_root_files/{self.events.metadata['dataset']}/{self.events.metadata['filename'].split('/')[-1].replace('.root','')}_{int(self.events.metadata['entrystart'])}_{int(self.events.metadata['entrystop'])}_vars.parquet")
