import fastjet as fj
import awkward as ak
import vector
import uproot
import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import mplhep as hep
hep.style.use("CMS")
plt.rcParams.update({'font.size': 14})


def deltaPhi(phi1, phi2):
    dphi = phi1 - phi2
    dphi = np.where(dphi > np.pi, dphi - 2 * np.pi, dphi)
    dphi = np.where(dphi < -np.pi, dphi + 2 * np.pi, dphi)
    return abs(dphi)

mass = sys.argv[1] #in GeV           
rinv = sys.argv[2]
ptcut = int(sys.argv[3]) # in GeV
jet = [1,2,3,4]
data = []

events = uproot.open("/ceph/mgais/SVJ_std2_UL2018_scouting_truth_study_with_isr_NEW/SVJ_mMed-%sGeV_mDark-20GeV_rinv-%s_alpha-peak_13TeV/SVJ_mMed-%s_mDark-20_rinv-%s_alpha-peak.root"%(mass,rinv,mass,rinv))['mmtree/Events']

vector.register_awkward()
PFcands = ak.zip({
    "pt": events['PFCands_pt'].array(),
    "eta": events['PFCands_eta'].array(),
    "phi": events['PFCands_phi'].array(),
    "mass": events['PFCands_mass'].array(),
}, with_name="Momentum4D")


jetdef = fj.JetDefinition(fj.antikt_algorithm, 0.8)
cluster = fj.ClusterSequence(PFcands, jetdef)
ak8 = cluster.inclusive_jets(ptcut)

# manually sort the jets based on pt because python fastjet does not do so reliably
sorted_indices = ak.argsort(ak8.pt, ascending=False)
ak8 = ak8[sorted_indices]

preselection = abs(ak8.eta) < 2.4
nFatJet = (ak.count(ak8.pt[preselection], axis=1) >= 2) & (events['scouting_trig'].array() == 1)
# separate pt_cut on the first two jets, this is to potentially allow the third jet to have lower pt
pt_cut = ak8.pt[preselection][nFatJet][:,1] >= 150

jet_pt = ak8.pt[preselection][nFatJet][pt_cut]
jet_eta = ak8.eta[preselection][nFatJet][pt_cut]
jet_phi = ak8.phi[preselection][nFatJet][pt_cut]

pfcand_pt = events['PFCands_pt'].array()[nFatJet][pt_cut]
pfcand_eta = events['PFCands_eta'].array()[nFatJet][pt_cut]
pfcand_phi = events['PFCands_phi'].array()[nFatJet][pt_cut]

const = cluster.constituents(ptcut)
const = const[sorted_indices]
const_pt = const.pt[preselection][nFatJet][pt_cut]
const_eta = const.eta[preselection][nFatJet][pt_cut]
const_phi = const.phi[preselection][nFatJet][pt_cut]

dq_pt = events['MatrixElementGenParticle_pt'].array()[nFatJet][pt_cut]
dq_eta = events['MatrixElementGenParticle_eta'].array()[nFatJet][pt_cut]
dq_phi = events['MatrixElementGenParticle_phi'].array()[nFatJet][pt_cut]

isr_pt = events['ISRGluonGenParticle_pt'].array()[nFatJet][pt_cut]
isr_eta = events['ISRGluonGenParticle_eta'].array()[nFatJet][pt_cut]
isr_phi = events['ISRGluonGenParticle_phi'].array()[nFatJet][pt_cut]

print("signal efficiency: ", len(jet_pt)/len(events['run'].array()))
for n in range(100):
    matched = 0
    for ndark in range(len(dq_eta[n])):
        for njet in range(len(jet_pt[n])):
            dR = np.sqrt(abs(dq_eta[n,ndark]-jet_eta[n,njet])**2 + deltaPhi(dq_phi[n,ndark],jet_phi[n,njet])**2)
            if dR <= 0.8:
                matched += 1
                break
                
    if matched < 2:
        fig,ax = plt.subplots(figsize=(7,6))
        sizes = pfcand_pt[n] * 10
        # draw all PF candidates in grey
        plt.scatter(pfcand_eta[n], pfcand_phi[n], marker='.', s=sizes, label=f"PF candidates", color="grey")
        # then redraw those that are jet constituents in color
        print("event: ", n)
        print("DQ: ", dq_eta[n], dq_phi[n])
        print("ISR: ", isr_eta[n], isr_phi[n])
        for njet in range(len(jet_pt[n])):
            print("constituents: ", const_eta[n,njet], const_phi[n,njet])
            sizes=const_pt[n,njet] * 10
            plt.scatter(const_eta[n,njet], const_phi[n,njet], marker='.', s=sizes, label=f"Jet {njet+1}")
        # draw the dark quarks
        plt.scatter(dq_eta[n], dq_phi[n], marker='x', label=f"dark quark", color="red")
        plt.scatter(isr_eta[n], isr_phi[n], marker='x', label=f"ISR gluon", color="green")
        plt.xlim(-2.4,2.4)
        plt.ylim(-3.14,3.14)
        plt.xlabel(r"$\eta$")
        plt.ylabel(r"$\phi$")
        hep.cms.text("Simulation Private Work", fontsize=15)
        plt.title(r"m$_{\mathrm{Z'}}$ = $%s\,$GeV, $r_\mathrm{inv}=%s$"%(mass,rinv), fontsize=12, loc='right')
        plt.legend(loc=1, frameon=True)
        fig.savefig(f'/web/mgais/public_html/scouting_truth_study/pfcand_eta_phi_with_ISR/eta_phi_mZ-{mass}_rinv-{rinv}_ptcut-{ptcut}_event-{n}.png',bbox_inches='tight',dpi=300)
        fig.savefig(f'/web/mgais/public_html/scouting_truth_study/pfcand_eta_phi_with_ISR/eta_phi_mZ-{mass}_rinv-{rinv}_ptcut-{ptcut}_event-{n}.pdf',bbox_inches='tight')
        plt.close()
