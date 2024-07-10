# plot delta phi, delta eta, delta R between the two matched jets
# the goal is to check if there is an easy way to identify the matched jets when we loosen the requirements on the 3rd jet
# make sure you have the correct fastjet version that works on awkward arrays

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
ptcut = int(sys.argv[3]) # pt cut for 3rd jet onwards in GeV

events = uproot.open("/ceph/mgais/SVJ_std2_UL2018_scouting_truth_study_with_isr_NEW/SVJ_mMed-%sGeV_mDark-20GeV_rinv-%s_alpha-peak_13TeV/SVJ_mMed-%s_mDark-20_rinv-%s_alpha-peak.root"%(mass,rinv,mass,rinv))['mmtree/Events']

print(len(events['run'].array()))

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

# python fastjet jet sorting does not work consistently
sorted_indices = ak.argsort(ak8.pt, ascending=False)
ak8 = ak8[sorted_indices]

# jet and event selections
preselection = abs(ak8.eta) < 2.4
nFatJet = (ak.count(ak8.pt[preselection], axis=1) >= 2) & (events['scouting_trig'].array() == 1)
# require the first two jets to be above 150 GeV, 3rd jet onward are set by input ptcut
pt_cut = ak8.pt[preselection][nFatJet][:,1] >= 150

jet_pt = ak8.pt[preselection][nFatJet][pt_cut]
jet_eta = ak8.eta[preselection][nFatJet][pt_cut]
jet_phi = ak8.phi[preselection][nFatJet][pt_cut]

dq_pt = events['MatrixElementGenParticle_pt'].array()[nFatJet][pt_cut]
dq_eta = events['MatrixElementGenParticle_eta'].array()[nFatJet][pt_cut]
dq_phi = events['MatrixElementGenParticle_phi'].array()[nFatJet][pt_cut]


print("signal efficiency: ", len(jet_pt)/len(events['run'].array()))
delta_phi = []
delta_eta = []
delta_R = []
for n in range(len(jet_pt)):
    matched_eta = []
    matched_phi = []
    for dq in range(len(dq_eta[n])):
        for jet in range(len(jet_pt[n])):
            dR = np.sqrt(abs(dq_eta[n,dq]-jet_eta[n,jet])**2 + deltaPhi(dq_phi[n,dq],jet_phi[n,jet])**2)
            if dR <= 0.8:
                matched_eta.append(jet_eta[n,jet])
                matched_phi.append(jet_phi[n,jet])
                break
                
    if len(matched_eta) >= 2:
        #calculate deltas
        delta_phi.append(deltaPhi(matched_phi[0], matched_phi[1]))
        delta_eta.append(abs(matched_eta[0] - matched_eta[1]))
        delta_R.append(np.sqrt(abs(matched_eta[0] - matched_eta[1])**2 + deltaPhi(matched_phi[0], matched_phi[1])**2))

#plot deltas
"""
fig,ax = plt.subplots(figsize=(7,6))
bins = np.linspace(0,5,30)
plt.hist(delta_eta, histtype='step', bins=bins)
plt.xlabel(r"matched jet $\Delta\eta$")
plt.ylabel(r"arbitrary units")
hep.cms.text("Simulation Private Work", fontsize=15)
plt.figtext(0.65, 0.6, r"m$_{\mathrm{Z'}}$ = $%s\,$GeV"%mass, fontsize=13)
plt.figtext(0.65, 0.56, r"r$_{\mathrm{inv}}$ = %s"%rinv, fontsize=13)
#plt.legend(loc=1, frameon=True)
fig.savefig(f'/web/mgais/public_html/scouting_truth_study/matched_jet_plots/matched_jet_delta_eta_mZ-{mass}_rinv-{rinv}_ptcut-{ptcut}.png',bbox_inches='tight',dpi=300)
fig.savefig(f'/web/mgais/public_html/scouting_truth_study/matched_jet_plots/matched_jet_delta_eta_mZ-{mass}_rinv-{rinv}_ptcut-{ptcut}.pdf',bbox_inches='tight')
plt.close()
"""

fig,ax = plt.subplots(figsize=(7,6))
bins = np.linspace(0,3.14,30)
plt.hist(delta_phi, histtype='step', bins=bins)
plt.xlabel(r"matched jet $\Delta\phi$")
plt.ylabel(r"arbitrary units")
hep.cms.text("Simulation Private Work", fontsize=14)
plt.title(r"m$_{\mathrm{Z'}}$ = $%s\,$GeV, $r_\mathrm{inv}=%s$"%(mass,rinv), fontsize=11, loc='right')
#plt.legend(loc=1, frameon=True)
fig.savefig(f'/web/mgais/public_html/scouting_truth_study/matched_jet_plots/matched_jet_delta_phi_mZ-{mass}_rinv-{rinv}_ptcut-{ptcut}.png',bbox_inches='tight',dpi=300)
fig.savefig(f'/web/mgais/public_html/scouting_truth_study/matched_jet_plots/matched_jet_delta_phi_mZ-{mass}_rinv-{rinv}_ptcut-{ptcut}.pdf',bbox_inches='tight')
plt.close()

fig,ax = plt.subplots(figsize=(7,6))
bins = np.linspace(0,5,30)
plt.hist(delta_R, histtype='step', bins=bins)
plt.xlabel(r"matched jet $\Delta R$")
plt.ylabel(r"arbitrary units")
hep.cms.text("Simulation Private Work", fontsize=14)
plt.title(r"m$_{\mathrm{Z'}}$ = $%s\,$GeV, $r_\mathrm{inv}=%s$"%(mass,rinv), fontsize=11, loc='right')
#plt.legend(loc=1, frameon=True)
fig.savefig(f'/web/mgais/public_html/scouting_truth_study/matched_jet_plots/matched_jet_delta_R_mZ-{mass}_rinv-{rinv}_ptcut-{ptcut}.png',bbox_inches='tight',dpi=300)
fig.savefig(f'/web/mgais/public_html/scouting_truth_study/matched_jet_plots/matched_jet_delta_R_mZ-{mass}_rinv-{rinv}_ptcut-{ptcut}.pdf',bbox_inches='tight')
plt.close()
