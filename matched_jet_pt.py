#plot pt distribution of the subleading jets matched to a dark quark to estimate how much signal is lost depending on the jet pt cut

import sys
import uproot
import awkward as ak 
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
ptcut = int(sys.argv[2]) # in GeV
rinv = [0.7, 0.5, 0.3]

fig,ax = plt.subplots(figsize=(7,6))

for r in rinv:
    #open root files
    events = uproot.open("/ceph/mgais/SVJ_std2_UL2018_scouting_truth_study_with_isr_NEW/SVJ_mMed-%sGeV_mDark-20GeV_rinv-%s_alpha-peak_13TeV/SVJ_mMed-%s_mDark-20_rinv-%s_alpha-peak.root"%(mass,r,mass,r))['mmtree/Events']
    
    #ensure at least 2 Fatjets, add pt, eta and trigger cuts
    preselection = (events['FatJet_pt'].array() > ptcut) & (abs(events['FatJet_eta'].array()) < 2.4)
    nFatJet = (ak.count(events['FatJet_pt'].array()[preselection], axis=1) >= 2) & (events['scouting_trig'].array() == 1)

    jet_pt = events['FatJet_pt'].array()[preselection][nFatJet]
    print("jet_pt: ",jet_pt)
    jet_eta = events['FatJet_eta'].array()[preselection][nFatJet]
    print("jet_eta: ",jet_eta)
    jet_phi = events['FatJet_phi'].array()[preselection][nFatJet]    

    dq_pt = events['MatrixElementGenParticle_pt'].array()[nFatJet]
    dq_eta = events['MatrixElementGenParticle_eta'].array()[nFatJet]
    dq_phi = events['MatrixElementGenParticle_phi'].array()[nFatJet]


    print("signal efficiency: ",len(jet_pt)/len(events['run'].array()))
    matched_leading_jet_pt = []
    matched_subleading_jet_pt = []

    for n in range(len(jet_pt)):
        is_matched = 0
        for ndark in range(len(dq_eta[n])):
            #for njet in range(len(jet_pt[n])):
            for njet in range(2):
                dR = np.sqrt(abs(dq_eta[n,ndark]-jet_eta[n,njet])**2 + deltaPhi(dq_phi[n,ndark],jet_phi[n,njet])**2)
                if dR <= 0.8:
                    # save pt of leading and subleading jets if matched
                    if njet == 0:
                        matched_leading_jet_pt.append(jet_pt[n,njet])
                    if njet == 1:
                        matched_subleading_jet_pt.append(jet_pt[n,njet])

    plt.hist(matched_subleading_jet_pt, histtype='step', density=True, bins=np.linspace(0,700,50), label=r'$r_\mathrm{inv}=%s$'%r)
    #plt.hist(jet_pt[:,1], histtype='step', density=True, bins=np.linspace(0,700,50), label=r'with unmatched, $r_\mathrm{inv}=%s$'%r)
    


hep.cms.text("Simulation Private Work", fontsize=15)
plt.figtext(0.65, 0.7, r"m$_{\mathrm{Z'}}$ = $%s\,$GeV"%mass, fontsize=12)
plt.figtext(0.65, 0.66, r"$p_{\mathrm{T}} > %s\,$GeV, $|\eta| < 2.4$"%ptcut, fontsize=12)
plt.figtext(0.65, 0.62, r"$\geq$ 2 AK8 jets", fontsize=12)
plt.xlabel(r"matched subleading AK8 jet $p_\mathrm{T}$ [GeV]")
plt.ylabel("norm. to unit area")
#plt.ylabel("arbitrary units")
plt.legend(loc=1)
fig.savefig(f'/web/mgais/public_html/scouting_truth_study/FatJet_plots/matched_subleading_AK8_pt_mZ-{mass}_ptcut-{ptcut}.png',bbox_inches='tight',dpi=300)
fig.savefig(f'/web/mgais/public_html/scouting_truth_study/FatJet_plots/matched_subleading_AK8_pt_mZ-{mass}_ptcut-{ptcut}.pdf',bbox_inches='tight')



