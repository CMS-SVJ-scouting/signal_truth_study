# plots 2d histograms of the matching efficiency of dark quarks depending on their pt and eta

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
    """Compute delta Phi between two objects, such that it falls into [-pi,pi]"""
    dphi = phi1 - phi2
    dphi = np.where(dphi > np.pi, dphi - 2 * np.pi, dphi)
    dphi = np.where(dphi < -np.pi, dphi + 2 * np.pi, dphi)
    return abs(dphi)

def deltaR(eta1, phi1, eta2, phi2):
    """Compute delta R between two objects."""
    delta_eta = eta1 - eta2
    delta_phi = deltaPhi(phi1, phi2)
    return np.sqrt(delta_eta**2 + delta_phi**2)


mass = sys.argv[1] # mediator mass in GeV
ptcut = 150 # AK8 jet pt cut in GeV
#ptcut = int(sys.argv[2]) 
rinv = [0.7, 0.5, 0.3]

for r in rinv:
    # open root files
    events = uproot.open("/ceph/mgais/SVJ_std2_UL2018_scouting_truth_study_with_isr_NEW/SVJ_mMed-%sGeV_mDark-20GeV_rinv-%s_alpha-peak_13TeV/SVJ_mMed-%s_mDark-20_rinv-%s_alpha-peak.root"%(mass,r,mass,r))['mmtree/Events']

    # AK8 jet requirements
    preselection = (events['FatJet_pt'].array() > ptcut) & (abs(events['FatJet_eta'].array()) < 2.4)
    # event level cuts
    nFatJet = (ak.count(events['FatJet_pt'].array()[preselection], axis=1) >= 2) & (events['scouting_trig'].array() == 1)

    jet_pt = events['FatJet_pt'].array()[preselection][nFatJet]
    jet_eta = events['FatJet_eta'].array()[preselection][nFatJet]
    jet_phi = events['FatJet_phi'].array()[preselection][nFatJet]    

    dq_pt = events['MatrixElementGenParticle_pt'].array()[nFatJet]
    dq_eta = events['MatrixElementGenParticle_eta'].array()[nFatJet]
    dq_phi = events['MatrixElementGenParticle_phi'].array()[nFatJet]

    pt_bins = np.linspace(0,1000,21)
    eta_bins = np.linspace(0,4,21)
    efficiency_hist = np.zeros((len(pt_bins)-1,len(eta_bins)-1))
    total_hist = np.zeros((len(pt_bins)-1,len(eta_bins)-1))

    # match the dark quarks to jets and save matched dark quarks separately
    for n in range(len(jet_pt)):
        matched_dq_pt = []
        matched_dq_eta = []
        total_dq_pt = []
        total_dq_eta = []
        for dq in range(len(dq_pt[n])):
            total_dq_pt.append(dq_pt[n,dq])
            total_dq_eta.append(dq_eta[n,dq])
            for jet in range(len(jet_pt[n])):
                dR = deltaR(dq_eta[n,dq], dq_phi[n,dq], jet_eta[n,jet], jet_phi[n,jet])
                if dR <= 0.8:
                    matched_dq_pt.append(dq_pt[n,dq])
                    matched_dq_eta.append(dq_eta[n,dq])
                    break # to avoid double matching
                    
        # Fill one histogram with matched dark quarks and another one with all dark quarks
        hist, _, _ = np.histogram2d(matched_dq_pt, np.abs(matched_dq_eta), bins=[pt_bins, eta_bins])
        efficiency_hist += hist
        total_dqs, _, _ = np.histogram2d(total_dq_pt, np.abs(total_dq_eta), bins=[pt_bins, eta_bins])
        total_hist += total_dqs

    # Divide histograms to calculate matching efficiency for each bin
    out = np.full(efficiency_hist.shape, np.nan) # fill the output array with nans to differentiate between bins with low efficiency and empty bins
    efficiency_hist = np.divide(efficiency_hist, total_hist, out=out, where=total_hist!=0)
    #efficiency_hist = np.divide(efficiency_hist, total_hist, out=np.zeros_like(efficiency_hist), where=total_hist!=0)
    #print(efficiency_hist)

    # Plot
    fig,ax = plt.subplots(figsize=(8,6))
    plt.imshow(efficiency_hist.T, extent=[pt_bins[0], pt_bins[-1], eta_bins[0], eta_bins[-1]], aspect='auto', origin='lower')
    plt.colorbar(label='matching efficiency')
    plt.xlabel('dark quark $p_T$ [GeV]')
    plt.ylabel('dark quark $|\eta|$')
    hep.cms.text("Simulation Private Work", fontsize=14)
    plt.title(r"m$_{\mathrm{Z'}}$ = $%s\,$GeV, $r_\mathrm{inv}=%s$"%(mass,r), fontsize=11, loc='right') 
    r_inv = str(r).replace('.','')
    fig.savefig(f'/web/mgais/public_html/scouting_truth_study/2D_pt_eta/2D_eff_mZ-{mass}_rinv-{r_inv}.png',bbox_inches='tight',dpi=300)
    fig.savefig(f'/web/mgais/public_html/scouting_truth_study/2D_pt_eta/2D_eff_mZ-{mass}_rinv-{r_inv}.pdf',bbox_inches='tight')



