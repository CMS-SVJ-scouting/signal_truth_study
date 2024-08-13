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


#mode = sys.argv[1] #dq or isr
mass = sys.argv[1] #in GeV
rinv = [0.7, 0.5, 0.3]

bkg = uproot.open("/ceph/mgais/svj_final/QCD_tagger/final_QCD_tagger1100/combined/QCD_weighted.root")['events']
bkg_nFatJet = (bkg['n_fatjet'].array() >= 2) & (bkg['scouting_trig'].array() == 1)
bkg_pt = bkg["FatJet_pt"].array()[bkg_nFatJet]
bkg_eta = bkg["FatJet_eta"].array()[bkg_nFatJet]
bkg_phi = bkg["FatJet_phi"].array()[bkg_nFatJet]
bkg_delta_phi = deltaPhi(bkg_phi[:,0], bkg_phi[:,1])
bkg_weights = bkg["evtweight"].array()[bkg_nFatJet]

for var in ['pt','eta','phi', 'delta_phi']:
    fig,ax = plt.subplots(figsize=(7,6))
    #binnings
    if var == 'pt':
        bins = np.linspace(0,1000,50)
        plt.hist(bkg_pt[:,0],histtype='step', bins=bins, weights=bkg_weights, density=True, label=r"QCD")
    if var == 'eta':
        bins = np.linspace(-4,4,50)
        plt.hist(bkg_eta[:,0],histtype='step', bins=bins, weights=bkg_weights, density=True, label=r"QCD")
    if var == 'phi':
        bins = np.linspace(-3.5,3.5,50)
        plt.hist(bkg_phi[:,0],histtype='step', bins=bins, weights=bkg_weights, density=True, label=r"QCD")
    if var == 'delta_phi':
        bins = np.linspace(0,3.5,50)
        plt.hist(bkg_delta_phi,histtype='step', bins=bins, weights=bkg_weights, density=True, label=r"QCD")

    for r in rinv:
        #open root files
        sig = uproot.open("/ceph/mgais/SVJ_std2_UL2018_scouting_truth_study_with_isr_NEW/SVJ_mMed-%sGeV_mDark-20GeV_rinv-%s_alpha-peak_13TeV/SVJ_mMed-%s_mDark-20_rinv-%s_alpha-peak.root"%(mass,r,mass,r))['mmtree/Events']

        #ensure at least 2 Fatjets
        nFatJet = (sig['nFatJet'].array() >= 2) & (sig['scouting_trig'].array() == 1)

        if var == 'delta_phi':
            branch = f"FatJet_phi"
            plt.hist(deltaPhi(sig[branch].array()[nFatJet][:,0], sig[branch].array()[nFatJet][:,1]),histtype='step', bins=bins, density=True, label=r"$signal, r_\mathrm{inv} = $%s"%r)
        else:
            branch = f"FatJet_{var}"
            #plt.hist(ak.flatten(sig[branch].array()),histtype='step', bins=bins, density=True, label=r"$r_\mathrm{inv} = $%s"%r)
            plt.hist(sig[branch].array()[nFatJet][:,0],histtype='step', bins=bins, density=True, label=r"$signal, r_\mathrm{inv} = $%s"%r)

    hep.cms.text("Simulation Private Work", fontsize=20)
    plt.figtext(0.65, 0.6, r"m$_{\mathrm{Z'}}$ = $%s\,$GeV"%mass, fontsize=13)
    plt.xlabel(f"leading jet {var}")
    plt.ylabel("norm. to unit area")
    plt.legend(loc=1)
    fig.savefig(f'/web/mgais/public_html/scouting_truth_study/delta_phi_jets/{var}_mZprime-{mass}.png',bbox_inches='tight',dpi=300)
    fig.savefig(f'/web/mgais/public_html/scouting_truth_study/delta_phi_jets/{var}_mZprime-{mass}.pdf',bbox_inches='tight')
