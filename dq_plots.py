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


for var in ['pt','eta','phi', 'delta_phi']:
    fig,ax = plt.subplots(figsize=(7,6))
    for r in rinv:
        #open root files
        events = uproot.open("/ceph/mgais/SVJ_std2_UL2018_scouting_truth_study_with_isr_NEW/SVJ_mMed-%sGeV_mDark-20GeV_rinv-%s_alpha-peak_13TeV/SVJ_mMed-%s_mDark-20_rinv-%s_alpha-peak.root"%(mass,r,mass,r))['mmtree/Events']

        #ensure at least 2 Fatjets
        nFatJet = (events['nFatJet'].array() >= 2) & (events['scouting_trig'].array() == 1)

        if var == 'delta_phi':
            branch = f"MatrixElementGenParticle_phi"
            bins = np.linspace(-3.5,3.5,50)
            plt.hist(deltaPhi(events[branch].array()[nFatJet][:,0], events[branch].array()[nFatJet][:,1]),histtype='step', bins=bins, density=True, label=r"$r_\mathrm{inv} = $%s"%r)
        else:
            branch = f"MatrixElementGenParticle_{var}"

            #binnings
            if var == 'pt':
                bins = np.linspace(0,1000,50)
            if var == 'eta':
                bins = np.linspace(-4,4,50)
            if var == 'phi':
                bins = np.linspace(-3.5,3.5,50)

            #plt.hist(ak.flatten(events[branch].array()),histtype='step', bins=bins, density=True, label=r"$r_\mathrm{inv} = $%s"%r)
            plt.hist(events[branch].array()[nFatJet][:,0],histtype='step', bins=bins, density=True, label=r"$r_\mathrm{inv} = $%s"%r)

    hep.cms.text("Simulation Private Work", fontsize=20)
    plt.figtext(0.65, 0.6, r"m$_{\mathrm{Z'}}$ = $%s\,$GeV"%mass, fontsize=13)
    plt.xlabel(f"dark quark {var}")
    plt.ylabel("norm. to unit area")
    plt.legend(loc=1)
    fig.savefig(f'/web/mgais/public_html/scouting_truth_study/dark_quark_plots/dark_quark_{var}_mZprime-{mass}.png',bbox_inches='tight',dpi=300)
    fig.savefig(f'/web/mgais/public_html/scouting_truth_study/dark_quark_plots/dark_quark_{var}_mZprime-{mass}.pdf',bbox_inches='tight')
