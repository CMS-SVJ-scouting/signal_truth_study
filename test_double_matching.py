import sys
import uproot
import awkward as ak 
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import mplhep as hep
hep.style.use("CMS")
plt.rcParams.update({'font.size': 14})
from collections import Counter

def deltaPhi(phi1, phi2):
    dphi = phi1 - phi2
    dphi = np.where(dphi > np.pi, dphi - 2 * np.pi, dphi)
    dphi = np.where(dphi < -np.pi, dphi + 2 * np.pi, dphi)
    return abs(dphi)


mode = sys.argv[1] #dq or isr
mass = sys.argv[2] #in GeV
ptcut = int(sys.argv[3]) # in GeV
jet = [1,2,3]
rinv = [0.7, 0.5, 0.3]
data = []


for r in rinv:
    #open root files
    events = uproot.open("/ceph/mgais/SVJ_std2_UL2018_scouting_truth_study_inv_parts/SVJ_mMed-%sGeV_mDark-20GeV_rinv-%s_alpha-peak_13TeV/SVJ_mMed-%s_mDark-20_rinv-%s_alpha-peak.root"%(mass,r,mass,r))['mmtree/Events']

    #ensure at least 2 Fatjets, add pt and eta cuts
    preselection = (events['FatJet_pt'].array() > ptcut) & (abs(events['FatJet_eta'].array()) < 2.4)
    nFatJet = (ak.count(events['FatJet_pt'].array()[preselection], axis=1) >= 2) & (events['scouting_trig'].array() == 1)
    jet_pt = events['FatJet_pt'].array()[preselection][nFatJet]
    print("jet_pt: ",jet_pt)
    jet_eta = events['FatJet_eta'].array()[preselection][nFatJet]
    print("jet_eta: ",jet_eta)
    jet_phi = events['FatJet_phi'].array()[preselection][nFatJet]    

    #to match dark quarks
    if mode == 'dq':
        dq_pt = events['MatrixElementGenParticle_pt'].array()[nFatJet]
        dq_eta = events['MatrixElementGenParticle_eta'].array()[nFatJet]
        dq_phi = events['MatrixElementGenParticle_phi'].array()[nFatJet]

    #to match ISR 
    if mode == 'isr':
        dq_pt = events['ISRGluonGenParticle_pt'].array()[nFatJet]
        dq_eta = events['ISRGluonGenParticle_eta'].array()[nFatJet]
        dq_phi = events['ISRGluonGenParticle_phi'].array()[nFatJet]


    hist = np.zeros(len(jet))
    both_matched = 0
    one_matched = 0
    none_matched = 0
    matched_list = []
    dq_matchlist = []
    print("signal efficiency: ", len(jet_pt)/len(events['run'].array()))

    for n in range(len(jet_pt)):
        matched = 0
        for ndark in range(len(dq_eta[n])):
            dq_matches = 0
            for njet in range(len(jet_pt[n])):
                dR = np.sqrt(abs(dq_eta[n,ndark]-jet_eta[n,njet])**2 + deltaPhi(dq_phi[n,ndark],jet_phi[n,njet])**2)
                if dR <= 0.8:
                    try:
                        hist[njet] += 1
                    except:
                        print("matched with jet 4 or higher")
                    finally:
                        matched += 1
                        dq_matches += 1
                        #break
            if dq_matches > 1:
                print(f"Dark quark matched {dq_matches} times")
            dq_matchlist.append(dq_matches)
        matched_list.append(matched)
    matched_counter = Counter(matched_list)
    dq_counter = Counter(dq_matchlist)
    print(matched_counter)
    print(dq_counter)
    print("both matched: ", matched_counter[2])
    print("one matched: ", matched_counter[1])
    print("none matched: ", matched_counter[0])
    if mode == "dq":
        hist = np.append(hist,matched_counter[2])
        hist = np.append(hist,matched_counter[1])
    hist = np.append(hist,matched_counter[0])
    hist = hist/len(jet_pt)
    hist = list(np.around(np.array(hist),3))
    print(hist)
    data.append(hist)



