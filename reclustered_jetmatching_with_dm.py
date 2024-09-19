import fastjet as fj
import vector
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
jet_radius = float(sys.argv[4])
jet_type = "AK" + str(int(10*jet_radius))
jet = [1,2,3]
rinv = [0.7, 0.5, 0.3]
data = []


for r in rinv:
    #open root files
    events = uproot.open("/ceph/mgais/SVJ_std2_UL2018_scouting_truth_study_inv_parts/SVJ_mMed-%sGeV_mDark-20GeV_rinv-%s_alpha-peak_13TeV/SVJ_mMed-%s_mDark-20_rinv-%s_alpha-peak.root"%(mass,r,mass,r))['mmtree/Events']


    vector.register_awkward()
    PFcands = ak.zip({
        "pt": ak.concatenate([events['PFCands_pt'].array(),events['DarkMatterParticles_pt'].array()], axis=1),
        "eta": ak.concatenate([events['PFCands_eta'].array(),events['DarkMatterParticles_eta'].array()], axis=1),
        "phi": ak.concatenate([events['PFCands_phi'].array(),events['DarkMatterParticles_phi'].array()], axis=1),
        "mass": ak.concatenate([events['PFCands_mass'].array(),events['DarkMatterParticles_mass'].array()], axis=1),
    }, with_name="Momentum4D")

    print(events['PFCands_pt'].array())
    print(events['DarkMatterParticles_pt'].array())
    print(ak.concatenate([events['PFCands_pt'].array(),events['DarkMatterParticles_pt'].array()], axis=1))
    
    jetdef = fj.JetDefinition(fj.antikt_algorithm, jet_radius)
    cluster = fj.ClusterSequence(PFcands, jetdef)
    jets = cluster.inclusive_jets(ptcut)

    sorted_indices = ak.argsort(jets.pt, ascending=False)
    jets = jets[sorted_indices]
    print(jets)
    print(jets.pt)
    print(jets.eta)
        
    preselection = abs(jets.eta) < 2.4
    nFatJet = (ak.count(jets.pt[preselection], axis=1) >= 2) #& (events['scouting_trig'].array() == 1)

    #to match dark quarks
    if mode == 'dq':
        dq_sel = abs(events['MatrixElementGenParticle_eta'].array()) < 2.4
        n_dq = ak.count(events['MatrixElementGenParticle_pt'].array()[dq_sel], axis=1) == 2
        dq_pt = events['MatrixElementGenParticle_pt'].array()[nFatJet & n_dq]
        dq_eta = events['MatrixElementGenParticle_eta'].array()[nFatJet & n_dq]
        dq_phi = events['MatrixElementGenParticle_phi'].array()[nFatJet & n_dq]

    #to match ISR 
    if mode == 'isr':
        dq_pt = events['ISRGluonGenParticle_pt'].array()[nFatJet]
        dq_eta = events['ISRGluonGenParticle_eta'].array()[nFatJet]
        dq_phi = events['ISRGluonGenParticle_phi'].array()[nFatJet]

    jet_pt = jets.pt[preselection][nFatJet & n_dq]
    jet_eta = jets.eta[preselection][nFatJet & n_dq]
    jet_phi = jets.phi[preselection][nFatJet & n_dq]


    hist = np.zeros(len(jet))
    both_matched = 0
    one_matched = 0
    none_matched = 0
    matched_list = []
    print("signal efficiency: ", len(jet_pt)/len(events['run'].array()))
    print("number of events: ", len(jet_pt))
    
    for n in range(len(jet_pt)):
        matched = 0
        for ndark in range(len(dq_eta[n])):
            for njet in range(len(jet_pt[n])):
                dR = np.sqrt(abs(dq_eta[n,ndark]-jet_eta[n,njet])**2 + deltaPhi(dq_phi[n,ndark],jet_phi[n,njet])**2)
                if dR <= jet_radius:
                    try:
                        hist[njet] += 1
                    except:
                        print("matched with jet 4 or higher")
                    finally:
                        matched += 1
                        break
        matched_list.append(matched)
    matched_counter = Counter(matched_list)
    print(matched_counter)
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



fig,ax = plt.subplots(figsize=(7,6))
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)

im = ax.imshow(data,vmin=0,vmax=1)

# Show all ticks and label them with the respective list entries
if mode == "dq":
    jet.append(r"both $q_D$"+"\n matched")
    jet.append(r"one $q_D$"+"\n unmatched")
    jet.append(r"both $q_D$"+"\n unmatched")
if mode == "isr":
    jet.append(r"unmatched")
    
ax.set_xticks(np.arange(len(jet)))
ax.set_yticks(np.arange(len(rinv)))
ax.set_xticklabels(jet)
ax.set_yticklabels(rinv)
ax.tick_params(which='minor',bottom = False,left=False,right=False,top=False)
ax.tick_params(which='major',right=False,top=False,left=False,bottom=False)
ax.set_xlabel(f"{jet_type} jet number")
ax.set_ylabel(r"$r_\mathrm{inv}$")

# Create colorbar
cbar = ax.figure.colorbar(im, cax=cax)
if mode == 'dq':
    cbar.ax.set_ylabel("dark quark matching probability",loc='center', fontsize=13)
if mode == 'isr':
    cbar.ax.set_ylabel("ISR matching probability",loc='center', fontsize=13)

# Loop over data dimensions and create text annotations.
for i in range(len(rinv)):
    for j in range(len(jet)):
        text = ax.text(j, i, data[i][j],
                       ha="center", va="center", color="w")

ax.set_title(r"$m_\mathrm{Z'}=%s\,$GeV"%mass,loc='right',fontsize=13)
fig.tight_layout()

hep.cms.text("Simulation Private Work", fontsize=17, ax=ax)
if mode == "dq":
    fig.savefig(f'/web/mgais/public_html/scouting_truth_study/{jet_type}_with_DM/{mode}_{mass}_{jet_type}_matching_pt{ptcut}_without_trigger_with_dq_eta.png',bbox_inches='tight',dpi=300)
    fig.savefig(f'/web/mgais/public_html/scouting_truth_study/{jet_type}_with_DM/{mode}_{mass}_{jet_type}_matching_pt{ptcut}_without_trigger_with_dq_eta.pdf',bbox_inches='tight')
if mode == "isr":
    fig.savefig(f'/web/mgais/public_html/scouting_truth_study/{jet_type}_with_DM/{mode}_{mass}_{jet_type}_matching_pt{ptcut}_dq_cut.png',bbox_inches='tight',dpi=300)
    fig.savefig(f'/web/mgais/public_html/scouting_truth_study/{jet_type}_with_DM/{mode}_{mass}_{jet_type}_matching_pt{ptcut}_dq_cut.pdf',bbox_inches='tight')
