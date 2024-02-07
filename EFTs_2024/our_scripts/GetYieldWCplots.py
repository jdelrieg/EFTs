import gzip, pickle, json, os, uproot, sys
from coffea import hist
from topcoffea.plotter.plotter import plotter
import matplotlib.pyplot as plt
import numpy as np

from topcoffea.plotter.plotter import plotter, DrawComp, Rebin, loadHistos, GetHisto

# Path to histograms
path ='/nfs/fanae/user/jriego/TFGs/EFTs_2024/new_tfg/pruebattH.pkl.gz'

# Categories to integrate
integrate_categories = ['sample']

cat_2lss = ['2lss_4t_m_4j', '2lss_4t_m_5j', '2lss_4t_m_6j', '2lss_4t_m_7j', '2lss_4t_p_4j', '2lss_4t_p_5j', '2lss_4t_p_6j', '2lss_4t_p_7j', '2lss_CR_1j', '2lss_CR_2j', '2lss_m_4j', '2lss_m_5j', '2lss_m_6j', '2lss_m_7j', '2lss_p_4j', '2lss_p_5j', '2lss_p_6j', '2lss_p_7j']
cat_3l = ['3l_m_offZ_1b_2j', '3l_m_offZ_1b_3j', '3l_m_offZ_1b_4j', '3l_m_offZ_1b_5j', '3l_m_offZ_2b_2j', '3l_m_offZ_2b_3j', '3l_m_offZ_2b_4j', '3l_m_offZ_2b_5j', '3l_onZ_1b_2j', '3l_onZ_1b_3j', '3l_onZ_1b_4j', '3l_onZ_1b_5j', '3l_onZ_2b_2j', '3l_onZ_2b_3j', '3l_onZ_2b_4j', '3l_onZ_2b_5j', '3l_p_offZ_1b_2j', '3l_p_offZ_1b_3j', '3l_p_offZ_1b_4j', '3l_p_offZ_1b_5j', '3l_p_offZ_2b_2j', '3l_p_offZ_2b_3j', '3l_p_offZ_2b_4j', '3l_p_offZ_2b_5j']
cat_4l = ['4l_2j', '4l_3j', '4l_4j']
categories = {
  'systematic' : 'nominal',
  'appl' : ['isSR_3l'],
  'channel' : cat_3l,
}
catname = '3l'

# Select the observable
observable = 'ht'

# Rebin the histogram
rebin = {observable : 8}

# Finally, get the histogram
h = GetHisto(path, observable, categories=categories, integrate=integrate_categories, rebin=rebin)

def GetYieldForWCval(h, WC, val):
  WCs = np.array(h._wcnames)
  vals = np.zeros_like(WCs, dtype=float)
  vals += np.where(WCs==WC, val, vals)
  h.set_wilson_coeff_from_array(vals)
  y = h.values()[list(h.values().keys())[0]].sum()
  return y

SM = GetYieldForWCval(h, 'cpt', 0)

WCs = 'cpt, ctp, cptb, cQlMi, cQq81, cQq11, cQl3i, ctq8, ctlTi, ctq1, ctli, cQq13, cbW, cpQM, cpQ3, ctei, cQei, ctW, ctlSi, cQq83, ctZ, ctG'.replace(' ', '').split(',')
x = np.linspace(-3, 3, 1000)

fig, ax = plt.subplots(1, 1, figsize=(7,7))

for wc in WCs[0:10]:
  yields = [ GetYieldForWCval(h, wc, ix)/SM for ix in x]
  ax.plot(x, yields, '-', label=wc)

ax.set_xlabel('WC value')
ax.set_ylabel('Yield/SM')
legend = ax.legend(loc='upper center')
CMS  = plt.text(0., 1., r"$\bf{CMS}$ Preliminary", fontsize=16, horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes)
CMS  = plt.text(0.8, 1., catname, fontsize=16, horizontalalignment='right', verticalalignment='bottom', transform=ax.transAxes)

fig.savefig('yieldsWC.png')
