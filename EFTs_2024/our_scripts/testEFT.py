import gzip, pickle, json, os, uproot, sys
from coffea import hist
from topcoffea.plotter.plotter import plotter
import matplotlib.pyplot as plt
import numpy as np

from topcoffea.plotter.plotter import plotter, DrawComp, Rebin, loadHistos, GetHisto

# Path to histograms
path = '/nfs/fanae/user/jriego/TFGs/EFTs_2024/new_tfg/pruebattH.pkl.gz'

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

# Select the observable
observable = 'ht'

# Rebin the histogram
rebin = {observable : 8}

# Finally, get the histogram
h = GetHisto(path, observable, categories=categories, integrate=integrate_categories, rebin=rebin)

### Now, draw for some WC values

# The first is SM, the rest have some modified WCs
WClist = [{}, {'cQq83':1.}, {'ctW':1.}, {'ctp':1.}, {'cQq83':-1.}, {'ctW':-1.}, {'ctp':-1.}]

colors = ['k', 'green', 'red', 'blue', 'green', 'red', 'blue']

# Labels for the legend
labels = ['SM', '$c_{Qq83} =\pm 1$', '$c_{tW}=\pm 1$', '$c_{tp}=\pm1$', None, None, None]


# Let's do it
DrawComp(h, colors=colors, doFill=1,labels=labels,EFT_WCdict=WClist, outname='temp.png')

