from collections import defaultdict, OrderedDict
import gzip
import pickle
import json
import os
import uproot
import matplotlib.pyplot as plt
import numpy as np
from coffea import hist, processor
from coffea.hist import plot
from cycler import cycler
from topcoffea.plotter.OutText import OutText
import mplhep as hep

#lumi=??


path2 = '/nfs/fanae/user/jriego/TFGs/EFTs_2024/prueba_picoZ/DY50_UL16.pkl.gz'#'histos_ttbar_noPUGood/TTTo2L2Nu.pkl.gz'

pathdat1 = '/nfs/fanae/user/jriego/TFGs/EFTs_2024/prueba_picoZ/SingleMuon.pkl.gz'#'histos_ttbar_PUGood/TTTo2L2Nu.pkl.gz'
pathdat2 = '/nfs/fanae/user/jriego/TFGs/EFTs_2024/prueba_picoZ/SingleElectron.pkl.gz'
#pathdat3 = '/nfs/fanae/user/jriego/TFGs/EFTs_2024/prueba_picoZ/DoubleMuon.pkl.gz'
#pathdat4 = '/nfs/fanae/user/jriego/TFGs/EFTs_2024/prueba_picoZ/MuonEG.pkl.gz'
#pathdat5 = '/nfs/fanae/user/jriego/TFGs/EFTs_2024/prueba_picoZ/DoubleEG.pkl.gz'

data=[pathdat1,pathdat2]

# Select variable, channel and cuts
var = 'invmass'


cut = '2j1b'#'base'
outname = 'kk.png'#channel=['emOS']


channel=['2los_ee_CRZ_0j']#,'2los_mm_CRZ_0j']
appl="isSR_2lOS"#"isSR_2lOS"


histogramas={}                             #Defino un diccionario donde voy a almacenar los histogramas de ese conjunto
for d in range(len(data)):
	print('Opening path: ', data[d])
	h={}
	hists={}                              #Necesitamos otro diccionario para cada .pkl.gz
	with gzip.open(data[d]) as fin:
		hin=pickle.load(fin)
		print('>> looking for histograms...')
		for k in hin.keys():
			if k in hists: hists[k]+=hin[k]
			else:          hists[k]=hin[k]
	histogramas[d]=hists[var]             #En cada key del diccionario principal almacenamos el diccionario del .pkl.gz en cuestion
	histogramas[d] = histogramas[d].integrate('channel', channel)
	histogramas[d] = histogramas[d].integrate('appl', appl)
	histogramas[d] = histogramas[d].integrate('systematic','nominal')
	histogramas[d] = histogramas[d].sum('sample')	                     #Hacemos las integraciones necesarias, de forma que en cada key del diccionario principal tengamos ya el histograma preparado para plotear solo con llamar a esa key
	
	
	

hdat=histogramas[0]
for d in range(1,len(data)):
	hdat=hdat+histogramas[d]   



hists2 = {}
with gzip.open(path2) as fin2:
  hin2 = pickle.load(fin2)
  print(' >> looking for histograms...')
  for k in hin2.keys():
    if k in hists2: hists2[k]+=hin2[k]
    else:          hists2[k]=hin2[k]



# Select the histogram var, channel and cut
h2 = hists2[var]
h2 = h2.integrate('channel', channel)
h2 = h2.integrate('appl', appl)
#h2 = h2.integrate('sumcharge', ['ch+','ch-'])
h2 = h2.integrate('systematic','nominal')
h2 = h2.sum('sample')
#h2.scale(lumi)
# Integrate over samples

from matplotlib import gridspec
fig=plt.figure(figsize=(14,12))
gs = gridspec.GridSpec(2, 1, height_ratios=[2,1],hspace=0.1)
ax0 = plt.subplot(gs[0])

#fig, axs = plt.subplots(2, sharex=True, gridspec_kw={'height_ratios': [3, 1]},figsize=(14,12)) 






hist.plot1d(h2,clear=False, line_opts={'color':'black','linewidth':3})
hist.plot1d(hdat, clear=False,line_opts={'color':'orange','linewidth':3})
#hist.plot1d(h_u,clear=False,line_opts={'color':'red','linewidth':3})

line_sf, = plt.plot([-1, -2, -3], label='Without syst', color='black')
line_nosf, = plt.plot([-3, -2, -1], label='With syst', color='orange')
line_down, = plt.plot([-3, -2, -1], label='PU DOWN', color='blue')
line_up, = plt.plot([-3, -2, -1], label='PU UP', color='red')

plt.ylabel('',fontsize=20)
plt.xlabel('',fontsize=20)
plt.xticks(fontsize=0)
plt.yticks(fontsize=22)
plt.legend(handles=[line_sf, line_nosf,line_up, line_down],fontsize=22)
ax1 = plt.subplot(gs[1])
hist.plotratio(num=hdat,denom=h2,ax=ax1,clear=False,unc='num', error_opts={'color':'red','linewidth':3,'linestyle': '-'})
#hist.plotratio(num=h_d,denom=h,ax=ax1,clear=False,unc='num', error_opts={'color':'blue','linewidth':3,'linestyle': '-'})
plt.ylabel('Syst./Nom.',fontsize=22)
plt.xlabel('PU',fontsize=22)
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)
fig.savefig(outname)
print('Output histogram saved in %s'%outname)

exit()




hep.style.use("CMS")


hep.histplot([hdat.to_hist(),h2.to_hist()],stack=False,histtype="fill",binticks=True,label=['DAta','MC'],color=['orchid','blue'])#,ax=axs[0])

hep.cms.label(llabel='Academic', data=True, lumi=40.4,year=2016,com=13)#,ax=axs[0])

#plt.yticks(fontsize=10)
#plt.tick_params(axis='x', which='both', labelsize=10)#, size=20)

plt.legend(frameon=True,framealpha=1,loc='center',fancybox=False,edgecolor='white')

plt.savefig(outname)
print('Output histogram saved in %s'%outname)

exit()
