
#wcnames=["cpt", "ctp", "cptb", "cQlMi", "cQq81", "cQq11", "cQl3i", "ctq8", "ctlTi", "ctq1", "ctli", "cQq13", "cbW", "cpQM", "cpQ3", "ctei", "cQei", "ctW", "ctlSi", "cQq83", "ctZ",  "ctG"]


#ESTE ES EL SCRIPT MODELO PARA CREAR HISTOGRAMAS CON VARIACIONES DE 6 wcs Y UNA TABLA DE YIELDS


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


#wcnames=["cpt", "ctp", "cptb", "cQlMi", "cQq81", "cQq11", "cQl3i", "ctq8", "ctlTi", "ctq1", "ctli", "cQq13", "cbW", "cpQM", "cpQ3", "ctei", "cQei", "ctW", "ctlSi", "cQq83", "ctZ",  "ctG"]
wcnames=["ctp","cQq13", "cbW", "cpQ3", "ctW", "cQq83"]

#Seleccion de canales, variable, region de aplicacion...:
cat_2lss = ['2lss_4t_m_4j', '2lss_4t_m_5j', '2lss_4t_m_6j', '2lss_4t_m_7j', '2lss_4t_p_4j', '2lss_4t_p_5j', '2lss_4t_p_6j', '2lss_4t_p_7j', '2lss_CR_1j', '2lss_CR_2j', '2lss_m_4j', '2lss_m_5j', '2lss_m_6j', '2lss_m_7j', '2lss_p_4j', '2lss_p_5j', '2lss_p_6j', '2lss_p_7j']
cat_3l = ['3l_m_offZ_1b_2j', '3l_m_offZ_1b_3j', '3l_m_offZ_1b_4j', '3l_m_offZ_1b_5j', '3l_m_offZ_2b_2j', '3l_m_offZ_2b_3j', '3l_m_offZ_2b_4j', '3l_m_offZ_2b_5j', '3l_onZ_1b_2j', '3l_onZ_1b_3j', '3l_onZ_1b_4j', '3l_onZ_1b_5j', '3l_onZ_2b_2j', '3l_onZ_2b_3j', '3l_onZ_2b_4j', '3l_onZ_2b_5j', '3l_p_offZ_1b_2j', '3l_p_offZ_1b_3j', '3l_p_offZ_1b_4j', '3l_p_offZ_1b_5j', '3l_p_offZ_2b_2j', '3l_p_offZ_2b_3j', '3l_p_offZ_2b_4j', '3l_p_offZ_2b_5j']
cat_4l = ['4l_2j', '4l_3j', '4l_4j']


var = 'Rlj_rebin'
channel = cat_3l #De momento solo exactamente 4j, luego ampliaremos para coger tambien los de mas jets
					  #Y OJO que no es nada intuitivo como coger esto. 											hay que coger 'lep_chan_lst'_'njets'

						#Sobre estas categorias, entiendo que p es carga positiva de los leptones, m negativa, y luego sin el 4t significa privados de 4 tops, y con 4t 'enriquecido' en ellos
						#COMO PONER LAS DE 4t, porque tal y como estan escritas ahi arriba no me las coge
						
appl = 'isSR_3l'#"isSR_2lSS" #'isSR_4l',

path1='histos_signal_rebinned/signal_ttH.pkl.gz'  #Cargo los datos    #aqui va ttH
#path1np='histos_signal/signal_ttH_np.pkl.gz' 
peso1=38564.22247337481   #Modificar para cada señal (ver con vi en el .json que valor coge cada uno)

path2='histos_signal_rebinned/tHq_6WC_bis.pkl.gz'
#path2np='histos_signal/signal_tHq_np.pkl.gz'
peso2=2090.4472420636366

path3='histos_signal_rebinned/signal_ttZ.pkl.gz'
#path3np='histos_signal/signal_ttZ_np.pkl.gz'
peso3=20993.713027504105

path4='histos_signal_rebinned/signal_tZq.pkl.gz'
#path4np='histos_signal_rebinned/signal_tZq_np.pkl.gz'
peso4=1933.3616499403668

path5='histos_signal_rebinned/signal_ttW.pkl.gz'
#path5np='histos_signal/signal_ttW_np.pkl.gz'
peso5=11485.9609381048

#path6='histos_signal/signal_tttt.pkl.gz'
#path6np='histos_signal/signal_tttt_np.pkl.gz'
#peso6=187.55443212752408


signals=[path3]#,path2,path3,path4,path5]#,path1np,path2,path2np]#,path2,path3,path4,path5]#th4]#,path2,path3,path4,path5]#,path6]
pesos=[peso3]#,peso1,peso2,peso2]  #OJO que ahora hay que poner 2 veces cada peso para tener en cuenta los np


histogramas={}                             #Defino un diccionario donde voy a almacenar los histogramas de ese conjunto
for d in range(len(signals)):
	print('Opening path: ', signals[d])
	h={}
	hists={}                              #Necesitamos otro diccionario para cada .pkl.gz
	with gzip.open(signals[d]) as fin:
		hin=pickle.load(fin)
		print('>> looking for histograms...')
		for k in hin.keys():
			if k in hists: hists[k]+=hin[k]
			else:          hists[k]=hin[k]
	histogramas[d]=hists[var]             #En cada key del diccionario principal almacenamos el diccionario del .pkl.gz en cuestion
	histogramas[d] = histogramas[d].integrate('channel', channel)
	#if (d % 2)==0:
	histogramas[d] = histogramas[d].integrate('appl', appl)
	histogramas[d] = histogramas[d].integrate('systematic','nominal')
	histogramas[d] = histogramas[d].sum('sample')	                     #Hacemos las integraciones necesarias, de forma que en cada key del diccionario principal tengamos ya el histograma preparado para plotear solo con llamar a esa key
	histogramas[d].scale(1./pesos[d])
	


histograma=histogramas[0]
for d in range(1,len(signals)):
	histograma=histograma+histogramas[d]       #Esto es una trampita que hay que hacer para sumar los elementos de cada grupo (en este caso con hacer h1+h2 ya directamente suma los bines, no es como el tt5TeV que habia que tirar del diccionario)



minneg, =plt.plot([-1],label='-1',color='limegreen')
sm, =plt.plot([-1],label='ME',color='k')
maxpos, =plt.plot([-1],label='1',color='darkorange')




yield_minneg=np.zeros(6)
yield_maxpos=np.zeros(6)
yield_sm=np.zeros(6)


fig1, _axs1=plt.subplots(2,3,gridspec_kw={"height_ratios": (3, 1)},sharex=True) 
fig1.subplots_adjust(hspace=.17,wspace=.417,right=0.98,left=0.095)  #Hacer mayor wspace para dejar margenes entre histogramas. Para ver cuanto usar plt.subplot_tool() (space) plt.show() (DESPUES DEL FOR)
axs1 = _axs1.flatten()


fig2, _axs2=plt.subplots(2,3,gridspec_kw={"height_ratios": (3, 1)},sharex=True) 
fig2.subplots_adjust(hspace=.17,wspace=.417,right=0.98,left=0.095)
axs2 = _axs2.flatten()


	
vector=histograma.values()
yield_sm=vector[list(vector.keys())[0]].sum()

hsm=histograma.copy() #HAGO ESTO PARA NO HACER UNA ASIGNACION DOBLE




i=-1
for j in range(len(wcnames)):
	initwc=np.zeros(22)
	if j!=0 and j != 1 and j != 2 and j != 3 and j != 4 and j !=5 :   #Esta linea es pa coger solo los 6 coefs que me interesan  CLAVE  
	
			continue													
	else:
		i=i+1
		if i<3:

			
			initwc[j]=1
			histograma.set_wilson_coeff_from_array(initwc)
			vector=histograma.values()
			yield_maxpos=vector[list(vector.keys())[0]].sum()
			histograma.scale(yield_sm/yield_maxpos)
			hist.plot1d(histograma,clear=False,line_opts={'color':'darkorange'},ax=axs1[i])
			
			


			hist.plotratio(num=histograma,denom=hsm,ax=axs1[i+3],clear=False, error_opts={'color':'darkorange','marker':'.','elinewidth':1.5},denom_fill_opts={},guide_opts={},unc='num')
			histograma.scale(yield_maxpos/yield_sm) #recuperaria el histograma con wc=1


			initwc[j]=-1
			histograma.set_wilson_coeff_from_array(initwc)
			vector=histograma.values()
			yield_minneg=vector[list(vector.keys())[0]].sum()
			histograma.scale(yield_sm/yield_minneg)
			hist.plot1d(histograma,clear=False,line_opts={'color':'limegreen'},ax=axs1[i])
			
			
			
			hist.plotratio(num=histograma,denom=hsm,ax=axs1[i+3],clear=False, error_opts={'color':'limegreen','marker':'.','elinewidth':1},denom_fill_opts={},guide_opts={},unc='num')
			histograma.scale(yield_minneg/yield_sm)#recuperaria histo con wc=-1


			hist.plot1d(hsm,clear=False,line_opts={'color':'k'},ax=axs1[i])

			axs1[i].legend(handles=[minneg,sm,maxpos],fontsize=8,title=r"$\bf{"+ wcnames[j] + "}$")  
			axs1[i].get_legend().shadow = True
			axs1[i].set_xlim(-0.5,3)
			#axs1[i].set_ylim(0,4e-9)
			axs1[i+3].set_ylim(0,3) 
			axs1[i].set_ylabel("Nº Sucesos")
			axs1[i+3].set_ylabel("EFT/ME")
			
		else:

			
			initwc[j]=1
			histograma.set_wilson_coeff_from_array(initwc)
			vector=histograma.values()
			yield_maxpos=vector[list(vector.keys())[0]].sum()
			histograma.scale(yield_sm/yield_maxpos)

			hist.plot1d(histograma,clear=False,line_opts={'color':'darkorange'},ax=axs2[i-3])
			


			hist.plotratio(num=histograma,denom=hsm,ax=axs2[i],clear=False, error_opts={'color':'darkorange','marker':'.','elinewidth':1.5},denom_fill_opts={},guide_opts={},unc='num')
			
			histograma.scale(yield_maxpos/yield_sm) #recuperaria el histograma con wc=1


			initwc[j]=-1
			histograma.set_wilson_coeff_from_array(initwc)
			vector=histograma.values()
			yield_minneg=vector[list(vector.keys())[0]].sum()
			histograma.scale(yield_sm/yield_minneg)
			hist.plot1d(histograma,clear=False,line_opts={'color':'limegreen'},ax=axs2[i-3])
			
			
			
			hist.plotratio(num=histograma,denom=hsm,ax=axs2[i],clear=False, error_opts={'color':'limegreen','marker':'.','elinewidth':0.7},denom_fill_opts={},guide_opts={},unc='num')
			histograma.scale(yield_minneg/yield_sm)#recuperaria histo con wc=-1

			hist.plot1d(hsm,clear=False,line_opts={'color':'k'},ax=axs2[i-3])

			axs2[i-3].legend(handles=[minneg,sm,maxpos],fontsize=8,title=r"$\bf{"+ wcnames[j] + "}$")  
			axs2[i-3].get_legend().shadow = True
			axs2[i].set_xlim(-0.5,3)
			#axs2[i-3].set_ylim(0,4e-9)
			axs2[i].set_ylim(0,3)
			axs2[i].set_ylabel("EFT/ME")
			axs2[i-3].set_ylabel("Nº Sucesos")

#axs1[0].set_ylim(0,3e-7)
#axs1[4].set_ylim(0,2.5)
#axs2[4].set_ylim(0.5,6)
#axs2[4].set_ylim(0,2.5)

fig1.suptitle('$ \Delta \eta_{lep-jet}$ en categoría 3l (señal tHq)',fontstyle='italic',color='navy')
fig1.savefig("Renormalizaciones/3l_tHq_incertidumbres/deltaetalj_3l_tHq_1.png") 


fig2.suptitle('$\Delta \eta_{lep-jet}$ en categoría 3l (señal tHq)',fontstyle='italic',color='navy')
fig2.savefig("Renormalizaciones/3l_tHq_incertidumbres/deltaetalj_3l_tHq_2.png") 
 
			
exit()	
	
