import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import glob
import os
import re
import sys
from matplotlib.pyplot import figure
from scipy.signal import find_peaks
from scipy import integrate
import natsort
from matplotlib import cm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes



def main():
 
    #dirHere = os.getcwd()
    #print(dirHere)

    #run parameters. TODO: Fix for kR,kS -> bigK.
    sze = 120.0
    r1 = 54.5 
    r2 = 58.5

    rStrt = 109.5
    kGrp = 30.0
    tEnd = 50.0

    #vG = -2.0*np.sin(2.0*np.pi*(sze-kS)/(1.*sze))
    #tMe = np.abs(8.0*( ( (sze-1) + (r1+0.5-s0) )/ (1.*vG) ))
    #MeEdge = ( ( (sze-1) - s0 ) / (1.*vG) )

    #graphics settings.
    alfa=0.75
  

    
    #for dirW in sorted(glob.glob('e0.10E*K30.00')):
    for dirW in sorted(glob.glob('e0.00E*K30.00')):

        print(dirW)
        eNval = float(dirW[1:5])
        print(eNval)
        bigKval = float(dirW[-5:-1])
        print(bigKval)
        typVal = dirW[5]
        print(typVal)


        #Go through the output data, one directory at a time.
        #Go into the subdirectory.
        os.chdir(dirW)
        #print(dirW)


        #Manipulations of the probabilities for visualization.
        os.chdir("probs")
        print()
        print("Using probs from Seed #0")
        print()
        #dirHere = os.getcwd()
        #print(dirHere)

        #Add a time column to the prob files.
        for pFileR in sorted(glob.glob('prob_*.dat')):
            if pFileR.endswith('wTime.dat'):
                break
            else: 
                #print(pFileR)
                tme = float(pFileR[5:-4])
                #print(tme)
                #filepath=dirCreate+'/'+pFileR
                colnames=['indxR','probR']
                datProbR = pd.read_csv(pFileR, delimiter='\t',names=colnames)
                datProbR['time'] = tme
                #print(datProbR)

                #datProbR['rNew'] = datProbR.indxR
                sortd_dfR = datProbR.sort_values(by=['indxR'], ascending=True)

                prefix = pFileR[:-4]
                suffix = ".dat"
                addendum = "_wTime"
                filepath = prefix + addendum + suffix
                #print(filepath)
                #ADDS THE TIME COLUMN. Should do this in main.cpp
                sortd_dfR.to_csv(filepath, sep='\t', index=False, header=None)
        probRfiles = []

        for pFileR in sorted(glob.glob('prob_*wTime.dat')):
          #  print(pFileR)
            probRfiles.append(pFileR)
        probRfilesSorted = natsort.natsorted(probRfiles)
        #print(len(probRfilesSorted))

        xRarrAll = []
        tarrAll = []
        PsiArrAll = []
        colnames=['indxR','probR','time']
        for pFileR in probRfilesSorted:
            #print(pFileR)
            datProbR = pd.read_csv(pFileR, delimiter='\t',names=colnames)
            xRarr = np.array(datProbR.indxR)
            xRarrAll.append(xRarr)
            tarr = np.array(datProbR.time)
            tarrAll.append(tarr)
            psiarr = np.array(datProbR.probR)
            #psiarr = psiarr[:400]
            PsiArrAll.append(psiarr)
        xRr=np.vstack(xRarrAll)
        yTr=np.vstack(tarrAll)
        zPsiR=np.vstack(PsiArrAll)
        tMe = np.abs(yTr).max()

        plt.rcParams["font.family"] = "serif"
        plt.rcParams["font.serif"] = ["Times New Roman"] + plt.rcParams["font.serif"]


        fig, ax1 = plt.subplots(ncols=1, nrows=1, figsize=(7,6)) #, subplot_kw={'projection': '3d'})

        #print(np.shape(xRr),np.shape(yTr),np.shape(zPsiR))
        #exit()


        y = xRr #arrays of dye indices.
        x = yTr #arrays of times.
        z = zPsiR
        #print(len(z),len(x),len(y))
        #print(z)

        # x and y are bounds, so z should be the value *inside* those bounds.
        # Therefore, remove the last value from the z array.
        #z = z[:-1, :-1]
        z_min, z_max = 0., np.abs(z).max()

        x2 = [r1, r1]
        y2 = [0, tMe]
        ax1.plot(x2, y2, color="white", linewidth=0.5)
        x2 = [r2, r2]
        y2 = [0, tMe]
        ax1.plot(x2, y2, color="white", linewidth=0.5)

        #rEnd = rStrt-2.0*np.sin(2*np.pi*((1.0*kGrp)/sze))*tEnd
        rEnd = rStrt-2.0*tEnd
        x3 = [rStrt, rEnd]
        y3 = [0, tEnd]
        #print(rEnd)
     #   ax1.plot(x3, y3, color="black", linewidth=0.75)

        c1 = ax1.pcolormesh(y, x, z, cmap='gist_ncar', alpha=alfa)

        fontsze=14
        ax1.set_xlabel('Dye index r', linespacing=4, fontsize=fontsze)
        ax1.set_ylabel(r'Time $\tau = \omega_J t$',fontsize=fontsze)

        # x and y are bounds, so z should be the value *inside* those bounds.
        # Therefore, remove the last value from the z array.
        z = z[:-1, :-1]
        z_min, z_max = 0., np.abs(z).max()

        c = ax1.pcolormesh(y, x, z, cmap='gist_ncar', alpha=alfa)

        ax1.tick_params(axis='x', labelsize=12)
        ax1.tick_params(axis='y', labelsize=12)

        #labelW21 = "delta = {}".format(eNval)
        #labelWK = "K = {}".format(1.*bigKval)
        #labelT = "Sim Typ = {}".format(typVal)
        #ax1.text(9, 15, labelW21, fontweight="bold", fontsize=17,color="white")
        #ax1.text(9, 12.5, labelWK, fontweight="bold", fontsize=17,color="white")
        #ax1.text(9, 17.5, labelT, fontweight="bold", fontsize=17,color="white")



        ax1.set_xlim(0,sze-1)
        ax1.set_ylim(0,tMe)


        fg_color = 'white'
        bg_color = 'black'

        #cbaxes = inset_axes(ax1, width="4%", height="30%", loc=3) 
        #cb1=fig.colorbar(c1, cax=cbaxes, ticks=[0.,0.01,0.02,0.03,0.04], orientation='vertical')

        # set colorbar tick color
        #cb1.ax.yaxis.set_tick_params(color=fg_color)

        # set colorbar edgecolor 
        #cb1.outline.set_edgecolor(fg_color)

        # set colorbar ticklabels
        #plt.setp(plt.getp(cb1.ax.axes, 'yticklabels'), color=fg_color)


        fig.tight_layout()
        plt.show()
        #plt.savefig("../probsRandSbound_w0.70.png")
    

        os.chdir("../../")

    #exit()
                 

def sums_and_such():

    exit()

    return


if __name__ == '__main__':
    main()
