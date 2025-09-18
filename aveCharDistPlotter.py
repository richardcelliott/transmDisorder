import matplotlib.pyplot as plt
import numpy as np
np.bool=np.bool_
np.int = np.int32
np.float = np.float64
import pandas as pd
from matplotlib.ticker import FuncFormatter, MultipleLocator

def main():

    fig, [ax1,ax2] = plt.subplots(1,2,figsize=(12,8))

    #colnames0=["wval","aveVal","stdev"]
    df_all = pd.read_csv('noiseEstats.csv',sep="\t") #,names=colnames0)
    #print(df_all)

    #ax1.plot([0.0, 1.01], [1.0, 1.0], color='k', linestyle=':', linewidth=1)
    #ax1.plot([0.0, 1.01], [0.0, 0.0], color='k', linestyle=':', linewidth=1)
    ax1.errorbar(df_all['wVal'],df_all['charLength0'], yerr=df_all['errCharLen'], label='Ave. characteristic Length $\mathregular{d_C}$',fmt='.', markersize='8', ecolor='dodgerblue',capsize=4, elinewidth=4)
    ax2.errorbar(df_all['aveVal'],df_all['charLength0'], yerr=df_all['errCharLen'], label=r'Ave. char, Length $\mathregular{d_C}$ vs. $T$',fmt='.', markersize='8', ecolor='darkgreen',capsize=4, elinewidth=4, color='darkblue')
    #ax1.plot(df_all['wVal'],1.0-(df_all['wVal']*df_all['wVal']), label='$\mathregular{1-\lambda^2}$', color='green', linewidth=2)

    #df_all = pd.read_csv('noiseJstats.csv',sep="\t") #,names=colnames0)
    #plt.errorbar(df_all['wVal'],df_all['aveVal'], yerr=df_all['stdev'], label='Noise in $\mathregular{J_r}$',fmt='.', markersize='10', ecolor='darkorange',capsize=4, elinewidth=1.5,color='darkorange')
    #plt.plot(df_all['wVal'],1.0-(3.14*df_all['wVal']*df_all['wVal']), label='$\mathregular{1-\pi\lambda^2}$', color='red', linewidth=2)

    ax1.tick_params(axis='both', which='major', labelsize=14)
    ax2.tick_params(axis='both', which='major', labelsize=14)

    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax2.set_yscale("log")

    #ax1.gca().set_xlabel(r'$\lambda$')
    ax1.set_xlabel(r'$\log(\sigma_N)$',fontsize=18)
    ax1.set_ylabel(r'$\log(d_C)$',fontsize=18, weight='bold')
    ax2.set_xlabel(r'$T$',fontsize=18)
    ax2.set_ylabel(r'$\log(d_C)$',fontsize=18, weight='bold')

    ax1.legend(loc='upper right')
    ax2.legend(loc='upper left')



    ax1.set_xlim([0.0, 0.68])
    ax1.set_ylim([0.0, 22050.0])



    #ax2.yaxis.set_major_locator(plt.MultipleLocator(np.pi / 6))
    #ax2.yaxis.set_minor_locator(plt.MultipleLocator(np.pi / 12))
    #ax2.yaxis.set_major_formatter(plt.FuncFormatter(multiple_formatter()))

    ax2.set_xlim([0, 1.01])
    ax2.set_ylim([0.0, 5400.0])

    ax1.text(0.04,15,'One line, noise in site E, n=10',weight='bold')

    #adjust y-axis label position
    #x2.yaxis.set_label_coords(-.13, .5)

    #adjust x-axis label position 
    #ax2.xaxis.set_label_coords(.5, -.14)


    plt.show()
    #plt.savefig('siteEnTransm_oneLine.png')


def multiple_formatter(denominator=6, number=np.pi, latex='\pi'):
    def gcd(a, b):
        while b:
            a, b = b, a%b
        return a
    def _multiple_formatter(x, pos):
        den = denominator
        num = np.int(np.rint(den*x/number))
        com = gcd(num,den)
        (num,den) = (int(num/com),int(den/com))
        if den==1:
            if num==0:
                return r'$0$'
            if num==1:
                return r'$%s$'%latex
            elif num==-1:
                return r'$-%s$'%latex
            else:
                return r'$%s%s$'%(num,latex)
        else:
            if num==1:
                return r'$\frac{%s}{%s}$'%(latex,den)
            elif num==-1:
                return r'$\frac{-%s}{%s}$'%(latex,den)
            else:
                return r'$\frac{%s%s}{%s}$'%(num,latex,den)
    return _multiple_formatter

class Multiple:
    def __init__(self, denominator=6, number=np.pi, latex='\pi'):
        self.denominator = denominator
        self.number = number
        self.latex = latex

    def locator(self):
        return plt.MultipleLocator(self.number / self.denominator)

    def formatter(self):
        return plt.FuncFormatter(multiple_formatter(self.denominator, self.number, self.latex))

if __name__ == '__main__':
    main()

