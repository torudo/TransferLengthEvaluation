"""Author: Toby Rudolph
Date: 09/07/2021
To create linear Fits of TLM measurements
You need these modules installed
You can install a module by using "pip install modulename"
"""


import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator #for ticks
import matplotlib.font_manager as font_manager
import matplotlib.font_manager
from scipy import optimize
from scipy import stats

newparams = {
#from pylab import plot, show, savefig, xlim, figure, hold, ylim, legend, boxplot, setp, axes
      'font.size': 20,
      'figure.figsize': (8, 7),
      'lines.linewidth': 2.5,
      'lines.markersize' : 10
      }

plt.rcParams["font.family"] = "Arial"
matplotlib.rcParams['font.serif'] = 'Palatino'
plt.rcParams['text.latex.preamble'] = [r'\usepackage{siunitx}']
plt.rcParams.update(newparams)

def Set_Stuff():
    ax = plt.gca()
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.tick_params(which='major', direction='in',bottom=True, top=False, left=True, right=True, labelsize = 20, length=8, width=2)
    ax.tick_params(which='minor', direction='in',bottom=False, top=False, left=True, right=True, length=4.5, width=1.25)#,color = '#0f0f0f50')

def linfunc(x, m, b):
    return m * x + b

"""Aktuelle Datenausfgabe """
def Plots_single(filepath, file,ausgabepath, Treat = "BLS", Type="n type", wannaprint = 0):
    #print(file)
    plotfile = filepath+r'\\'+file # filename_dir[1]

    mydir = ausgabepath+"/"+Type
    CHECK_FOLDER = os.path.isdir(mydir)
    if not CHECK_FOLDER:
        os.makedirs(mydir)
        print("created folder : ", mydir)
    else:
        print(mydir, "\nfolder already exists.")

    df = pd.read_csv(plotfile, delimiter = r'#',skiprows=[0,23,24,25,26,27,28,29,30] )
    cols = df.columns
    slopes = []
    intercept = []

    fig, axes = plt.subplots(ncols=2, nrows=4,sharex=True, sharey=True,constrained_layout=True,figsize=(12,9))

    for i in range(0,len(axes)):
        for j in range(0,len(axes[0])):
            #print(axes[i,j])
            axes[i,j].set_axis_off()

    j = 0
    k = 0
    for count, i in enumerate(cols[1:]):
        ax = axes[k,j]
        ax.set_axis_on()
        x1,y1 = df['Messung'], df[i]
        ax.plot(x1,y1, marker="x", linewidth = 0)
        x = np.linspace(-0.04,0.04,100)
        popt, pcov = optimize.curve_fit(linfunc,x1,y1, p0=[1, 1.])
        ax.plot(x, linfunc(x, popt[0], popt[1]), label='fit')
        pcov = np.sqrt(np.diag(pcov))
        ax.set_title(i)
    #    print("Slope = {0:1.3f} +-{1:1.3f}".format(popt[0],pcov[0]))
        slopes.append(popt[0])
        intercept.append(popt[1])
        #print(k,j)
        if k<=2:
            k+=1
        else:
            k = 0
            j +=1
        if (k == 2& j ==1):
            break
        #plt.tight_layout()
        plt.title(file)
        image_format = 'png' # e.g .png, .svg, etc.
        image_name = '.png'
#        fig.savefig(ausgabepath+"/"+file[:-4]+image_name, format=image_format)#, dpi=1200, transparent=True)
        #plt.show()
#        plt.close()


    fig, ax = plt.subplots(figsize=(12,9))
    #print(slopes)
    x1 = [476,594,688,801,910,1013,1110]
    slopes = np.array(slopes)
    intercept = np.array(intercept)
    y1 = 1/slopes
    #print(y1)

    plt.plot(x1,y1,linewidth=0,marker='x')
    x = np.linspace(500,1100,100)
    popt, pcov = optimize.curve_fit(linfunc,x1,y1, p0=[0.001, 1.])
    #print(popt, pcov)
    plt.plot(x, linfunc(x, popt[0], popt[1]), label='fit')
    pcov = np.sqrt(np.diag(pcov))
    # print(pcov)
    #print("Slope = {0:1.3f} +-{1:1.3f}".format(popt[0]*1000,pcov[0]))
    plt.title(file[:-4])
    plt.xlabel('x-distance')
    plt.ylabel('1/R - from slopes')

    slope_r_x = popt[0]*1000
    intercept_r_x = popt[1]
    Rc = intercept_r_x / 2
    Lt = intercept_r_x / (2*slope_r_x)*1000
    rho = Rc*Lt*1/10000

    rhoso = intercept_r_x**2/(slope_r_x*40)
    m = slope_r_x
    deltam = pcov[0]
    deltab = pcov[1]
    b = intercept_r_x
    rhoerror = np.sqrt(( deltab*2*b/(40*m) )**2+(deltam*b**2/(40*m**2))**2)

    if wannaprint == 1:
        # print('inverse slopes')
        # print(y1)
        # print(slope_r_x)
        # print(intercept_r_x)
        # print(Rc)
        # print(Lt)
        print(rhoso)
        print(rhoerror)
        print("rho={0:1.4f} +/- {1:1.7f}".format(rhoso,rhoerror))
        #print(rho)

    #print("Steigung, Achenabschnitt:",beta, betastd, chiq)
    result = "y=m*x+b \n m=({0:1.4f})\n y={1:1.6f} pm ".format(popt[0],popt[1])
    #''+'m=('+str(np.round(popt[0],3))+r'$\pm$'+str(np.round(pcov[0],3))+') \n'+'b=('+str(np.round(popt[1],2))+r'$\pm$'+str(np.round(pcov[1],2))+') \n')
    plt.text(0.1, 0.9,result,bbox=dict(boxstyle="square",ec='k',fc='w'),fontsize=13, ha='center', va='center', transform=ax.transAxes)

    file1 = open(ausgabepath+r"\\"+Treat+"_"+Type+"_rho_data.txt","a")
    file1.write(file[:-4]+"\t {0:1.5f}\n".format(rho))
    file1.close()

    plt.tight_layout()
    image_format = 'png' # e.g .png, .svg, etc.
    image_name = '.png'
#    fig.savefig(mydir+"/"+Treat+"_R_over_X_"+file[:-4]+image_name, format=image_format)#, dpi=1200, transparent=True)
    #plt.show()
#    plt.close()

# For all Data
def Plots_All(filepath,ausgabepath, Treat = "BLS", Type="n type", wannaprint = 0):
    filename_dir = os.listdir(filepath)
    i = 1

    mydir = ausgabepath+"/"+Type
    CHECK_FOLDER = os.path.isdir(mydir)
    if not CHECK_FOLDER:
        os.makedirs(mydir)
        print("created folder : ", mydir)
    else:
        print(mydir, "\nfolder already exists.")

    for i in range(0,len(filename_dir)):
    #for i in range(0,1):
        file = filename_dir[i]
        #print(file)
        plotfile = filepath+r"\\"+file # filename_dir[1]
        #print(filename_dir)
        #print(plotfile)
        df = pd.read_csv(plotfile, delimiter = r'#',skiprows=[0,23,24,25,26,27,28,29,30] )
        cols = df.columns
        slopes = []
        intercept = []

        fig, axes = plt.subplots(ncols=2, nrows=4,sharex=True, sharey=True,constrained_layout=True,figsize=(12,9))

        for i in range(0,len(axes)):
            for j in range(0,len(axes[0])):
                #print(axes[i,j])
                axes[i,j].set_axis_off()

        j = 0
        k = 0
        for count, i in enumerate(cols[1:]):
            ax = axes[k,j]
            ax.set_axis_on()
            x1,y1 = df['Messung'], df[i]
            ax.plot(x1,y1, marker="x", linewidth = 0)
            x = np.linspace(-0.04,0.04,100)
            popt, pcov = optimize.curve_fit(linfunc,x1,y1, p0=[1, 1.])
            ax.plot(x, linfunc(x, popt[0], popt[1]), label='fit')
            pcov = np.sqrt(np.diag(pcov))
            ax.set_title(i)
        #    print("Slope = {0:1.3f} +-{1:1.3f}".format(popt[0],pcov[0]))
            slopes.append(popt[0])
            intercept.append(popt[1])
            #print(k,j)
            if k<=2:
                k+=1
            else:
                k = 0
                j +=1
            if (k == 2& j ==1):
                break
        #plt.tight_layout()
        plt.title(file)
        image_format = 'png' # e.g .png, .svg, etc.
        image_name = '.png'
        fig.savefig(ausgabepath+"/"+file[:-4]+image_name, format=image_format)#, dpi=1200, transparent=True)
        mydir+"/"+Treat+"_R_over_X_"+file[:-4]+image_name
        #plt.show()
        plt.close()

        # From here the evaluation of the slopes begins
        fig, ax = plt.subplots(figsize=(12,9))
        #print(slopes)
        x1 = [476,594,688,801,910,1013,1110]
        slopes = np.array(slopes)
        intercept = np.array(intercept)
        y1 = 1/slopes

        plt.plot(x1,y1,linewidth=0,marker='x')
        x = np.linspace(500,1100,100)
        popt, pcov = optimize.curve_fit(linfunc,x1,y1, p0=[0.001, 1.])
        plt.plot(x, linfunc(x, popt[0], popt[1]), label='fit')
        pcov = np.sqrt(np.diag(pcov))
        #print("Slope = {0:1.3f} +-{1:1.3f}".format(popt[0]*1000,pcov[0]))
        plt.title(file[:-4])
        plt.xlabel('x-distance')
        plt.ylabel('1/R - from slopes')

        slope_r_x = popt[0]*1000
        intercept_r_x = popt[1]
        Rc = intercept_r_x / 2
        Lt = intercept_r_x / (2*slope_r_x)*1000
        rho = Rc*Lt*1/10000

        if wannaprint == 1:
            print('inverse slopes')
            print(y1)
            print(slope_r_x)
            print(intercept_r_x)
            print(Rc)
            print(Lt)
            print(rho)
        result = "y=m*x+b \n m=({0:1.3f})\n y={1:1.3f}".format(popt[0],popt[1])
        plt.text(0.1, 0.9,result,bbox=dict(boxstyle="square",ec='k',fc='w'),fontsize=13, ha='center', va='center', transform=ax.transAxes)

        file1 = open(ausgabepath+r"\\"+Treat+"_"+Type+"_rho_data.txt","a")
        file1.write(file[:-4]+"\t {0:1.5f}\n".format(rho))
        file1.close()

        plt.tight_layout()
        image_format = 'png' # e.g .png, .svg, etc.
        image_name = '.png'
        fig.savefig(mydir+"/"+Treat+"_R_over_X_"+file[:-4]+image_name, format=image_format)#, dpi=1200, transparent=True)
        #plt.show()
        plt.close()



"""Make changes here :) """

ausgabepath = r"..."
filepath = r"..."
treatment = "filepath.split("\\")[-1]"
Plots_All(filepath, ausgabepath, Treat = treatment, Type=treatment+"Plots", wannaprint = 0)
