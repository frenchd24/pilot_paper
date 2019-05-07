 #!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)


$Id: pilot_make_histograms.py v 1.2 8/11/16

Comes from: filament_make_histograms.py, v 1.0 08/11/2015

but updated for the entire data set for the pilot paper. Makes photometry plots for the 
full galaxy table. (11/17/2015)

v1.1 Make nicer looking Texify'd plots (4/27/16)

v1.2: More formatting updates (8/11/16) DEPRECIATED.... SEE 'pilot_make_Lstar_histograms2.py'

'''

import sys
# import os
import csv
from pylab import *
# import atpy
import math
import getpass
import scipy.optimize as optimization
import pickle
import itertools
from utilities import *
from matplotlib import rc

fontScale = 18
rc('text', usetex=True)
rc('font', size=18, family='serif', weight='normal')
rc('xtick.major',size=8,width=0.6)
rc('xtick.minor',size=5,width=0.6)
rc('ytick.major',size=8,width=0.6)
rc('ytick.minor',size=5,width=0.6)
rc('xtick',labelsize = fontScale)
rc('ytick',labelsize = fontScale)
rc('axes',labelsize = fontScale)
rc('xtick', labelsize = fontScale)
rc('ytick',labelsize = fontScale)
# rc('font', weight = 450)
# rc('axes',labelweight = 'bold')
rc('axes',linewidth = 1)


def schechter(m,phi,mstar,alpha):
    # construct a Schechter luminosity function and return the associated density function
    
#     s = 0.4*log(10)*phi*(10**(0.4*(mstar-m)))**(alpha +1) * exp(-10**(0.4*(mstar-m)))
#     s = 0.4*log(10)*phi*(10**(0.4*(m-mstar)*(alpha +1))) * exp(-10**(0.4*(m-mstar)))

    s = 0.4 * math.log10(10) * phi * (10**(0.4*(mstar-m)))**(alpha +1) * exp(-10**(0.4*(mstar-m)))

    return s
    
    
def apparentMag(M,d):
    # convert M to apparent magnitude at distance d (input in Mpc)
    
    return 5*log10(d*10**6)-5+ M
    
    
    
def absoluteMag_noExtinc(m,d):
    # m is apparent magnitude, d is distance in Mpc
    M = float(m) - 5*math.log10((float(d)*10**6)/10)
    return M


def absoluteMag(m,d,e):
    # m is apparent magnitude, d is distance in Mpc, e is extinction E(B-V)
    M = float(m) - 5*math.log10((float(d)*10**6)/10) - 3.1*float(e)
    return M
    

# def extinction(m,e):
#     # m is apparent magnitude, e is extinction E(B-V)
#     
#     M = float(m) - 3.1*float(e)
#     return M
    

def lstarValue(mstar,m):
    # calculate and return L/Lstar
    lratio = 10**(-0.4*(m-mstar))
    return lratio



##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################



def make_histogram_band_apparent(d1,vlo,vhi,label,bins,saveDirectory,save):
    # make a histogram of apparent magnitudes for the data set d1
    #
    # d_gt = {'Bmedian':gt_BmedianList,'BLstar':gt_BLstarList,'ra':gt_raList,'dec':gt_decList,'dist':gt_distList,'vcorr':gt_vcorrList}
    
    
    
    if label == 'gt':
        phot = d1['Bmedian']
        err_up = array(phot)*0.01
        err_lo = array(phot)*0.01
        ra = d1['ra']
        dec = d1['dec']
        dist = d1['dist']
        vcorr = d1['vcorr']
    
    else:
        phot = d1['phot']
        err_up = d1['err_up']
        err_lo = d1['err_lo']
        ra = d1['ra']
        dec = d1['dec']
        dist = d1['dist']
        vcorr = d1['vcorr']
    
    cutPhot = []
    cutErrUp = []
    cutErrLo = []
    cutRA = []
    cutDec = []
    cutDist = []
    cutVcorr = []
    
    print 'len(phot): ',len(phot)
    print 'err_up: ',len(err_up)
    print 'err_lo: ',len(err_lo)
    print 'ra: ',len(ra)
    print 'dec: ',len(dec)
    print 'dist: ',len(dist)
    print 'vcorr: ',len(vcorr)
    
    
    for p,eup,elo,r,d,dis,v in zip(phot,err_up,err_lo,ra,dec,dist,vcorr):
        if float(v) >=vlo and float(v)<=vhi:
            print 'ra, dec: ',r,' ',d
            
            # convert to absolute mag if using Bmedian values
            if label == 'gt':
                absPhot = float(p)
                cutPhot.append(absPhot)
                cutErrUp.append(eup)
                cutErrLo.append(elo)            
                cutRA.append(float(r))
                cutDec.append(float(d))
                cutDist.append(float(dis))
                cutVcorr.append(float(v))
                
            else:
                cutPhot.append(apparentMag(float(p),float(dis)))
                cutErrUp.append(eup)
                cutErrLo.append(elo)            
                cutRA.append(float(r))
                cutDec.append(float(d))
                cutDist.append(float(dis))
                cutVcorr.append(float(v))         
            
    print 'cutVcorr: ',len(cutVcorr)

    fig = figure()
    ax = gca()
    
    errorsUp = array(cutErrUp) - array(cutPhot)
    errorsDown = array(cutPhot) - array(cutErrLo)
    
            
#         counts,bins = histogram(cutBLstar,bins=10)
#         step = abs(bins[0]-bins[1])
#         plot1 = ax.plot(bins[:-1],counts,lw=0,marker='.',ms=10)


    # do it without error bars
#         counts,returnbins = histogram(cutBLstar,bins=bins)
#         maxHeight = max(counts)+2
# 
#         plot1 = ax.hist(cutPhot,bins,histtype='bar')
#         ylim(0,maxHeight)
#         
#         title('{0}-band photometry'.format(label))
#         ax.invert_xaxis()
#         xlabel('Absolute magnitude')
#         ylabel('Number')


##########################################################################################
    # try to add error bars
    
    setbins = bins
    counts,bin_edges = histogram(cutPhot,bins=setbins)
    
    plot1 = ax.hist(cutPhot,setbins,histtype='bar')
#     ax.invert_xaxis()
    
    ylabel("Number")
    xlabel("Apparent Magnitude")
    title("Photometry in {0} band".format(label))
    
#     print 'bins: ',counts
#     maxHeight = max(counts)+2
#     ylim(0,maxHeight)
    
    
#         counts,binEdges = histogram(cutPhot,bins=bins)
#         countsUp,binEdgesUp = histogram(cutErrUp,bins=bins)
#         bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
#         width = 0.05
#         
#         maxHeight = max(counts)+2
# 
# #         plot1 = ax.hist(cutPhot,bins,histtype='bar')
# 
#         print 'bincenters: ',len(bincenters)
#         print 'cutErrUp: ',len(cutErrUp)
#         print 'binEdges: ',len(binEdges)
#         
#         
#         plot1 = ax.bar(bincenters,counts,width=width,color='b',yerr=countsUp)
#         ylim(0,maxHeight)
#         
#         title('{0}-band photometry'.format(label))
#         ax.invert_xaxis()
#         xlabel('Absolute magnitude')
#         ylabel('Number')


#     plot1 = scatter(cutPhot,errorsDown)
#     title('{0}-band photometry and errors'.format(label))
#     ax.invert_xaxis()
#     xlabel('Absolute magnitude')
#     ylabel('Errors')
    

#         ax.set_xscale('log')
    
#         norm = cm.colors.Normalize(vmax=max(cutBLstar), vmin=min(cutBLstar))
#         cmap = cm.gnuplot2
#         
#         plot1 = scatter(cutRA,cutDec, c = cutPhot, cmap = cmap,s = 30)
# 
#         
# #         plot1 = ax.plot(bins[:-1],counts,lw=0,marker='.',ms=10)
# #         ax.set_yscale('log')
# # 
#         cbar = plt.colorbar(plot1,cmap=cmap,orientation='vertical')
#         cbar.set_label('Galaxy Lstar ratios')
        
    if save:
        savefig('{0}{1}_histogram_apparent_{2}-{3}.pdf'.format(saveDirectory,label,vlo,vhi),format='pdf')
    else:
        show()



##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################


def make_histogram_band(d1,vlo,vhi,label,bins,saveDirectory,save):
    # unpack data
    # d_gt = {'Bmedian':gt_BmedianList,'BLstar':gt_BLstarList,'ra':gt_raList,'dec':gt_decList,'dist':gt_distList,'vcorr':gt_vcorrList}
    
    
    if label == 'gt':
        phot = d1['Bmedian']
        err_up = array(phot)*0.01
        err_lo = array(phot)*0.01
        ra = d1['ra']
        dec = d1['dec']
        dist = d1['dist']
        vcorr = d1['vcorr']
    
    else:
        phot = d1['phot']
        err_up = d1['err_up']
        err_lo = d1['err_lo']
        ra = d1['ra']
        dec = d1['dec']
        dist = d1['dist']
        vcorr = d1['vcorr']
    
    cutPhot = []
    cutErrUp = []
    cutErrLo = []
    cutRA = []
    cutDec = []
    cutDist = []
    cutVcorr = []
    
    print 'len(phot): ',len(phot)
    print 'err_up: ',len(err_up)
    print 'err_lo: ',len(err_lo)
    print 'ra: ',len(ra)
    print 'dec: ',len(dec)
    print 'dist: ',len(dist)
    print 'vcorr: ',len(vcorr)
    
    
    for p,eup,elo,r,d,dis,v in zip(phot,err_up,err_lo,ra,dec,dist,vcorr):
        if float(v) >=vlo and float(v)<=vhi:
        
            # convert to absolute mag if using Bmedian values
            if label == 'gt':
                extinc = 0.0126
                absPhot = absoluteMag(float(p),float(dis),extinc)
                cutPhot.append(absPhot)
                cutErrUp.append(absPhot*0.05)
                cutErrLo.append(absPhot*0.05)            
                cutRA.append(float(r))
                cutDec.append(float(d))
                cutDist.append(float(dis))
                cutVcorr.append(float(v))
                
            else:
                cutPhot.append(float(p))
                cutErrUp.append(eup)
                cutErrLo.append(elo)            
                cutRA.append(float(r))
                cutDec.append(float(d))
                cutDist.append(float(dis))
                cutVcorr.append(float(v))         
            
    print 'cutVcorr: ',len(cutVcorr)

    fig = figure()
    ax = gca()
    
    errorsUp = array(cutErrUp) - array(cutPhot)
    errorsDown = array(cutPhot) - array(cutErrLo)
    
            
#         counts,bins = histogram(cutBLstar,bins=10)
#         step = abs(bins[0]-bins[1])
#         plot1 = ax.plot(bins[:-1],counts,lw=0,marker='.',ms=10)


    # do it without error bars
#         counts,returnbins = histogram(cutBLstar,bins=bins)
#         maxHeight = max(counts)+2
# 
#         plot1 = ax.hist(cutPhot,bins,histtype='bar')
#         ylim(0,maxHeight)
#         
#         title('{0}-band photometry'.format(label))
#         ax.invert_xaxis()
#         xlabel('Absolute magnitude')
#         ylabel('Number')


##########################################################################################
    # try to add error bars
    
    setbins = bins
    counts,bin_edges = histogram(cutPhot,bins=setbins)
    
    plot1 = ax.hist(cutPhot,setbins,histtype='bar')
    ax.invert_xaxis()
    
    ylabel("Number")
    xlabel("Absolute Magnitude")
    title("Photometry in {0} band".format(label))
    
#     print 'bins: ',counts
#     maxHeight = max(counts)+2
#     ylim(0,maxHeight)
    
    
#         counts,binEdges = histogram(cutPhot,bins=bins)
#         countsUp,binEdgesUp = histogram(cutErrUp,bins=bins)
#         bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
#         width = 0.05
#         
#         maxHeight = max(counts)+2
# 
# #         plot1 = ax.hist(cutPhot,bins,histtype='bar')
# 
#         print 'bincenters: ',len(bincenters)
#         print 'cutErrUp: ',len(cutErrUp)
#         print 'binEdges: ',len(binEdges)
#         
#         
#         plot1 = ax.bar(bincenters,counts,width=width,color='b',yerr=countsUp)
#         ylim(0,maxHeight)
#         
#         title('{0}-band photometry'.format(label))
#         ax.invert_xaxis()
#         xlabel('Absolute magnitude')
#         ylabel('Number')


#     plot1 = scatter(cutPhot,errorsDown)
#     title('{0}-band photometry and errors'.format(label))
#     ax.invert_xaxis()
#     xlabel('Absolute magnitude')
#     ylabel('Errors')
    

#         ax.set_xscale('log')
    
#         norm = cm.colors.Normalize(vmax=max(cutBLstar), vmin=min(cutBLstar))
#         cmap = cm.gnuplot2
#         
#         plot1 = scatter(cutRA,cutDec, c = cutPhot, cmap = cmap,s = 30)
# 
#         
# #         plot1 = ax.plot(bins[:-1],counts,lw=0,marker='.',ms=10)
# #         ax.set_yscale('log')
# # 
#         cbar = plt.colorbar(plot1,cmap=cmap,orientation='vertical')
#         cbar.set_label('Galaxy Lstar ratios')

        
    if save:
        savefig('{0}{1}_histogram_{2}-{3}.pdf'.format(saveDirectory,label,vlo,vhi),format='pdf')
    else:
        show()


##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################


def make_histogram_lstar(d1,vlo,vhi,saveDirectory,save):
    # unpack data
    # d_gt = {'Bmedian':gt_BmedianList,'BLstar':gt_BLstarList,'ra':gt_raList,'dec':gt_decList,'dist':gt_distList,'vcorr':gt_vcorrList}
    
    phot = d1['Bmedian']
    BLstar = d1['BLstar']            
    ra = d1['ra']
    dec = d1['dec']
    dist = d1['dist']
    vcorr = d1['vcorr']
    
    
    cutPhot = []
    cutBLstar = []
    cutRA = []
    cutDec = []
    cutDist = []
    cutVcorr = []
    
    for p,b,r,d,dis,v in zip(phot,BLstar,ra,dec,dist,vcorr):
#             print 'v:', v
        if float(v) >=vlo and float(v)<=vhi:
            print 'v: ',v
            cutPhot.append(float(p))
            cutBLstar.append(log10(float(b)))
            cutRA.append(float(r))
            cutDec.append(float(d))
            cutDist.append(float(dis))
            cutVcorr.append(float(v))
            

    fig = figure(figsize=(14,10))
    ax = gca()

    print 'cutBLstar: ',cutBLstar
    print
    print 'max: ',max(cutBLstar)
    print 'min: ',min(cutBLstar)
    print 'med: ',median(cutBLstar)
    print
    
#         step = abs(bins[0]-bins[1])
#         plot1 = ax.plot(bins[:-1],counts,lw=0,marker='.',ms=10)

    setbins = arange(-3,2,0.1)
    counts,bins = histogram(cutBLstar,bins=setbins)
    
    plot1 = ax.hist(cutBLstar,setbins,histtype='bar',color='grey')
    print 'bins: ',counts
    maxHeight = max(counts)+(max(counts)*0.25)
    ylim(0,maxHeight)
    
#     ylabel(r'\rm Number',fontsize=14)
    ylabel(r'$\rm Number$')
#     xlabel(r'$\rm log_{10} (L/L_*)$',fontsize=14)
    xlabel(r'$\rm log_{10} (L/L_*)$')

#     title(r'$L_*$ Histogram for vcorr range = ({0}, {1})'.format(vlo,vhi))
    
    # x coordinate adjustment for annotations
    xAn = -0.1
    
    # label color
    lcolor = 'black'
    
    # label size
    lsize = 14
    
    # draw and annotate Lstar = 1 line
    axvline(x=log10(1),linewidth=1, color=lcolor)
    l1 = r'$L_*=1.0$'
    annotate(s=l1,xy=(log10(1)+xAn,maxHeight*0.92),size=lsize)
    
    # draw Lstar = 0.5 line
    axvline(x=log10(0.5),linewidth=1, color=lcolor)
    l_5 = r'$L_*=0.5$'
    annotate(s=l_5,xy=(log10(0.5)+xAn,maxHeight*0.87),size=lsize)
  
    # draw Lstar = 0.1 line
    axvline(x=log10(0.1),linewidth=1, color=lcolor)
    l_1 = r'$L_*=0.1$'
    annotate(s=l_1,xy=(log10(0.1)+xAn,maxHeight*0.92),size=lsize)
    
    # draw Lstar = 0.05 line
    axvline(x=log10(0.05),linewidth=1, color=lcolor)
    l_05 = r'$L_*=0.05$'
    annotate(s=l_05,xy=(log10(0.05)+xAn,maxHeight*0.87),size=lsize)
    
    # draw Lstar = 0.01 line
    axvline(x=log10(0.01),linewidth=1, color=lcolor)
    l_01 = r'$L_*=0.01$'
    annotate(s=l_01,xy=(log10(0.01)+xAn,maxHeight*0.92),size=lsize)

    
    
#         ax.set_xscale('log')

    

#         norm = cm.colors.Normalize(vmax=max(cutBLstar), vmin=min(cutBLstar))
#         cmap = cm.gnuplot2
#         
#         plot1 = scatter(cutRA,cutDec, c = cutPhot, cmap = cmap,s = 30)
#         ax.invert_xaxis()
#         
# #         hist2d(ra,dec,bins=70)
# #         plot1 = ax.scatter(cutRA,cutDec)
# 
#         ylabel('dec')
#         xlabel('ra')
#         
# #         plot1 = ax.plot(bins[:-1],counts,lw=0,marker='.',ms=10)
# #         ax.set_yscale('log')
# 
#         title('filament region galaxies')
# 
#         cbar = plt.colorbar(plot1,cmap=cmap,orientation='vertical')
#         cbar.set_label('Galaxy Lstar ratios')

    if save:
        savefig('{0}Lstar_histogram_{1}-{2}.pdf'.format(saveDirectory,vlo,vhi),format='pdf')
    else:
        show()


##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################


def make_colormap(d1,vlo,vhi,saveDirectory,save):
    # make colormap of filament region, with colors indicating brightest galaxies
    
    phot = d1['phot']
    ra = d1['ra']
    dec = d1['dec']
    dist = d1['dist']
    vcorr = d1['vcorr']
    
    cutPhot = []
    cutRA = []
    cutDec = []
    cutDist = []
    cutVcorr = []
    
    for p,r,d,dis,v in zip(phot,ra,dec,dist,vcorr):
        if float(v) >=vlo and float(v)<=vhi:
            cutPhot.append(p)
            cutRA.append(r)
            cutDec.append(d)
            cutDist.append(dis)
            cutVcorr.append(v)
            
    
    # bin up all the counts
    counts,bins = histogram(phot,bins=50)
    step = abs(bins[0]-bins[1])


    fig = figure()
    ax = gca()

    norm = cm.colors.Normalize(vmax=max(cutPhot), vmin=min(cutPhot))
    cmap = cm.RdBu
    
    plot1 = scatter(cutRA,cutDec, c = cutPhot, cmap = cmap,s = 30)
    ax.invert_xaxis()
    
#         hist2d(ra,dec,bins=70)
#         plot1 = ax.scatter(cutRA,cutDec)

    ylabel('dec')
    xlabel('ra')
    
#         plot1 = ax.plot(bins[:-1],counts,lw=0,marker='.',ms=10)
#         ax.set_yscale('log')

    title('Filament region galaxies')

    cbar = plt.colorbar(plot1,cmap=cmap,orientation='vertical')
    cbar.set_label('Galaxy {0}-band absolute magnitude'.format(dataset))

    if save:
        savefig('{0}colormap_{1}-{2}.pdf'.format(saveDirectory,vlo,vhi),format='pdf')
    else:
        show()
        
        
        
def main():
        
    '''
    The data set looks like this:
    
    d_b_j = {'phot':b_j_phot,'err_up':b_j_up_phot,'err_lo':b_j_lo_phot,'ra':b_j_raList,'dec':b_j_decList,'dist':b_j_distList,'vcorr':b_j_vcorrList}
    d_b_all = {'phot':b_all_phot,'err_up':b_all_up_phot,'err_lo':b_j_lo_phot,'ra':b_all_raList,'dec':b_all_decList,'dist':b_all_distList,'vcorr':b_all_vcorrList}
    d_j = {'phot':j_phot,'err_up':j_up_phot,'err_lo':j_lo_phot,'ra':j_raList,'dec':j_decList,'dist':j_distList,'vcorr':j_vcorrList}
    d_k = {'phot':k_phot,'err_up':k_up_phot,'err_lo':k_lo_phot,'ra':k_raList,'dec':k_decList,'dist':k_distList,'vcorr':k_vcorrList}
    d_h = {'phot':h_phot,'err_up':h_up_phot,'err_lo':h_lo_phot,'ra':h_raList,'dec':h_decList,'dist':h_distList,'vcorr':h_vcorrList}
    d_g = {'phot':g_phot,'err_up':g_up_phot,'err_lo':g_lo_phot,'ra':g_raList,'dec':g_decList,'dist':g_distList,'vcorr':g_vcorrList}
    d_i = {'phot':i_phot,'err_up':i_up_phot,'err_lo':i_lo_phot,'ra':i_raList,'dec':i_decList,'dist':i_distList,'vcorr':i_vcorrList}
    d_r = {'phot':r_phot,'err_up':r_up_phot,'err_lo':r_lo_phot,'ra':r_raList,'dec':r_decList,'dist':r_distList,'vcorr':r_vcorrList}
    d_v = {'phot':v_phot,'err_up':v_up_phot,'err_lo':v_lo_phot,'ra':v_raList,'dec':v_decList,'dist':v_distList,'vcorr':v_vcorrList}
    d_gt = {'Bmedian':gt_BmedianList,'BLstar':gt_BLstarList,'ra':gt_raList,'dec':gt_decList,'dist':gt_distList,'vcorr':gt_vcorrList}
    d_all = {'ra':all_raList,'dec':all_decList,'dist':all_distList,'vcorr':all_vcorrList}
    d_full = {'phot':full_phot,'err_up':full_errUp,'err_lo':full_errLo,'name':full_name,'source':full_source,'ra':full_raList,'dec':full_decList,'dist':full_distList,'vcorr':full_vcorrList}
                
    # add them all to the final dictionary
    d = {'d_b_j':d_b_j,\
    'd_b_all':d_b_all,\
    'd_j':d_j,\
    'd_k':d_k,\
    'd_h':d_h,\
    'd_g':d_g,\
    'd_i':d_i,\
    'd_r':d_r,\
    'd_v':d_v,\
    'd_gt':d_gt,\
    'd_all':d_all,\
    'd_full':d_full}  

    '''
    
#     raRange = (208,276)
#     decRange = (41,72)
#     make_phot_dataset(raRange,decRange)
    
    # open and read the galaxy table    
    if getpass.getuser() == 'David':
        photFilename = '/Users/David/Research_Documents/inclination/fullPhotPickle.p'
        galaxyFilename = '/Users/David/Research_Documents/gt/NewGalaxyTable5.csv'
        saveDirectory = '/Users/David/Research_Documents/inclination/git_inclination/pilot_paper_code/plots3/'
        
    elif getpass.getuser() == 'frenchd':
        photFilename = '/usr/users/frenchd/inclination/fullPhotPickle.p'
        galaxyFilename = '/usr/users/frenchd/gt/NewGalaxyTable5.csv'
        saveDirectory = '/usr/users/frenchd/inclination/git_inclination/pilot_paper_code/plots3/'
    else:
        print 'Could not determine username. Exiting.'
        sys.exit()
    
    # open the galaxy file
    theFile = open(galaxyFilename,'rU')
    reader = csv.DictReader(theFile)
    
    # open the pickle file to get at the compiled photometry file
    photFile = open(photFilename,'rU')
    d = pickle.load(photFile)
    photFile.close()
    
    
    # what to make?
    makeHistogramBand = False
    makeHistogramBandApparent = False
    makeHistogramLstar = True
    makeColormap = False
    

##########################################################################################
    if makeHistogramLstar:
        # plot a CDF for Lstar values
        dataset = 'd_gt'
        d1 = d[dataset]
        
        # velocity limits
        vlo = 0
        vhi = 10000
        
        save = True

        # make histogram plot of Lstar values
        make_histogram_lstar(d1,vlo,vhi,saveDirectory,save)
    
    
##########################################################################################
    if makeHistogramBand:
        # make histogram_band
        dataset = 'd_gt'
        d1 = d[dataset]
        
        # velocity limits
        vlo = 0
        vhi = 10000
        
        print "lstar ratio here: ",lstarValue(-19.57,-19)
    
        # make histogram plot photometry values
        label = 'gt'
        bins = 40
        
        save = False
        make_histogram_band(d1,vlo,vhi,label,bins,saveDirectory,save)
    
    
##########################################################################################
    if makeHistogramBandApparent:
        # make histogram_band
        dataset = 'd_gt'
        d1 = d[dataset]
        
        # velocity limits
        vlo = 0
        vhi = 10000
    
        # make histogram plot photometry values
        label = 'gt'
        bins = 40
        
        save = False
        make_histogram_band_apparent(d1,vlo,vhi,label,bins,saveDirectory,save)
    

##########################################################################################
    if makeColormap:
        # make colormap of photometry values
        dataset = 'd_b_j'
        d1 = d[dataset]
    
        label = 'Photometry values'
        bins = 25

        # velocity limits
        vlo = 0
        vhi = 10000
        
        save = False
        
        make_colormap(d1,vlo,vhi,saveDirectory,save)
    

    # close the galaxy table
    theFile.close()



if __name__=="__main__":
    main()
    
    