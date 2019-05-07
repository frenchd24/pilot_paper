#!/usr/bin/python

''' plots N(HI) image for ATCA HCO+ absorption sightlines
#originally written by EZBC, modified by CEM (10/14)
'''

from math import *
import pyfits as pf
import numpy as np

''' Plotting Functions
'''

def plot_image(image=None, header=None, title=None, 
        savedir='./', sources=None, names=None, colormap=None, filename=None, show=True):

    # Import external modules
    import matplotlib.pyplot as plt
    import matplotlib
    import numpy as np
    from mpl_toolkits.axes_grid1 import ImageGrid
    import pyfits as pf
    import matplotlib.pyplot as plt
    import pywcsgrid2 as wcs
    import pywcs
    from pylab import cm # colormaps
    from matplotlib.patches import Polygon

    # Set up plot aesthetics
    plt.clf()
    plt.rcdefaults()

    #color_cycle = [colormap(i) for i in np.linspace(0, 0.9, len(flux_list))]
    fontScale = 15
    params = {#'backend': .pdf',
              'axes.labelsize': fontScale*5/6,
              'axes.titlesize': fontScale,
              'text.fontsize': fontScale,
              'legend.fontsize': fontScale*1/2,
              'xtick.labelsize': fontScale,
              'ytick.labelsize': fontScale,
              'font.weight': 500,
              'axes.labelweight': 500,
              'text.usetex': False,
              'figure.figsize': (8,7),
              'figure.titlesize': fontScale
              #'axes.color_cycle': color_cycle # colors of different plots
             }
    plt.rcParams.update(params)

    # Create figure instance
    fig = plt.figure()

    nrows_ncols=(1,1)
    ngrids=1

    imagegrid = ImageGrid(fig, (1,1,1),
                 nrows_ncols=nrows_ncols,
                 ngrids=ngrids,
                 axes_pad=-10,
                 axes_class=(wcs.Axes,
                             dict(header=header)),
                 aspect=True,
                 label_mode='L',
                 share_all=True)

    # create axes
    ax = imagegrid[0]

    im = ax.imshow(image,
            interpolation='nearest',origin='lower',
            cmap=colormap,
	    norm=matplotlib.colors.LogNorm()
            )

    # Asthetics
    ax.set_display_coord_system("gal")
    ax.set_ticklabel_type("dms","dms")
    ax.set_xlabel('l (degrees)',)
    ax.set_ylabel('b (degrees)',)
    ax.locator_params(axis="x", nbins=20)
    ax.locator_params(axis="y", nbins=10)
    #ax.yaxis.set_ticks([-20,-30,-40,-50,-60,-70,-80])
    ax.swap_tick_coord()
    #ax.minorticks_on()

    ax.set_xlim(0,196)
    ax.set_ylim(0,140)

    cmap=cm.Greys
    cmap.set_bad(color='white')

    ax.grid(color='black',linewidth=0.6)
    wcs_header = pywcs.WCS(header)

    #--young stellar objects, from Casetti-Dinescu et al. 2014
    yso_lon=[275.3, 277.7, 288.8, 290.8, 300.7, 309.4, 298.8]
    yso_lat=[10.5, 13.1, 12.2, 9.6, -11.8, -14.1, -13.9]
    yso_names=['B02','B03','B14','B15','A15','A19','A08']
    for i in enumerate(yso_names):
	val=i[0]
	yso_name=i[1]
        source = wcs_header.wcs_sky2pix([[yso_lon[val],yso_lat[val], 0]], 0)[0]
        #casetti_dinescu = ax.scatter(source[0],source[1], s=80, c='green', marker='+')
        #ax.annotate(yso_name, xy=(source[0]+1,source[1]), size=10, color='green')

    #---HI absorption non-detections from Matthews et al. 2009
    matthews_lon=[301.878598, 300.077311, 294.192187]
    matthews_lat=[-67.888957, -49.554906, -48.034732]
    matthews_names=['J0053-4914', 'J0110-6727', 'J0154-6800']
    for i in enumerate(matthews_names):
	val=i[0]
	matthews_name=i[1]
        source = wcs_header.wcs_sky2pix([[matthews_lon[val],matthews_lat[val], 0]], 0)[0]
        matthews = ax.scatter(source[0],source[1], s=90, facecolor='cyan', edgecolor='none', marker='v')
    #---and detections
    matthews_lon=[ 298.953217]
    matthews_lat=[-48.749733]
    matthews_names=['J0119-6809']
    for i in enumerate(matthews_names):
	val=i[0]
	matthews_name=i[1]
        source = wcs_header.wcs_sky2pix([[matthews_lon[val],matthews_lat[val], 0]], 0)[0]
        ax.scatter(source[0],source[1], s=100, facecolor='cyan', edgecolor='black', marker='v')
        #ax.annotate(matthews_name, xy=(source[0]+1,source[1]), size=10, color='blue')

    #---HI absorption non-detections from Kobulnicky & Dickey 1999
    kd_lon=[316.976119, 310.5204811,  291.809271,  270.5504581, 286.3683817, 279.4727650]
    kd_lat=[-45.952796, -47.9689541,  -42.351407,  -36.0720093, -27.1582861, -20.1301744]
    kd_names=['B2300-683', 'B2353-686', 'B0242-729', 'B0506-612', 'B0637-752', 'B0743-673']
    for i in enumerate(kd_names):
	val=i[0]
	kd_name=i[1]
        source = wcs_header.wcs_sky2pix([[kd_lon[val],kd_lat[val], 0]], 0)[0]
        kob_dickey = ax.scatter(source[0],source[1], s=90, facecolor='cyan', edgecolor='none', marker='v')
    #---and detections:
    kd_lon=[297.5486000, 293.4395878]
    kd_lat=[-40.0450765, -37.5524697]
    kd_names=['B0202-765', 'B0312-770']
    for i in enumerate(kd_names):
	val=i[0]
	kd_name=i[1]
        source = wcs_header.wcs_sky2pix([[kd_lon[val],kd_lat[val], 0]], 0)[0]
        ax.scatter(source[0],source[1], s=100, facecolor='cyan', edgecolor='black', marker='v')
        #ax.annotate(kd_name, xy=(source[0]+1,source[1]), size=10, color='cyan')


    #---HI absorption detections from Dickey et al. 2000 (SMC) *** There is also a list of non-detections
    d_lon=[304.126, 303.088, 302.839, 302.701, 302.812, 302.379, 302.365, 301.915, 300.958, 300.997, 300.689, 300.649, 300.437]
    d_lat=[-42.7280, -44.5292, -44.6758, -44.5983, -43.9154, -46.0042, -42.9093, -45.4698, -44.6100, -43.8142, -43.6281, -43.6637, -43.7625]
    d_names=['J003824', 'J004956', 'J005218', 'J005238', 'J005337', 'J005611', 'J005732', 'J010029', 'J011005', 'J011049', 'J011412', 'J011432', 'J011628']
    for i in enumerate(d_names):
	val=i[0]
	d_name=i[1]
        source = wcs_header.wcs_sky2pix([[d_lon[val],d_lat[val], 0]], 0)[0]
        #dickey_smc = ax.scatter(source[0],source[1], s=30, c='m', marker='s')
        #ax.annotate(d_name, xy=(source[0]+1,source[1]), size=10, color='m')

    #---CO clouds from Mizuno et al. 2006
    mizuno_lon=[298.273, 297.818, 297.043, 296.816, 296.452, 295.594, 295.375, 295.335]
    mizuno_lat=[-42.1125, -41.9165, -42.0779, -42.0071, -41.4252, -41.5586, -41.7980, -41.7837]
    mizuno_names=['A','B','C','D','E','F','G','H']
    for i in enumerate(mizuno_names):
	val=i[0]
	mizuno_name=i[1]
        source = wcs_header.wcs_sky2pix([[mizuno_lon[val],mizuno_lat[val], 0]], 0)[0]
        mizuno_co = ax.scatter(source[0],source[1], s=60, edgecolor='magenta', facecolor='none', marker='x')
        #ax.annotate(mizuno_name, xy=(source[0]+1,source[1]), size=10, color='m')


    #---CO cloud from Muller et al. 2003
    muller_lon=297.04417
    muller_lat=-42.078172
    source = wcs_header.wcs_sky2pix([[muller_lon,muller_lat, 0]], 0)[0]
    muller_co = ax.scatter(source[0],source[1], s=60, edgecolor='magenta', facecolor='none', marker='x')
    #ax.annotate('CO Clouds', xy=(source[0]+1,source[1]), size=10, color='yellow')

    #---HI shells from Muller et al. 2003
    #import pyregion
    #import pyfits
    #reg = pyregion.open('/d/leffe2/cmurray/data/muller_2003_shells.reg') # ShapeList
    #print reg[0].coord_list
    #reg2 = reg[0].as_imagecoord(header)
    #msk = reg[0].get_mask(shape=im.shape)
    #ax.imshow(msk, origin="lower", cmap="gray")

    #---H2 in dense filament: Richter et al. 2003 (logN(H2)=16.6\pm0.5)
    richter_lon=279.33
    richter_lat=-32.79
    richter_name='HD 36521'
    source = wcs_header.wcs_sky2pix([[richter_lon,richter_lat, 0]], 0)[0]
    #richter_03_h2 = ax.scatter(source[0],source[1], s=80, c='green', marker='h')
    #ax.annotate(richter_name, xy=(source[0]+1,source[1]), size=10, color='green')

    #---H2 in Leading Arm: Sembach et al. 2000 (logN(H2)=16.8\pm0.1) ***off the plot
    richter_lon=287.46
    richter_lat=22.95
    richter_name='HVC 287.5 +22.95'
    source = wcs_header.wcs_sky2pix([[richter_lon,richter_lat, 0]], 0)[0]
    sembach_h2 = ax.scatter(source[0],source[1], s=80, facecolor='green', edgecolor='none', marker='H')
    #ax.annotate(richter_name, xy=(source[0]+1,source[1]), size=10, color='green')

    #---H2 in towards Fairall 9: Richter et al. 2001 (logN(H2)=16.4\pm0.5) 
    richter_lon=295.1
    richter_lat=-57.8
    richter_name='Fairall 9'
    source = wcs_header.wcs_sky2pix([[richter_lon,richter_lat, 0]], 0)[0]
    richter_01_h2 = ax.scatter(source[0],source[1], s=80, c='green', marker='H')
    #ax.annotate(richter_name, xy=(source[0]+1,source[1]), size=10, color='green')

    #---High-res HI bridge footprint from Muller et al. 2003
    mul_hi_lon=[289.329, 293.552, 296.460, 293.445]
    mul_hi_lat=[-42.569, -44.575, -40.060, -38.664]
    for i in enumerate(mul_hi_lon):
	val=i[0]
	source = wcs_header.wcs_sky2pix([[mul_hi_lon[val],mul_hi_lat[val], 0]], 0)[0]
	#muller_03_hi = ax.scatter(source[0],source[1], s=20, edgecolor='purple', facecolor='purple', marker='s')  

    #---Shapley 1940 young blue stars
    shapley_lon=295.34676
    shapley_lat=-41.874962
    source = wcs_header.wcs_sky2pix([[shapley_lon,shapley_lat, 0]], 0)[0]
    #shapley = ax.scatter(source[0],source[1], s=120, c='blue',  marker='+')  
    #ax.annotate('Shapley blue', xy=(source[0],source[1]), size=10, color='blue')

    #---HCO+ ATCA sightlines
    for i in enumerate(names):
	val=i[0]
	name=i[1]
        source = wcs_header.wcs_sky2pix([[sources[val][0],sources[val][1], 0]], 0)[0]
        atca_hcop = ax.scatter(source[0],source[1], s=90, facecolor='red', edgecolor='none', marker='^', linewidth=1)
        #ax.annotate(name, xy=(source[0],source[1]), size=10, color='red')

    detect_lon=293.8509748
    detect_lat=-31.3708847
    source = wcs_header.wcs_sky2pix([[detect_lon,detect_lat, 0]], 0)[0]
    ax.scatter(source[0],source[1], s=90, facecolor='red', edgecolor='black', marker='^', linewidth=1)
    #ax.annotate('Shapley blue', xy=(source[0],source[1]), size=10, color='blue')


    ax.legend([atca_hcop, matthews, mizuno_co, sembach_h2], ['HCO+ absorption (this work)',  'HI absorption (Kobulnickey+ 1999, Matthews+ 2009)', 'CO emission (Muller+ 2003, Mizuno+ 2006)', 'H2 absorption (Richter+ 2001)'], loc=4, scatterpoints=1)

    cbaxes= fig.add_axes([0.9,0.18,0.015,0.64])
    cb = fig.colorbar(im, cax=cbaxes)
    #cb.ax.minorticks_on()
    #cmin, cmax = cb.get_clim()
    #ticks = np.logspace(cmin,cmax,num=200)
    #cb.set_ticks(ticks)
    #cb.set_ticklabels(['%.3g' %10**t for t in ticks])
    cb.set_label(r'N(HI) ($10^{18}$ cm$^{-2}$), Putman et al. 2003')

    plt.savefig(savedir + filename)

def load_fits(filename,return_header=True):
    ''' Loads a fits file.
    '''

    import pyfits as pf

    f = pf.open(filename)

    if return_header:
        return f[0].data,f[0].header
    else:
        return f[0].data

def ga2equ(ga):

    l = radians(ga[0])
    b = radians(ga[1])

    # North galactic pole (J2000) -- according to Wikipedia
    pole_ra = radians(192.859508)
    pole_dec = radians(27.128336)
    posangle = radians(122.932-90.0)

    ra = atan2( (cos(b)*cos(l-posangle)), (sin(b)*cos(pole_dec) - cos(b)*sin(pole_dec)*sin(l-posangle)) ) + pole_ra
    dec = asin( cos(b)*cos(pole_dec)*sin(l-posangle) + sin(b)*sin(pole_dec) )
    return np.array([degrees(ra), degrees(dec)])

def main():

    import numpy as np
    from os import system,path
    import csv

    # define directory locations
    output_dir = '/d/leffe2/cmurray/python/'
    figure_dir = '/d/leffe2/cmurray/python/figures/'
    im_dir = '/d/leffe2/cmurray/python/'
    target_dir='/d/leffe2/cmurray/data/ATCA/'

    # load image
    im_data, im_header = load_fits(im_dir + \
		'fig5data.fits',
            return_header=True)

    im_data=np.squeeze(im_data)

    with open(target_dir+'targets.csv') as csvfile:
        sources = [(float(x),float(y)) for x,y in csv.reader(csvfile, delimiter=',')]

    names=open(target_dir+'names.txt')

    # Plot HI image
    plot_image(image=im_data, header=im_header,
        savedir=figure_dir, sources=sources,
        names=names, colormap='Greys',
        filename='test2.ps',show=0)



if __name__ == '__main__':
    main()

