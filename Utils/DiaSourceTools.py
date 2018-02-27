
import lsst.afw.geom as afwGeom
import lsst.afw.display.utils as displayUtils
import matplotlib.pyplot as plt
import numpy as np
import lsst.afw.table as afwTable
from lsst.meas.algorithms.detection import SourceDetectionTask
import lsst.meas.algorithms.detection as sDet
from astropy.table import Table

import astropy.coordinates as coord
import astropy.units as u
import re

def get_time_mosaic(butler, dataid_list, source, frame=1, equalize=False, title="time_mosaic"):
    mosaic = displayUtils.Mosaic(gutter=5, background=3, mode="x")
    
    
    
    for dataid in dataid_list:
        
        mosaic_temp = displayUtils.Mosaic(gutter=0, background=0, mode="y")
        
        diffExp = butler.get("deepDiff_differenceExp", dataid)
        sciExp  = butler.get("calexp", dataid)
        tmpExp = butler.get("deepDiff_warpedExp", dataid)
        bgConf = sDet.BackgroundConfig()
        background,tmpExp = sDet.estimateBackground(tmpExp,bgConf,True)
        
        s1 = get_stamp(source, sciExp)
        s2 = get_stamp(source, tmpExp)
        s3 = get_stamp(source, diffExp)
        
        if equalize == True:
            s1 = equalize_image(s1)
            s2 = equalize_image(s2)
            s3 = equalize_image(s3)
        
        mosaic_temp.append(s1.getMaskedImage())
        mosaic_temp.append(s2.getMaskedImage())
        mosaic_temp.append(s3.getMaskedImage())
        m = mosaic_temp.makeMosaic(frame=None, display=None).clone()
        mosaic.append(m)
        
    mosaic.makeMosaic(frame=frame, title=title)
    
def equalize_image(u_stamp):
    stamp = u_stamp.clone()
    image = stamp.getMaskedImage().getImage()
    imarr = image.getArray()
    max_im = np.max(imarr) 
    min_im = np.min(imarr)
    
    tot = max_im - min_im
    
    
    for x in range(image.getWidth()):
        for y in range(image.getHeight()):
            image.set(x,y, ((image.get(x,y)-min_im)/tot) *255.)
            
    return stamp

def get_stamp(source, exposure, offset=10):

    
        bbox = source.getFootprint().getBBox()

       
	sourceRa = source.getRa()
	sourceDec = source.getDec()


	wcs = exposure.getWcs()


	

	Center = afwGeom.Point2I(wcs.skyToPixel(sourceRa, sourceDec))
	
	
	height= bbox.getHeight()/2
	width= bbox.getWidth()/2

	height=9
	width=9
	
	centerX= (bbox.getEndX()+bbox.getBeginX())/2
	centerY= (bbox.getEndY()+bbox.getBeginY())/2

	Begin = afwGeom.Point2D(centerX - height, centerY - width)
	Begin = afwGeom.Point2I(Begin)
	
	End = afwGeom.Point2D(centerX + height+1, centerY + width+1)
	End = afwGeom.Point2I(End)

	print centerX, centerY
	print Begin, End
	print exposure.getX0(), exposure.getY0()

        ExpOrig = afwGeom.Point2I(exposure.getX0(), exposure.getY0())



        correctedBegin = Begin- ExpOrig
        correctedEnd = End - ExpOrig

    
        correctedBegin= afwGeom.Point2I(correctedBegin.getX()-offset,correctedBegin.getY()-offset )
        correctedEnd = afwGeom.Point2I(correctedEnd.getX()+offset,correctedEnd.getY()+offset )
	
		
	bboxT = afwGeom.Box2I(correctedBegin,correctedEnd) 
	
	
	
	#bboxT = bbox
        print bboxT.toString()
        print exposure.getBBox().toString()
	return exposure.Factory(exposure,bboxT, deep=True)

def get_fluxes_and_sigmas(source_list, flux_variable):
	sigmas = []
	fluxes = []
	for source in source_list:
		flux = source.get(flux_variable)
		sigma = source.get(flux_variable+"Sigma")
		if np.isnan(flux) == False and np.isnan(sigma) == False:
		
			fluxes.append(flux)
			sigmas.append(sigma)
	return fluxes, sigmas

def get_sources_over_sigma(source_list, sigma_t, flux_variable):

	sources = []
	for source in source_list:
		flux = source.get(flux_variable)
		if np.abs(flux) > sigma_t:
			sources.append(source)
	
	return sources

def get_naive_dipole_probability(source):

    
    pos = source.get("ip_diffim_NaiveDipoleFlux_pos_flux")
    neg = np.abs(source.get("ip_diffim_NaiveDipoleFlux_neg_flux"))
    tot = pos+neg
    pos_per = pos/tot
    neg_per = neg/tot
    
    #print pos_per, neg_per
    
    if pos_per < 0.65 and neg_per < 0.65:
        return 1.0
    else:
        return 0.0

def detect_diasources(exposure, doSmooth=False, threshold=5.5):

    schema = afwTable.SourceTable.makeMinimalSchema()
    table = afwTable.SourceTable.make(schema)

    config = SourceDetectionTask.ConfigClass()
    config.thresholdPolarity = "both"
    config.thresholdValue = threshold
    config.reEstimateBackground = False
    config.thresholdType = "pixel_stdev"

    detectionTask = SourceDetectionTask(config=config, schema=schema)
    table = afwTable.SourceTable.make(schema)
    results = detectionTask.makeSourceCatalog(table=table, exposure=exposure , doSmooth=doSmooth )

    fpSet = results.fpSets.positive
    fpSet.merge(results.fpSets.negative, 2, 2, False)
    diaSources = afwTable.SourceCatalog(table)
    fpSet.makeSources(diaSources)


    results.sources=diaSources
    results.fpSets.positive = fpSet
    return results

def get_cumulative_flux(stamp, plane_mask="DETECTED", positive=True):

	

    mi = stamp.getMaskedImage()
    m = mi.getMask()
    values = []
    for x in range(stamp.getWidth()):
        for y in range(stamp.getHeight()):

            planeb_mask = m.getPlaneBitMask(plane_mask)

            val = mi.getImage().get(x,y)

            if planeb_mask & m[x,y].get(0,0) != 0:
                #print mi.get(x,y)

                if positive:
                    if val > 0:
                        values.append(np.abs(val))
                else:
                    if val < 0:
                        values.append(np.abs(val))


    if len(values) > 0:
        values.sort()



        base = [i for i in range(len(values))]


        cumulative = np.cumsum(values)
        cumulative = cumulative/cumulative[-1]

	return cumulative
    else:
	return None



def plot_cumulative_flux(stamp, plane_mask="DETECTED", positive=True):

    mi = stamp.getMaskedImage()
    m = mi.getMask()
    values = []
    for x in range(stamp.getWidth()):
        for y in range(stamp.getHeight()):
            
            planeb_mask = m.getPlaneBitMask(plane_mask)
            
            val = mi.getImage().get(x,y)
         
            if planeb_mask & m[x,y].get(0,0) != 0:
                #print mi.get(x,y)
                
                if positive:
                    if val > 0:
                        values.append(np.abs(val))
                else:
                    if val < 0:
                        values.append(np.abs(val))
                        
        
    if len(values) > 0:           
        values.sort()
      


        base = [i for i in range(len(values))]
        
        
        cumulative = np.cumsum(values)
        cumulative = cumulative/cumulative[-1]
         
        first = False
        f = 0
        second = False
        s = 0
        third = False
        t = 0
        forth = False
        fo = 0
        for val, b in zip(cumulative, base):
	    #print b, float(b+1) / float(base[-1]+1), 1.0- val
            if not forth and float(b+1) / float(base[-1]+1) > 0.90:
                print "10% of sources contribute " +str((1.0-val)*100)+ "% of total flux"
                forth = True
                fo=b 
            if not first and float(b+1) / float(base[-1]+1) > 0.25:
                print float(b+1) / float(base[-1]+1), b
                print "75% of sources contribute " +str((1.0-val)*100)+ "% of total flux"
                first = True
                f=b
            if not second and float(b+1) / float(base[-1]+1) > 0.50:
                print "50% of sources contribute " +str((1.0-val)*100)+ "% of total flux"
                second = True
                s=b
            if not third and float(b+1) / float(base[-1]+1) > 0.75:
                print "25% of sources contribute " +str((1.0-val)*100)+ "% of total flux"
                third = True
                t=b
            
        plt.step(base, cumulative)
        axes = plt.axes()
        axes.set_xlim(xmax = base[-1])
        
        plt.axvline(f, color='r', linestyle='dashed', linewidth=2)
        plt.axvline(s, color='g', linestyle='dashed', linewidth=2)
        plt.axvline(t, color='b', linestyle='dashed', linewidth=2)
        plt.axvline(fo, color='y', linestyle='dashed', linewidth=2)
        
        plt.show()

	return cumulative
    else:
	print "No values for mask "+ plane_mask
	return None

def source_distance(src1, src2):
    ra, dec = src1['ra'], src1['dec']
    ra2, dec2 = src2['ra'], src2['dec']
            
    return np.sqrt((float(ra)-float(ra2))**2+(float(dec)-float(dec2))**2)/3.14159*180*3600


def threshold_light_curves(light_curves, threshold):
    t_light_curves = [lc for lc in light_curves if len(lc) >= threshold]
    return t_light_curves

def build_light_curve_from_snls_file(data, coord):


    lightcurve = {}
    lightcurve['bandpass'] = []
    lightcurve['mjd'] = []
    lightcurve['ra'] = []
    lightcurve['dec'] = []
    lightcurve['flux'] = []
    lightcurve['flux_error'] = []
    lightcurve['zp'] = []
    lightcurve['zpsys'] = []

    

    for bandpass, mjd, flux, error in data:

        #print 'yep',visit
        lightcurve['bandpass'].append(str('sdss' + bandpass))
        lightcurve['mjd'].append(float(mjd))
        lightcurve['ra'].append(coord.ra.radian)
        lightcurve['dec'].append(coord.dec.radian)
        lightcurve['flux'].append(float(flux))
        lightcurve['flux_error'].append(float(error))
        #lightcurve['flux'].append(src['base_CircularApertureFlux_12_0_flux'])
        #lightcurve['flux_error'].append(src['base_CircularApertureFlux_12_0_fluxSigma'])
        lightcurve['zp'].append(25.0)
        lightcurve['zpsys'].append('ab')

    lc = Table(data=lightcurve)
    return lc

def build_light_curve_from_snls_file_2(data, coord, id, z):


    lightcurve = {}
    lightcurve['id'] = []
    lightcurve['z'] = []
    lightcurve['bandpass'] = []
    lightcurve['mjd'] = []
    lightcurve['ra'] = []
    lightcurve['dec'] = []
    lightcurve['flux'] = []
    lightcurve['flux_error'] = []
    lightcurve['zp'] = []
    lightcurve['zpsys'] = []

    

    for bandpass, mjd, flux, error in data:

        #print 'yep',visit
        lightcurve['id'].append(id)
        lightcurve['z'].append(z)
        lightcurve['bandpass'].append(str('sdss' + bandpass))
        lightcurve['mjd'].append(float(mjd))
        lightcurve['ra'].append(coord.ra.radian)
        lightcurve['dec'].append(coord.dec.radian)
        lightcurve['flux'].append(float(flux))
        lightcurve['flux_error'].append(float(error))
        #lightcurve['flux'].append(src['base_CircularApertureFlux_12_0_flux'])
        #lightcurve['flux_error'].append(src['base_CircularApertureFlux_12_0_fluxSigma'])
        lightcurve['zp'].append(25.0)
        lightcurve['zpsys'].append('ab')

    lc = Table(data=lightcurve)
    return lc


def build_lightcurve(source_list):
    """
    Assemble a light curve data table from available files.
    """

    bandpasses = ['r']


    lightcurve = {}
    lightcurve['classification'] = []
    lightcurve['bandpass'] = []
    lightcurve['mjd'] = []
    lightcurve['ra'] = []
    lightcurve['dec'] = []
    lightcurve['flux'] = []
    lightcurve['flux_error'] = []
    lightcurve['zp'] = []
    lightcurve['zpsys'] = []


    for visit, src in source_list:

        #print 'yep',visit
        lightcurve['classification'].append(src['classification_dipole'])
        lightcurve['bandpass'].append(str('sdss' + bandpasses[0]))
        
        lightcurve['mjd'].append(src['visit'])
        lightcurve['ra'].append(src['coord_ra'])
        lightcurve['dec'].append(src['coord_dec'])
        lightcurve['flux'].append(src['base_CircularApertureFlux_3_0_flux'])
        lightcurve['flux_error'].append(src['base_CircularApertureFlux_3_0_fluxSigma'])
        #lightcurve['flux'].append(src['base_CircularApertureFlux_12_0_flux'])
        #lightcurve['flux_error'].append(src['base_CircularApertureFlux_12_0_fluxSigma'])
        lightcurve['zp'].append(25.0)
        lightcurve['zpsys'].append('ab')
    lightcurve = Table(data=lightcurve)
    return lightcurve

def match_with_lc(validation_array, candidate_array):
    matches = []
    for lc in validation_array:
        #print "light curve"
        val = {"ra":lc[0]["ra"], "dec":lc[0]["dec"]}
        #print np.rad2deg(lc[0]["ra"]), np.rad2deg(lc[0]["dec"])
        for i, slc in enumerate(candidate_array):

            comp = {"ra":np.mean(slc["ra"]), "dec":np.mean(slc["dec"])}
            if source_distance(val, comp)<1:
                print i
                matches.append((lc,slc))


    print len(matches)
    return matches


def load_SNLS_SN():
    
    f = open('/renoir_data_02/jpreyes/lsst_data/sn_control/J_A+A_523_A7_table9.dat.txt','r')
    data_elems = f.read()
    elems = re.findall('^(.*?D3.*?)\\|(.*?)\\|(.*?)\\|(.*?)\\|(.*?)\\|(.*?)\\|(.*?)\\|', data_elems, re.MULTILINE)
    f.close()

    f = open('/renoir_data_02/jpreyes/lsst_data/sn_control/J_A+A_523_A7_table10.dat.txt','r')
    data = f.read()
    f.close()
    
    snls_array = []
    for sn in elems:

            c = coord.SkyCoord(sn[1], unit=(u.hourangle, u.deg))

            m = re.findall('^'+str(sn[0])+'\\|(r|g|z|i)\\|(.*?)\\|(.*?)\\|(.*?)$', data, re.MULTILINE)

            snls_lc = build_light_curve_from_snls_file_2(m, c, sn[0], sn[-1])

            if len(m)>0:
                #print sn[0], c.ra.deg, c.dec.deg

                #plt.errorbar(snls_lc['mjd'], snls_lc['flux'], yerr=snls_lc['flux_error'], fmt='.', color='blue')
                #

                snls_array.append(snls_lc)
    return snls_array


def build_lightcurve6(source_list, flux_parameter, filters):
    """
    Assemble a light curve data table from available files.
    """

    bandpasses = filters


    lightcurve = {}
    lightcurve['classification'] = []
    lightcurve['bandpass'] = []
    lightcurve['mjd'] = []
    lightcurve['ra'] = []
    lightcurve['dec'] = []
    lightcurve['flux'] = []
    lightcurve['flux_error'] = []
    lightcurve['zp'] = []
    lightcurve['zpsys'] = []
    #lightcurve['ccd'] = []


    for src in source_list:

        #print 'yep',visit
        lightcurve['classification'].append(src['ip_diffim_ClassificationDipole_value'])
        lightcurve['bandpass'].append(str('sdss' + bandpasses[int(src['filter'])]))

        lightcurve['mjd'].append(src['mjd'])
        lightcurve['ra'].append(src['coord_ra'])
        lightcurve['dec'].append(src['coord_dec'])
        lightcurve['flux'].append(src[flux_parameter])
        lightcurve['flux_error'].append(src[flux_parameter+"Sigma"])
        #lightcurve['flux'].append(src['base_CircularApertureFlux_12_0_flux'])
        #lightcurve['flux_error'].append(src['base_CircularApertureFlux_12_0_fluxSigma'])
        lightcurve['zp'].append(25.0)
        lightcurve['zpsys'].append('ab')
        #lightcurve['ccd'].append(src['ccd'])
    lightcurve = Table(data=lightcurve)
    return lightcurve


def build_lightcurve_tuple(source_list, flux_parameter):
    """
    Assemble a light curve data table from available files.
    """

    bandpasses = ['r']


    lightcurve = {}
    lightcurve['classification'] = []
    lightcurve['bandpass'] = []
    lightcurve['mjd'] = []
    lightcurve['ra'] = []
    lightcurve['dec'] = []
    lightcurve['flux'] = []
    lightcurve['flux_error'] = []
    lightcurve['zp'] = []
    lightcurve['zpsys'] = []


    for visit, src in source_list:

        #print 'yep',visit
        lightcurve['classification'].append(src['classification_dipole'])
        lightcurve['bandpass'].append(str('sdss' + bandpasses[0]))
        
        lightcurve['mjd'].append(src['visit'])
        lightcurve['ra'].append(src['coord_ra'])
        lightcurve['dec'].append(src['coord_dec'])
        lightcurve['flux'].append(src[flux_parameter])
        lightcurve['flux_error'].append(src[flux_parameter+"Sigma"])
        #lightcurve['flux'].append(src['base_CircularApertureFlux_12_0_flux'])
        #lightcurve['flux_error'].append(src['base_CircularApertureFlux_12_0_fluxSigma'])
        lightcurve['zp'].append(25.0)
        lightcurve['zpsys'].append('ab')
    lightcurve = Table(data=lightcurve)
    return lightcurve


def build_lightcurve(source_list, flux_parameter):
    """
    Assemble a light curve data table from available files.
    """

    bandpasses = ['r']


    lightcurve = {}
    lightcurve['classification'] = []
    lightcurve['bandpass'] = []
    lightcurve['mjd'] = []
    lightcurve['ra'] = []
    lightcurve['dec'] = []
    lightcurve['flux'] = []
    lightcurve['flux_error'] = []
    lightcurve['zp'] = []
    lightcurve['zpsys'] = []


    for src in source_list:

        #print 'yep',visit
        lightcurve['classification'].append(src['classification_dipole'])
        lightcurve['bandpass'].append(str('sdss' + bandpasses[0]))
        
        lightcurve['mjd'].append(src['visit'])
        lightcurve['ra'].append(src['coord_ra'])
        lightcurve['dec'].append(src['coord_dec'])
        lightcurve['flux'].append(src[flux_parameter])
        lightcurve['flux_error'].append(src[flux_parameter+"Sigma"])
        #lightcurve['flux'].append(src['base_CircularApertureFlux_12_0_flux'])
        #lightcurve['flux_error'].append(src['base_CircularApertureFlux_12_0_fluxSigma'])
        lightcurve['zp'].append(25.0)
        lightcurve['zpsys'].append('ab')
    lightcurve = Table(data=lightcurve)
    return lightcurve

def build_lightcurve2(source_list, flux_parameter):
    """
    Assemble a light curve data table from available files.
    """

    bandpasses = ['r']


    lightcurve = {}
    lightcurve['classification'] = []
    lightcurve['bandpass'] = []
    lightcurve['mjd'] = []
    lightcurve['ra'] = []
    lightcurve['dec'] = []
    lightcurve['flux'] = []
    lightcurve['flux_error'] = []
    lightcurve['zp'] = []
    lightcurve['zpsys'] = []
    lightcurve['ccd'] = []


    for src in source_list:

        #print 'yep',visit
        lightcurve['classification'].append(src['classification_dipole'])
        lightcurve['bandpass'].append(str('sdss' + bandpasses[0]))
        
        lightcurve['mjd'].append(src['visit'])
        lightcurve['ra'].append(src['coord_ra'])
        lightcurve['dec'].append(src['coord_dec'])
        lightcurve['flux'].append(src[flux_parameter])
        lightcurve['flux_error'].append(src[flux_parameter+"Sigma"])
        #lightcurve['flux'].append(src['base_CircularApertureFlux_12_0_flux'])
        #lightcurve['flux_error'].append(src['base_CircularApertureFlux_12_0_fluxSigma'])
        lightcurve['zp'].append(25.0)
        lightcurve['zpsys'].append('ab')
        lightcurve['ccd'].append(src['ccd'])
    lightcurve = Table(data=lightcurve)
    return lightcurve


def build_lightcurve3(source_list, flux_parameter):
    """
    Assemble a light curve data table from available files.
    """

    bandpasses = ['r']


    lightcurve = {}
    lightcurve['classification'] = []
    lightcurve['bandpass'] = []
    lightcurve['mjd'] = []
    lightcurve['ra'] = []
    lightcurve['dec'] = []
    lightcurve['flux'] = []
    lightcurve['flux_error'] = []
    lightcurve['zp'] = []
    lightcurve['zpsys'] = []
    #lightcurve['ccd'] = []


    for src in source_list:

        #print 'yep',visit
        lightcurve['classification'].append(src['classification_dipole'])
        lightcurve['bandpass'].append(str('sdss' + bandpasses[0]))
        
        lightcurve['mjd'].append(src['mjd'])
        lightcurve['ra'].append(src['coord_ra'])
        lightcurve['dec'].append(src['coord_dec'])
        lightcurve['flux'].append(src[flux_parameter])
        lightcurve['flux_error'].append(src[flux_parameter+"Sigma"])
        #lightcurve['flux'].append(src['base_CircularApertureFlux_12_0_flux'])
        #lightcurve['flux_error'].append(src['base_CircularApertureFlux_12_0_fluxSigma'])
        lightcurve['zp'].append(25.0)
        lightcurve['zpsys'].append('ab')
        #lightcurve['ccd'].append(src['ccd'])
    lightcurve = Table(data=lightcurve)
    return lightcurve

def build_lightcurve4(source_list, flux_parameter, filter):
    """
    Assemble a light curve data table from available files.
    """

    bandpasses = [filter]


    lightcurve = {}
    lightcurve['classification'] = []
    lightcurve['bandpass'] = []
    lightcurve['mjd'] = []
    lightcurve['ra'] = []
    lightcurve['dec'] = []
    lightcurve['flux'] = []
    lightcurve['flux_error'] = []
    lightcurve['zp'] = []
    lightcurve['zpsys'] = []
    #lightcurve['ccd'] = []


    for src in source_list:

        #print 'yep',visit
        lightcurve['classification'].append(src['classification_dipole'])
        lightcurve['bandpass'].append(str('sdss' + bandpasses[0]))
        
        lightcurve['mjd'].append(src['mjd'])
        lightcurve['ra'].append(src['coord_ra'])
        lightcurve['dec'].append(src['coord_dec'])
        lightcurve['flux'].append(src[flux_parameter])
        lightcurve['flux_error'].append(src[flux_parameter+"Sigma"])
        #lightcurve['flux'].append(src['base_CircularApertureFlux_12_0_flux'])
        #lightcurve['flux_error'].append(src['base_CircularApertureFlux_12_0_fluxSigma'])
        lightcurve['zp'].append(25.0)
        lightcurve['zpsys'].append('ab')
        #lightcurve['ccd'].append(src['ccd'])
    lightcurve = Table(data=lightcurve)
    return lightcurve

def build_lightcurve5(source_list, flux_parameter, filter):
    """
    Assemble a light curve data table from available files.
    """

    bandpasses = [filter]


    lightcurve = {}
    lightcurve['classification'] = []
    lightcurve['bandpass'] = []
    lightcurve['mjd'] = []
    lightcurve['ra'] = []
    lightcurve['dec'] = []
    lightcurve['flux'] = []
    lightcurve['flux_error'] = []
    lightcurve['zp'] = []
    lightcurve['zpsys'] = []
    #lightcurve['ccd'] = []


    for src in source_list:

        #print 'yep',visit
        lightcurve['classification'].append(src['ip_diffim_ClassificationDipole_value'])
        lightcurve['bandpass'].append(str('sdss' + bandpasses[0]))
        
        lightcurve['mjd'].append(src['mjd'])
        lightcurve['ra'].append(src['coord_ra'])
        lightcurve['dec'].append(src['coord_dec'])
        lightcurve['flux'].append(src[flux_parameter])
        lightcurve['flux_error'].append(src[flux_parameter+"Sigma"])
        #lightcurve['flux'].append(src['base_CircularApertureFlux_12_0_flux'])
        #lightcurve['flux_error'].append(src['base_CircularApertureFlux_12_0_fluxSigma'])
        lightcurve['zp'].append(25.0)
        lightcurve['zpsys'].append('ab')
        #lightcurve['ccd'].append(src['ccd'])
    lightcurve = Table(data=lightcurve)
    return lightcurve

def multi_match_catalogs(catalog_list, data_refs):
    multi_matches = None

    for catalog, data_ref in zip(catalog_list, data_refs):

            if multi_matches is None and len(catalog)>0:
                multi_matches = afwTable.MultiMatch(catalog[0].schema, {'visit':int, 'ccd':int}, radius=afwGeom.Angle(1./3600., afwGeom.degrees))
            if multi_matches is not None:
                multi_matches.add(catalog, {'visit':data_ref["visit"], 'ccd':data_ref["ccd"] })

    results = multi_matches.finish(removeAmbiguous=False)  
    return results

def get_light_curves_from_multimatch_results(results):
    light_curves = []
    i = 0
    current = -1
    while i < len(results):
        result = results[i]
        if current == -1 or current != result['object']:
            lc = [(result['visit'],result)]
            light_curves.append(lc)
            current = result['object']
        else:
            light_curves[-1].append((result['visit'],result))
        i+=1
    return light_curves

def get_light_curves_from_multimatch_results2(results):
    light_curves = []
    i = 0
    current = -1
    while i < len(results):
        result = results[i]
        if current == -1 or current != result['object']:
            lc = [result]
            light_curves.append(lc)
            current = result['object']
        else:
            light_curves[-1].append(result)
        i+=1
    return light_curves

def match_control_group(lcs, validation_lcs):
    matches = []
    for v_lc in validation_lcs:
        for lc in lcs:
            if source_distance(lc[0], v_lc[0]) < 1:
                matches.append((lc,v_lc))
                break
                
    return matches
    
