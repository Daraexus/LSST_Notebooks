
import lsst.afw.geom as afwGeom
import lsst.afw.display.utils as displayUtils
import matplotlib.pyplot as plt
import numpy as np
import lsst.afw.table as afwTable
from lsst.meas.algorithms.detection import SourceDetectionTask



def get_stamp(source, exposure, offset=10):

    
        bbox = source.getFootprint().getBBox()

        mos = displayUtils.Mosaic()
        Begin = afwGeom.Point2I(bbox.getBeginX(), bbox.getBeginY())
        End = afwGeom.Point2I(bbox.getEndX(), bbox.getEndY())

        ExpOrig = afwGeom.Point2I(exposure.getX0()-1, exposure.getY0()-1)



        correctedBegin = bbox.getBegin()- ExpOrig
        correctedEnd = bbox.getEnd() - ExpOrig

        correctedBegin= afwGeom.Point2I(correctedBegin.getX()-offset,correctedBegin.getY()-offset )
        correctedEnd = afwGeom.Point2I(correctedEnd.getX()+offset,correctedEnd.getY()+offset )
	bboxT = afwGeom.Box2I(correctedBegin,correctedEnd) 

        #print bboxT.toString
        return exposure.Factory(exposure,bboxT, True)

def get_fluxes_and_sigmas(source_list, flux_variable):
	sigmas = []
	fluxes = []
	for source in source_list:
		flux = source.get(flux_variable)
		if np.isnan(flux) == False:
			fluxes.append(flux)
			sigmas.append(source.get(flux_variable+"Sigma"))

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

def detect_diasources(exposure, doSmooth=False):

    schema = afwTable.SourceTable.makeMinimalSchema()
    table = afwTable.SourceTable.make(schema)

    config = SourceDetectionTask.ConfigClass()
    config.thresholdPolarity = "both"
    config.thresholdValue = 5.5
    config.reEstimateBackground = False
    config.thresholdType = "pixel_stdev"

    detectionTask = SourceDetectionTask(config=config, schema=schema)
    table = afwTable.SourceTable.make(schema)
    results = detectionTask.makeSourceCatalog(table=table, exposure=exposure , doSmooth=doSmooth )

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
