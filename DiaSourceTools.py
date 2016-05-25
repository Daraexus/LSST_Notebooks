
import lsst.afw.geom as afwGeom
import lsst.afw.display.utils as displayUtils
import numpy as np

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
