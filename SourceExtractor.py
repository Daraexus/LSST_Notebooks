import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
import lsst.afw.coord as afwCoord
import lsst.afw.math as afwMath
import lsst.afw.detection as afwDetection
import lsst.afw.display.utils as displayUtils
import lsst.afw.table as afwTable
import lsst.meas.algorithms as measAlg
import lsst.afw.display.ds9 as ds9
import numpy as np
import os

def Show_image_and_footprints(exposure, initial_frame):


    frame = initial_frame
    #exposure = afwImage.ExposureF(image_path)
    mi = exposure.getMaskedImage()
    im = mi.getImage()
		    
    mask = mi.getMask()
    mask &= ~(mask.getPlaneBitMask("DETECTED")) 
			        
				    
    print mask.getMaskPlaneDict().keys()
					    
    #
    #Subtract background
    #
    back_size = 64
    bctrl = afwMath.BackgroundControl(im.getWidth()//back_size + 1, im.getHeight()//back_size + 1)
    backobj = afwMath.makeBackground(im, bctrl)
    im -= backobj.getImageF("LINEAR")
					          
    nsigma = 5
    stats = afwMath.makeStatistics(im, afwMath.MEANCLIP | afwMath.STDEVCLIP)
    threshold5 = stats.getValue(afwMath.MEANCLIP) + nsigma*stats.getValue(afwMath.STDEVCLIP)
    print "Mean of the image is %.3f counts, stdev is %.3f counts" % (stats.getValue(afwMath.MEANCLIP), stats.getValue(afwMath.STDEVCLIP))
    print "The %d-sigma detection threshold is at %.2f counts" % (nsigma, threshold5)
																          
																	      
    #
    #Smooth image with a
    # 1 2 1
    # 2 4 2
    # 1 2 1
    # filter
    oneD = afwMath.PolynomialFunction1D([2, 0, -1])
    kernel = afwMath.SeparableKernel(3, 3, oneD, oneD)
    smoothedIm = im.Factory(im.getDimensions())
    afwMath.convolve(smoothedIm, im, kernel)
																	      
    npixMin=1
    threshold = afwDetection.Threshold(threshold5)
    fs0 = afwDetection.FootprintSet(smoothedIm, threshold, npixMin)
    #for footprint in fs0.getFootprints():

        #print footprint.getCentroid()
        #displayUtils.drawFootprint(footprint, frame=frame)
    
    
    #ds9.mtv(mi, frame=frame, title="Measured sources")
    grow =1
    if grow > 0:
        isotropic = False
        fs = afwDetection.FootprintSet(fs0, grow, isotropic)
    #for footprint in fs.getFootprints():
    #    displayUtils.drawFootprint(footprint, frame=frame, ctype=ds9.BLUE)
        
    fs.setMask(mi.getMask(), "DETECTED")

   
    
    apRad=3
    # Define the measurements we want to make
    ctrlCentroid = measAlg.SdssCentroidControl()
    ctrlAperture = measAlg.SincFluxControl()
    ctrlAperture.radius2 = apRad
    ctrlShape = measAlg.SdssShapeControl()
    ctrlFlux = measAlg.GaussianFluxControl()
    ctrlPsf = measAlg.PsfFluxControl()
    schema = afwTable.SourceTable.makeMinimalSchema()
    algorithms = [
        measAlg.MeasureSourcesBuilder().addAlgorithm(ctrlCentroid).build(schema),
        measAlg.MeasureSourcesBuilder().addAlgorithm(ctrlAperture).build(schema),
    	measAlg.MeasureSourcesBuilder().addAlgorithm(ctrlShape).build(schema),
    	measAlg.MeasureSourcesBuilder().addAlgorithm(ctrlFlux).build(schema),
    	measAlg.MeasureSourcesBuilder().addAlgorithm(ctrlPsf).build(schema)
    ]
    cat = afwTable.SourceCatalog(schema)
    
    table = cat.table
    table.defineCentroid("centroid.sdss")
    table.defineApFlux("flux.sinc")
    table.defineShape("shape.sdss")
    table.defineModelFlux("flux.gaussian")
    table.definePsfFlux("flux.psf")
    # Measure sources
    fs.makeSources(cat)
    print "Measuring %d objects" % (len(cat))
    
    for source in cat:
        for alg in algorithms:
            alg.apply(source, exposure)
            
 #We dont want ds9. Do we?   
 #   ds9.setMaskTransparency(80)
 #   with ds9.Buffering():
 #       for source in cat:
 #           xs = source.getCentroid().getX()-exposure.getX0()
 #           ys = source.getCentroid().getY()-exposure.getY0()
            #print source.getCentroid().getX(),exposure.getX0()
#            print source.getApFlux()
#            ds9.dot("+", xs, ys, frame=frame, size=10, ctype = ds9.YELLOW)
            
    frame += 1    
    return cat, exposure


from lsst.meas.algorithms.detection import SourceDetectionTask
from lsst.meas.algorithms.measurement import SourceMeasurementTask
from lsst.pipe.tasks.measurePsf import MeasurePsfTask
def Get_and_Measure_sources(exposure):

    schema = afwTable.SourceTable.makeMinimalSchema()

    config = SourceDetectionTask.ConfigClass()
    config.reEstimateBackground = False
    detectionTask = SourceDetectionTask(config=config, schema=schema)

    table = afwTable.SourceTable.make(schema)
    sources = detectionTask.run(table, exposure, sigma=2).sources

    #schema =  sources.schema
    config = SourceMeasurementTask.ConfigClass()
    config.slots.psfFlux = "flux.sinc" # use of the psf flux is hardcoded in secondMomentStarSelector
    measureTask = SourceMeasurementTask(schema, config=config)

    config = MeasurePsfTask.ConfigClass()
    starSelector = config.starSelector.apply()
    starSelector.config.badFlags = ["flags.pixel.edge", "flags.pixel.cr.center", "flags.pixel.interpolated.center", "flags.pixel.saturated.center"]
    psfDeterminer = config.psfDeterminer.apply()
    psfDeterminer.config.sizeCellX = 128
    psfDeterminer.config.sizeCellY = 128
    psfDeterminer.config.spatialOrder = 1
    psfDeterminer.config.nEigenComponents = 3
    measurePsfTask = MeasurePsfTask(config=config, schema=schema)

    table = afwTable.SourceTable.make(schema)
    sources = detectionTask.run(table, exposure, sigma=2).sources
    measureTask.measure(exposure, sources)


    return sources, exposure






