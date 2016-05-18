import lsst.pipe.tasks.ingest as Ingester
import lsst.pipe.tasks.processCcd as Processer
import lsst.pipe.base.argumentParser as ArgumentParser
import lsst.daf.persistence as dafPersist
import lsst.pipe.tasks.makeSkyMap as SkyMapper

import lsst.pipe.tasks.reportPatches as PatchReporter
import lsst.pipe.tasks.makeCoaddTempExp as TempexpCoadder
import lsst.pipe.tasks.assembleCoadd as Assembler

import lsst.afw.image as afwImage
import lsst.meas.algorithms as measAlg
import lsst.ip.diffim as ipDiffim
import lsst.ip.diffim.diffimTools as diffimTools


import lsst.meas.algorithms.detection as sDet


def SubtractImages(fwhm=None, templateExp=None, scienceExp=None):
    config = ipDiffim.ImagePsfMatchTask.ConfigClass()
    #Kernel name is Alard-Lupton by default
    config.kernel.name = "AL"
    subconfig = config.kernel.active
    #Fwhm will be calculated if none is provide
    defFwhm         = fwhm
    fwhmS = defFwhm
    fwhmT = defFwhm
    psfmatch = ipDiffim.ImagePsfMatchTask(config)
	        
    results  = psfmatch.subtractExposures(templateExp, scienceExp,
    templateFwhmPix = fwhmT, scienceFwhmPix = fwhmS, doWarping=True)
    differenceExposure = results.subtractedExposure
									            
    return differenceExposure

def SubtractBackground(ImagePath=None):
  
    Exp = afwImage.ExposureF(ImagePath)
    bgConf = sDet.BackgroundConfig()
    background,Exp0 = sDet.estimateBackground(Exp,bgConf,True)

    return Exp0

def SubractImagesWithBkgRemoval(TemplatePath, SciencePath):

    templateExp0 = SubtractBackground(ImagePath=TemplatePath)
    scienceExp0 = SubtractBackground(ImagePath=SciencePath)
    #templateExp0 = afwImage.ExposureF(TemplatePath)
    #scienceExp0 = afwImage.ExposureF(SciencePath)
					        
    dif0= SubtractImages(fwhm=2.5,templateExp=templateExp0, scienceExp=scienceExp0)
						       
    #dif0 = scienceExp0.getMaskedImage()
    #dif0 -= templateExp0.getMaskedImage()
							           
    return afwImage.ExposureF(dif0)
