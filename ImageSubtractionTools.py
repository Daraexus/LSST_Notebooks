import random
import math
import lsst.afw.math as afwMath
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.table as afwTable

import lsst.ip.diffim as ipDiffim
import lsst.ip.diffim.diffimTools as diffimTools
from lsst.ip.diffim.makeKernelBasisList import makeKernelBasisList
from lsst.ip.diffim import utils as diUtils
from lsst.ip.diffim import diffimLib

import lsst.pipe.tasks.ingest as Ingester
import lsst.pipe.tasks.processCcd as Processer

from lsst.pipe.tasks.registerImage import RegisterTask
from lsst.meas.algorithms import SourceDetectionTask, SourceMeasurementTask, starSelectorRegistry, PsfAttributes, SingleGaussianPsf
from lsst.meas.deblender import SourceDeblendTask
from lsst.ip.diffim import ImagePsfMatchTask, DipoleMeasurementTask, DipoleAnalysis, DiaCatalogSourceSelectorConfig, DiaCatalogSourceSelector, KernelCandidateF, SourceFlagChecker, KernelCandidateF, cast_KernelCandidateF, makeKernelBasisList, KernelCandidateQa, DiaCatalogSourceSelector, DiaCatalogSourceSelectorConfig

import lsst.meas.astrom as measAstrom
import lsst.daf.base as dafBase
import lsst.pipe.base as pipeBase
import lsst.daf.persistence as dafPersist
import lsst.pipe.tasks.imageDifference as ImageDifferencer
import lsst.meas.algorithms.detection as sDet

#Hardcoded initialization parameters. Not sure this is very pythonic.
config = ipDiffim.ImagePsfMatchTask.ConfigClass()
config.kernel.name = "AL"
subconfig = config.kernel.active
defFwhm         = None
scienceFwhmPix = defFwhm
templateFwhmPix = defFwhm
psfmatch = ipDiffim.ImagePsfMatchTask(config)
IDTask = ImageDifferencer.ImageDifferenceTask
IDConfig = IDTask.ConfigClass()
IDifTask = IDTask()
IDifTask.__init__()

FwhmPerSigma = 2 * math.sqrt(2 * math.log(2))
IqrToSigma = 0.741


def SubtractBackground(exposure):
        
        bgConf = sDet.BackgroundConfig()
	background,Exp0 = sDet.estimateBackground(exposure,bgConf,True)
	return Exp0

def ComputeInitialParameters(exposure):
	ctr = afwGeom.Box2D(exposure.getBBox(afwImage.PARENT)).getCenter()
	psfAttr = PsfAttributes(exposure.getPsf(), afwGeom.Point2I(ctr))
	sigmaOrig = psfAttr.computeGaussianWidth(psfAttr.ADAPTIVE_MOMENT)
	return ctr, psfAttr, sigmaOrig

def Convolve(dest_im, src_im, preConvPsf):
	convControl = afwMath.ConvolutionControl()
	afwMath.convolve(dest_im, src_im, preConvPsf.getLocalKernel(), convControl)
	print dest_im
	return dest_im

#Uses a double Gaussian for convolution
def SimplifiedPreConvolve(exposure, sigmaOrig, sciencePsf):
	mi = exposure.getMaskedImage()
	dest_mi = mi.Factory(mi.getDimensions())
	kWidth, kHeight =  sciencePsf.getLocalKernel().getDimensions()
        preConvPsf = SingleGaussianPsf(kWidth, kHeight, sigmaOrig)
	
	return Convolve(dest_mi, mi, preConvPsf), sigmaOrig * math.sqrt(2)

#Uses the PSF model of the image for convolution
def PreConvolve(exposure, sigmaOrig, sciencePsf):
	convControl = afwMath.ConvolutionControl()
	mi = exposure.getMaskedImage()
        dest_mi = mi.Factory(mi.getDimensions())
        psf = exposure.getPsf()
	preConvPsf = psf
	return Convolve(dest_mi, mi, preConvPsf), sigmaOrig * math.sqrt(2)

def SelectSources(exposure, sigma):
	return psfmatch.getSelectSources(exposure, sigma, doSmooth = not False, idFactory=None)

def CalculateNumberOfBasisFunctions(scienceSigma, templateSigma):
	return len(makeKernelBasisList(IDifTask.subtract.config.kernel.active, referenceFwhmPix=scienceSigma * FwhmPerSigma, targetFwhmPix=templateSigma * FwhmPerSigma))

def SelectKernelSources(exposure, selectSources):

	astrometer = measAstrom.Astrometry(measAstrom.MeasAstromConfig())
	astromRet = astrometer.useKnownWcs(selectSources, exposure=exposure)
	matches = astromRet.matches
	kernelSources = selectSources
	random.shuffle(kernelSources, random.random)
	controlSources = kernelSources[IDConfig.controlStepSize]
	kernelSources = [k for i,k in enumerate(kernelSources) if i % IDConfig.controlStepSize]
	return kernelSources


def AlignImages(templateSources, templateExposure, selectSources, scienceExposure):
	wcsResults = IDifTask.fitAstrometry(templateSources, templateExposure, selectSources)
	warpedExp = IDifTask.register.warpExposure(templateExposure, wcsResults.wcs, scienceExposure.getWcs(), scienceExposure.getBBox(afwImage.PARENT))
	return warpedExposure

def SubtractImages(templateExposure, templateSigma,  scienceExposure, scienceSigma, kernelSources):


	subtractRes = IDifTask.subtract.subtractExposures(
                templateExposure=templateExposure,
                scienceExposure=scienceExposure,
                scienceFwhmPix=scienceSigma * FwhmPerSigma,
                templateFwhmPix=templateSigma * FwhmPerSigma,
                candidateList=kernelSources,
                convolveTemplate=IDConfig.convolveTemplate,
                doWarping=not IDConfig.doUseRegister
            )
	return subtractRes.subtractedExposure



