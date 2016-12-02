#
# LSST Data Management System
# Copyright 2012 LSST Corporation.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#
import math
import random
import numpy

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.daf.base as dafBase
import lsst.afw.geom as afwGeom
import lsst.afw.math as afwMath
import lsst.afw.table as afwTable
import lsst.meas.astrom as measAstrom
from lsst.pipe.tasks.registerImage import RegisterTask
from lsst.meas.algorithms import SourceDetectionTask, PsfAttributes, SingleGaussianPsf, \
     SecondMomentStarSelectorTask, SecondMomentStarSelectorConfig
from lsst.ip.diffim import ImagePsfMatchTask, DipoleMeasurementTask, DipoleAnalysis, \
    SourceFlagChecker, KernelCandidateF, cast_KernelCandidateF, makeKernelBasisList, \
    KernelCandidateQa, DiaCatalogSourceSelectorTask, DiaCatalogSourceSelectorConfig, \
    GetCoaddAsTemplateTask, GetCalexpAsTemplateTask
import lsst.ip.diffim.diffimTools as diffimTools
import lsst.ip.diffim.utils as diUtils

import DiaSourceTools as DSTools
import DCTools

FwhmPerSigma = 2 * math.sqrt(2 * math.log(2))
IqrToSigma = 0.741

class SourceClassificationConfig(pexConfig.Config):
    """Config for ImageDifferenceTask
    """
    coaddName = pexConfig.Field(
        doc="coadd name: typically one of deep or goodSeeing",
        dtype=str,
        default="deep",
    )

    print "here"

class SourceClassificationTaskRunner(pipeBase.TaskRunner):
    @staticmethod
    def getTargetList(parsedCmd, **kwargs):
        return pipeBase.TaskRunner.getTargetList(parsedCmd,
                                                 **kwargs)


class SourceClassificationTask(pipeBase.CmdLineTask):
    """Subtract an image from a template and measure the result
    """
    ConfigClass = SourceClassificationConfig
    RunnerClass = SourceClassificationTaskRunner
    _DefaultName = "SourceClassification"

    def __init__(self, **kwargs):
        #do nothing
	pipeBase.CmdLineTask.__init__(self, **kwargs)
	print "initialization 2"

    @pipeBase.timeMethod
    def run(self, sensorRef, templateIdList=None):
        #self.log.info("Processing %s" % (sensorRef.dataId))
        
        visit = sensorRef.dataId['visit']
        ccd = sensorRef.dataId['ccd']
        sigma = 6
        
     	algMetadata = dafBase.PropertyList()
        schema = afwTable.SourceTable.makeMinimalSchema()
        dipoleMeasurement = DipoleMeasurementTask

        config = dipoleMeasurement.ConfigClass()
        dipoleMeasurement = dipoleMeasurement(schema, algMetadata=algMetadata)
        
        source_catalog = None
        classification = []
        
        diffExp = sensorRef.get("deepDiff_differenceExp")   
        
        
        results = DSTools.detect_diasources(diffExp, doSmooth=False, threshold=sigma)
        try:
            diaSources =  DCTools.merge_sources(results, schema, algMetadata)
            dipoleMeasurement.run(diaSources, diffExp)
            #print diaSources
            source_catalog = DCTools.remove_nan_sources(diaSources)           
        except Exception, e: 
             print "exception"
             print e
        
        
        if source_catalog is None:
            classification =  [str(visit)+"-"+str(ccd)]
            
        img_label = str(visit)+"-"+str(ccd)
        
        for source in source_catalog:
            
            bbox = diffExp.getBBox()
            wcs = diffExp.getWcs()
            
            center = afwGeom.Point2I(wcs.skyToPixel(source.getRa(), source.getDec()))
            
            try:    
                stamp = DSTools.get_stamp(source, diffExp, offset=10)

                if stamp is not None:
                    stamp_clone = stamp.clone()
                    results = DSTools.detect_diasources(stamp_clone, doSmooth=False,threshold=sigma)


                    if len(results.fpSets.positive.getFootprints()) == 0 and len(results.fpSets.negative.getFootprints()) > 0:  
                        label = 'Negative'
                        classification.append((label, source.get("id"), img_label))
                    elif len(results.fpSets.positive.getFootprints()) > 0 and len(results.fpSets.negative.getFootprints()) == 0:
                        label = DCTools.classify_point_source(stamp_clone, center)
                        classification.append((label, source.get("id"), img_label))
                    else:


                            stamp_clone = stamp.clone()
                            results = DSTools.detect_diasources(stamp_clone, doSmooth=False,threshold=6)
                            photo, geom = DCTools.get_dipole_lobes_metrics(stamp_clone)

                            p_dipole, g_dipole = DCTools.photo_dipole(photo), DCTools.geom_dipole(geom)

                            if p_dipole and g_dipole:
                                label = "Dipole type I"
                            elif p_dipole and not g_dipole:
                                label = "Dipole type II"
                            elif not p_dipole and g_dipole:
                                label = "Fringe"
                            else:
                                label = "Artifact"
                            classification.append((label, source.get("id"), img_label))
                        
            except Exception, e:

                    print e
                    label = "Artifact"
                    classification.append((label, source.get("id"), img_label))
                            
                            
            
            n = 0
            if label == "Dipole type I":
                n= 0
            elif label == "Dipole type II":
                n=1
            elif label == "Fringe":
                n = 3
            elif label == "Artifact":
                n = 4
            elif label == "Negative":
                n = 5
            elif label == "Positive":
                n = 6
            elif label == "Point positive":
                n = 7

            source["classification_dipole"] = n
            
        source_catalog.writeFits("/renoir_data_00/jpreyes/stacks/notebook_files/catalogs/3sigma/"+img_label+".fits")
        
        
        return pipeBase.Struct(
                   classification = classification
        )



    def _getConfigName(self):
        """Return the name of the config dataset
        """
        return "%sDiff_config" % (self.config.coaddName,)

    def _getMetadataName(self):
        """Return the name of the metadata dataset
        """
        return "%sDiff_metadata" % (self.config.coaddName,)

    #def getSchemaCatalogs(self):
    #    """Return a dict of empty catalogs for each catalog dataset produced by this task."""
    #    diaSrc = afwTable.SourceCatalog(self.schema)
    #    diaSrc.getTable().setMetadata(self.algMetadata)
    #    return {self.config.coaddName + "Diff_diaSrc": diaSrc}
	

    @classmethod
    def _makeArgumentParser(cls):
        """Create an argument parser
        """
        parser = pipeBase.ArgumentParser(name=cls._DefaultName)
        parser.add_id_argument("--id", "calexp", help="data ID, e.g. --id visit=12345 ccd=1,2")
        #parser.add_id_argument("--templateId", "calexp", doMakeDataRefList=False,
        #                       help="Optional template data ID (visit only), e.g. --templateId visit=6789")
        return parser

