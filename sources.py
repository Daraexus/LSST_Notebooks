import lsst.meas.astrom as measAstrom
import lsst.daf.base as dafBase
import lsst.pipe.base as pipeBase
import lsst.daf.persistence as dafPersist

import lsst.afw.table as afwTable
import lsst.afw.table
#import lsst.meas.algorithms. as sourceDetector
schema = afwTable.SourceTable.makeMinimalSchema()
table = afwTable.SourceTable.make(afwTable.SourceTable.makeMinimalSchema())

sCatalog1 = afwTable.SourceCatalog.readFits(" /renoir_data_02/jpreyes/lsst_data/CFHTLS_v11/output/deepDiff/05AL01/D3/2005-06-29/r/diaSrc-800720-14.fits")
sCatalog2 = afwTable.SourceCatalog.readFits(" /renoir_data_02/jpreyes/lsst_data/CFHTLS_v11/output/deepDiff/05AL01/D3/2005-06-29/r/diaSrc-800721-14.fits")
sCatalog3 = afwTable.SourceCatalog.readFits(" /renoir_data_02/jpreyes/lsst_data/CFHTLS_v11/output/deepDiff/05AL01/D3/2005-06-29/r/diaSrc-800718-14.fits")
sCatalog4 = afwTable.SourceCatalog.readFits(" /renoir_data_02/jpreyes/lsst_data/CFHTLS_v11/output/deepDiff/05AL01/D3/2005-06-29/r/diaSrc-800719-14.fits")


print sCatalog1.schema
