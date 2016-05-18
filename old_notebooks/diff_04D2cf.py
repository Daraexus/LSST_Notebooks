
# First test of image subtraction
# ============================

# Start with  source loadlsst, the setup command we need (here obs_cfht), then run ipython notebook
# 
# Next lines correspond to import for this notebook 

# In[ ]:

import numpy as num
import lsst.afw.image as afwImage
import lsst.meas.algorithms as measAlg
import lsst.ip.diffim as ipDiffim
import lsst.ip.diffim.diffimTools as diffimTools


# Now we defined the location of the two input image and the output image

# In[3]:

sciencePath     = "/renoir_data_00/fouchez/lsst/data/CFHTLS/output_hide/coadd_r_dir/deepCoadd/r/0/0,11.fits"
templatePath    = "/renoir_data_00/fouchez/lsst/data/CFHTLS/output/coadd_r_dir_05/deepCoadd/r/0/0,11.fits"
outputPath      = "/renoir_data_00/fouchez/lsst/data/CFHTLS/output/diff/0,11.fits"


# Now read the images with afwImage.ExposureF

# In[5]:

templateExposure = afwImage.ExposureF(templatePath)
scienceExposure  = afwImage.ExposureF(sciencePath)


# Some configuration (no clear idea of what it is)

# In[7]:

config = ipDiffim.ImagePsfMatchTask.ConfigClass()
config.kernel.name = "AL"
subconfig = config.kernel.active


# * Interlude ...  Here a test of latex : $x^2$  ... 
# 
# Now we run the subtraction, via the ipDiffim.ImagePsfMatchTask class 

# In[10]:

defFwhm         = 3.5
fwhmS = defFwhm
fwhmT = defFwhm
psfmatch = ipDiffim.ImagePsfMatchTask(config)
results  = psfmatch.subtractExposures(templateExposure, scienceExposure,
                            templateFwhmPix = fwhmT, scienceFwhmPix = fwhmS)


# Finalize by writing the output file of the image difference

# In[13]:

differenceExposure = results.subtractedExposure
differenceExposure.writeFits(outputPath)

