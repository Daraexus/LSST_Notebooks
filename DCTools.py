import pickle
import sncosmo
import numpy as np
import lsst.daf.base as dafBase
import lsst.daf.persistence as dafPersist

import re

import lsst.afw.display.ds9 as ds9

import matplotlib.pyplot as plt

from lsst.ip.diffim import ImagePsfMatchTask, DipoleMeasurementTask
from lsst.meas.base import SingleFrameMeasurementConfig
import lsst.afw.table as afwTable

import lsst.daf.base as dafBase
import DiaSourceTools as DSTools

import lsst.meas.algorithms.detection as sDet
import lsst.afw.image as afwImage

from lsst.meas.algorithms.detection import SourceDetectionTask

import skimage
from skimage import measure
from skimage import data
from skimage import filters
from scipy import ndimage
from skimage.measure import regionprops
from skimage.measure import moments
from skimage.measure import moments_central


import lsst.afw.geom as afwGeom
import random

def SubtractBackground(exp):
  
    bgConf = sDet.BackgroundConfig()
    background,Exp0 = sDet.estimateBackground(exp,bgConf,True)

    return background, Exp0

def merge_sources(results, schema, algMetadata):
    
    table = afwTable.SourceTable.make(schema)
    table.setMetadata(algMetadata)
    fpSet = results.fpSets.positive
    fpSet.merge(results.fpSets.negative, 1, 1, False)
    diaSources = afwTable.SourceCatalog(table)
    fpSet.makeSources(diaSources)
    return diaSources

def remove_nan_sources(diaSources):
    c_diaSources = diaSources.copy()
    c_diaSources.clear()
    for diaSrc in diaSources:
        val = diaSrc.get("coord_ra").asArcseconds()
        if np.isnan(val) == False:
            c_diaSources.append(diaSrc)
        
    return c_diaSources

def get_source_stamp(src, visit, filter, ccds, offset=0):
    
    for ccd in ccds:

        if butler.datasetExists("deepDiff_differenceExp", {'visit': visit , 'filter':filter , 'ccd':ccd}):

            diffExp = butler.get("deepDiff_differenceExp", {'visit': visit , 'filter':filter , 'ccd':ccd})
            bbox = diffExp.getBBox()
            wcs = diffExp.getWcs()
            
            c = afwGeom.Point2I(wcs.skyToPixel(src.getRa(), src.getDec()))
            
            if bbox.contains(c):
                psf = diffExp.getPsf()
                shape = psf.computeShape()
                sigma = shape.getDeterminantRadius()
                #print sigma
                
                return DSTools.get_stamp(src, diffExp, offset=offset), c
            
    return None, None  

def get_dipole_lobes_metrics(stamp):
    
    w,h = stamp.getWidth(), stamp.getHeight()
    pos_fp = stamp.getMaskedImage().getMask().getPlaneBitMask("DETECTED")
    neg_fp = stamp.getMaskedImage().getMask().getPlaneBitMask("DETECTED_NEGATIVE")

    img_arr, mask_arr, var_arr = stamp.getMaskedImage().getArrays()
    
    values = (mask_arr & pos_fp == pos_fp)
    values2 = values & (img_arr > 0)
    pos_pixels = len(img_arr[values])
    pos_flux=  np.sum(abs(img_arr[values2]))

    values = (mask_arr & neg_fp == neg_fp)
    values2 = values & (img_arr < 0)
    neg_flux= np.sum(abs(img_arr[values2]))
    neg_pixels = len(img_arr[values])

    values = mask_arr & (neg_fp|pos_fp) == (neg_fp|pos_fp)
    inter_flux = np.sum(abs(img_arr[values]))
    inter_pixels = len(img_arr[values])

    values = mask_arr & (neg_fp|pos_fp) != 0
    total_flux = np.sum(abs(img_arr[values]))
    mask_pixels = len(img_arr[values])
    
    fluxes = [inter_flux/total_flux, pos_flux/total_flux, neg_flux/total_flux]
    geom = [float(inter_pixels)/mask_pixels, float(pos_pixels)/mask_pixels, float(neg_pixels)/mask_pixels]
    

    return fluxes, geom
    

def show_stamps(stamp, title='Stamp'):

    img_arr, mask_arr, var_arr = stamp.getMaskedImage().getArrays()
    w,h = stamp.getWidth(), stamp.getHeight()
    plt.figure(figsize=(16,5))
    
    plt.suptitle(title)
    
    plt.subplot(1,3,1)
    plt.imshow(img_arr, origin='lower', vmin=img_arr.min(), vmax=img_arr.max(), cmap='gray', extent=(0,w-1, 0, h-1), interpolation='none')


    plt.subplot(1,3,2)
    plt.imshow(mask_arr, origin='lower',  vmin=0, vmax=mask_arr.max(), cmap=cm, extent=(0,w-1, 0, h-1), interpolation='none')


    plt.subplot(1,3,3)
    layer_mask = mask_arr.copy()
    layer_mask = np.ma.masked_where(layer_mask==0, layer_mask)
    plt.imshow(img_arr, origin='lower', vmin=img_arr.min(), vmax=img_arr.max(), cmap='gray', extent=(0,w-1, 0, h-1), interpolation='none')
    plt.imshow(layer_mask, origin='lower', alpha=0.3, vmin=0, vmax=mask_arr.max(), cmap=cm, extent=(0,w-1, 0, h-1), interpolation='none')
    
def visualize_dipoles_and_planes(stamp, source):
    
    plCak = source.getFootprint().getPeaks()
    
    w,h = stamp.getWidth(), stamp.getHeight()
    img_arr, mask_arr, var_arr = stamp.getMaskedImage().getArrays()

    plt.figure(figsize=(16,5))

    plt.subplot(1,2,1)
    plt.imshow(img_arr, origin='lower', vmin=img_arr.min(), vmax=img_arr.max(), cmap='gray', extent=(0,w-1, 0, h-1), interpolation='none')

    for pk in pkCat:
        if pk.getPeakValue() < 0:
            plt.plot(pk.getIx()-stamp.getX0(), pk.getIy()-stamp.getY0(), 'bo')
        else:
            plt.plot(pk.getIx()-stamp.getX0(), pk.getIy()-stamp.getY0(), 'ro')

    plt.subplot(1,2,2)
    im = plt.imshow(mask_arr, origin='lower', vmin=0, vmax=100, cmap=cm, extent=(0,w-1, 0, h-1), interpolation='none')
    
def photo_dipole(photo_prop):
    return photo_prop[1] < 0.66 and photo_prop[2] < 0.66

def geom_dipole(geom_prop):
    return geom_prop[0]<0.33

def detect_point_source(stamp, center, alpha=1):
    
    img_arr, mask_arr, var_arr = stamp.getMaskedImage().getArrays()
    
    mean = np.mean(img_arr)
    std = np.std(img_arr)
    T_p = mean + alpha*std
    T_n = mean - alpha*std

    #thresholding_n = img_arr > T_n
    thresholding_n = img_arr < T_p
    thresholding_p = img_arr < T_p
    thresholding = thresholding_p & thresholding_n
    
    w, h = stamp.getWidth(), stamp.getHeight()
    
    t_img_arr = img_arr.copy()
    t_img_arr[thresholding] = 0

    t_img_arr = t_img_arr.astype(float)

    m_c = moments_central(t_img_arr, center.getX()-stamp.getX0(), center.getY()-stamp.getY0(), order=2)
    

    cov = [[m_c[2,0]/m_c[0,0], m_c[1,1]/m_c[0,0]], [m_c[1,1]/m_c[0,0], m_c[0,2]/m_c[0,0]]]
    #print(np.matrix(cov))
    vals, vectors = np.linalg.eig(cov)
    
    return ( 2*np.sqrt(vals[0]), 2*np.sqrt(vals[1])), t_img_arr 

def classify_point_source(stamp, center):
        try:
            axes, im = detect_point_source(stamp, center)

            h_m, w_m = stamp.getWidth()/2, stamp.getHeight()/2

            if np.isnan(axes[0])==True or np.isnan(axes[1])==True or axes[0] > h_m or axes[1] > w_m:


                return "Point positive"
            else:
                return "Positive"

        except np.linalg.linalg.LinAlgError as e:
                return "Point positive"
