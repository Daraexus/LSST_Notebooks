import lsst.daf.persistence as dafPersist
import matplotlib.pyplot as plt
import numpy as np

import lsst.afw.table as afwTable
import lsst.afw.geom as afwGeom
import Utils.DiaSourceTools as DSTools
from astropy.time import Time

import lsst.afw.display.ds9 as ds9

import re

from multiprocessing import Pool
from functools import partial

from astropy.table import Column

from astropy.io import ascii


def get_light_curves_per_patch(butler, filter, dates, patch):

    multi_matches = None
    results = []
    stack_lcs = []
    for date in dates:
        t = Time(date)
        date_catalog = None
        #print date
        if butler.datasetExists("dayDiff_class_diaSrc", dataId={"filter":filter,"tract":0, "patch":patch, 'date':date}):
            date_catalog = butler.get("dayDiff_class_diaSrc", dataId={"filter":filter,"tract":0, "patch":patch, 'date':date})
            if multi_matches is None:
                multi_matches= afwTable.MultiMatch(date_catalog.schema, {'mjd':'D'}, radius=afwGeom.Angle(1./3600., afwGeom.degrees))
            multi_matches.add(date_catalog, {'mjd':int(t.mjd)})


    if multi_matches is not None:
        results = multi_matches.finish(removeAmbiguous=False)
        lcs = DSTools.get_light_curves_from_multimatch_results2(results)
        t_lcs = DSTools.threshold_light_curves(lcs, 3)
        for t_lc in t_lcs:
            stack_lcs.append(DSTools.build_lightcurve4(t_lc, "base_CircularApertureFlux_4_5_flux", filter))
            
    return patch, stack_lcs

def function(args):
    return get_light_curves_per_patch(*args)

def build_multi_filter_lc(lcs):
    final_lc = []
 
    for lc in lcs:
        added = False
        
        for i, f_lc in enumerate(final_lc):
            s1 = {'ra':np.mean(lc["ra"]), 'dec':np.mean(lc["dec"])}
            s2 = {'ra':np.mean(f_lc["ra"]), 'dec':np.mean(f_lc["dec"])}
      
            if DSTools.source_distance(s1,s2)<1.0:             
                added = True
                final_lc[i] = vstack([lc, f_lc])
		break
                  #print final_lc[i]
        if added == False:
                  final_lc.append(lc)
                  
    return final_lc

def FluxToMagnitud( flux,  zp):
    return zp-(2.5*np.log10(flux))

def FluxErrorToMagnitud(flux, error):
    return (2.5/np.log(10))*(error/flux)
    

DATADIR="/datadec/cppm/jpreyes/CFHT_Complete"
directory = DATADIR+"/detect_testSN_2/"
butler = dafPersist.Butler(directory) 

patches = []
patches_file=open(DATADIR+"/patches.txt", "rb")
for line in patches_file:
    line = line.replace('\n','')
    text = re.split("=| ",line)
    patches.append(text[-1])
patches_file.close()

dates = []
days_file=open(DATADIR+"/days_total.txt", "rb")

for day in days_file:
    day = day.replace('\n','')
    dates.append(day)

days_file.close()

p = Pool(56)



params = [(butler, 'r', dates, patch) for patch in patches]
params2 = [(butler, 'g', dates, patch) for patch in patches]
params3 = [(butler, 'i', dates, patch) for patch in patches]
params4 = [(butler, 'z', dates, patch) for patch in patches]


params.extend(params2)
params.extend(params3)
params.extend(params4)

res = p.map(function, params)

from astropy.table import Table, vstack
lcs = {}
for r in res:
    if len(r[1])>0:
        if lcs.has_key(r[0]):
            lcs[r[0]].extend(r[1])
        else:
            lcs[r[0]] = r[1]

lcs2 = lcs
lcs = lcs.values()
print len(lcs)
p = Pool(100)
mf_lcs = p.map(build_multi_filter_lc, lcs)

tot_lcs = []
for p_lc in mf_lcs:
    tot_lcs.extend(p_lc)
    
#print tot_lcs



for i, lc in enumerate(tot_lcs):
    l = lc.copy()
    ms = []
    m_es = []
	
    l['ra']=l['ra'].astype(float)
    l['dec']=l['dec'].astype(float)

    for row in lc:
        
        row['zp'] = 30.0
        
        if np.isnan(row['flux']) == False and row['flux'] > 0: 
            f =  FluxToMagnitud(row['flux'], row['zp'])
            e = FluxErrorToMagnitud(row['flux'],row['flux_error'])
            ms.append(f)
            m_es.append(e)
        else:
            ms.append(np.nan)
            m_es.append(np.nan)
    c_ms = Column(ms, name='magnitude')
    c_m_es = Column(m_es, name='magnitude_error')
    l.add_column(c_ms)
    l.add_column(c_m_es)
        
    ascii.write(l, '/datadec/cppm/jpreyes/light_curves/'+str(i)+'.dat')
