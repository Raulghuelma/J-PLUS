#!/usr/bin/env python
# coding: utf-8

# In[9]:


from astropy.io import fits
import numpy as np
import os.path
import sys
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy import units as u
import csv
from astropy.table import Table
import pandas as pd
from astropy.nddata import Cutout2D


# In[ ]:





# In[54]:


fits1 = fits.util.get_testdata_filepath('1000001-JPLUS-00421-v2_J0660_swp.fits.fz')
hdul = fits.open(fits1)
hdr= hdul[1]


# In[1]:


#hdr.header


# In[55]:


Pix_ref_RA = hdr.header['CRPIX1']
Pix_ref_DEC = hdr.header['CRPIX2']
WCS_PR_RA = hdr.header['CRVAL1']
WCS_PR_DEC = hdr.header['CRVAL2']
PROJ_RA = hdr.header['CD1_1']
PROJ_DEC = hdr.header['CD2_2']
NPIX1 = hdr.header['NAXIS1']
NPIX2 = hdr.header['NAXIS2']
print(WCS_PR_RA)
print(WCS_PR_DEC)


# In[8]:


f = open('Coordenadas.txt', 'r')
lines=f.readlines()
DEC=lines[1]
RA=lines[0]
#print(RA)
#print(DEC)    
f.close()


# In[25]:


delta_RA = float(WCS_PR_RA) - float(RA)
delta_DEC = float(WCS_PR_DEC) - float(DEC)
print(delta_RA)
print(delta_DEC)


# In[26]:


deltapx_RA = -delta_RA/float(PROJ_RA)
deltapx_DEC = -delta_DEC/float(PROJ_DEC)
print(deltapx_RA)
print(deltapx_DEC)


# In[27]:


POS_RA = Pix_ref_RA + deltapx_RA
POS_DEC = Pix_ref_DEC + deltapx_DEC
print(POS_RA)
print(POS_DEC)


# In[61]:


position = (POS_RA, POS_DEC)
if (POS_RA>NPIX1-300):
    size = (NPIX1-POS_RA,300)
if (POS_DEC>NPIX2-300):
    size = (300, NPIX2-POS_DEC)
if ((POS_RA>NPIX1-300)and(POS_DEC>NPIX2-300)):
    size = (NPIX1-POS_RA, NPIX2-POS_DEC)
if (POS_RA<300):
    size = (POS_RA, 300)
if (POS_DEC<300):
    size = (300, POS_DEC)
if ((POS_RA<300)and(POS_DEC<300)):
    size = (POS_RA, POS_DEC)
if ((POS_RA>NPIX1-300)and(POS_DEC<300)):
    size = (NPIX1-POS_RA, POS_DEC)
if ((POS_RA<300)and(POS_DEC>NPIX2-300)):
    size = (POS_RA, NPIX2-POS_DEC)
else:
    size = (300, 300)
cutout = Cutout2D(hdr.data, position, size)


# In[43]:


hdr.data = cutout.data
wcs = WCS(hdr.header)
hdr.header.update(cutout.wcs)
cutout_Filename = 'SN2008az.fits'
hdr.writeto(cutout_Filename, overwrite=True)


# In[56]:


fits2=fits.open('SN2008az.fits')
hdul = fits2
hdu= hdul[1]
#hdu.header


# In[ ]:




