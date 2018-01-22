#!/usr/bin/env python 

import sys,os,re, math
#import numpy as np
from astropy.io import fits

i = 0
radtodeg = 180. / math.pi
degtorad = math.pi /180.
for arg in sys.argv:
  if (i>0):
    if (arg[0] != "-"):
      myimage = arg
      mylabel = myimage.replace("fit","lbl")
      myout = myimage.replace(".fit","_out.fit")
      if os.path.exists(myimage):
        hdu = fits.open(myimage)
        hdr = hdu[1].header
        dim1 = hdr['NAXIS1']
        dim2 = hdr['NAXIS2']
        flipflag = hdr['P_FLPROT']
        crval1 = hdr['S_RA']
        crval2 = hdr['S_DEC']
        crpix1 = 0.5 + (dim1 / 2.)
        crpix2 = 0.5 + (dim2 / 2.)
        blra = hdr['S_RA1']
        bldec = hdr['S_DEC1']
        brra = hdr['S_RA2']
        brdec = hdr['S_DEC2']
        tlra = hdr['S_RA3']
        tldec = hdr['S_DEC3']
        trra = hdr['S_RA4']
        trdec = hdr['S_DEC4']
        blpx = 0.5
        blpy = 0.5
        brpx = dim1 + 0.5
        brpy = 0.5
        tlpx = 0.5
        tlpy = dim2 + 0.5
        trpx = dim1 + 0.5
        trpy = dim2 + 0.5

        modxy = (trpx - blpx)*(trpx - blpx) + (trpy - blpy)*(trpy - blpy)
        sinxy = (trpy - blpy) / math.sqrt(modxy)

        # from spherical right-angle triangles
        cosc = math.cos(degtorad*(trra-brra)) * math.cos(degtorad*(trdec-brdec))
        cosa = math.cos(degtorad*(trra-blra)) * math.cos(degtorad*(trdec-bldec))
        cosb = math.cos(degtorad*(brra-blra)) * math.cos(degtorad*(brdec-bldec))
        sina = math.sin(math.acos(cosa))
        sinb = math.sin(math.acos(cosb))
        cosad = (cosc - cosa * cosb) / (sina*sinb)
        
        theta = math.asin(sinxy) - math.acos(cosad)

        if os.path.exists(mylabel):
          label = open(mylabel,'r')
          for line in label:
            if (line == "END\r\n" or line == "END\n"):
              break
            temp = line.split('=')
            if (temp[0].strip(' ') == "VCO:SPHERICAL_RADIUS"):
              radiuss = (temp[1].strip(' ')).split(' ')
              radius = float(radiuss[0])*1000.
            if (temp[0].strip(' ') == "HORIZONTAL_PIXEL_FOV"):
              hfovrads = (temp[1].strip(' ')).split(' ')
              hfovrad = float(hfovrads[0])
            if (temp[0].strip(' ') == "VERTICAL_PIXEL_FOV"):
              vfovrads = (temp[1].strip(' ')).split(' ')
              vfovrad = float(vfovrads[0])

          cdelt1 = hfovrad * radtodeg # from horizontal pixel field of view in the label (in degrees)
          cdelt2 = vfovrad * radtodeg # from vertical pixel field of view in the label (in degrees)
          pc11 = math.cos(theta)
          pc12 = -math.sin(theta)
          pc21 = math.sin(theta)
          pc22 = math.cos(theta)

          if flipflag == 10:
            cdelt2 = -cdelt2 # images are y-flipped

          if flipflag == 3:
            pc11 = pc11
            pc12 = -pc12
            pc21 = pc21
            pc22 = -pc22
                                        # images are rotated CCW 270 deg 
          hdunew = fits.PrimaryHDU(hdu[1].data)
          hdtuple = (repr(hdr[7:])).split('\n')
          hdunew.header['CDELT1'] = cdelt1
          hdunew.header['CDELT2'] = cdelt2
          hdunew.header['PC1_1'] = pc11
          hdunew.header['PC1_2'] = pc12
          hdunew.header['PC2_1'] = pc21
          hdunew.header['PC2_2'] = pc22
          hdunew.header['CD1_1'] = pc11 * cdelt1
          hdunew.header['CD1_2'] = pc12 * cdelt1
          hdunew.header['CD2_1'] = pc21 * cdelt2
          hdunew.header['CD2_2'] = pc22 * cdelt2
        else:
          print 'WARNING!! Cannot find label %s!\n' %(mylabel)

        hdunew.header['CRVAL1'] = crval1
        hdunew.header['CRVAL2'] = crval2
        hdunew.header['CRPIX1'] = crpix1
        hdunew.header['CRPIX2'] = crpix2
        hdunew.header['CTYPE1'] = 'RA---TAN'
        hdunew.header['CTYPE2'] = 'DEC--TAN'
        hdunew.header['A_RADIUS'] = radius
        hdunew.header['B_RADIUS'] = radius
        hdunew.header['C_RADIUS'] = radius
        hdunew.writeto(myout)
        hdu.close()
      else:
        print 'Cannot find image %s!\n' %(myimage)
        exit(1)
  i = i+1

