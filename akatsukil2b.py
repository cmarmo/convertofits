#!/usr/bin/env python 

import sys,os,re, math
import numpy as np
from astropy.io import fits

i = 0
radtodeg = 180. / math.pi
degtorad = math.pi /180.
for arg in sys.argv:
  if (i>0):
    if (arg[0] != "-"):
      myimage = arg
      mylabel = myimage.replace("fit","lbl")
      mygeom = myimage.replace("_l2b_","_geo_")
      mygeom4 = myimage.replace("_l2b_","_geo4_")
      myout = myimage.replace(".fit","_out.fit")
      if os.path.exists(myimage):
        hdu = fits.open(myimage)
        hdr = hdu[1].header
        dim1 = hdr['NAXIS1']
        dim2 = hdr['NAXIS2']

        if os.path.exists(mylabel):
          label = open(mylabel,'r')
          for line in label:
            if (line == "END\r\n" or line == "END\n"):
              break
            temp = line.split('=')
            if (temp[0].strip(' ') == "VCO:SPHERICAL_RADIUS"):
              radiuss = (temp[1].strip(' ')).split(' ')
              radius = float(radiuss[0])*1000.

        else:
          print('WARNING!! Cannot find label {0:s}!\n'.format(mylabel))

        hdr['A_RADIUS'] = radius
        hdr['B_RADIUS'] = radius
        hdr['C_RADIUS'] = radius

        if os.path.exists(mygeom):
          hdgeom = fits.open(mygeom)
          if os.path.exists(mygeom4):
            hdgeom4 = fits.open(mygeom4)

            ## Sparse description of coordinates

            hdr['CRPIX1'] = 0.5
            hdr['CRVAL1'] = 1.0
            hdr['CRPIX2'] = 0.5
            hdr['CRVAL2'] = 1.0

            hdr['CD1_1'] = 2.
            hdr['CD2_2'] = 2.

            hdr['PC1_1'] = 1.0
            hdr['PC2_2'] = 1.0

            hdr['CDELT1'] = 2.
            hdr['CDELT2'] = 2.


            hdr['CTYPE1'] = 'VELN-TAB'
            hdr['PS1_0 '] = 'WCS-TAB '
            hdr['PS1_1 '] = 'COORDS  '
            hdr['PS1_2 '] = 'LNINDEX '
            hdr['PV1_1 '] = 1
            hdr['PV1_2 '] = 1
            hdr['PV1_3 '] = 1

            hdr['CTYPE2'] = 'VELT-TAB'
            hdr['PS2_0 '] = 'WCS-TAB '
            hdr['PS2_1 '] = 'COORDS  '
            hdr['PS2_2 '] = 'LTINDEX '
            hdr['PV2_1 '] = 1
            hdr['PV2_2 '] = 1
            hdr['PV2_3 '] = 2

            lon = hdgeom['Longitude'].data
            lat = hdgeom['Latitude'].data
            lt = hdgeom['Local time'].data
            pa = hdgeom['Phase angle'].data
            ia = hdgeom['Incidence angle'].data
            ea = hdgeom['Emission angle'].data
            aa = hdgeom['Azimuthal angle'].data
            [m,n] = lon.shape
            #newpa = np.array([])
            #newia = np.array([])
            #newea = np.array([])
            #newaa = np.array([])
            LLlon = hdgeom4['LL Longitude'].data
            LLlat = hdgeom4['LL Latitude'].data
            LLlt = hdgeom4['LL Local time'].data
            LLpa = hdgeom4['LL Phase angle'].data
            LLia = hdgeom4['LL Incidence angle'].data
            LLea = hdgeom4['LL Emission angle'].data
            LLaa = hdgeom4['LL Azimuthal angle'].data
            [m4,n4] = LLlon.shape
            hdgeom4.close()

            # Building long index
            indexln = np.zeros((m+m4))
            for i in range(0,m4-1):
              lntest = LLlon[:,i]
              if len(lntest[~np.isnan(lntest)]):
                indexln[2*i] = (2*i)+1
            for i in range(0,m-1):
              lntest = lon[:,i]
              if len(lntest[~np.isnan(lntest)]):
                indexln[(2*i)+1] = 2*(i+1)
            newindexln = indexln[np.nonzero(indexln)]
            xmax = len(newindexln)

            # Building lat index
            indexlt = np.zeros((n+n4))
            for j in range(0,n4-1):
              lttest = LLlat[j,:]
              if len(lttest[~np.isnan(lttest)]):
                indexlt[2*j] = (2*j)+1
            for j in range(0,n-1):
              lttest = lat[j,:]
              if len(lttest[~np.isnan(lttest)]):
                indexlt[(2*j)+1] = 2*(j+1)
            newindexlt = indexlt[np.nonzero(indexlt)]
            ymax = len(newindexlt)

            # Building coordinate array
            coordsarr = np.zeros((ymax,xmax,2))

            for j in range(0,ymax):
              for i in range(0,xmax):
                  indj = newindexlt[j]
                  indi = newindexln[i]
                  if indj % 2 == 0:
                      myj = int(indj / 2) - 1
                      myi = int(indi / 2) - 1
                      #if ~np.isnan(lon[myj,myi]) and ~np.isnan(lat[myj,myi]):
                      coordsarr[j,i,0] = lon[myj,myi]
                      coordsarr[j,i,1] = lat[myj,myi]
                      #coordsarr[j,i,2] = ia[myj,myi]
                      #coordsarr[j,i,3] = ea[myj,myi]
                  else:
                      myj = int((indj - 1) / 2)
                      myi = int((indi - 1) / 2)
                      #if ~np.isnan(LLlon[myj,myi]) and ~np.isnan(LLlat[myj,myi]):
                      coordsarr[j,i,0] = LLlon[myj,myi]
                      coordsarr[j,i,1] = LLlat[myj,myi]
                      #coordsarr[j,i,2] = LLia[myj,myi]
                      #coordsarr[j,i,3] = LLea[myj,myi]

            coordsarr[~np.isnan(coordsarr)]

          else:
            print('Cannot find geom file {0:s}!\n'.format(mygeom4))

          hdgeom.close()

          ttype = ["COORDS", "LNINDEX", "LTINDEX"]
          tform= [str(2*xmax*ymax)+'D', str(xmax)+'D', str(ymax)+'D']
          tunit = ['deg','','']
          tdim = '(2,' + str(xmax) + ',' + str(ymax) + ')'

          coords = fits.Column(name=ttype[0], unit=tunit[0], format=tform[0], dim=tdim, array=[coordsarr])
          indln  = fits.Column(name=ttype[1], format=tform[1], array=[newindexln])
          indlt  = fits.Column(name=ttype[2], format=tform[2], array=[newindexlt])
          tbhdu = fits.BinTableHDU.from_columns([coords,indln,indlt])
          tbhdu.header.set('extname', 'WCS-TAB ')
          tbhdu.header.set('extver ', 1)
        else:
          print('Cannot find geom file {0:s}!\n'.format(mygeom))

        hdu.append(tbhdu)
        hdu.writeto(myout)
        hdu.close()
      else:
        print('Cannot find image {0:s}!\n'.format(myimage))
        exit(1)
  i = i+1

