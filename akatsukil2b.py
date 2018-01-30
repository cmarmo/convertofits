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

          # from spherical right-angle triangles

          coscp = math.cos(hfovrad*degtorad*(trpx-brpx)) * math.cos(vfovrad*degtorad*(trpy-brpy))
          cosap = math.cos(hfovrad*degtorad*(trpx-blpx)) * math.cos(vfovrad*degtorad*(trpy-blpy))
          cosbp = math.cos(hfovrad*degtorad*(brpx-blpx)) * math.cos(vfovrad*degtorad*(brpy-blpy))
          sinap = math.sin(math.acos(cosap))
          sinbp = math.sin(math.acos(cosbp))
          cosxy = (coscp - cosap * cosbp) / (sinap*sinbp)

          cosc = math.cos(degtorad*(trra-brra)) * math.cos(degtorad*(trdec-brdec))
          cosa = math.cos(degtorad*(trra-blra)) * math.cos(degtorad*(trdec-bldec))
          cosb = math.cos(degtorad*(brra-blra)) * math.cos(degtorad*(brdec-bldec))
          sina = math.sin(math.acos(cosa))
          sinb = math.sin(math.acos(cosb))
          cosad = (cosc - cosa * cosb) / (sina*sinb)

          theta = math.acos(cosxy) - math.acos(cosad)


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
            pc22 = -pc22 # images are rotated CCW 270 deg 
          hdr['CDELT1'] = cdelt1
          hdr['CDELT2'] = cdelt2
          hdr['PC1_1'] = pc11
          hdr['PC1_2'] = pc12
          hdr['PC2_1'] = pc21
          hdr['PC2_2'] = pc22
          hdr['CD1_1'] = pc11 * cdelt1
          hdr['CD1_2'] = pc12 * cdelt1
          hdr['CD2_1'] = pc21 * cdelt2
          hdr['CD2_2'] = pc22 * cdelt2
        else:
          print('WARNING!! Cannot find label {0:s}!\n'.format(mylabel))

        hdr['CRVAL1'] = crval1
        hdr['CRVAL2'] = crval2
        hdr['CRPIX1'] = crpix1
        hdr['CRPIX2'] = crpix2
        hdr['CTYPE1'] = 'RA---TAN'
        hdr['CTYPE2'] = 'DEC--TAN'
        hdr['A_RADIUS'] = radius
        hdr['B_RADIUS'] = radius
        hdr['C_RADIUS'] = radius

        if os.path.exists(mygeom):
          hdgeom = fits.open(mygeom)
          if os.path.exists(mygeom4):
            hdgeom4 = fits.open(mygeom4)

            ## Sparse description of coordinates
            hdr['CRVAL1a'] = 1.0
            hdr['CRVAL2a'] = 1.0
            hdr['CRPIX1a'] = 0.5
            hdr['CRPIX2a'] = 0.5
            hdr['CTYPE1a'] = 'VELN-TAB'
            hdr['CTYPE2a'] = 'VELT-TAB'
            hdr['CD1_1a'] = 0.5
            hdr['CD1_2a'] = 0.0
            hdr['CD2_1a'] = 0.0
            hdr['CD2_2a'] = 0.5
            hdr['PC1_1a'] = 0.5
            hdr['PC1_2a'] = 0.0
            hdr['PC2_1a'] = 0.0
            hdr['PC2_2a'] = 0.5
            hdr['CDELT1a'] = 0.5
            hdr['CDELT2a'] = 0.5
            hdr['PS1_0a '] = 'WCS-TAB '
            hdr['PS1_1a '] = 'COORDS  '
            hdr['PS1_2a '] = 'INDEXLN '
            hdr['PS2_0a '] = 'WCS-TAB '
            hdr['PS2_1a '] = 'COORDS  '
            hdr['PS2_2a '] = 'INDEXLT '
            hdr['PV1_1a '] = 1
            hdr['PV2_1a '] = 1
            hdr['PV1_2a '] = 1
            hdr['PV2_2a '] = 1
            hdr['PV1_3a '] = 1
            hdr['PV2_3a '] = 2
            ## sparse description of local time
            #hdr['CRVAL1b'] = 1.0
            #hdr['CRPIX1b'] = 0.5
            #hdr['CTYPE1b'] = 'VET--TAB'
            #hdr['CD1_1b'] = 0.5
            #hdr['PS1_0b '] = 'WCS-TAB '
            #hdr['PS1_1b '] = 'LOCTIME '
            #hdr['PS1_2b '] = 'INDEX   '
            #hdr['PV1_1b '] = 1.0
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
            LLpa = hdgeom['Phase angle'].data
            LLia = hdgeom['Incidence angle'].data
            LLea = hdgeom['Emission angle'].data
            LLaa = hdgeom['Azimuthal angle'].data
            [m4,n4] = LLlon.shape
            hdgeom4.close()

            # Building long index
            indexln = np.zeros([m+m4])
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
            indexlt = np.zeros([n+n4])
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
            coordsarr = np.zeros([ymax,xmax,2])

            for j in range(0,ymax):
              for i in range(0,xmax):
                  indj = newindexlt[j]
                  indi = newindexln[i]
                  if indj % 2 == 0:
                      myj = int(indj / 2) - 1
                      myi = int(indi / 2) - 1
                      coordsarr[j,i,0] = lon[myj,myi]
                      coordsarr[j,i,1] = lat[myj,myi]
                  else:
                      myj = int((indj - 1) / 2)
                      myi = int((indi - 1) / 2)
                      coordsarr[j,i,0] = LLlon[myj,myi]
                      coordsarr[j,i,1] = LLlat[myj,myi]
                  #tlt.append(LLlt[i,j])
                  #newpa = np.append(newpa,LLpa[i,j])
                  #newia = np.append(newia,LLia[i,j])
                  #newea = np.append(newea,LLea[i,j])
                  #newaa = np.append(newaa,LLaa[i,j])

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
          tbhdu.header.set('extname', 'WCS-TAB')
        else:
          print('Cannot find geom file {0:s}!\n'.format(mygeom))

        hdu.append(tbhdu)
        hdu.writeto(myout)
        hdu.close()
      else:
        print('Cannot find image {0:s}!\n'.format(myimage))
        exit(1)
  i = i+1

