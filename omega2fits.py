#!/usr/bin/env python 

import sys,os,re
import numpy as np
from astropy.io import fits
from struct import *

myversion = 1.0

# HISTORY
# 10/3/2017 first version 

i = 0
for arg in sys.argv:
  if (i>0):
    if (arg[0] != "-"):
      myimage = arg

      ### Each omega cube is a three channel cube ###
      myfitsIR1 = myimage.replace("qub","IR1.fits")
      myfitsIR2 = myimage.replace("qub","IR2.fits")
      myfitsV = myimage.replace("qub","V.fits")
      mygeom = myimage.replace("qub","nav")

      if os.path.exists(myimage):
        toparse = myimage
        ### parsing PDS label of .qub file ### 
        myfile = open(toparse,'r')
        myskip = 0
        rbytes = 1
        lastlabrec = 0
        myarch = 'littleendian'
        lablines = 0
        for line in myfile:
          lablines = lablines + 1
          if (line == "END\r\n" or line == "END\n"):
            break
          temp = line.split('=')
          if (temp[0].strip(' ') == "RECORD_BYTES"):
            rbytes = int(temp[1].strip(' '))
          if (temp[0].strip(' ') == "^QUBE"):
            mystr = temp[1].strip(' ')
            tests = re.match('["a-zA-Z]',mystr)
            if (tests == None):
              lastlabrec = int(mystr.strip('\r\n')) - 1
          if (temp[0].strip(' ') == "CORE_ITEMS"):
            mydim = (temp[1].strip(' \r\n')).split(',')
            mysamples = int(re.sub(r'[^\w]', ' ', mydim[0]))
            mybands = int(re.sub(r'[^\w]', ' ', mydim[1]))
            mylines = int(re.sub(r'[^\w]', ' ', mydim[2]))
          if (temp[0].strip(' ') == "CORE_ITEM_TYPE"):
            tbyte = (temp[1].strip(' ')).split('_')
            if (tbyte[0] == "LSB" or tbyte[0] == "PC" or tbyte[0] == "VAX"):
                myarch = 'littleendian'
            if (tbyte[0] == "MSB" or tbyte[0] == "MAC" or tbyte[0] == "SUN"):
                myarch = 'bigendian'
          if (temp[0].strip(' ') == "CORE_ITEM_BYTES"):
            imbytes = int(temp[1].strip(' '))
            sbit = 8 * imbytes
          if (temp[0].strip(' ') == "SAMPLE_SUFFIX_ITEM_BYTES"):
            sidebytes = int(temp[1].strip(' '))
          if (temp[0].strip(' ') == "SUFFIX_ITEMS"):
            mydim = (temp[1].strip(' \r\n')).split(',')
            ssamples = int(re.sub(r'[^\w]', ' ', mydim[0]))
            sbands = int(re.sub(r'[^\w]', ' ', mydim[1]))
          if (temp[0].strip(' ') == "BAND_SUFFIX_ITEM_BYTES"):
            backbytes = int(temp[1].strip(' '))
          continue

        IRbands = 128
        dimIR = (IRbands, mylines, mysamples)
        dimV = (mybands-(2*IRbands), mylines, mysamples)
        if (sbit == 8):
          mytype = np.int8
          myform = ""
        if (sbit == 16):
          if (tbyte[1] == "UNSIGNED"):
            mytype = np.uint16
            myform = "H" 
          else:
            mytype = np.int16
            myform = "h"
        if (sbit == 32):
          if (len(tbyte)>2 and tbyte[2] == "INTEGER"):
            mytype = np.int32
            myform = "i"
          else:
            mytype = np.float32
            myform = "f"
        if (sbit == 64):
          if (len(tbyte)>2 and tbyte[2] == "INTEGER"):
            mytype = np.long
            myform = "l"
          else:
            mytype = np.float64
            myform = "d"
        bsqIR1 = np.zeros(dimIR,dtype=mytype)
        bsqIR2 = np.zeros(dimIR,dtype=mytype)
        bsqV = np.zeros(dimV,dtype=mytype)
        form = str(mysamples) + myform
        if ((rbytes != 0) and (lastlabrec != 0)):
          myskip = rbytes * lastlabrec
        myfile.seek(myskip)
        bytex = mysamples * imbytes
        if not sidebytes:
          sidebytes = 0
        if not backbytes:
          backbytes = 0

        ### reading data ###
        for y in range(0, mylines):
          for z in range(0, mybands):
            contents=myfile.read(bytex)
            if (z<IRbands):
              bsqIR1[z,y,:]=unpack_from(form, contents)
            elif (z<2*IRbands):
              bsqIR2[z-IRbands,y,:]=unpack_from(form, contents)
            else:
              bsqV[z-(2*IRbands),y,:]=unpack_from(form, contents)
            offset=myfile.tell()+backbytes
            myfile.seek(offset)

          offset=myfile.tell()+(sidebytes*mysamples*sbands)
          myfile.seek(offset)

        myfile.close()

        if os.path.exists(mygeom):
          toparse = mygeom
          ### parsing PDS label from geometry cube ### 
          myfile = open(toparse,'r')
          myskip = 0
          rbytes = 1
          lastlabrec = 0
          myarch = 'littleendian'
          lablines = 0
          for line in myfile:
            lablines = lablines + 1
            if (line == "END\r\n" or line == "END\n"):
              break
            temp = line.split('=')
            if (temp[0].strip(' ') == "RECORD_BYTES"):
              rbytes = int(temp[1].strip(' '))
            if (temp[0].strip(' ') == "^QUBE"):
              mystr = temp[1].strip(' ')
              tests = re.match('["a-zA-Z]',mystr)
              if (tests == None):
                lastlabrec = int(mystr.strip('\r\n')) - 1
            if (temp[0].strip(' ') == "CORE_ITEM_BYTES"):
              imbytes = int(temp[1].strip(' '))

          ### Reading coordinates from geometry cube ###
          if ((rbytes != 0) and (lastlabrec != 0)):
            myskip = rbytes * lastlabrec
          myfile.seek(myskip)
          dim = (1,mysamples,mylines)
          longIR1 = np.zeros(dim,dtype=np.int32)
          latIR1 = np.zeros(dim,dtype=np.int32)
          incIR1 = np.zeros((mylines,mysamples),dtype=np.int32)
          emIR1 = np.zeros((mylines,mysamples),dtype=np.int32)
          paIR1 = np.zeros((mylines,mysamples),dtype=np.int32)
          longIR2 = np.zeros(dim,dtype=np.int32)
          latIR2 = np.zeros(dim,dtype=np.int32)
          incIR2 = np.zeros((mylines,mysamples),dtype=np.int32)
          emIR2 = np.zeros((mylines,mysamples),dtype=np.int32)
          paIR2 = np.zeros((mylines,mysamples),dtype=np.int32)
          longV = np.zeros(dim,dtype=np.int32)
          latV = np.zeros(dim,dtype=np.int32)
          incV = np.zeros((mylines,mysamples),dtype=np.int32)
          emV = np.zeros((mylines,mysamples),dtype=np.int32)
          paV = np.zeros((mylines,mysamples),dtype=np.int32)
          bytex = mysamples * imbytes
          form = str(mysamples) + 'i'
          for y in range(0, mylines):
            offset=myfile.tell()+(6*bytex)
            myfile.seek(offset)
            contents=myfile.read(bytex)
            longIR1[0,:,y]=unpack_from(form, contents)
            contents=myfile.read(bytex)
            latIR1[0,:,y]=unpack_from(form, contents)
            contents=myfile.read(bytex)
            incIR1[y,:]=unpack_from(form, contents)
            contents=myfile.read(bytex)
            emIR1[y,:]=unpack_from(form, contents)
            contents=myfile.read(bytex)
            paIR1[y,:]=unpack_from(form, contents)
            offset=myfile.tell()+(10*bytex)
            myfile.seek(offset)
            contents=myfile.read(bytex)
            longIR2[0,:,y]=unpack_from(form, contents)
            contents=myfile.read(bytex)
            latIR2[0,:,y]=unpack_from(form, contents)
            contents=myfile.read(bytex)
            incIR2[y,:]=unpack_from(form, contents)
            contents=myfile.read(bytex)
            emIR2[y,:]=unpack_from(form, contents)
            contents=myfile.read(bytex)
            paIR2[y,:]=unpack_from(form, contents)
            offset=myfile.tell()+(10*bytex)
            myfile.seek(offset)
            contents=myfile.read(bytex)
            longV[0,:,y]=unpack_from(form, contents)
            contents=myfile.read(bytex)
            latV[0,:,y]=unpack_from(form, contents)
            contents=myfile.read(bytex)
            incV[y,:]=unpack_from(form, contents)
            contents=myfile.read(bytex)
            emV[y,:]=unpack_from(form, contents)
            contents=myfile.read(bytex)
            paV[y,:]=unpack_from(form, contents)
            offset=myfile.tell()+(10*bytex)
            myfile.seek(offset)
          factor = 0.0001
          newlongIR1 = factor * longIR1
          newlatIR1 = factor * latIR1
          newincIR1 = factor * incIR1
          newemIR1 = factor * emIR1
          newpaIR1 = factor * paIR1
          newlongIR2 = factor * longIR2
          newlatIR2 = factor * latIR2
          newincIR2 = factor * incIR2
          newemIR2 = factor * emIR2
          newpaIR2 = factor * paIR2
          newlongV = factor * longV
          newlatV = factor * latV
          newincV = factor * incV
          newemV = factor * emV
          newpaV = factor * paV
        else:
          print 'Cannot find geometry for image %s!\n' %(myimage)
          continue

        ### Writing FITS files ###
        tdim = '(1,' + str(mysamples) + ',' + str(mylines)+ ')'

        # SWIR-C channel #
        hduIR1 = fits.PrimaryHDU(bsqIR1)
        tbhduIR1 = fits.BinTableHDU.from_columns([fits.Column(name='COORDS1', unit='deg', dim=tdim, format=str(mylines*mysamples)+'D', array=newlongIR1), fits.Column(name='COORDS2', unit='deg', dim=tdim, format=str(mylines*mysamples)+'D', array=newlatIR1)])
        hduincIR1 = fits.ImageHDU(newincIR1)
        hduemIR1 = fits.ImageHDU(newemIR1)
        hdupaIR1 = fits.ImageHDU(newpaIR1)
        listIR1 = fits.HDUList([hduIR1,tbhduIR1,hduincIR1,hduemIR1,hdupaIR1])
        IR1hdr = listIR1[0].header
        IR1hdr.set('ctype1', 'RA---TAB')
        IR1hdr.set('ctype2', 'DEC--TAB')
        IR1hdr.set('CRPIX1', 1.0)
        IR1hdr.set('CRPIX2', 1.0)
        IR1hdr.set('CRVAL1', 1.0)
        IR1hdr.set('CRVAL2', 1.0)
        IR1hdr.set('CD1_1', 1.0)
        IR1hdr.set('CD2_2', 1.0)
        IR1hdr.set('CD1_2', 0.0)
        IR1hdr.set('CD2_1', 0.0)
        IR1hdr.set('PS1_0 ', 'WCS-TAB ')
        IR1hdr.set('PS1_1 ', 'COORDS1 ')
        IR1hdr.set('PS2_0 ', 'WCS-TAB ')
        IR1hdr.set('PS2_1 ', 'COORDS2 ')
        IR1hdr.set('PV1_3 ', 1.0)
        IR1hdr.set('PV2_3 ', 2.0)

        IR1thdr = listIR1[1].header
        IR1thdr.set('extname', 'WCS-TAB')
        IR1ihdr = listIR1[2].header
        IR1ihdr.set('extname', 'INCIDENCE')
        IR1ehdr = listIR1[3].header
        IR1ehdr.set('extname', 'EMISSION')
        IR1pahdr = listIR1[4].header
        IR1pahdr.set('extname', 'PHASE-ANGLE')
        listIR1.writeto(myfitsIR1)

        # SWIR-L channel #
        hduIR2 = fits.PrimaryHDU(bsqIR2)
        tbhduIR2 = fits.BinTableHDU.from_columns([fits.Column(name='COORDS1', unit='deg', dim=tdim, format=str(mysamples*mylines)+'D', array=newlongIR2), fits.Column(name='COORDS2', unit='deg', dim=tdim, format=str(mysamples*mylines)+'D', array=newlatIR2)])
        listIR2 = fits.HDUList([hduIR2,tbhduIR2])
        IR2hdr = listIR2[0].header
        hduIR2.writeto(myfitsIR2)

        # VNIR channel #
        hduV = fits.PrimaryHDU(bsqV)
        tbhduV = fits.BinTableHDU.from_columns([fits.Column(name='COORDS1', unit='deg', dim=tdim, format=str(mysamples*mylines)+'D', array=newlongV), fits.Column(name='COORDS2', unit='deg', dim=tdim, format=str(mysamples*mylines)+'D', array=newlatV)])
        listV = fits.HDUList([hduV,tbhduV])
        Vhdr = listV[0].header
        hduV.writeto(myfitsV)

      else:
        print 'Cannot find image %s!\n' %(myimage)
        exit(1)
  i = i+1
