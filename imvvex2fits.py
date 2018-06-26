#!/usr/bin/env python 

import sys,os,re
import numpy as np
from astropy.io import fits
from struct import *

myversion = 1.1

# HISTORY
# 8/6/2018 It works with raw data: cubes are just translated in FITS 
# 22/3/2018 first version 

i = 0
for arg in sys.argv:
  if (i>0):
    if (arg[0] != "-"):
      myimage = arg
      ext = re.compile("QUB")
      test = ext.match(myimage)
      if (test == None):
        myfits = myimage.replace("CAL","fits")
        mygeom = myimage.replace("CAL","GEO")
      else:
        myfits = myimage.replace("QUB","fits")

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
          if (temp[0].strip(' ') == "PRODUCER_ID"):
            author = (temp[1].strip(' ')).strip('\r\n')
          if (temp[0].strip(' ') == "INSTRUMENT_HOST_NAME"):
            telescop = (temp[1].strip(' ')).strip('\r\n')
          if (temp[0].strip(' ') == "INSTRUMENT_ID"):
            instrume = (temp[1].strip(' ')).strip('\r\n')
          if (temp[0].strip(' ') == "INSTRUMENT_MODE_ID"):
            instrmode = int((temp[1].strip(' ')).strip('\r\n'))
          if (temp[0].strip(' ') == "TARGET_NAME"):
            target = (temp[1].strip(' ')).strip('\r\n')
          if (temp[0].strip(' ') == "RECORD_BYTES"):
            rbytes = int(temp[1].strip(' '))
          if (temp[0].strip(' ') == "^QUBE"):
            mystr = temp[1].strip(' ')
            tests = re.match('["a-zA-Z]',mystr)
            if (tests == None):
              lastlabrec = int(mystr.strip('\r\n')) - 1
          if (temp[0].strip(' ') == "CORE_ITEMS"):
            mydim = (temp[1].strip(' \r\n')).split(',')
            mybands = int(re.sub(r'[^\w]', ' ', mydim[0]))
            mysamples = int(re.sub(r'[^\w]', ' ', mydim[1]))
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
          if (temp[0].strip(' ') == "SUFFIX_BYTES"):
            sbytes = int(temp[1].strip(' '))
          if (temp[0].strip(' ') == "SUFFIX_ITEMS"):
            mydim = (temp[1].strip(' \r\n')).split(',')
            sbands = int(re.sub(r'[^\w]', ' ', mydim[0]))
            ssamples = int(re.sub(r'[^\w]', ' ', mydim[1]))
            slines = int(re.sub(r'[^\w]', ' ', mydim[2]))
          continue

        dim = (mylines, mysamples, mybands)
        if (sbit == 8):
          mytype = np.int8
          myform = ""
        if (sbit == 16):
          if (tbyte[1] == "UNSIGNED"):
            mytype = np.uint16
            myform = "H"
            bzero = 32768.0
            bscale = 1.
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
        form = str(mybands) + myform
        if ((rbytes != 0) and (lastlabrec != 0)):
          myskip = rbytes * lastlabrec
        myfile.seek(myskip)
        bytex = mybands * imbytes
        try: sbytes
        except:
          sbytes = 0

        ### reading data ###
        if (instrmode==7 and mybands==3456):
          dim7 = (1, 1, mybands)
          bsq = np.zeros(dim7,dtype=mytype)
          bsdark = np.zeros(dim7,dtype=mytype)
          contents=myfile.read(bytex)
          bsq[0,0,:]=unpack_from(form, contents)
          offset=myfile.tell()+ (sbands*sbytes)
          myfile.seek(offset)
          bsdark[0,0,:]=unpack_from(form, contents)        
        else :
          bsq = np.zeros(dim,dtype=mytype)
          for y in range(0, mylines):
            for x in range(0, mysamples):
              contents=myfile.read(bytex)
              bsq[y,x,:]=unpack_from(form, contents)
              offset=myfile.tell()+ (sbands*sbytes)
              myfile.seek(offset)

            offset=myfile.tell()+ sbytes*ssamples*(mybands+sbands)
            myfile.seek(offset)

        myfile.close()

        ### writing FITS ###
        hdu = fits.PrimaryHDU(bsq)
        try:
          hdudark = fits.ImageHDU(bsdark)
          listfits = fits.HDUList([hdu,hdudark])
          hdrdark = hdudark.header
          hdrdark['EXTNAME'] = 'DARK'
        except:
          hdr = hdu.header
        hdr['TELESCOP'] = telescop
        hdr['AUTHOR'] = author
        hdr['OBJECT'] = target
        hdr['INSTRUME'] = instrume
        try:
          listfits.writeto(myfits)
        except:
          hdu.writeto(myfits)

      else:
        print 'Cannot find image %s!\n' %(myimage)
        exit(1)

  i = i+1
