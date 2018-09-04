#!/usr/bin/env python 

import sys,os,re
import numpy as np
from astropy.io import fits
from struct import *

myversion = 1.0

# HISTORY
# 04/09/2018 ext number as argument
# 16/04/2017 first version 

i = 0
for arg in sys.argv:
  if (i>0):
    if (arg[0] == "-"):
      if (arg[1] == "n"):
        i = i + 1
        nexp = int(sys.argv[i])
    else :
      if (len(arg) > 2):
        myfits = arg + ".fits"
        ### Each hirise fits will contain extensions ###
        hdul = fits.HDUList()
        hdul.append(fits.PrimaryHDU())
        myext = 0
        while (myext<nexp):
          myimg0 = arg + str(myext) + "_0.IMG"
          myimg1 = arg + str(myext) + "_1.IMG"
          if os.path.exists(myimg0) and os.path.exists(myimg1):
            print myimg0
            print myimg1
            ### parsing PDS label ### 
            myfile0 = open(myimg0,'r')
            myfile1 = open(myimg1,'r')
            for line in myfile0:
              if (line == "END\r\n" or line == "END\n"):
                break
              temp = line.split('=')
              if (temp[0].strip(' ') == "^IMAGE"):
                rbytes = int(((temp[1].strip(' ')).split(" "))[0])
              if (temp[0].strip(' ') == "LINES"):
                mylines = int(temp[1].strip(' '))
              if (temp[0].strip(' ') == "LINE_SAMPLES"):
                mysamples = int(temp[1].strip(' '))
              if (temp[0].strip(' ') == "SAMPLE_BITS"):
                sbit = int(temp[1].strip(' '))
              if (temp[0].strip(' ') == "SAMPLE_TYPE"):
                tbyte = (temp[1].strip(' ')).split('_')
                if (tbyte[0] == "LSB" or tbyte[0] == "PC" or tbyte[0] == "VAX"):
                  myarch = 'littleendian'
                if (tbyte[0] == "MSB" or tbyte[0] == "MAC" or tbyte[0] == "SUN"):
                  myarch = 'bigendian'
              if (temp[0].strip(' ') == "LINE_PREFIX_BYTES"):
                prefix = int(temp[1].strip(' '))
              if (temp[0].strip(' ') == "LINE_SUFFIX_BYTES"):
                suffix = int(temp[1].strip(' '))
              continue

            dim = (mylines, mysamples*2)
            if (sbit == 8):
              if (tbyte[1] == "UNSIGNED"):
                mytype = np.uint8
                myform = "B"
              else:
                mytype = np.int8
                myform = "b"
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

          bsq = np.zeros(dim,dtype=mytype)
          form = str(mysamples) + myform
          if myarch=="bigendian":
            form = ">" + str(mysamples) + myform
          if (rbytes != 0):
            myskip = rbytes + prefix
            myfile0.seek(myskip)
            myfile1.seek(myskip)

            bytex = mysamples * sbit / 8

          ### reading data ###
          for y in range(0, mylines):
            contents=myfile0.read(bytex)
            bsq[y,mysamples:mysamples*2]=unpack_from(form, contents)
            contents=myfile1.read(bytex)
            bsq[y,0:mysamples]=unpack_from(form, contents)
            offset=myfile0.tell()+suffix+prefix
            myfile0.seek(offset)
            myfile1.seek(offset)

          myfile0.close()
          myfile1.close()

          myext = myext + 1
          hdul.append(fits.ImageHDU(data=bsq,do_not_scale_image_data=True,uint=True))
          hdr = hdul[myext].header
          hdr.set('ctype1', 'MALG-AZP')
          hdr.set('ctype2', 'MALN-AZP')
          hdr.set('CRPIX1', 1.0)
          hdr.set('CRPIX2', 1.0)
          hdr.set('CRVAL1', (myext-1) * mysamples * 2.)
          hdr.set('CRVAL2', 1.0)
          hdr.set('CD1_1', 1.0)
          hdr.set('CD2_2', 1.0)
          hdr.set('CD1_2', 0.0)
          hdr.set('CD2_1', 0.0)

        ### Writing FITS files ###
        hdul.writeto(myfits)
      else:
        i = i+1
  i = i+1
