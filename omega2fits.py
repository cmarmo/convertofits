#!/usr/bin/env python

import sys
import os
import re

import numpy as np
import codecs as cd

from astropy.io import fits
from struct import unpack_from

myversion = 2.0

# HISTORY
# 25/3/2020 python3 compatibility
# 10/3/2017 first version

i = 0
for arg in sys.argv:
    if (i > 0):
        if (arg[0] != "-"):
            myimage = arg

            # Each omega cube is a three channel cube #
            myfitsIR1 = myimage.replace("qub", "IR1.fits")
            myfitsIR2 = myimage.replace("qub", "IR2.fits")
            myfitsV = myimage.replace("qub", "V.fits")
            mygeom = myimage.replace("qub", "nav")

            if os.path.exists(myfitsIR1):
                print(f"Fits image {myfitsIR1} already exists! \
                    Please remove it.")
                exit(1)
            if os.path.exists(myfitsIR2):
                print(f"Fits image {myfitsIR2} already exists! \
                    Please remove it.")
                exit(1)
            if os.path.exists(myfitsV):
                print(f"Fits image {myfitsV} already exists! \
                    Please remove it.")
                exit(1)

            if os.path.exists(myimage):
                toparse = myimage
                # parsing PDS label of .qub file #
                myskip = 0
                rbytes = 1
                lastlabrec = 0
                myarch = 'littleendian'
                lablines = 0
                with open(toparse, 'rb') as myfile:
                    header = ''
                    bt = myfile.read(1)
                    while (bt != b'\r\n'):
                        try:
                            header += cd.decode(bt)
                        except Exception:
                            print("end of the cube header.")
                            break
                        bt = myfile.read(1)

                    lines = header.split('\r\n')
                    for line in lines:
                        lablines = lablines + 1
                        temp = line.split('=')
                        if (temp[0].strip(' ') == "RECORD_BYTES"):
                            rbytes = int(temp[1].strip(' '))
                        if (temp[0].strip(' ') == "^QUBE"):
                            mystr = temp[1].strip(' ')
                            tests = re.match('["a-zA-Z]', mystr)
                            if (tests is None):
                                lastlabrec = int(mystr.strip('\r\n')) - 1
                        if (temp[0].strip(' ') == "CORE_ITEMS"):
                            mydim = (temp[1].strip(' \r\n')).split(',')
                            mysamples = int(re.sub(r'[^\w]', ' ', mydim[0]))
                            mybands = int(re.sub(r'[^\w]', ' ', mydim[1]))
                            mylines = int(re.sub(r'[^\w]', ' ', mydim[2]))
                        if (temp[0].strip(' ') == "CORE_ITEM_TYPE"):
                            tbyte = (temp[1].strip(' ')).split('_')
                            if (tbyte[0] == "LSB" or tbyte[0] == "PC" or
                                    tbyte[0] == "VAX"):
                                myarch = 'littleendian'
                            if (tbyte[0] == "MSB" or tbyte[0] == "MAC" or
                                    tbyte[0] == "SUN"):
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
                        if (len(tbyte) > 2 and tbyte[2] == "INTEGER"):
                            mytype = np.int32
                            myform = "i"
                        else:
                            mytype = np.float32
                            myform = "f"
                    if (sbit == 64):
                        if (len(tbyte) > 2 and tbyte[2] == "INTEGER"):
                            mytype = np.long
                            myform = "l"
                        else:
                            mytype = np.float64
                            myform = "d"

                    bsqIR1 = np.zeros(dimIR, dtype=mytype)
                    bsqIR2 = np.zeros(dimIR, dtype=mytype)
                    bsqV = np.zeros(dimV, dtype=mytype)
                    form = str(mysamples) + myform
                    if ((rbytes != 0) and (lastlabrec != 0)):
                        myskip = rbytes * lastlabrec
                    myfile.seek(myskip)
                    bytex = mysamples * imbytes
                    if not sidebytes:
                        sidebytes = 0
                    if not backbytes:
                        backbytes = 0

                    # reading data #
                    for y in range(0, mylines):
                        for z in range(0, mybands):
                            contents = myfile.read(bytex)
                            if (z < IRbands):
                                bsqIR1[z, y, :] = unpack_from(form, contents)
                            elif (z < 2*IRbands):
                                bsqIR2[z-IRbands, y, :] = unpack_from(
                                    form, contents)
                            else:
                                bsqV[z-(2*IRbands), y, :] = unpack_from(
                                    form, contents)
                            offset = myfile.tell()+backbytes
                            myfile.seek(offset)

                        offset = myfile.tell()+(sidebytes*mysamples*sbands)
                        myfile.seek(offset)

                if os.path.exists(mygeom):
                    toparse = mygeom
                    # parsing PDS label from geometry cube #
                    myskip = 0
                    rbytes = 1
                    lastlabrec = 0
                    myarch = 'littleendian'
                    lablines = 0
                    with open(toparse, 'rb') as myfile:
                        header = ''
                        bt = myfile.read(1)
                        while (bt != b'\r\n'):
                            try:
                                header += cd.decode(bt)
                            except Exception:
                                print("end of the geometry header.")
                                break
                            bt = myfile.read(1)

                        lines = header.split('\r\n')
                        for line in lines:
                            lablines = lablines + 1
                            temp = line.split('=')
                            if (temp[0].strip(' ') == "RECORD_BYTES"):
                                rbytes = int(temp[1].strip(' '))
                            if (temp[0].strip(' ') == "^QUBE"):
                                mystr = temp[1].strip(' ')
                                tests = re.match('["a-zA-Z]', mystr)
                                if (tests is None):
                                    lastlabrec = int(mystr.strip('\r\n')) - 1
                            if (temp[0].strip(' ') == "CORE_ITEM_BYTES"):
                                imbytes = int(temp[1].strip(' '))

                        # Reading coordinates from geometry cube #
                        if ((rbytes != 0) and (lastlabrec != 0)):
                            myskip = rbytes * lastlabrec
                        myfile.seek(myskip)
                        dim = (mylines, mysamples, 2)
                        coordsIR1 = np.zeros(dim, dtype=np.int32)
                        incIR1 = np.zeros((mylines, mysamples), dtype=np.int32)
                        emIR1 = np.zeros((mylines, mysamples), dtype=np.int32)
                        paIR1 = np.zeros((mylines, mysamples), dtype=np.int32)
                        coordsIR2 = np.zeros(dim, dtype=np.int32)
                        incIR2 = np.zeros((mylines, mysamples), dtype=np.int32)
                        emIR2 = np.zeros((mylines, mysamples), dtype=np.int32)
                        paIR2 = np.zeros((mylines, mysamples), dtype=np.int32)
                        coordsV = np.zeros(dim, dtype=np.int32)
                        incV = np.zeros((mylines, mysamples), dtype=np.int32)
                        emV = np.zeros((mylines, mysamples), dtype=np.int32)
                        paV = np.zeros((mylines, mysamples), dtype=np.int32)
                        bytex = mysamples * imbytes
                        form = str(mysamples) + 'i'
                        for y in range(0, mylines):
                            offset = myfile.tell()+(6*bytex)
                            myfile.seek(offset)
                            contents = myfile.read(bytex)
                            coordsIR1[y, :, 0] = unpack_from(form, contents)
                            contents = myfile.read(bytex)
                            coordsIR1[y, :, 1] = unpack_from(form, contents)
                            contents = myfile.read(bytex)
                            incIR1[y, :] = unpack_from(form, contents)
                            contents = myfile.read(bytex)
                            emIR1[y, :] = unpack_from(form, contents)
                            contents = myfile.read(bytex)
                            paIR1[y, :] = unpack_from(form, contents)
                            offset = myfile.tell()+(10*bytex)
                            myfile.seek(offset)
                            contents = myfile.read(bytex)
                            coordsIR2[y, :, 0] = unpack_from(form, contents)
                            contents = myfile.read(bytex)
                            coordsIR2[y, :, 1] = unpack_from(form, contents)
                            contents = myfile.read(bytex)
                            incIR2[y, :] = unpack_from(form, contents)
                            contents = myfile.read(bytex)
                            emIR2[y, :] = unpack_from(form, contents)
                            contents = myfile.read(bytex)
                            paIR2[y, :] = unpack_from(form, contents)
                            offset = myfile.tell()+(10*bytex)
                            myfile.seek(offset)
                            contents = myfile.read(bytex)
                            coordsV[y, :, 0] = unpack_from(form, contents)
                            contents = myfile.read(bytex)
                            coordsV[y, :, 1] = unpack_from(form, contents)
                            contents = myfile.read(bytex)
                            incV[y, :] = unpack_from(form, contents)
                            contents = myfile.read(bytex)
                            emV[y, :] = unpack_from(form, contents)
                            contents = myfile.read(bytex)
                            paV[y, :] = unpack_from(form, contents)
                            offset = myfile.tell()+(10*bytex)
                            myfile.seek(offset)
                            factor = 0.0001
                else:
                    print(f'Cannot find geometry for {myimage}!\n')
                    continue

                # Writing FITS files #

                # TAB formatting #
                form = str(2*mysamples*mylines)
                ttype = "COORDS"
                tform = form+'J'
                tunit = 'deg'
                tdim = '(2,' + str(mysamples) + ',' + str(mylines) + ')'

                # WCS-TAB header info #
                head = fits.Header()
                head.set('ctype1', 'MALN-TAB')
                head.set('ctype2', 'MALT-TAB')
                head.set('CRPIX1', 1.0)
                head.set('CRPIX2', 1.0)
                head.set('CRVAL1', 1.0)
                head.set('CRVAL2', 1.0)
                head.set('CD1_1', 1.0)
                head.set('CD2_2', 1.0)
                head.set('CD1_2', 0.0)
                head.set('CD2_1', 0.0)
                head.set('PC1_1', 1.0)
                head.set('PC2_2', 1.0)
                head.set('PC1_2', 0.0)
                head.set('PC2_1', 0.0)
                head.set('CDELT1', 1.0)
                head.set('CDELT2', 1.0)
                head.set('PS1_0 ', 'WCS-TAB ')
                head.set('PS1_1 ', 'COORDS ')
                head.set('PS2_0 ', 'WCS-TAB ')
                head.set('PS2_1 ', 'COORDS ')
                head.set('PV1_1 ', 1)
                head.set('PV2_1 ', 1)
                head.set('PV1_2 ', 1)
                head.set('PV2_2 ', 1)
                head.set('PV1_3 ', 1)
                head.set('PV2_3 ', 2)

                # SWIR-C channel #
                hduIR1 = fits.PrimaryHDU(bsqIR1)
                tbhduIR1 = fits.BinTableHDU.from_columns(
                    [fits.Column(name=ttype, unit=tunit, dim=tdim,
                                 format=tform, array=[coordsIR1])])
                hduincIR1 = fits.ImageHDU(incIR1)
                hduemIR1 = fits.ImageHDU(emIR1)
                hdupaIR1 = fits.ImageHDU(paIR1)
                listIR1 = fits.HDUList(
                    [hduIR1, tbhduIR1, hduincIR1, hduemIR1, hdupaIR1])
                IR1hdr = listIR1[0].header
                IR1hdr += head
                IR1thdr = listIR1[1].header
                IR1thdr.set('extname', 'WCS-TAB')
                IR1thdr.set('TSCAL1', factor)
                IR1ihdr = listIR1[2].header
                IR1ihdr.set('extname', 'INCIDENCE')
                IR1ihdr.set('BSCALE', factor)
                IR1ehdr = listIR1[3].header
                IR1ehdr.set('extname', 'EMISSION')
                IR1ehdr.set('BSCALE', factor)
                IR1pahdr = listIR1[4].header
                IR1pahdr.set('extname', 'PHASE-ANGLE')
                IR1pahdr.set('BSCALE', factor)
                listIR1.writeto(myfitsIR1)

                # SWIR-L channel #
                hduIR2 = fits.PrimaryHDU(bsqIR2)
                tbhduIR2 = fits.BinTableHDU.from_columns(
                    [fits.Column(name=ttype, unit=tunit, dim=tdim,
                                 format=tform, array=[coordsIR2])])
                hduincIR2 = fits.ImageHDU(incIR2)
                hduemIR2 = fits.ImageHDU(emIR2)
                hdupaIR2 = fits.ImageHDU(paIR2)
                listIR2 = fits.HDUList(
                    [hduIR2, tbhduIR2, hduincIR2, hduemIR2, hdupaIR2])
                IR2hdr = listIR2[0].header
                IR2hdr += head
                IR2thdr = listIR2[1].header
                IR2thdr.set('extname', 'WCS-TAB')
                IR2thdr.set('TSCAL1', factor)
                IR2ihdr = listIR2[2].header
                IR2ihdr.set('extname', 'INCIDENCE')
                IR2ihdr.set('BSCALE', factor)
                IR2ehdr = listIR2[3].header
                IR2ehdr.set('extname', 'EMISSION')
                IR2ehdr.set('BSCALE', factor)
                IR2pahdr = listIR2[4].header
                IR2pahdr.set('extname', 'PHASE-ANGLE')
                IR2pahdr.set('BSCALE', factor)
                listIR2.writeto(myfitsIR2)

                # VNIR channel #
                hduV = fits.PrimaryHDU(bsqV)
                tbhduV = fits.BinTableHDU.from_columns(
                    [fits.Column(name=ttype, unit=tunit, dim=tdim,
                                 format=tform, array=[coordsV])])
                hduincV = fits.ImageHDU(incV)
                hduemV = fits.ImageHDU(emV)
                hdupaV = fits.ImageHDU(paV)
                listV = fits.HDUList([hduV, tbhduV, hduincV,
                                      hduemV, hdupaV])
                Vhdr = listV[0].header
                Vhdr += head
                Vthdr = listV[1].header
                Vthdr.set('extname', 'WCS-TAB')
                Vthdr.set('TSCAL1', factor)
                Vihdr = listV[2].header
                Vihdr.set('extname', 'INCIDENCE')
                Vihdr.set('BSCALE', factor)
                Vehdr = listV[3].header
                Vehdr.set('extname', 'EMISSION')
                Vehdr.set('BSCALE', factor)
                Vpahdr = listV[4].header
                Vpahdr.set('extname', 'PHASE-ANGLE')
                Vpahdr.set('BSCALE', factor)
                listV.writeto(myfitsV)
        else:
            print(f'Cannot find image {myimage}!\n')
            exit(1)
    i = i+1
