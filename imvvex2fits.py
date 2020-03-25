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
# 25/3/2020 Python3 compatibility
# 4/8/2018 It works with calibrated imaging data
# 22/3/2018 first version

i = 0
for arg in sys.argv:
    if (i > 0):
        if (arg[0] != "-"):
            myimage = arg
            ext = re.compile("QUB")
            test = ext.match(myimage)
            if (test is None):
                myfits = myimage.replace("CAL", "fits")
                mygeom = myimage.replace("CAL", "GEO")
            else:
                myfits = myimage.replace("QUB", "fits")

            if os.path.exists(myfits):
                print(f"Fits image {myfits} already exists! Please remove it.")
                exit(1)
            if os.path.exists(myimage):
                toparse = myimage
                # parsing PDS label of .qub file #
                myskip = 0
                rbytes = 1
                lastlabrec = 0
                myarch = '>'
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
                        if (temp[0].strip(' ') == "PRODUCER_ID"):
                            author = (temp[1].strip(' ')).strip('\r\n')
                        if (temp[0].strip(' ') == "INSTRUMENT_HOST_NAME"):
                            telescop = (
                                temp[1].strip(' ')).strip('\r\n').strip('"')
                        if (temp[0].strip(' ') == "INSTRUMENT_ID"):
                            instrume = (temp[1].strip(' ')).strip('\r\n')
                        if (temp[0].strip(' ') == "INSTRUMENT_MODE_ID"):
                            instrmode = int((temp[1].strip(' ')).strip('\r\n'))
                        if (temp[0].strip(' ') == "TARGET_NAME"):
                            target = (
                                temp[1].strip(' ')).strip('\r\n').strip('"')
                        if (temp[0].strip(' ') == "RECORD_BYTES"):
                            rbytes = int(temp[1].strip(' '))
                        if (temp[0].strip(' ') == "^QUBE"):
                            mystr = temp[1].strip(' ')
                            tests = re.match('["a-zA-Z]', mystr)
                            if (tests is None):
                                lastlabrec = int(mystr.strip('\r\n')) - 1
                        if (temp[0].strip(' ') == "CORE_ITEMS"):
                            mydim = (temp[1].strip(' \r\n')).split(',')
                            mybands = int(re.sub(r'[^\w]', ' ', mydim[0]))
                            mysamples = int(re.sub(r'[^\w]', ' ', mydim[1]))
                            mylines = int(re.sub(r'[^\w]', ' ', mydim[2]))
                        if (temp[0].strip(' ') == "CORE_ITEM_TYPE"):
                            tbyte = (temp[1].strip(' ')).split('_')
                            if (tbyte[0] == "LSB" or tbyte[0] == "PC" or
                                    tbyte[0] == "VAX"):
                                myarch = '<'
                            if (tbyte[0] == "MSB" or tbyte[0] == "MAC" or
                                    tbyte[0] == "SUN"):
                                myarch = '>'
                        if (temp[0].strip(' ') == "CORE_ITEM_BYTES"):
                            imbytes = int(temp[1].strip(' '))
                            sbit = 8 * imbytes
                        if (temp[0].strip(' ') == "CORE_NULL"):
                            blank = float(temp[1].strip(' '))
                        if (temp[0].strip(' ') == "CORE_VALID_MINIMUM"):
                            datamin = float(temp[1].strip(' '))
                        if (temp[0].strip(' ') == "CORE_UNIT"):
                            bunit = temp[1].strip(' \r\n').strip('"')
                        if (temp[0].strip(' ') == "SUFFIX_BYTES"):
                            sbytes = int(temp[1].strip(' '))
                        if (temp[0].strip(' ') == "SUFFIX_ITEMS"):
                            mydim = (temp[1].strip(' \r\n')).split(',')
                            sbands = int(re.sub(r'[^\w]', ' ', mydim[0]))
                            ssamples = int(re.sub(r'[^\w]', ' ', mydim[1]))
                            slines = int(re.sub(r'[^\w]', ' ', mydim[2]))
                        if (temp[0].strip(' ') == "LINE_SUFFIX_NAME"):
                            names = (temp[1].strip(' \r\n')).split(',')
                        if (temp[0].strip(' ') == "LINE_SUFFIX_UNIT"):
                            units = (temp[1].strip(' \r\n')).split(',')

                    dim = (mybands, mylines, mysamples)
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
                    form = myarch + str(mybands) + myform
                    if ((rbytes != 0) and (lastlabrec != 0)):
                        myskip = rbytes * lastlabrec
                    myfile.seek(myskip)
                    bytex = mybands * imbytes
                    try:
                        sbytes
                    except Exception:
                        sbytes = 0

                    # reading data #
                    bsq = np.zeros(dim, dtype=mytype)
                    wl = np.zeros(mybands, dtype=mytype)
                    fwhm = np.zeros(mybands, dtype=mytype)
                    uncrtnt = np.zeros(mybands, dtype=mytype)
                    for y in range(0, mylines):
                        for x in range(0, mysamples):
                            contents = myfile.read(bytex)
                            bsq[:, y, x] = unpack_from(form, contents)
                            offset = myfile.tell() + (sbands*sbytes)
                            myfile.seek(offset)

                    contents = myfile.read(bytex)
                    wl = unpack_from(form, contents)
                    contents = myfile.read(bytex)
                    fwhm = unpack_from(form, contents)
                    contents = myfile.read(bytex)
                    uncrtnt = unpack_from(form, contents)

                    if os.path.exists(mygeom):
                        toparse = mygeom
                        # parsing PDS label from geometry cube #
                        myskip = 0
                        rbytes = 1
                        lastlabrec = 0
                        myarch = '>'
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
                                        lastlabrec = int(mystr.strip(
                                            '\r\n')) - 1
                                if (temp[0].strip(' ') == "CORE_ITEMS"):
                                    mydim = (temp[1].strip(' \r\n')).split(',')
                                    mybands = int(re.sub(r'[^\w]', ' ',
                                                         mydim[0]))
                                    mysamples = int(re.sub(r'[^\w]', ' ',
                                                           mydim[1]))
                                    mylines = int(re.sub(r'[^\w]', ' ',
                                                         mydim[2]))
                                if (temp[0].strip(' ') == "CORE_ITEM_TYPE"):
                                    tbyte = (temp[1].strip(' ')).split('_')
                                    if (tbyte[0] == "LSB" or
                                            tbyte[0] == "PC" or
                                            tbyte[0] == "VAX"):
                                        myarch = '<'
                                    if (tbyte[0] == "MSB" or
                                            tbyte[0] == "MAC" or
                                            tbyte[0] == "SUN"):
                                        myarch = '>'
                                if (temp[0].strip(' ') == "CORE_ITEM_BYTES"):
                                    imbytes = int(temp[1].strip(' '))

                            # Reading coordinates from geometry cube #
                            if ((rbytes != 0) and (lastlabrec != 0)):
                                myskip = rbytes * lastlabrec
                            myfile.seek(myskip)
                            dim = (mylines, mysamples, 2)
                            data = np.zeros((mybands, mysamples, mylines),
                                            dtype=np.int32)
                            coords = np.zeros(dim, dtype=np.int32)
                            coordsa = np.zeros(dim, dtype=np.int32)
                            coordsb = np.zeros(dim, dtype=np.int32)
                            coordsul = np.zeros(dim, dtype=np.int32)
                            coordsaul = np.zeros(dim, dtype=np.int32)
                            coordsur = np.zeros(dim, dtype=np.int32)
                            coordsaur = np.zeros(dim, dtype=np.int32)
                            coordsbr = np.zeros(dim, dtype=np.int32)
                            coordsabr = np.zeros(dim, dtype=np.int32)
                            coordsbl = np.zeros(dim, dtype=np.int32)
                            coordsabl = np.zeros(dim, dtype=np.int32)
                            inc = np.zeros((mylines, mysamples),
                                           dtype=np.int32)
                            em = np.zeros((mylines, mysamples),
                                          dtype=np.int32)
                            pa = np.zeros((mylines, mysamples),
                                          dtype=np.int32)
                            incc = np.zeros((mylines, mysamples),
                                            dtype=np.int32)
                            emc = np.zeros((mylines, mysamples),
                                           dtype=np.int32)
                            pac = np.zeros((mylines, mysamples),
                                           dtype=np.int32)
                            topo = np.zeros((mylines, mysamples),
                                            dtype=np.int32)
                            topoc = np.zeros((mylines, mysamples),
                                             dtype=np.int32)
                            slant = np.zeros((mylines, mysamples),
                                             dtype=np.int32)
                            lt = np.zeros((mylines, mysamples),
                                          dtype=np.int32)
                            bytex = mybands * imbytes
                            form = myarch + str(mybands) + 'i'
                            for y in range(0, mylines):
                                for x in range(0, mysamples):
                                    contents = myfile.read(bytex)
                                    data[:, x, y] = unpack_from(form, contents)
                                coordsul[y, :, 0] = data[0, :, y]
                                coordsur[y, :, 0] = data[1, :, y]
                                coordsbr[y, :, 0] = data[2, :, y]
                                coordsbl[y, :, 0] = data[3, :, y]
                                coordsul[y, :, 1] = data[4, :, y]
                                coordsur[y, :, 1] = data[5, :, y]
                                coordsbr[y, :, 1] = data[6, :, y]
                                coordsbl[y, :, 1] = data[7, :, y]
                                coords[y, :, 0] = data[8, :, y]
                                coords[y, :, 1] = data[9, :, y]
                                inc[y, :] = data[10, :, y]
                                em[y, :] = data[11, :, y]
                                pa[y, :] = data[12, :, y]
                                topo[y, :] = data[13, :, y]
                                slant[y, :] = data[14, :, y]
                                lt[y, :] = data[15, :, y]
                                coordsaul[y, :, 0] = data[16, :, y]
                                coordsaur[y, :, 0] = data[17, :, y]
                                coordsabr[y, :, 0] = data[18, :, y]
                                coordsabl[y, :, 0] = data[19, :, y]
                                coordsaul[y, :, 1] = data[20, :, y]
                                coordsaur[y, :, 1] = data[21, :, y]
                                coordsabr[y, :, 1] = data[22, :, y]
                                coordsabl[y, :, 1] = data[23, :, y]
                                coordsa[y, :, 0] = data[24, :, y]
                                coordsa[y, :, 1] = data[25, :, y]
                                incc[y, :] = data[26, :, y]
                                emc[y, :] = data[27, :, y]
                                pac[y, :] = data[28, :, y]
                                topoc[y, :] = data[29, :, y]
                                coordsb[y, :, 0] = data[30, :, y]
                                coordsb[y, :, 1] = data[31, :, y]
                            factor = 0.0001
                    else:
                        print(f'Cannot find geometry for image {myimage}!\n')
                        continue
                # writing FITS #

                # TAB formatting #
                form = str(2*mysamples*mylines)
                ttype = "COORDS"
                tform = form+'J'
                tunit = 'deg'
                tdim = '(2,' + str(mysamples) + ',' + str(mylines) + ')'

                # WCS-TAB header info #
                head = fits.Header()

                # Lat Long IAU2000
                head.set('wcsname', 'Body-fixed planetocentric')
                head.set('ctype1', 'VELN-TAB')
                head.set('ctype2', 'VELT-TAB')
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

                # Lat Long on 60Km cloud layer reference
                head.set('wcsnamea', '60Km Venus cloud layer reference')
                head.set('ctype1a', 'VELN-TAB')
                head.set('ctype2a', 'VELT-TAB')
                head.set('CRPIX1a', 1.0)
                head.set('CRPIX2a', 1.0)
                head.set('CRVAL1a', 1.0)
                head.set('CRVAL2a', 1.0)
                head.set('CD1_1a', 1.0)
                head.set('CD2_2a', 1.0)
                head.set('CD1_2a', 0.0)
                head.set('CD2_1a', 0.0)
                head.set('PC1_1a', 1.0)
                head.set('PC2_2a', 1.0)
                head.set('PC1_2a', 0.0)
                head.set('PC2_1a', 0.0)
                head.set('CDELT1a', 1.0)
                head.set('CDELT2a', 1.0)
                head.set('PS1_0a ', 'WCS-CTAB')
                head.set('PS1_1a ', 'COORDS ')
                head.set('PS2_0a ', 'WCS-CTAB')
                head.set('PS2_1a ', 'COORDS ')
                head.set('PV1_1a ', 1)
                head.set('PV2_1a ', 1)
                head.set('PV1_2a ', 1)
                head.set('PV2_2a ', 1)
                head.set('PV1_3a ', 1)
                head.set('PV2_3a ', 2)

                # RA DEC J2000
                head.set('wcsnameb', 'RA DEC J2000')
                head.set('ctype1b', 'RA---TAB')
                head.set('ctype2b', 'DEC--TAB')
                head.set('CRPIX1b', 1.0)
                head.set('CRPIX2b', 1.0)
                head.set('CRVAL1b', 1.0)
                head.set('CRVAL2b', 1.0)
                head.set('CD1_1b', 1.0)
                head.set('CD2_2b', 1.0)
                head.set('CD1_2b', 0.0)
                head.set('CD2_1b', 0.0)
                head.set('PC1_1b', 1.0)
                head.set('PC2_2b', 1.0)
                head.set('PC1_2b', 0.0)
                head.set('PC2_1b', 0.0)
                head.set('CDELT1b', 1.0)
                head.set('CDELT2b', 1.0)
                head.set('PS1_0b ', 'WCS-WTAB')
                head.set('PS1_1b ', 'COORDS ')
                head.set('PS2_0b ', 'WCS-WTAB')
                head.set('PS2_1b ', 'COORDS ')
                head.set('PV1_1b ', 1)
                head.set('PV2_1b ', 1)
                head.set('PV1_2b ', 1)
                head.set('PV2_2b ', 1)
                head.set('PV1_3b ', 1)
                head.set('PV2_3b ', 2)

                hdu = fits.PrimaryHDU(bsq)
                tbhdu = fits.BinTableHDU.from_columns(
                             [fits.Column(name=ttype, unit=tunit, dim=tdim,
                                          format=tform, array=[coords])])
                tbahdu = fits.BinTableHDU.from_columns(
                             [fits.Column(name=ttype, unit=tunit, dim=tdim,
                                          format=tform, array=[coordsa])])
                tbbhdu = fits.BinTableHDU.from_columns(
                             [fits.Column(name=ttype, unit=tunit, dim=tdim,
                                          format=tform, array=[coordsb])])
                # Wavelength TAB formatting #
                wtform = '1E'
                tbwhdu = fits.BinTableHDU.from_columns(
                             [fits.Column(name=names[0].strip('("'),
                                          unit=units[0].strip('("'),
                                          format=wtform, array=wl),
                              fits.Column(name=names[1].strip('"'),
                                          unit=units[1].strip('"'),
                                          format=wtform, array=fwhm),
                              fits.Column(name=names[2].strip('")'),
                                          unit=units[2].strip('")'),
                                          format=wtform, array=uncrtnt)])
                hduinc = fits.ImageHDU(inc)
                hduem = fits.ImageHDU(em)
                hdupa = fits.ImageHDU(pa)
                hduincc = fits.ImageHDU(incc)
                hduemc = fits.ImageHDU(emc)
                hdupac = fits.ImageHDU(pac)
                hdutopo = fits.ImageHDU(topo)
                hdutopoc = fits.ImageHDU(topoc)
                hduslant = fits.ImageHDU(slant)
                hdult = fits.ImageHDU(lt)
                listhdu = fits.HDUList([hdu, tbhdu, tbahdu, tbbhdu, tbwhdu,
                                        hduinc, hduem, hdupa, hdutopo,
                                        hduslant, hdult, hduincc, hduemc,
                                        hdupac, hdutopoc])

                hdr = listhdu[0].header
                if (mytype == np.float32 or mytype == np.float64):
                    bsq[bsq == blank] = np.nan
                else:
                    hdr['BLANK'] = blank
                hdr['DATAMIN'] = datamin
                hdr['BUNIT'] = bunit
                hdr['TELESCOP'] = telescop
                hdr['AUTHOR'] = author
                hdr['OBJECT'] = target
                hdr['INSTRUME'] = instrume
                hdr += head

                thdr = listhdu[1].header
                tahdr = listhdu[2].header
                tbhdr = listhdu[3].header
                tbwhdr = listhdu[4].header
                thdr.set('extname', 'WCS-TAB')
                tahdr.set('extname', 'WCS-CTAB')
                tbhdr.set('extname', 'WCS-WTAB')
                tbwhdr.set('extname', 'WAVELENGTH')
                thdr.set('TSCAL1', factor)
                tahdr.set('TSCAL1', factor)
                tbhdr.set('TSCAL1', factor)

                listhdu[5].header.set('extname', 'INCIDENCE')
                listhdu[5].header.set('BSCALE', factor)
                listhdu[6].header.set('extname', 'EMERGENCE')
                listhdu[6].header.set('BSCALE', factor)
                listhdu[7].header.set('extname', 'PHASE')
                listhdu[7].header.set('BSCALE', factor)
                listhdu[8].header.set('extname', 'TOPOGRAPHY')
                listhdu[8].header.set('BSCALE', factor)
                listhdu[9].header.set('extname', 'SLANT ANGLE')
                listhdu[9].header.set('BSCALE', factor)
                listhdu[10].header.set('extname', 'LOCAL TIME')
                listhdu[10].header.set('BSCALE', factor)
                listhdu[11].header.set('extname', 'CLOUD INCIDENCE')
                listhdu[11].header.set('BSCALE', factor)
                listhdu[12].header.set('extname', 'CLOUD EMERGENCE')
                listhdu[12].header.set('BSCALE', factor)
                listhdu[13].header.set('extname', 'CLOUD PHASE')
                listhdu[13].header.set('BSCALE', factor)
                listhdu[14].header.set('extname', 'CLOUD TOPOGRAPHY')
                listhdu[14].header.set('BSCALE', factor)

                listhdu.writeto(myfits)

        else:
            print(f'Cannot find image {myimage}!\n')
            exit(1)

    i = i+1
