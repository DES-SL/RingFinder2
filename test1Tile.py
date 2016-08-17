import pylab as plt
from StandAloneRingFinder import RingFinder as RingFinder
import colorImage
import fitsio
import os
import sys
import exceptions
import string
import getopt
from BlueRings import BlueRings 
plt.ion()
 
class RingFindforTile(exceptions.Exception):
    def __init__(self, ):
        self.MakeTag = 0
        self.inpar = {}
        self.headers = {}
    " Read fits file, extract header parameters and images "    
    def produce(self,infile):
        fitsfile = os.path.normpath(infile)
        fits = fitsio.FITS(fitsfile,'rw')

        prihdr = fits[0].read_header()
        self.ncol = prihdr["NAXIS1"]
        self.nrow = prihdr["NAXIS2"]
        self.objects = []
        self.images = {}
        self.objext = {}
        
        " loop on all hdus and create a list of objects and its extensions"
        for hdu in fits:
            header = hdu.read_header()
            objectN = header["Object"]
            if objectN not in self.objects:
                self.objects.append(objectN)
                self.objext[objectN] = {'g':0,'r':0,'i':0,'z':0}
        fits.close()    

#        print self.objects
        " Loop on objects and select g,r,i images for each "
        fits = fitsio.FITS(fitsfile,'rw')
        ext = 0
        for hdu in fits:
            header = hdu.read_header()
            curObj = header["Object"]
            curBand = header["band"]
            typeE = header["TYPE"] # possible 'IMAGE' 'WEIGHT' 'PSF'  

            if curObj in self.objects:
                if typeE.find("IMAGE") >=0:
                    objind = self.objext.get(curObj)
                    if curBand.find('g') >= 0: objind['g'] = ext
                    if curBand.find('r') >= 0: objind['r'] = ext
                    if curBand.find('i') >= 0: objind['i'] = ext
                    if curBand.find('z') >= 0: objind['z'] = ext
                    self.objext[curObj] = objind

                if typeE.find("WEIGHT") >=0:
                    objind = self.objext.get(curObj)
                    if curBand.find('g') >= 0: objind['gw'] = ext
                    if curBand.find('r') >= 0: objind['rw'] = ext
                    if curBand.find('i') >= 0: objind['iw'] = ext     
                    if curBand.find('z') >= 0: objind['zw'] = ext     
                    self.objext[curObj] = objind

                if typeE.find("PSF") >=0:
                    objind = self.objext.get(curObj)
                    if curBand.find('g') >= 0: objind['gp'] = ext
                    if curBand.find('r') >= 0: objind['rp'] = ext
                    if curBand.find('i') >= 0: objind['ip'] = ext     
                    if curBand.find('z') >= 0: objind['zp'] = ext     
                    self.objext[curObj] = objind


            ext+=1        
        fits.close()
        " loop on objects and create PNG images "
        fits = fitsio.FITS(fitsfile,'rw')
        interesting=[ "3017263885", "3017249430", "3017265746", "3017263357", "3017255895", "3017261463", "3017263527"]

        for objectN in self.objects:
            objext = self.objext[objectN]
            outDir = "./gallery/"
            if not os.path.exists(outDir):
                os.makedirs(outDir)
            imdict = {}
            sigdict= {}
            psfdict= {}

            for s in ['g','r','i','z']:
                sw=s+"w"
                sp=s+"p"
                a=90

                ext=objext[s]
                header = fits[ext].read_header()
                imdict[s] = fits[ext].read()[a:-a,a:-a]

                ext=objext[sw]
                header = fits[ext].read_header()
                sigdict[s] = (fits[ext].read()[a:-a,a:-a])**-2

                ext=objext[sp]
                header = fits[ext].read_header()
                psfdict[s] = fits[ext].read()
            print objectN
            #["3017249430"]

            #BR=BlueRings(imdict,sigdict,psfdict)
            
            #if objectN in interesting or int(objectN)>3017265667:
            if 1==1:
                BR=BlueRings(imdict,sigdict,psfdict)
                if BR.residualAnalyse(.2):
                    inter=BR.plot()
                else:
                    inter=""
                #if inter!="":interesting.append(objectN)
                #print "[",
                #for inter in interesting:
                #    print "\"%s\","%inter,
                #print"]"
                if inter !="": 
                    BR2=BlueRings(imdict,sigdict,psfdict,psfmode='match')
                    BR2.plot()
            """
            plt.subplot(131)
            color = colorImage.ColorImage()
            colorimage = color.createModel(imdict['g'],imdict['r'],imdict['i'])
            plt.imshow(colorimage,interpolation="none")
            plt.subplot(132)

            #psfmode="crossconvolve"
            #RF=RingFinder(imdict['g'],sigdict['i']**-2,sigdict['g']**-2,imdict['i'],psfdict['g'],psfdict['i'],0.265,1e12,1e12,visualize=False,psfmode=psfmode)
            #plt.imshow(RF.Dshow,interpolation="none")
            #plt.show()

            BR=BlueRings(imdict,sigdict,psfdict)

            colorimage = color.colorize(BR.subdict['g'],BR.subdict['r'],BR.subdict['i'])
            plt.imshow(colorimage,interpolation="none")

            plt.subplot(133)
            BR.subdict['g'][BR.subdict['g']<0]=0
            plt.imshow(BR.subdict['g'],interpolation="none")


            plt.draw()
            raw_input()
            """

if __name__ == "__main__":
    dir_input = "/home/ttemp/BlueRings/data/"
    file_input = "DES0005-0041_cutouts.fits"
    fitsfile = dir_input+file_input

    tile=RingFindforTile()
    tile.produce(fitsfile)
exit() 


a=40
for j in range(100):
    imgg=hdulist[0+j*15].data[a:-a,a:-a]
    sigg=hdulist[1+j*15].data[a:-a,a:-a]
    psfg=hdulist[2+j*15].data
    imgr=hdulist[3+j*15].data[a:-a,a:-a]
    sigr=hdulist[4+j*15].data[a:-a,a:-a]
    psfr=hdulist[5+j*15].data
    imgi=hdulist[6+j*15].data[a:-a,a:-a]
    sigi=hdulist[7+j*15].data[a:-a,a:-a]
    psfi=hdulist[8+j*15].data
    imgz=hdulist[9+j*15].data[a:-a,a:-a]
    sigz=hdulist[10+j*15].data[a:-a,a:-a]
    psfz=hdulist[11+j*15].data

    color = colorImage.ColorImage()
    colorimage = color.createModel(imgg,imgr,imgi)
  #plt.imshow(colorimage,interpolation="none")
  #plt.show()
    psfmode="none"
    RF=RingFinder(imgg,imgi,sigg,sigi,psfg,psfi,0.265,1e12,1e12,visualize=False,psfmode=psfmode)
    RFgz=RingFinder(imgg,imgz,sigg,sigi,psfg,psfz,0.265,1e12,1e12,visualize=False,psfmode=psfmode)
    RFrz=RingFinder(imgr,imgz,sigr,sigz,psfr,psfz,0.265,1e12,1e12,visualize=False,psfmode=psfmode)
    RFiz=RingFinder(imgi,imgz,sigi,sigz,psfi,psfz,0.265,1e12,1e12,visualize=False,psfmode=psfmode)

  #RFres=RF.ringfind(vb=True)
    
    ax=plt.subplot(1,1,1)
    size=colorimage.shape[0]
    
  #RF.Dshow[RF.fitmask]=numpy.nan
    plt.subplot(131)
    color.bMinusr = 0.8
    color.bMinusg = 0.4
    color.nonlin = 1.
    colored=color.createModel(imgg,imgr,imgi)
    plt.imshow(colored,interpolation="none")
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    plt.subplot(132)
    colorresid=color.colorize(RFgz.Dshow,RFrz.Dshow,RFiz.Dshow)
    d=0
  #plt.imshow(RF.Dshow,interpolation="none")
    plt.imshow(colorresid,interpolation="none")
  #plt.xlim([0,size-20])
  #plt.ylim([0,size-20])
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    
    plt.subplot(133)
    plt.imshow(RF.Dshow,interpolation="none")
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    
    plt.show()
    plt.cla()
