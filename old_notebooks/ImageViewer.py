from IPython.display import Image, display_png

import lsst.afw.display.ds9 as ds9
import os

def saveRemoteImage(path):
    ds9.ds9Cmd('saveimage ' +path)
    os.system("scp jpreyes@157.253.204.26:/home/jpreyes/"+  path +" /renoir_data_00/jpreyes/stacks/notebook_files/")

def displayFigure(path):
    i =  Image(path)
    display_png(i)

def displayImage(Image,frame=0,title="NoName",path="fig"):
    ds9.mtv(Image, frame=frame, title=title)
    ds9.ds9Cmd("zscale")
    ds9.ds9Cmd("zoom to fit")
    ds9.ds9Cmd("raise")
    ds9.ds9Cmd('saveimage '+path)
    print title
    displayFigure(path)

def showFrame(frame=0,title="NoTitle",path="fig"):
    ds9.ds9Cmd("zscale")
    # ds9.ds9Cmd("zoom to fit")
    ds9.ds9Cmd("raise")
    ds9.ds9Cmd('saveimage '+path)
    print title
    displayFigure(path)
		           
def displayImageWithSources(Image,frame=0,title="NoName",path="fig", sources=[]):
    ds9.mtv(Image, frame=frame, title=title)
    ds9.ds9Cmd("zscale")
    ds9.ds9Cmd("zoom to fit")
    ds9.ds9Cmd("raise")
					       
				           
    for source in sources:
         ds9.dot("+", source.getX()-Image.getX0()  , source.getY()-Image.getY0() , frame=frame, size = 25, ctype = ds9.RED)
 
    ds9.ds9Cmd('saveimage '+path)
    print title
    displayFigure(path)

	    
	   
