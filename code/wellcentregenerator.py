
##########################

# WELL CENTRE COORDINATE GENERATOR #

# See: github.com/LondonBiofoundry/ColonyPicker

# Author:	mara.donora@gmail.com

##########################

from skimage import io, draw
from skimage.filters import threshold_otsu, threshold_local, rank
from skimage.morphology import closing, disk, erosion, dilation
from skimage.measure import label, regionprops
from skimage.feature import match_template

import numpy as np
import sys
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import time
import csv
import math

import cv2



####
# VARIABLES
####

# Pixel to mm conversion factor (from camera data)
px2mm = 0.042311

# First Otsu threshold correction: default = 35 (wrt. 0-255 value range)
# This may be increased or decreased if the image you are working with has particularly high or low contrast globally
otsuCorrect1 = 35

# Second Otsu threshold correction: default = -0.03 (wrt. 0-1 value range)
# This may be increased or decreased if the image you are working with has particularly high or low contrast in the colony sites
otsuCorrect2 = -0.03

# First closing radius: default = 3 px
# For higher resolution images/more agressive cleaning during plate finding, increase this number
closing1 = 3

# Second closing radius: default = 3 px
# For higher resolution images/more agressive cleaning during colony finding, increase this number
closing2 = 3

# If your SELECT head has a y-offset like mine, use this value to offset the well zone boundaries
# (used to define which well zone, e.g. F5, the colony is in)
selectheadYcorrection = 20

# There is a function to exclude stab marks (at the centre of the wells) from the pick selection
# If there is a (hardware) offset between the device used to place the colonies and the SELECT head,
# use the offsets below to correct for it:
# Centre offset of stabs relative to what the select head expects:
stabx = 0.2 # in mm
staby = -1
# Radius of stab exclusion zone. Set to 0 if no stab exclusion is desired.
# stab tolerance radius
stabtol = 1.2

####
# LOCATE 96-WELL PLATE
####

# Load image
img = cv2.imread(sys.argv[1])

# Convert to gray
gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)

# Clean image of speckle
selem = disk(closing1)
gray = closing(gray, selem)

# Invert image
gray = 255-gray

# Local otsu thresholds calculated for whole image
selem = disk(40)
local_otsu = rank.otsu(gray, selem)+ otsuCorrect1

# Create thresholded image
thresh = np.where(gray>(local_otsu), 0, 255)

# Set up image for contour finding (convert from uint32 to uint8 - for some reason this is necessary)
thresh = cv2.convertScaleAbs(thresh)

# Re-invert
thresh = 255-thresh

# Find outer contour
cntrs, hierarchy = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)

# Join contours into a single numpy array
join_cntrs = np.concatenate(cntrs)

# Find minimum area rectangle enclosing the full set of contours
rotrect = cv2.minAreaRect(join_cntrs) # Yields ((centre x, y), (width, height), angle) or ...(height, width)... depending on rotation

# Optional: draw rectangle on image and show result
box = cv2.boxPoints(rotrect)
box = np.int0(box)
result = img.copy()
cv2.drawContours(result,[box],0,(0,0,255),4)
io.imshow(result)
plt.show()

# Get angle from rotated rectangle and center, width and height
theta = rotrect[-1]
if theta > 45:
	theta = theta-90
	center = rotrect[0]
	width = rotrect[1][1]
	height = rotrect[1][0]
else:
	center = rotrect[0]
	width = rotrect[1][0]
	height = rotrect[1][1]

# Info
print('Angle correction: ', theta,"deg")

# cv2.warpAffine expects shape in (length, height)
shape = (img.shape[1], img.shape[0])

# Define rotation matrix and rotate image
matrix = cv2.getRotationMatrix2D( center=center, angle=theta, scale=1 )
img = cv2.warpAffine( src=img, M=matrix, dsize=shape )

# Get crop size
SCALEx = int( center[0] - width/2  )
SCALEy = int( center[1] - height/2 )

# Crop image
img = img[SCALEy:int(SCALEy+height), SCALEx:int(SCALEx+width)]

# Convert to greyscale for colony finding
image = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)/255

####

plateH, plateW = image.shape


ys = []
for i in range(8):
	ys.append(int(plateH*(0.1045*i+0.1566)))

xs = []
for i in range(12):
	xs.append(int(plateW*(0.0697*i+0.1158)))




# print(wellcentres)

wellcentres = np.zeros((96,2))
c = 0
for x in xs:
	for y in reversed(ys):
		wellcentres[c,0] = x
		wellcentres[c,1] = y
		c +=1


with open('wellcentres.csv', mode='w', newline='') as coordcsv:
	writer = csv.writer(coordcsv, delimiter=',')
	#header
	writer.writerow(['WCx', 'WCy', 'Row', 'Column'])
	for i in range(len(wellcentres)):
		col = 1 + (math.floor(i/8))%12
		destinationrownum = 1 + (i%8)

		if destinationrownum == 8:
			row = 'H'
		elif destinationrownum == 7:
			row = 'G'
		elif destinationrownum == 6:
			row = 'F'
		elif destinationrownum == 5:
			row = 'E'
		elif destinationrownum == 4:
			row = 'D'
		elif destinationrownum == 3:
			row = 'C'
		elif destinationrownum == 2:
			row = 'B'
		elif destinationrownum == 1:
			row = 'A'

		writer.writerow([wellcentres[i,0],wellcentres[i,1],row,col])