
##########################

# WELL CENTRE COORDINATE CALIBRATOR #

# See: github.com/LondonBiofoundry/ColonyPicker

# (Requires an image of a calibration plate, i.e. an agar plate which has been stabbed at (0,0) coordinates by the SELECT head)

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

#####
# OPTIONAL #
# Display diagnositc images?
showimages = False

####
# Load image and do all the things:

# load image as HSV and select saturation
img = cv2.imread(sys.argv[1])
hh, ww, cc = img.shape

# convert to gray
gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)

selem = disk(3) # Increase to ignore larger speckles/more aggressively clean image. Default = 3
gray = closing(gray, selem)
gray = 255-gray

# threshold the grayscale image
# print(np.average(gray))
# ret, thresh = cv2.threshold(gray,np.average(gray)-60,255,0)

selem = disk(40)
local_otsu = rank.otsu(gray, selem)+35
thresh = np.where(gray>(local_otsu), 0, 255)

thresh = cv2.convertScaleAbs(thresh)
thresh = 255-thresh


# find outer contour
cntrs, hierarchy = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)

# cntrs = cntrs[0] if len(cntrs) == 2 else cntrs[1]

# print(cntrs)
out = np.zeros_like(img)
cv2.drawContours(out, cntrs, -1, 255, 3)

join_cnts = np.concatenate(cntrs)

rotrect = cv2.minAreaRect(join_cnts)
box = cv2.boxPoints(rotrect)
box = np.int0(box)

center = rotrect[0]
width = rotrect[1][0]
height = rotrect[1][1]

# draw rotated rectangle on copy of img as result
result = img.copy()
cv2.drawContours(result,[box],0,(0,0,255),2)

# get angle from rotated rectangle
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

print(theta,"deg")

image = result.copy()
shape = ( result.shape[1], result.shape[0] ) # cv2.warpAffine expects shape in (length, height)

matrix = cv2.getRotationMatrix2D( center=center, angle=theta, scale=1 )
image = cv2.warpAffine( src=image, M=matrix, dsize=shape )

SCALEx = int( center[0] - width/2  )
SCALEy = int( center[1] - height/2 )

image = image[SCALEy:int(SCALEy+height), SCALEx:int(SCALEx+width)]
image = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)/255
####

# Read well centre coordinates from file
# Coords are defined during calibration using stab marks from Felix SELECT head
# Calibration is required for each new SELECT head and/or imaging setup
wellcenters = np.genfromtxt('wellcentres.csv', delimiter=',', names=True, dtype=None, encoding='utf8')

# Generate row/column bounds
wellbounds = np.stack((wellcenters['WCx'], wellcenters['WCy']), axis=1)
wellbounds = wellbounds + 100

# Info
print('Filename: ', sys.argv[1])

# Mask plate edges
# Calculate edge widths based on standard proportions of a rectangular plate:
plateH, plateW = image.shape
print('Plate size in pixels: ' + str(plateH) + 'x' + str(plateW))
ratioY1 = 0.09
ratioY2 = 0.9
ratioX1 = 0.065
ratioX2 = 0.92

mask = np.zeros(shape=image.shape[0:2], dtype="bool")
mask[0:int(plateH*ratioY1), :] = True
mask[:, int(plateW*ratioX2):-1] = True
mask[:, 0:int(plateW*ratioX1)] = True
mask[int(plateH*ratioY2):-1, :] = True
# Mask to white-grey
image[mask] = np.average(image[int(plateH*ratioY2), int(plateW*ratioX1):int(plateW*ratioX2)])

# Create mask and threshold each well wrt itself 
# (rather than globally across the plate - this neutralises cross plate gradients)
mask = np.ones(shape=image.shape[0:2], dtype="bool")
for i in range(len(wellcenters)):
	x = wellcenters[i][0]
	y = wellcenters[i][1] - 20 # Correction for the fact that my SELECT head skews north, such that plate edge is a problem at row A
	rr, cc = draw.disk((y, x), radius=85, shape=image.shape[0:2]) # Increase radius if needed; be wary of plate edge
	mask[rr,cc] = False

	# Threshold well:
	thresh = threshold_otsu(image[rr,cc])-0.03
	image[rr,cc] = np.where(image[rr,cc]< thresh, 0.2, 1)

# Mask to white
image[mask] = 1

# Closing function to clean up image
selem = disk(3) # Increase to ignore larger speckles/more aggressively clean image. Default = 3
image2 = closing(image, selem)
# Then grow black regions to form central masses
image2 = erosion(image, selem)

# OPTIONAL: Display image after masking, thresholding and closing:
if showimages:
	io.imshow(image2)
	plt.show()

# Get regions with region properties in image
label_image = label(image2, background=1, return_num=False, connectivity=1)

# Get colony properties - size and centroid coordinates
props = []
for region in regionprops(label_image=label_image):
	cy = region.centroid[0]
	cx = region.centroid[1]
	size = region.area
	props.append((cy, cx, size))

# Set all colonies to same greyscale value for visualisation
visual_label_image = label_image
visual_label_image[visual_label_image > 0] = 1

# Set up plot for later
plt.imshow(visual_label_image, cmap='Blues')

# Sort colonies by size
props.sort(key=lambda tup: tup[2])

# Unzip properties tuple
y,x,size = zip(*props)

# Total colonies
totCol = len(x)

# Set size segments to all
lowerbound = 0
upperbound = len(x)

# 'OP' array will contain Outputs for .csv file
# Columns: Row, Column, Xoffset, Yoffset, Size
OP = []

# List of wells from which a pick has already been identified
prevWell = []
picksX = []
picksY = []

# Get well columns and rows and bin into wells - A1, B1, etc
# Goes by reverse size order so to pick the largest from each well, if single picks is turned on
for i in reversed(range(len(x[lowerbound:upperbound]))):
		j = i + lowerbound
		WellR = 'NA'
		WellC = 0

		# columns
		if x[j] < wellbounds[1][0]:
			WellC = 1
		elif x[j] < wellbounds[9][0]:
			WellC = 2
		elif x[j] < wellbounds[17][0]:
			WellC = 3
		elif x[j] < wellbounds[25][0]:
			WellC = 4
		elif x[j] < wellbounds[33][0]:
			WellC = 5
		elif x[j] < wellbounds[41][0]:
			WellC = 6
		elif x[j] < wellbounds[49][0]:
			WellC = 7
		elif x[j] < wellbounds[57][0]:
			WellC = 8
		elif x[j] < wellbounds[65][0]:
			WellC = 9
		elif x[j] < wellbounds[73][0]:
			WellC = 10
		elif x[j] < wellbounds[81][0]:
			WellC = 11
		elif x[j] < wellbounds[90][0]:
			WellC = 12
		else:
			WellC = 0

		if y[j] < wellbounds[7][1]:
			WellR = 'H'
		elif y[j] < wellbounds[6][1]:
			WellR = 'G'
		elif y[j] < wellbounds[5][1]:
			WellR = 'F'
		elif y[j] < wellbounds[4][1]:
			WellR = 'E'
		elif y[j] < wellbounds[3][1]:
			WellR = 'D'
		elif y[j] < wellbounds[2][1]:
			WellR = 'C'
		elif y[j] < wellbounds[1][1]:
			WellR = 'B'
		elif y[j] < wellbounds[0][1]:
			WellR = 'A'
		else:
			WellR = 'NA'

		# Define picks and append to output array
		WellCode = WellR + str(WellC)
		if WellCode not in prevWell: 
			OP.append((WellR, WellC, x[j], y[j], size[j]))
			prevWell.append(WellCode)
			picksX.append(x[j])
			picksY.append(y[j])

# Optional: View locations of picks
if showimages:
	plt.scatter(picksX, picksY, marker='.', color='red')
	plt.show()

# Sort by row:
OP.sort(key=lambda tup: tup[0])
# And by column:
OP.sort(key=lambda tup: tup[1])

# Optional:
# print(OP)

# unzip tuple
ro,col,xcoord,ycoord,size = zip(*OP)

print('Number of colonies found in the plate: ', str(len(xcoord)), '  <-- this should be 96. Check diagnostic images if not.')

# Write to .csv file
with open('wellcentres.csv', mode='w', newline='') as coordcsv:
	writer = csv.writer(coordcsv, delimiter=',')
	#header
	writer.writerow(['WCx', 'WCy', 'WCr','WCc'])
	
	for i in range(len(ro)):


		writer.writerow([xcoord[i],ycoord[i],ro[i],col[i]])