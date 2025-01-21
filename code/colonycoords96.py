##########################

# COLONY IDENTIFICATION ALGORITHM #

# Finds the coordinates of colonies relative to well center coordinates of a 96-well layout
# A Felix liquid handling robot can then use these coordinates to pick and place each colony

# 	Uses a contour finding algorithm to locate a 96-well plate in an image
#	[Global-local rank otsu thresholds -> contour finding -> minimum area rectangle fit -> rotate and crop]
#	Original image is thus auto-modified to a non-rotated plate with no border
# 	Well centers are interpreted from pixel distances (assumes a succesful plate find)
# 	Colonies are found using (new) local Otsu thresholds and region seeking
#	Colonies are sorted by size and can be selected by size percentile

# OPTIONS:
#	[EITHER]	Sort by size percentiles (e.g. pick the smallest 20%, or the middle 20%-80%, or the largest 50%, or all colonies (0%-100%) etc)
#	[OR]		Pick X colonies from each well location (picks the largest colonies available)

# See README on Github #ADD GITHUB LINK

# Call from command prompt using syntax:
#	[EITHER]	python [this script] [image file] [lower bound colony size percentile] [upper bound colony size percentile]
#	[OR]		python [this script] [image file] [Pick X colonies from each well?] [Number of colonies to pick per well]
# 	e.g.	python colonycoords96.py colonies.tif 20 80
#	or		python colonycoords96.py colonies2.tif y 3
#				NB: affirmatives are 'y, Y, yes, Yes, YES'
#				If no arguments are specified it will default to picking all colonies

# Input:	Image must be of same pixel density (recommended: 600dpi) as imaging system used for calibration (i.e. pixel size of plate must be the same)
#			Rotations will be corrected (up to +-45deg or so)
#			Must have a small buffer region around the plate
#			Avoid dark regions in other parts of the image - crop/digital zoom to avoid dark-tone borders
#			Calibration to your imaging setup and SELECT head necessary (see GitHub)

# Output:	[filename]_coords.csv
# 			Output .csv is readable by Felix protocol (see Github)

# Author:	mara.donora@gmail.com

# TO DO:
#			Write up stab exclusion section here
#			Test with plates
#			Write up
#			Upload to github + link here

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

######

####
# VARIABLES
####

# Diagnostic images - switch to True to view the image as it goes through the process
showimages = False

# Pixel to mm conversion factor (from camera data)
px2mm = 0.042311

# First Otsu threshold correction: default = 35 (wrt. 0-255 value range)
# This may be increased or decreased if the image you are working with has particularly high or low contrast globally
otsuCorrect1 = 35

# Second Otsu threshold correction: default = -0.03 (wrt. 0-1 value range)
# This may be increased or decreased if the image you are working with has particularly high or low contrast in the colony sites
otsuCorrect2 = -0.05

# First closing radius: default = 3 px
# For higher resolution images/more agressive cleaning during plate finding, increase this number
closing1 = 3

# Second closing radius: default = 4 px
# For higher resolution images/more agressive cleaning during colony finding, increase this number
closing2 = 4


# There is a function to exclude stab marks (at the centre of the wells) from the pick selection
# If there is a (hardware) offset between the device used to place the colonies and the SELECT head,
# use the offsets below to correct for it:
# Centre offset of stabs relative to what the select head expects:
stabx = 0.2 # in mm
staby = -2
# Radius of stab exclusion zone. Set to 0 if no stab exclusion is desired.
# stab tolerance radius
stabtol = 1

# If your SELECT head has an xy-offset like mine, use this value to offset the well zone boundaries
# (used to define which well zone, e.g. F5, the colony is in)
selectheadXcorrection = int(stabx/px2mm)
selectheadYcorrection = int(staby/px2mm)

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
if showimages:
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
# LOCATE WELLS AND THRESHOLD EACH
####

# Read well centre coordinates from file
# (Coords are defined during calibration using stab marks from Felix SELECT head
# Calibration is required for each new SELECT head and/or imaging setup)
# (There is also a well centre generator script to generate centres from your image size using ratios, if preferred)
wellcentres = np.genfromtxt('wellcentres.csv', delimiter=',', names=True, dtype=None, encoding='utf8')

# Generate row/column bounds
wellbounds = np.stack((wellcentres['WCx'], wellcentres['WCy']), axis=1)
wellbounds = wellbounds + 100

# Define whether multiple or single colonies are picked from each well location
uniqueWells = False
pickperwell = 1
default = False
try:
	if sys.argv[2] in ['Yes', 'yes', 'YES', 'y', 'Y']:
		uniqueWells = True
		try:
			pickperwell = int(sys.argv[3])
			print('Picking ', str(pickperwell), ' colonies from each well')
		except:
			print('Picking 1 colony from each well (no number of picks per well specified; default is 1 pick per well)')
	else:
		try:
			print('Size percentiles selected: ', sys.argv[2], '% -> ', sys.argv[3], '%')
		except:
			print('Cannot interpret input arguments - defaulting to picking all colonies (except central stab mark regions)')
			default = True
except:
	print('No input arguments specified - defaulting to picking all colonies (except central stab mark regions)')
	default = True

# Get filename as string for later
filename=(str(sys.argv[1])).split('.')[0]

# Optional: Display cropped image
if showimages:
	io.imshow(image)
	plt.show()

# Create mask and threshold each well wrt itself 

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

# Mask around wells
# And grey out central well to exclude stab marks
mask = np.ones(shape=image.shape[0:2], dtype="bool")
for i in range(len(wellcentres)):
	x = wellcentres[i][0] + selectheadXcorrection
	y = wellcentres[i][1] + selectheadYcorrection # Correction for the fact that my SELECT head skews north, such that plate edge is a problem at row A
	rr, cc = draw.disk((y, x), radius=100, shape=image.shape[0:2]) # Increase radius if needed; be wary of plate edge
	mask[rr,cc] = False
	meanval = np.average(image[rr,cc])
	rr2, cc2 = draw.disk((int(y), int(x)), radius=int(stabtol/px2mm), shape=image.shape[0:2]) # Increase radius if needed; be wary of plate edge
	image[rr2,cc2] = meanval

	# rrcc = np.stack((rr,cc), axis = 1)
	# rrcc2 = np.stack((rr2,cc2), axis = 1)


	# Threshold well:

	# OPTION 1: pure thresholding option (may pick up stab marks or noise in empty wells):
	#########

	thresh = threshold_otsu(image[rr,cc])+otsuCorrect2
	image[rr,cc] = np.where(image[rr,cc]< thresh, 0.2, 1)
	#########

	# OPTION 2: enforce minimum variance condition to ignore empty wells:
	#########
	# pxrange= max(image[rr,cc]) - min(image[rr,cc])
	# if pxrange > 0.25 and np.var(image[rr,cc]) > 0.002:
	# 	thresh = threshold_otsu(image[rr,cc])-0.03
	# 	image[rr,cc] = np.where(image[rr,cc]< thresh, 0.2, 1)
	# else:
	# 	image[rr,cc] = 1
	#########

# Optional: Display image after thresholding with mask visible
if showimages:
	image[mask] = 0
	io.imshow(image)
	plt.show()

# Mask to white
image[mask] = 1

# Closing function to clean up image
selem = disk(closing2) 
image2 = dilation(image, selem)
selem = disk(closing2-1)
image2 = erosion(image2, selem)

# OPTIONAL: Display image after masking, thresholding and closing:
if showimages:
	io.imshow(image2)
	plt.show()

####
# FIND COLONY REGIONS
####

# Get regions with region properties in image
label_image = label(image2, background=1, return_num=False, connectivity=1)

# Get colony properties - size and centroid coordinates
props = []
for region in regionprops(label_image=label_image):
	cy = region.centroid[0]
	cx = region.centroid[1]
	size = region.area
	props.append((cy, cx, size))

# Set up plot for later

# Either: use image with regions identified
# 	Set all colonies to same greyscale value for visualisation
# visual_label_image = label_image
# visual_label_image[visual_label_image > 0] = 1
# plt.imshow(visual_label_image, cmap='Blues')

# Or: use original image
plt.imshow(img)

# Sort colonies by size
props.sort(key=lambda tup: tup[2])

# Unzip properties tuple
y,x,size = zip(*props)

# Total colonies
totCol = len(x)

# Determine size segments to address based on user input
# Ignored if picks per well is chosen by user
if default or uniqueWells:
	lowerbound = 0
	upperbound = len(x)
else:
	lowerbound = math.floor(len(x)*(float(sys.argv[2])/100))
	upperbound = math.floor(len(x)*(float(sys.argv[3])/100))

####
# CHOOSE COLONIES TO PICK AND GET COORDINATES
####

# 'OP' array will contain Outputs for .csv file
# Columns: Row, Column, Xoffset, Yoffset, Size, xpx, ypx
OP = []

# Lists for defining pick coordinates
picksX = []
picksY = []

# Generate a dictionary of wells with the number of picks already assigned:
# (Used if picks per well is indicated by user)
l = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
n = [1,2,3,4,5,6,7,8,9,10,11,12]
keys = []
for letter in l:
	for num in n:
		keys.append(letter + str(num))
d = dict((key, 0) for key in keys)
# Yields a dictionary of the form {A1 : 0, A2 : 0, ... , H12 : 0}

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

		# Calculate relative offsets in millimeters
		relative_x = 0
		relative_y = 0
		for k in range(len(wellcentres)):
			if wellcentres[k][2] == WellR and wellcentres[k][3] == WellC:
				relative_x = px2mm*(x[j]-wellcentres[k][0])
				relative_y = px2mm*(y[j]-wellcentres[k][1])

				# Avoid small numbers (the exponent notation probably doesn't vibe with the felix code)
				if abs(relative_x) < 0.01:
					relative_x = 0
				if abs(relative_y) < 0.01:
					relative_y = 0

				stab_flag = False
				# Avoid picking up stab marks (ignores centre marks within a +- 1mm tolerance)
				deviation = math.sqrt((relative_x-stabx)**2 + (relative_y-staby)**2)
				if  deviation < stabtol:
				# if -1 < relative_x and relative_x < 1.5 and -2 < relative_y and relative_y < 0.5: #NB our select head has an inherent x/y-offet
					stab_flag = True
					

				# Optional: Diagnostic print blocks

				# print('WellR=', str(WellR),'; WellC=', str(WellC))
				# print('x[j]=', str(x[j]), '; y[j]=', str(y[j]))
				# print('wellcentres[k][0]=', str(wellcentres[k][0]), '; wellcentres[k][1]=', str(wellcentres[k][1]))
				# print('relative_x=', str(relative_x), '; relative_y=', str(relative_y))
				# print('\n')

		# Define picks and append to output array
		WellCode = WellR + str(WellC)
		if uniqueWells:
			if d[WellCode] < pickperwell and not stab_flag:
				OP.append((WellR, WellC, relative_x, relative_y, size[j], x[j], y[j]))
				picksX.append(x[j])
				picksY.append(y[j])
				d[WellCode] += 1
		else:
			if not stab_flag:
				OP.append((WellR, WellC, relative_x, relative_y, size[j], x[j], y[j]))
				picksX.append(x[j])
				picksY.append(y[j])



####
# GENERATE .CSV FILE OF PICK COORDINATES
####

# Sort by row:
OP.sort(key=lambda tup: tup[0])
# And by column:
OP.sort(key=lambda tup: tup[1])

# unzip tuple
ro,col,xcoord,ycoord,size,picksx,picksy = zip(*OP)

# Create pick map
plt.scatter(picksx, picksy, marker='.', color='red', s=1)
for i in range(len(OP)):
	plt.annotate(str(i+1), (picksx[i],picksy[i]), size=3)
plt.savefig(filename+'_picks.png', dpi = 1200)

print('Number of colonies found in the plate: ', str(totCol), '; Number of picks: ', str(len(ro)))
print(filename+'_picks.png generated')

# Optional: View locations of picks
if showimages:
	plt.scatter(picksx, picksy, marker='.', color='red', s=5)
	plt.show()

# Write to .csv file
with open(filename + '_coords.csv', mode='w', newline='') as coordcsv:
	writer = csv.writer(coordcsv, delimiter=',')
	#header
	writer.writerow(['pickID', 'ro', 'col', 'xcoord','ycoord', 'size', 'destination plate', 'destination row', 'destination column'])
	
	for i in range(len(ro)):

		# Define destination plate, row and column
		# (This is not read by the Felix, but should match where each colony was placed for future reference)
		destinationplate = 1 + math.floor(i/96)
		destinationcolumn = 1 + (math.floor(i/8))%12
		destinationrownum = 1 + (i%8)

		if destinationrownum == 8:
			destinationrow = 'H'
		elif destinationrownum == 7:
			destinationrow = 'G'
		elif destinationrownum == 6:
			destinationrow = 'F'
		elif destinationrownum == 5:
			destinationrow = 'E'
		elif destinationrownum == 4:
			destinationrow = 'D'
		elif destinationrownum == 3:
			destinationrow = 'C'
		elif destinationrownum == 2:
			destinationrow = 'B'
		elif destinationrownum == 1:
			destinationrow = 'A'
		else:
			destinationrow = 'NA'

		writer.writerow([i+1,ro[i],col[i],xcoord[i],ycoord[i],size[i],destinationplate, destinationcolumn, destinationrow])

print(filename+'_coords.csv generated')