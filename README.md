# ColonyPicker

An automated colony picker using the CyBio Felix liquid handler.

_Author: Mara Donora - https://github.com/donora_

# Table of contents

1. Introduction
2. Overview
   1. Quick how-to
3. Detailed methodology
   1. Imaging
   2. Algorithm
      1. User input
      2. Plate finding
      3. Image conditioning
      4. Colony sorting
   3. Automated picking
   4. Deleting picks
   5. Accuracy
   6. Future work
4. How to get started
   1. Download
   2. Setup
   3. Calibration
   4. Optimisation variables

# Introduction

Colony picking is a common laboratory task which can be very time consuming when performed by hand. Automated colony pickers are available but at a price point which may be unaffordable for many laboratories (around £100k for a colony picker which can perform the task quickly, i.e. around 500 colonies per hour). Many laboratories have general purpose liquid handling robots, which may be able to perform the same task to a good degree of accuracy for a fraction of the price.

This document details a colony picking algorithm which requires a gel dock (or good quality camera) and a CyBio Felix liquid handling robot. The algorithm analyses a photo of an agar plate, finds colonies and identifies the (two-dimensional) centre of mass of each, and translates the pixel locations into relative x-y coordinates. It produces a .csv file of colony coordinates which is read by the Felix liquid handling robot. A Felix protocol directs a SELECT head (which has 8 individually actuatable tips) to pick each colony in turn and place them into recipient 96 well plates. The protocol picks circa 520 colonies per hour, i.e. a performance level comparable to dedicated colony picking devices. The use of tips represents a small running cost, but given the light usage of these consumables it _may_ be possible to sterilize and reuse boxes of tips. Tests showed a picking accuracy of 99% (including potential failures in the incubation stage), which is comparable to the accuracy demonstrated by dedicated colony picking systems.

# Overview

This tool is used to find colonies within an image of an agar plate, determine their relative coordinates, and pick them using the Felix liquid handling robot.

The algorithm is written in Python using the NumPy, scikit-image, OpenCV and Matplotlib libraries and is accessed from the command line. When calling the script, the user points to the image to be analysed, and can optionally pass arguments indicating whether to choose colonies by size percentile or to pick a certain number of colonies per well. The image of the colony plate is analysed using basic computer vision analysis to find the plate and apply rotation and crop transformations to locate the plate centrally and without padding. The colonies are found by image region analysis and their size and 2-D centre of mass are determined. The pixel location of each colony is translated into relative millimetre offsets with respect to the well-centre locations of the 96 well plate. The colonies are sorted by size and position, and, depending on the input arguments, a selection is made based on size or well location.

<img src="/readme_images/overview.PNG" width="800">

A .csv file is produced containing the pickID, colony coordinates (relative to well centres), colony size, and destination well/plate. (Destination well/plate is not utilised by the Felix itself, but indicates where the Felix protocol will place the colonies for later reference). A pick map is also produced for reference. A Felix protocol is run in the Composer software, which reads this .csv file and directs the Felix to the locations of the colonies. The Felix can be loaded with three sets of tips and three destination plates at a time; after these have been filled, the protocol will pause and request the user to replace the tips and plates.

## Quick how-to

_[This assumes you have already followed the Setup instructions, below]_

1. Take an image of your plate (face-down, with H1 in the top left and A1 in the bottom left). Place the image in the ColonyPicker folder.
2. Call the script from the command line, using the image name and other arguments (see below).
3. Check the pick map ([*filename*]\_picks.png) – these are the colonies which will be picked.
4. Copy the output file ([*filename*]\_coords.csv) and place it in the same folder as your Felix protocols.
5. In Composer, point to the new .csv file in the Felix protocol by double-clicking on the _Operate Work List_ function and pointing _Database File_ to [*filename*]\_coords.csv.
6. Set up deck with tips and destination plates (see below). Run protocol. Change tips/plates when necessary, if number of colonies requires it (i.e. every 288 colonies).

# Detailed methodology

## Imaging

A top-down image of the colony plate should be taken on a lightbox with the plate face-down (such that H1 is in the top left corner and A1 in the bottom left). _NB: Resolution circa 600dpi is ideal; many aspects of the algorithm are optimised for this resolution (thresholding, speckle clean-up, closing functions, etc). Variables are located at the top of the code and may be changed by the user if a different resolution is to be used. Pixel-per-mm calibration must be known and changed here - see setup instructions below._

A gel dock is ideal, but home-made gel imagers will also work so long as the camera does not have a fish-eye effect, and the image contrast is sufficient. If pixel-per-mm is not known, camera metadata or imageJ may be used to determine this.

The image should ideally have a small white border around the plate, and dark regions in this border should be avoided. (The plate-finding algorithm performs a local Otsu threshold globally to attempt to ignore amorphous dark regions, but speckle at the border will interfere with the process.) Plate rotations will be corrected up to +- 45 degrees.

<img src="/readme_images/original_plate.png" width="600">


## Algorithm

### _User input_

The algorithm is called from the command line. The user specifies the script to be run, the image file, and two further arguments. These can either define the lower and upper bounds of the size percentiles to be selected, or they can specify a number of colonies to pick from each well.

Either: `python [this script] [image file] [lower bound colony size percentile] [upper bound colony size percentile]`

Or: `python [this script] [image file] [Pick X colonies from each well?] [Number of colonies to pick per well]`

e.g. `> python colonycoords96.py colonies.tif 20 80`

or `> python colonycoords96.py colonies2.tif y 3`

The first command will pick the central 20%-80% of colonies by size. The second will pick three colonies from each well location. _NB: affirmatives are 'y, Y, yes, Yes, YES'. If no arguments are specified the program will default to picking all colonies._

### _Plate finding_

Thresholding is performed with a series of ranked local Otsu thresholds across the whole image. This produces a white image with the edges of the plate visible in black. After inversion, region perimeter contours are found (which define the bounds of the plate), and a minimum-area rectangle is generated to enclose these contours. _NB: the plate-finding part of the algorithm assumes that the padding region around the plate is uniformly white. If there is too much contrast variation in the padding region, the plate-finding algorithm will not work as intended._

<img src="/readme_images/found_plate1.png" width="600">

<img src="/readme_images/found_plate2.png" width="600">

The position and rotation of the rectangle enclosing the plate is used to transform the image, such that an image of the plate with rotation corrected is produced.

<img src="/readme_images/cropped_plate.png" width="600">

_(NB:There is an option at the top of the code to turn on image viewing during the processing steps, which can be useful for troubleshooting.)_

### _Image conditioning_

The image (with padding removed and rotation corrected) is masked with 96 circles to isolate the well regions. A local Otsu threshold is performed in each well region to distinguish colonies from background. The centre of each well is greyed out (at the average value for the well area) to exclude the stab marks created when placing the colonies; the radius and x-y offset of this exclusion zone may be adjusted if necessary, and may be set to zero if no stab marks are present (for example, if the colonies were introduced to the plate by a swipe). An Otsu threshold is performed on the well area. _NB: in the code, an option is commented which sets a minimum variance condition for the each well - i.e. if there is insufficient value variation within the well area the well area will be ignored. By default this feature is not enabled._

An asymmetric closing function (dilation/erosion) is performed to smooth jagged colony edges, remove speckle, and to attempt to better separate close neighbours. Colonies (black pixels on white background) are identified using a region-finding algorithm, which yields their size (in pixels) and the pixel coordinates of their two-dimensional centre of mass. [#A weakness of this algorithm is that joined colonies will be interpreted as one region and a centre of mass identified roughly between the two. A circularity condition could be incorporated into the colony identification process to help separate overlapping colonies.]

<img src="/readme_images/found_colonies.png" width="600">


### _Colony sorting_

The colonies are sorted by size and by well locations (A1, F5 etc). Pixel values of the centre of each well are known from prior calibration (see below – the coordinates are read from the wellcentres.txt file) and are used to determine relative x and y offsets with respect to each colony’s closest well centre. Finally, pixel-per-mm calibration is used to yield a relative offset in mm of each colony centre (which is what the Felix protocol requires).

Depending on the user input, the colonies will either be added to the output file (or skipped) depending on one of two mechanisms: either a size percentile slice is defined (e.g. the smallest 50% of colonies on the plate, or the middle 20%-80% by size, etc); or, X number of colonies from each well location will be picked (in this mode, colonies are selected by largest first).

An image is produced with the locations of the picks overlayed onto the original image. pickID (which is equivalent to the order in which the colonies are picked) is indicated next to each pick. This can be used as a reference to manually remove entries from the .csv if need be (see below).

<img src="/readme_images/pickmap.png" width="600">

<img src="/readme_images/pickmap_zoom.png" width="600">


## Automated picking




https://github.com/user-attachments/assets/14ef3e5a-5bcf-4f64-b8d0-44c43ef72290


<!--
<figure class="video_container">
  <video controls="true" allowfullscreen="true" poster="/readme_images/select_pick.jpg">
    <source src="/readme_images/picking_video_small.mp4" type="video/mp4">
  </video>
</figure>
-->

A protocol written in the Composer program performs the automated picking routine. The .csv file output from the previous step must either be named consistently each time to be read by the protocol, or the program can be edited to select the correct .csv file. The SELECT head is used – this consists of a set of 8 vertically aligned individually actuated pipettes.

<img src="/readme_images/select_pick.jpg" width="600">


First, the SELECT head picks up a column of 8 tips. The first 8 colonies are picked one by one using tips 1-8. The column of tips places the picks into the first column of the first destination plate, before replacing the tips. The next column of tips is picked up, the next 8 colonies are picked. The picks are deposited in the next column of the destination plate, and the tips are replaced.

<img src="/readme_images/felix_deck.jpg" width="600">


Tip boxes are placed in positions 10, 11 and 12, paired with destination plates in positions 7, 8 and 9 respectively. The protocol will pick colonies until the three destination plates are filled (288 picks), at which point it will pause and ask the user to replace the destination plates and tip boxes. It will then continue for another three destination plates, pause and ask the user to replace, etc, until all the colonies designated in the .csv file have been picked.

The colony plate is placed on an adapter which spans positions 1 and 4 (see figure, above). Custom x, y and z offsets are defined in the protocol to create a new ‘position 13’. The plate must be held in between positions 1 and 4 so that the SELECT head can access the full range of positions with the full range of tips (i.e. tip 8 can access row A, and tip 1 can access row H).

The Z-height determines the height at which picking is performed. Different volumes of agar will result in different surface heights – it is important to determine the correct height for your plate.

To ensure the plate is in a consistent x-y position, snug the plate into the bottom left corner of the plate position at the beginning of the protocol.

## Deleting picks

It may be necessary to delete picks from the .csv file, for example if there is text on the plate which has been interpreted as a colony. The map of picks produced by the script may be used to identify the pickID of the colonies to delete. Open the .csv file in Excel and delete most of the row – leave the destination plate, row and column entries. Choose the ‘shift cells up’ option (see figure below). When you have deleted the necessary picks, clean up the end of the .csv by deleting the surplus destination plate, row and column records (figure, below). Make sure to save as a _.csv_, not an _Excel .csv_ file (which may add special characters to the file).

<img src="/readme_images/csv_excel1.PNG" width="600">

<img src="/readme_images/csv_excel2.PNG" width="600">


## Accuracy

A plate with 367 identified colonies was automatically picked and placed into media. OD measurements were taken from the four destination plates; OD measurements from the blank wells were averaged and a 3-sigma threshold value was used to determine successful colony growth. By this criterion 363 of the 367 colonies grew, indicating a success rate of 99%. An image of the plate was taken after picking, with stab marks visible where each pick occurred. A comparison of this to the plate map shows that there was some small variation in the accuracy of the pick location, with no obvious bias towards any one direction. In almost all cases, the variation was not so much that the colony was missed. (See figure below; NB: some colony growth occurred between the first and second images).

<img src="/readme_images/sidebyside.PNG" width="800">


Tests with the SELECT head indicated that the average absolute error after calibration across the full process (i.e. incorporating error from image analysis, hardware tolerance, agar shrinking, tip variation) was 0.22 mm (σ = 0.11 mm) at the (0,0) coordinate of the wells.

## Future work

The major limitation of the algorithm as it currently stands is that it does not distinguish overlapping colonies very well. Slight overlaps may be fixed by the closing function, since it is uneven and preferentially erodes the colonies. However, a circularity metric could be applied to seek colonies with circular morphologies, thus avoiding overlaps. A more involved approach would be a machine-learning algorithm which finds each colony by appearance and might be taught to identify multiple overlapping instances. A fair set of training data would of course be required to do this well.

A GUI with options to add, remove and move picks would be useful.

A bracket to hold other plate types (e.g. circular plates) could be designed to expand the range of usable plasticware.

# How to get started

## Download

Download the package from this repo; it will contain the following files:

_**colonycoords96.py**_ – this is the main script.

_**wellcentregenerator.py**_ – generates _expected_ well centres using ratios – just call from command line with your image as an argument to run. This should be run _during initial setup_, as it will generate _centrecoords.txt_, which the main script requires in order to run.

_**wellcentrecalibrator.py**_ – this is the calibration script used to produce a calibrated wellcentres.txt file – call from command line using an image which has been stabbed at (0,0) coordinates to generate a _wellcentres.txt_ file calibrated to your SELECT head. See Calibration section below.

_**testimage.tif**_ – an image to test the script with.

_**zerocoords.csv**_ – an output file with 96 ‘picks’ at (0,0) locations in each well. Used for calibration and testing.

_**colonypicker.bms**_ – this is the Felix protocol.

## Setup

1. Download files from _download_files_ in this repo.
2. Place your image (e.g. _colonies.tif_) in the same folder.
3. Install dependencies from _requirements.txt_ in the root directory of this repo.
4. Run _wellcentregenerator.py_ by calling `> python wellcentregenerator.py colonies.tif`; this will generate _wellcentres.txt_.
5. Ideally, follow calibration steps below.

## Calibration

For accuracy, it is best to calibrate the script to your particular SELECT head. The easiest way to do this is as follows:

Use zerocoords.csv to run the automated picking protocol with a blank agar plate. Drop the z-height by approx. 0.5mm if required – you are aiming to produce visible stab marks in the agar.

<img src="/readme_images/stabbed_agar.png" width="600">


Once this protocol has finished, use your imaging setup to image the plate, and then run the image through the _wellcentrecalibrator.py_ script. This will output a new _wellcentres.txt_ file (overwriting the previous), with the well centres located at the stab marks. Thus, if your SELECT head has any inherent offsets or peculiarities, they will be corrected for when producing the relative x-y pick coordinates.

Depending on the hardware you use to place colonies, and your SELECT head, you may need to adjust the x/y offsets in the variables at the top of the script. Turn 'show images' on, then run the script and observe the mask locations and exclusion zones for accuracy, changing the variables if need be.

## Optimisation variables

The code has been optimised for my imaging setup (at 600dpi) and SELECT head. A different imaging setup and new hardware may require re-optimisation of some variables within the code. These variables are located at the beginning of the script and are explained in comment.
