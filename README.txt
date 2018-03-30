
        	
        Cell Separation Image Ananlysis Pipeline
        ########################################

	Copyright 2018 INRA - CNRS
	File author(s): Stéphane Verger <stephane.verger@ens-lyon.fr>
	File contributor(s): Guillaume Cerutti <guillaume.cerutti@inria.fr>

	Distributed under the Cecill-C License.
	See accompanying file LICENSE.txt or copy at
	http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html
       
        Please cite the original publication:
        Verger et al. (2018). The tension-adhesion feedback loop in plant
        epidermis.

	Links:
        DOI:...
        Github: https://github.com/sverger/Cell_separation_analysis

Description:
============                                                                
# This python script allows the semi-automatic analysis of cell separations   #
# in 2D images, and can further perform comparaisons of mean cell separation  #
# area between two conditions or genotypes, as well as the analysis of cell   #
# separation orientation and anisotropy of multiple conditions or genotypes.  #
# This is a "semi-automatic" pipeline, because the analysis requires a        #
# preliminary manual step that need to be be performed for each image before  #
# running the script: an appropriate threshold properly separating the cell   #
# signal from the cell separation background has to be defined manually.      #
# For more details see Verger et al. (2018). The tension-adhesion feedback    #
# loop in plant epidermis.                                                    #
# This script was developed on Linux (Ubuntu 14.04) in the tissuelab          #
# environment of the OpenAleaLab platform (github.com/VirtualPlants/tissuelab,#
# Cerutti G et al., (2017). FrontPlantSci 8:353. doi:10.3389/fpls.2017.00353).#
# It is designed to run with python 2.7x and has been tested on Linux         #
# (Ubuntu 14.04), Mac (OSX...?) and windows (...?), using the "recommended"   #
# install with miniconda (see "install" below).                               #

Prerequist: (e.g. using ImageJ)
============================
# For each image:
# - From a raw confocal Z-stack, make a maximal intensity Z-projection to     # 
#   obtain a 2D image.                                                        #
# - If necessary enhance the contrasts.                                       #
# - Smooth the image with a median filter to remove noise (Radius ~ 2 pixels).# 
# - With the threshold tool in ImageJ, dertermine the appropriate threshold   #
#   that separates the best the cell separations from the cells. Change the   #
#   lut to "Grays" for easier visualization.                                  #
# - Depending on the quality of the images it may be difficult or impossible  #
#   to segment the cell separtations based on a threshold. If some of your    #
#   images are in this situation, you may exclude them from your analysis. If #
#   most of your images are in this situation, this image analysis pipeline   #
#   is not appropriate for your case.                                         #
# - Then save the images in 8 bit .tif or .jpg and add "_xxxthld" at the end  #
#   of the name (e.g. "sample_1_055thld.tif). The value before thld is the    #
#   threshold value that will be used for the image segmentation (It has to   #
#   be 3 digits).                                                             #
# - Then the file arborescence has to be organised as such: A "main"          #
#   directory (updir), containing subdirectories for each condition/mutant,   #
#   each containing all the images corresponding to the given                 #
#   condition/mutant.                                                         #


Install:
========
# You can either:
# - Directly run the script if your python environment already has all the    #
#   required dependencies (list below).                                       #
#                                                                             #
# - Install the missing dependencies and run the script.                      #
#                                                                             #
# - (Recommended) Install miniconda and create a conda environment to run the # 
#   script. For this, follow the instruction in the INSTALL.txt file.         #


Dependencies:
=============
# Designed to run on linux based operating systems.                           #
# - Python 2.7x  	https://www.python.org/                               #
# - Python modules:                                                           #
# 	- matplotlib 	https://matplotlib.org/                               #
#	- nose 		http://nose.readthedocs.io/                           #
#	- numpy		www.numpy.org/                                        #  
# 	- pandas	https://pandas.pydata.org/                            #
# 	- pillow 	https://pillow.readthedocs.io/                        #
#      	- pycircstat	https://github.com/circstat/pycircstat                #
#	- scipy		https://www.scipy.org/                                #


Settings:
=========
# - Download/open the script "Cell_separation_analysis.py"                    #
#                                                                             #
# - Before running the script, define the paramenters in the section called   #
#   "Parameters".                                                             #
# - Define a directory containing all the data to analyse and compare (updir).#
# - Define the pixel size of your image in micrometer (pixel_size).           #
# - Define min and max area of crack to eliminate areas that are too small    #
# (min), and/or the background (max).                                         #
# - Define the threshold type "min" or "max". "max" will detect and segment   #
# zones with lowest intensity of signal (black gaps between separated cells). #
# - Output global cracks size analysis? True or False. This is if you want to #
# run the analysis of gap size comparing two mutants or conditions in the     #
# updir folder.                                                               #
# - Output global crack orientation analysis? True or False. This is if you   #
# want to run the analysis of gap orientation.                                #
# - Save the file to save the new parameters.                                 # 


running the script:
===================
# In a terminal run $ python or $ ipython.                                    #
# In the python console, type $ %run Cell_separation_analysis.py              #
#                                                                             #
# This should run the script, write some output in the python console and     #
# create the output images and files.                                         #


Output:
=======
# - For each image, a .csv file is created containing for each segmented cell #
#   separation, the label number, center position, area in pixels and um      #
#   square, the principal angle (orientation of the separation), the          #
#   anisotropy, the eigen values and vectors and the standard deviations.     #
# - Also for each image, an inverted version of the image, overlaid with the  #
#   segmented areas as well as a representation of the anisotropy and         #
#   principal angle for each area, is saved as a vectorial .pdf.              #
# - Finally, for each image, a polar histogram representing the distribution  #
#   of cell separation orientation in the image is created and saved as a     #
#   vectorial .pdf.                                                           #
# - If running the "Global_Output_Size", statistical tests will be run on the #
#   compared samples, and the summary of these tests will be saved in a .txt  #
#   file in the "updir", as well as a box plot as a vectorial .pdf.           #
# - If running the "Global_Polarhist_Output", a global polar histogarm will   #
#   be created for each condition/mutant and saved as a vectorial .pdf and    #
#   the circular mean angle, the resultant vector length and the mean         #
#   anisotropy will be computed and saved in a summary .txt file.             #





