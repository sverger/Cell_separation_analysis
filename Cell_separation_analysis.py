# -*- coding: utf-8 -*-
# -*- python -*-
#
#       Cell Separation Image Ananlysis Pipeline
#
#       Copyright 2018 INRA - CNRS
#
#       File author(s): Stéphane Verger <stephane.verger@ens-lyon.fr>
#
#       File contributor(s): Guillaume Cerutti <guillaume.cerutti@inria.fr>
#
#       Distributed under the Cecill-C License.
#       See accompanying file LICENSE.txt or copy at
#           http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html
#       
#       Please cite the original publication:
#       Verger et al. (2018). The tension-adhesion feedback loop in plant
#       epidermis.
#       
# 	Links:
#	DOI:...
#       Github...

###############################################################################
# Description:                                                                #
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
# This script was developed and run in the tissuelab environment of the       #
# OpenAleaLab platform (github.com/VirtualPlants/tissuelab, Cerutti G et al., #
# (2017). Front Plant Sci 8:353. doi:10.3389/fpls.2017.00353).                #
# It is designed to run with python 2.7x and has been tested on Linux (14.04) # 
# and Mac (OSX...?)                                                           #
#                                                                             #
# For more details on how to you this script, see README.txt                  #
#	github...                                                             #
###############################################################################

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.patheffects as patheffects
import numpy as np
import os
import pandas as pd
import pycircstat
from scipy import stats
from scipy.ndimage.io import imread
import scipy.ndimage as nd

# Parameters==================================================================#
###############################################################################
# Before running the script, define all the parameters.                       #
# - Define a directory containing all the data to analyse and compare (updir).#
# - Define the pixel size in micrometer (pixel_size).                         #
# - Define min and max area of crack to eliminate areas that are too small    #
# (min), and/or the background (max).                                         #
# - Define the threshold type "min" or "max". "max" will detect and segment   #
# zones with lowest intensity of signal (black gaps between separated cells). #
# - Output global cracks size analysis? True or False. This is if you want to #
# run the analysis of gap size comparing two mutants or conditions in the     #
# updir folder.                                                               #
# - Output global crack orientation analysis? True or False. This is if you   #
# want to run the analysis of gap orientation.                                #
###############################################################################

updir = './Test_files/minthreshold/' 

pixel_size = 0.363636  # 0.363636um square is the pixel size
min_area_of_crack = 100  # Size in pixels
max_area_of_crack = 100000
thld_type = 'max'
Global_Output_Size = True
Global_Polarhist_Output = True
# Parameters==================================================================#


# Functions===================================================================#
def Segment_cell_separation(img, threshold=150, threshold_type=thld_type, 
                            min_area=min_area_of_crack, 
                            max_area=max_area_of_crack):
 
    if threshold_type == 'max':
        img_regions = nd.label(img < int(threshold))
    elif threshold_type == 'min':
        img_regions = nd.label(img > int(threshold))
        
    labelled_regions = img_regions[0]
    labels = np.arange(img_regions[1])+1
    
    region_areas = dict(zip(labels, 
                            nd.sum(np.ones_like(labelled_regions), 
                                   labelled_regions, 
                                   index=labels)))
    if max_area is not None:
        regions_to_remove = np.array([r for r in region_areas.keys() if region_areas[r] > max_area])
        for l in regions_to_remove:
            labelled_regions[labelled_regions == l] = 0
    
    regions_to_remove = np.array([r for r in region_areas.keys() if region_areas[r] < min_area])
    for l in regions_to_remove:
        labelled_regions[labelled_regions == l] = 0
        
    return labelled_regions


def labelled_region_dataframe(labelled_regions, px_size=pixel_size):
    
    region_labels = np.unique(labelled_regions)[1:]
    region_data = dict()
    for field in ['label', 'center', 'center_x', 'center_y', 'area_px', 
                  'area_um_sq', 'principal_angle', 'anisotropy', 
                  'eigen_values', 'eigen_vectors', 'standard_deviations', 
                  'standard_deviations0', 'standard_deviations1']:
        
        region_data[field] = []

    for l in region_labels:
        coords = np.transpose(np.where(labelled_regions == l))
        region_center = np.mean(coords, axis=0)
        region_covariance = np.cov(coords, rowvar=False)
        eigen_values, eigen_vectors = np.linalg.eig(region_covariance)
        
        eigen_vectors = eigen_vectors[np.argsort(-np.abs(eigen_values))]
        eigen_values = eigen_values[np.argsort(-np.abs(eigen_values))]
        
        region_data['label'] += [l]
        region_data['center'] += [region_center]
        region_data['center_x'] += [region_center[0]]
        region_data['center_y'] += [region_center[1]]
        region_data['area_px'] += [len(coords)]
        region_data['area_um_sq'] += [len(coords)*px_size*px_size]
        region_data['eigen_values'] += [eigen_values]
        region_data['standard_deviations'] += [np.sqrt(np.abs(eigen_values))]
        region_data['standard_deviations0'] += \
            [np.sqrt(np.abs(eigen_values[0]))]
        region_data['standard_deviations1'] += \
            [np.sqrt(np.abs(eigen_values[1]))]
        region_data['eigen_vectors'] += [eigen_vectors]
        region_data['anisotropy'] += [np.abs(eigen_values[0]/eigen_values[1])]
        region_data['principal_angle'] += \
            [((np.sign(eigen_vectors[0][0]) *
               np.arccos(eigen_vectors[0][1])*180./np.pi)) % 180] 
    
    return pd.DataFrame.from_dict(region_data)

""" Script:
From the designated directory (updir),
finds subdirectories containing the experiments/conditions to analyse.
"""
for dirname in sorted(os.listdir(updir)):
  drpath = updir + dirname
  if os.path.isdir(drpath):
    # print drpath
    print "Folder analyzed: " + dirname
    
    """ For each .tif image in each directory: """
    for img_file in sorted(os.listdir(drpath)):
      if ("thld.jpg" in img_file) or ("thld.tif" in img_file):
        imgpath = drpath + "/" + img_file
        # print 'image path: ' + imgpath
        img = imread(imgpath).astype(float)
	if img.shape[-1]==3:
            img = img.mean(axis=-1)
        print "image analyzed: " + img_file
        # Extracts the threshold value from the title
        thld1 = img_file[:-8] 
        thld = thld1[-3:]
        # Trims the 0 at the begining of the value if it is lower than 100  
        if thld[:1]is '0':
            thld = thld[1:]

        """ Runs the segmentation."""
        labelled_regions = Segment_cell_separation(img, threshold=thld)
        print "segmentation done"

        """ Labels of segmented regions, shuffled for visualization."""
        labels = np.unique(labelled_regions)[1:]
        np.random.shuffle(labels)
        shuffled_labels = dict(zip(np.unique(labelled_regions)[1:],labels))
        shuffled_labels[0] = 0
        
	shuffled_label_regions = np.zeros_like(labelled_regions)
	for l in shuffled_labels.keys():
            shuffled_label_regions[labelled_regions==l] = shuffled_labels[l]
        labelled_regions = shuffled_label_regions
  
        """ Makes and saves an inverted (pixel value) .tif image."""
        figure = plt.figure(img_file + "_inv")
        figure.clf()
        # figure.patch.set_facecolor('w')
        figure.gca().set_xlim(0, img.shape[1])
        figure.gca().set_ylim(img.shape[0], 0)
        figure.set_size_inches(5, 5)
        figure.gca().imshow(img, cmap='gray_r', vmin=0, vmax=255)
        figure.savefig(imgpath[:-4]+"_inv.pdf", dpi=300)
  
        """ Makes and saves an overlay of the inverted image with the detected 
        labelled regions."""
        figure.gca().imshow(np.ma.masked_where(labelled_regions == 0, labelled_regions), 
                            cmap='winter', alpha=1, interpolation='none')
        figure.savefig(imgpath[:-4]+"_inv_sep_labels.pdf", dpi=300)
                
        """ Makes the data frame."""
        data = labelled_region_dataframe(labelled_regions)
        
        """ Generates and plots the separation orientation axes."""
        tensor_factor = 2
        for center, vals, vecs, angle, area in zip(data['center'],
                                                   data['standard_deviations'],
                                                   data['eigen_vectors'],
                                                   data['principal_angle'],
                                                   data['area_px']):
            
            figure.gca().scatter(center[1], center[0], color='k', s=0)
            for val, vec in zip(vals, vecs):
                figure.gca().plot([center[1]+tensor_factor*val*vec[1], 
                                   center[1]-tensor_factor*val*vec[1]], 
                                  [center[0]-tensor_factor*val*vec[0], 
                                   center[0]+tensor_factor*val*vec[0]], 
                                  color='r', linewidth=2, 
                                  path_effects=[patheffects.withStroke(linewidth=4, foreground='w')])
                
        """ Saves an overlay of the inverted image with the detected labelled
        regions and the cell separation orientations."""
        figure.savefig(imgpath[:-4]+"_inv_sep_orientations.pdf", dpi=300)

        """ Makes and saves a polar histogram of cracks orientations."""
        figure = plt.figure(img_file + "_hist")
        figure.clf()
        figure.patch.set_facecolor('w')
        figure.set_size_inches(5, 5)
        ax = figure.add_subplot(111, polar=True)
        
        colormap = 'plasma'
        n_bins = 36
        
        # weights = data['anisotropy']*data['area_px']
        # weights = data['anisotropy']
        # weights = data['area_px']
        # weights = np.ones_like(data['area_px'])
        # weights = np.array([np.max(np.abs(v)) for v in data['eigen_values']])
        weights = np.array([np.max(s) for s in data['standard_deviations']])
        
        for offset in [0, np.pi]:
            histo, bins, patches = figure.gca().hist(offset+data['principal_angle']*np.pi/180., 
                                                     bins=offset+np.linspace(0, 180, n_bins/2 + 1)*np.pi/180., 
                                                     ec='k', weights=weights)
            norm = colors.Normalize(0, histo.max())
            for h, p in zip(histo, patches):
                p.set_facecolor(cm.get_cmap(colormap)(norm(h)))
        figure.gca().set_yticks([])
        figure.savefig(imgpath[:-4]+"_sep_hist.pdf", dpi=300)
        
        # world.add(data,'crack_data')
        data.to_csv(imgpath[:-4]+"_sep_log.csv")
print "Done analysing each image"


"""Global_Output_Size:
Makes and saves a boxplot to compare the sum of crack area per image in the 
different conditions
"""
if Global_Output_Size is True:
    print "Global cell separation size analysis"
    fdatasize = open(updir + '/Global cell separation size analysis_'
                     'Summary.txt', 'w')  
    sums = []
    plotlabels = []
    sum_ = {}
    normpop = []
    for dirname in sorted(os.listdir(updir)):
        drpath = updir + dirname
        if os.path.isdir(drpath):
            print dirname
            sum_[dirname] = []
            # print drpath
      
            for filename in os.listdir(drpath):
                if '_sep_log.csv' in filename:
                    file_data = pd.read_csv(drpath + '/' + filename)
                    # print (drpath+"/"+filename)
                    sum_area_per_img = np.sum(file_data['area_um_sq'])
                    sum_[dirname].append(sum_area_per_img)
            print sum_[dirname]
            mean_area_ALL_img = np.mean(sum_[dirname])
            std_area_ALL_img = np.std(sum_[dirname])
            print mean_area_ALL_img
            print std_area_ALL_img
            sums.append(sum_[dirname])
            print sum_[dirname]
            fdatasize.write("==> " + dirname + "\n" +
                            "mean area of separation (+-Std) --> " + 
                            str(mean_area_ALL_img) + "+-" + 
                            str(std_area_ALL_img) + "\n")
   
            """ Shapiro's test for population normality."""
            w, p = stats.shapiro(sum_[dirname])
            print p
            if p > 0.05:
                normpop.append(True)
                print "population is normal"
                fdatasize.write("--> population is normal --> "
                                "p-value (Shapiro's test) :" + str(p) + "\n\n")
            else:
                normpop.append(False)
                print "population is NOT normal"
                fdatasize.write("--> population is NOT normal --> "
                                "p-value (Shapiro's test) :" + str(p) + "\n\n")
            plotlabels.append(dirname)
    print normpop
    print sums
    print plotlabels
    if False in normpop:
        """ Non parametric Wilcoxon rank sum test."""
        print "At least one sample does Not have a normal distibution" \
            "--> wilcoxon rank sum test"
        fdatasize.write("At least one sample does Not have a normal "
                        "distibution --> wilcoxon rank sum test" + "\n")
        statrank, prank = stats.ranksums(*sums)
        if prank > 0.05:
            print "--> populations are NOT statiscically different " \
                "--> p-value is " + str(prank)
            fdatasize.write("--> populations are NOT statiscically different "
                            "--> p-value (wilcoxon rank sum test) : " + 
                            str(prank) + "\n\n")
        else:
            print "--> populations are statiscically different " \
                "--> p-value is " + str(prank)
            fdatasize.write("--> populations are statiscically different "
                            "--> p-value (wilcoxon rank sum test) : " + 
                            str(prank) + "\n\n")
    
    else:     
        """ Bartlett's test for equal variance."""
        t, p2 = stats.bartlett(*sums)
        print p2
        if p2 > 0.05:
            equalvar = True
            print "Groups have equal variances"
            fdatasize.write("--> The groups have equal variances --> "
                            "p-value (Bartlett's test) : " + str(p2) + "\n\n")
        else:
            equalvar = False
            print "Groups DO NOT have equal variances"
            fdatasize.write("--> The groups Do NOT have equal variances --> "
                            "p-value (Bartlett's test) : " + str(p2) + "\n\n")
  
        """t-test (Welch's if variances are unequal)."""
        t2, prob = stats.ttest_ind(*sums, equal_var=equalvar)
        print prob
        if prob > 0.05:
            print "populations are NOT statiscically different\np-value is " \
                + str(prob)
            fdatasize.write("--> populations are NOT statiscically different "
                            "--> p-value (t-test, Welch's if variances are "
                            "unequal) : " + str(p2) + "\n\n")
        else:
            print "populations are statiscically different\np-value is " + \
                str(prob)
            fdatasize.write("--> populations are statiscically different "
                            "--> p-value (t-test, Welch's if variances are "
                            "unequal) : " + str(p2) + "\n\n")
    fdatasize.close()
    
    """Boxplot."""
    figure = plt.figure("Global_crack_Size")
    figure.clf()
    medianprops = {'color': 'red', 'linewidth': 2}
    boxprops = {'color': 'black', 'linewidth': 3}
    whiskerprops = {'color': 'black', 'linewidth': 3}
    capprops = {'color': 'black', 'linewidth': 3}
    flierprops = {'color': 'black', 'marker': 'x'}

    ax = plt.boxplot(sums, labels=plotlabels, widths=0.7, 
                     medianprops=medianprops, boxprops=boxprops, 
                     whiskerprops=whiskerprops, capprops=capprops, 
                     flierprops=flierprops, patch_artist=True)
    plt.tight_layout()
    plt.show()
    figure.savefig(updir + "Boxplot.jpg")
    figure.savefig(updir + "Boxplot.pdf")


"""Global_Polarhist_Output
Makes and saves a global polar histogarm for each experiment and gives the 
circ_mean angle, the resultant vector length and the mean anisotropy.
"""
if Global_Polarhist_Output is True:
    print "Global orientation analysis"
    Crack_data_ = {}
    for dirname in sorted(os.listdir(updir)):
        drpath = updir + dirname
        if os.path.isdir(drpath):
            print dirname
            Crack_data_[dirname] = []
            # print drpath

            fdata = open(drpath + "/" + dirname + "_Summary.txt", 'w')
            fdata.write("Directory used: \n" + updir + dirname + 
                        "\nFiles analysed: \n")
            for filename in os.listdir(drpath):
                if "_sep_log.csv" in filename:
                    fdata.write(filename + "\n")
                    file_data = pd.read_csv(drpath + "/" + filename)
                    principal_angle = file_data['principal_angle']
                    Rad_angle = np.deg2rad(principal_angle)
                    Rad_angle2 = (np.deg2rad(principal_angle))*2
                    circ_mean = np.rad2deg(stats.circmean(Rad_angle, 
                                                          high=np.pi, low=0))
                    circ_std = np.rad2deg(stats.circstd(Rad_angle, 
                                                        high=np.pi, low=0))
                    RVL = pycircstat.resultant_vector_length(Rad_angle2)
                    ani_mean = np.mean(file_data['anisotropy'])
                    ani_std = np.std(file_data['anisotropy'])
                    fdata.write("Circ_Mean Angle to x (+-Std) --> " + 
                                str(circ_mean) +  "+-" + str(circ_std) + 
                                " degrees\nResultant_vector_length " + 
                                str(RVL) + "\nMean Anisotropy (+-Std) --> " + 
                                str(ani_mean) + "+-" + str(ani_std)+ "\n\n")
                    Crack_data_[dirname].append(file_data)

            data = pd.concat(Crack_data_[dirname], ignore_index=True)
            principal_angles = data['principal_angle']
            Rad_angles = np.deg2rad(principal_angles)
            Rad_angles2 = (np.deg2rad(principal_angles))*2
            circ_mean_all = np.rad2deg(stats.circmean(Rad_angles, 
                                                      high=np.pi, low=0))
            circ_std_all = np.rad2deg(stats.circstd(Rad_angles, 
                                                    high=np.pi, low=0))
            RVL_all = pycircstat.resultant_vector_length(Rad_angles2)
            ani_mean_all = np.mean(data['anisotropy'])
            ani_std_all = np.std(data['anisotropy'])
            pval, U, Uc = pycircstat.raospacing(Rad_angles2)
            print pval
            print U
            print Uc
            if pval > 0.05:
                print "The populations IS uniformly distributed\np-value is " \
                    + str(pval)
                fdata.write("Global analysis:\nCirc_Mean Angle to x (+-Std) "
                            "--> " + str(circ_mean_all) + "+-" + 
                            str(circ_std_all) + 
                            " degrees\nResultant_vector_length " + 
                            str(RVL_all) + 
                            "\nMean Anisotropy (+-Std) --> " + 
                            str(ani_mean_all) + "+-" + str(ani_std_all) + 
                            "\nThe populations IS uniformly distributed\n"
                            "p-value is " + str(pval) + "\n")
            else:
                print "The populations is NOT uniformly distributed\np-value" \
                    "is " + str(pval)
                fdata.write("Global analysis:\nCirc_Mean Angle to x (+-Std) "
                            "--> " + str(circ_mean_all) + "+-" + 
                            str(circ_std_all) + 
                            " degrees\nResultant_vector_length " + 
                            str(RVL_all) + 
                            "\nMean Anisotropy (+-Std) --> " + 
                            str(ani_mean_all) + "+-" + str(ani_std_all) + 
                            "\nThe populations IS NOT uniformly distributed\n"
                            "p-value is " + str(pval) + "\n")          
            fdata.close()
      
            figure = plt.figure(0)
            figure.clf()
            figure.patch.set_facecolor('w')
            ax = figure.add_subplot(111, polar=True)
      
            colormap = 'plasma'
            n_bins = 36
      
            # weights = data['anisotropy']*data['area']
            weights = data['anisotropy']
            # weights = data['area']
            # weights = np.ones_like(data['area'])
            # weights = np.array([np.max(np.abs(v)) for v in data['eigen_values']])
            # weights = np.array([np.max(s) for s in np.array(data['standard_deviations'])])
            # weights = data['standard_deviations0']
      
            for offset in [0, np.pi]:
                histo, bins, patches = figure.gca().hist(offset+data['principal_angle']*np.pi/180., 
                                                         bins=offset+np.linspace(0, 180, n_bins/2 + 1)*np.pi/180., 
                                                         ec='k', weights=weights)
                norm = colors.Normalize(0, histo.max())
                for h, p in zip(histo, patches):
                    p.set_facecolor(cm.get_cmap(colormap)(norm(h)))
            figure.gca().set_yticks([])
            figure.savefig(drpath + "/" + dirname + "_cracks_polarhist.jpg")
            figure.savefig(drpath + "/" + dirname + "_cracks_polarhist.pdf")
            data.to_csv(drpath + "/" + dirname + "_data_all.csv")
            print "image saved"
