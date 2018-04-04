Installation procedure for the Cell Separation Image Ananlysis Pipeline
#######################################################################

Following this guide will help you install all the dependencies required for the script and ensure that the script will run properly.  


Download and install Miniconda
==============================

- Download the miniconda installer from the official website repo.continuum.io
	LINUX: https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
	MAC: https://repo.continuum.io/miniconda/Miniconda2-latest-MacOSX-x86_64.sh
	Windows: https://repo.continuum.io/miniconda/Miniconda2-latest-Windows-x86_64.exe

	You can also use wget to perform this download from a shell (LINUX or MAC):
	LINUX: $ wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh

- Install miniconda by running the installer:
	LINUX: Open a new terminal window, go to the directory where you downloaded the installer and run:
		$ bash Miniconda2-latest-Linux-x86_64.sh
		$ rm Miniconda2-latest-Linux-x86_64.sh
	MAC: Open a new terminal window, go to the directory where you downloaded the installer and run:
		$ bash Miniconda2-latest-MacOSX-x86_64.sh
		$ rm Miniconda2-latest-MacOSX-x86_64.sh
	Windows: Execute the installer and follow the instructions
		During the installation you will be asked a number of choices
		You can set up the directory of your choice when asked, e.g. ~/.miniconda
		Make sure to answer YES when asked to add conda to your PATH

- You should now have miniconda properly installed; test your installation by running conda in a terminal to make sure the command is found

Conda tips:
Get out of the current conda environment : $ source deactivate
View the available environments : $ conda env list


Create and activate a conda environment to run the script
=========================================================

- Create or download (in your "home" for Linux and Mac) the file called cell-sep-env.yml containing the list of dependencies

	name: cell-sep-env
	channels:
	  - defaults
	dependencies:
	  - python=2.7
	  - ipython
	  - numpy
	  - scipy
	  - matplotlib
	  - pandas
	  - nose
	  - pillow
	  - pip:
	      - pycircstat

In a terminal:
- Define a new conda environment
	$ conda env create -f cell-sep-env.yml

- Activate the environment
	$ source activate cell-sep-env


Check your installation
=======================
You can run the test samples in python in the newly installed environment

- Download the test folder...

- In a terminal run $ ipython

This will launch ipython.

- In the python console, type $ %run Cell_separation_analysis.py
 
This should run the script, output and create images.
