%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            PSM README
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

---------------------
      Contents
---------------------
0) Introduction
1) Compatibility
2) Version History
3) Dependencies
4) PSM Process



-----------------
0) Introduction
-----------------
This project generates Phase Slope Magnitude (PSM) images from raw MRI datasets. PSM images are positive-contrast images that emphasize spatial magnetic field variations to help distinguish SPIO-induced image hypointensity from intrinsic contrast.

The first publication that describe PSM's function and exhibits PSM results:

Enhanced positive-contrast visualization of paramagnetic contrast agents using phase images. (http://tinyurl.com/9r7qfba OR http://tinyurl.com/8en459a)

PSM was developed in Ahrens Lab, Pittsburgh NMR Center for Biomedical Research, Carnegie Mellon University. Financial support was provided in part by the U.S. National Institutes of Health.




-----------------
1) Compatibility
-----------------
PSM runs in MATLAB and has been tested in:
	-Windows 7/Vista 64-bit (no phase-unwrapping, see note below)
	-Mac OS (10.4)
	-Linux (Ubuntu 10.02)

Note: PSM cannot unwrap MRI phase data on Windows platforms. The phase unwrapping software used,
PRELUDE, is a binary executable that runs only in Mac OS X and Linux.




-------------------
2) Version History
-------------------
PSM-8-27-2012.zip
Second posting of code, some dependencies included, sparse documentation.




------------------
3) Dependencies
------------------
PSM is dependent on the following:

-MATLAB 2010a or newer
-neighbor_differences.m (MATLAB file, included)
-prelude_linux (Binary executable - Linux, included)
-prelude_macosx (Binary executable - MacOS, included)




---------------------
4) PSM Process
---------------------
The process of creating Phase Slope Magnitude images can be explained in 4 steps below.
A tutorial that follows this process can be found in TUTORIAL.TXT.

1) Import your MRI phase image data into MATLAB

2) Phase unwrap your phase image using PRELUDE

3) (Optional) High-pass filter your phase-unwrapped data

4) Run PSM
