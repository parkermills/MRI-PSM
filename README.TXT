%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            PSM README
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


---------------------
      Contents
--------------------- 
1) Compatibility
2) Version History
3) Dependencies
4) PSM Process
5) PSM Data Structure



-----------------
1) Compatibility
-----------------
PSM runs in MATLAB and has been tested in:
	-Windows Vista 64-bit (no pre-processing, see note below)
	-Windows 7 64-bit (no pre-processing, see note below)
	-Mac OS X
	-Ubuntu 10.02

Note: PSM cannot complete MRI data pre-processing on Windows platforms. Its phase unwrapping pre-processing step depends on a binary executable that runs only in Mac OS X and Linux. 




-----------------
2) Version History
-----------------
PDQ-8-27-12
Original posting of code, no dependencies included, sparse documentation.




-----------------
3) Dependencies
-----------------
PSM is dependent on the following MATLAB functions created by 3rd parties:

XXXXXXXXXXXX




-------------------------------------
4) PSM Process
-------------------------------------

A) Import data
MRIdata = import_fid(filename, k_space_type, ordering, RARE_factor, read_dim, reco_dim, noise, NEX_scheme_y,NEX_scheme_z, bigEndian);

B) Process data
MRIdata = XXXXXXXXXXXXXXXXX




-------------------------------------
5) PSM Data Structure
-------------------------------------
          .PSM                 Results from a run of PSM - Phase Slope Magnitude
              .PSM.mag                 (3D Float)    Phase slope magnitude image   (How rapidly phase is changing in space)
              .PSM.x                   (3D Float)    X-direction Phase slope magnitude image
              .PSM.y                   (3D Float)    Y-direction Phase slope magnitude image
              .PSM.z                   (3D Float)    Z-direction Phase slope magnitude image
              .PSA                     (3D Float)    Phase slope angle image       (The angle (in radians) at which phase is changing)
              .AC                      (3D Float)    Phase angle coherence image   (Measure of how regionally consistent PSA image is)
              .Composite               (3D Float)    Composite image               (PSM * AC. Highlights commonly-oriented rapid phase changes!)



