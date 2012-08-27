%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            PSM README
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


---------------------
      Contents
---------------------
0) Dependencies
1) PSM Process
2) PSM Data Structure



-----------------
0) Dependencies
-----------------
PSM is dependent on the following MATLAB functions created by 3rd parties:

XXXXXXXXXXXX


-------------------------------------
1) PSM Process
-------------------------------------

A) Import data
MRIdata = import_fid(filename, k_space_type, ordering, RARE_factor, read_dim, reco_dim, noise, NEX_scheme_y,NEX_scheme_z, bigEndian);

B) Process data
MRIdata = XXXXXXXXXXXXXXXXX




-------------------------------------
2) PSM Data Structure
-------------------------------------
          .PSM                 Results from a run of PSM - Phase Slope Magnitude
              .PSM.mag                 (3D Float)    Phase slope magnitude image   (How rapidly phase is changing in space)
              .PSM.x                   (3D Float)    X-direction Phase slope magnitude image
              .PSM.y                   (3D Float)    Y-direction Phase slope magnitude image
              .PSM.z                   (3D Float)    Z-direction Phase slope magnitude image
              .PSA                     (3D Float)    Phase slope angle image       (The angle (in radians) at which phase is changing)
              .AC                      (3D Float)    Phase angle coherence image   (Measure of how regionally consistent PSA image is)
              .Composite               (3D Float)    Composite image               (PSM * AC. Highlights commonly-oriented rapid phase changes!)



