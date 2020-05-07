# open-slam-c
We have converted the OpenSLAM code (https://github.com/OpenSLAM-org/openslam_bailey-slam) to C and parallelized using OpenMP.

# Running slam
To run either ekfSLAM or fastSLAM, navigate to the corresponding directory and run either `slam-omp` or `slam-seq`. A help menu is available for both programs. Test data can be found in the test folder. A landmark file and waypoint file are needed for the program to run. 
