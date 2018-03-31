# RapidCFD
### CFD toolbox running on CUDA

OpenFOAM solvers ported to Nvidia CUDA. All the calculations are done on the GPU giving a huge speed-up.
Still in development stage, waiting for your contribution!

### Features:
* most incompressible and compressible solvers on static mesh are available
* all the calculations are done on the GPU
* no overhead for GPU-CPU memory copy
* can run in parallel on multiple GPUs

### Compilation Instructions for Ubuntu 16.04:
* ensure CUDA 7.5 is not installed from Ubuntu repsitories 
* ensure you are using a nVidia driver compatible with CUDA 8
  e.g., for the nVidia Tesla K20: http://www.nvidia.com/download/driverResults.aspx/118962/en-us
* download CUDA 8.0 here: https://developer.nvidia.com/cuda-80-ga2-download-archive
* install CUDA 8.0 per instructions provided at time of download
* to compile in parallel, export WM_NCOMPPROCS=10 (as appropriate for your computer)
* software should compile without errors. Result: RapidCFD for Ubuntu 16.04
* ThirdParty-dev is needed for multiple GPU's.
