# RapidCFD
### CFD toolbox running on CUDA

OpenFOAM solvers ported to Nvidia CUDA. All the calculations are done on the GPU giving a huge speed-up.
Still in development stage, waiting for your contribution!

### Features:
* most incompressible and compressible solvers on static mesh are available
* all the calculations are done on the GPU
* no overhead for GPU-CPU memory copy
* can run in parallel on multiple GPUs

### This Fork:
* updated selected files so RapidCFD will compile under Ubuntu 16.04
* fixed interfaceProperties so contactAngle boundary conditions do not cuase thrust::system errors
