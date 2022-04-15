# YADPF-II
Yet Another Dynamic Programming Function II (C++ Version)

YADPF (MATLAB version) relies on vectorization to speed up the computation. Consequently, YAPDF memory requirements grow exponentially as the number of states and input variables grows. For this reason, we develop a C++ version of YAPDF, and we name it YADPF-II. 

In YAPDF-II, we avoid flattening the state and input variables resulting in a heavily nested for-loop structure in its source code. To increase the computational speed, we utilize OpenMP to parallelize some of these for-loops. The result is quite promising since it is not far different than its vectorized version.

Currently, there is only Windows implementation (Visual Studio 2019). Please extract the following three ZIP-files resulting in three new folders with similar names.  

- armadillo-10.8.2  
- boost-1.78.0.zip  
- lapack.zip  

These three files can be downloaded from [this dropbox link.](https://www.dropbox.com/sh/2hwz8nuxwkazo3y/AACsZaaEV4bneR7n74EC5-rda?dl=0) 

Extract them in the root directory of YADPF-II. Some of those three libraries are precompiled using Visual Studio 2018, while others are headers only.  

There are three demo file provided. Detailed documentations are coming soon.

manurunga@yandex.com
