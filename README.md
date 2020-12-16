# QuaHOG
Quadrature for High-Order Geometries

QuaHOG is a repository of MATLAB functions for producing quadrature rules over high-order (polynomial or rational parametric) trimmed/intersected geometries in 2D. We plan to add 3D support as well.

This code accompanies the paper
<a id="1">[1]</a> 
Gunderman, D., Weiss, K., Evans, J.A.. 
Spectral mesh-free quadrature for planar regions bounded by rational parametric curves. 
Computer Aided Design, 130(102944)  (2021).
https://www.sciencedirect.com/science/article/abs/pii/S0010448520301378?via%3Dihub

The methods listed in the paper are in the SPECTRAL_quads.m (for spectrally convergent quadrature over 2D regions bounded by rational curves) and SPECTRALPE_quads.m (same behavior as SPECTRAL_quads, but exact for polynomials up to a pre-specified degree)

# License
Copyright (c) 2020, Lawrence Livermore National Security, LLC. Produced at the Lawrence Livermore National Laboratory.

Copyrights and patents in the QuaHOG project are retained by contributors. No copyright assignment is required to contribute to QuaHOG.

See LICENSE for details.

Unlimited Open Source - BSD 3-clause Distribution LLNL-CODE- 816589
