README

Version 1.0, 19-06-2014.

Author:  Kang Wang and Hui Ji

This package contains the MATLAB implementation of wavelet frame based image deconvolution in the presence of kernel error

==========================================
Use of this code is free for research purposes only.

================================================================================================
Reference:
H. Ji and K. Wang, Robust image deconvolution with an inaccurate blur kernel, IEEE Transactions on Image Processing, 21(4), 1624-1634, Apr. 2012 

=================================================================
This implementation has been tested with MATLAB 7.12.0 (R2011a) on a PC with 64-bit Windows 7.

=================================================================
Installation and additional toolbox requirement

1. Unpack the contents of the compressed file to a new directory.

2.  Wavelet toolbox and  image processing toolbox of MATLAB are required to run the code

========================================================
Demo

Run demo.m for a demonstration of robust image deblurring
or run demo_details.m for the demonstration with the details of the change of PSNR values for each iteration

==========================================================
Main routines 

The directory contains the routines of framelet transforms, including decomposition, reconstruction and others.

APG3_gray.m solves the L1 norm related optimization problem arising from image deconvolution using the accelerated proximal gradient  method




 

