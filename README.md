# Adaptive Curvelet Thresholding for Filtering of Additive Stationary White/Colored Gaussian Noise

This package contains the [`denoising filter`](./ACT_filter.m) for additive stationary
white/colored(a.k.a. correlated) Gaussian noise based on the adaptive
curvelet thresholding (ACT) method presented in: 

  N. Eslahi and A. Aghagolzadeh, "Compressive Sensing Image Restoration
  Using Adaptive Curvelet Thresholding and Nonlocal sparse Regularization",
  IEEE Trans. Image Process., vol.25, no.7, pp. 3126-3140, Jul. 2016
  https://doi.org/10.1109/TIP.2016.2562563




Author:                Nasser Eslahi
                [(Tampere University, Finland)](https://www.tuni.fi/en)



## Contents
### The ACT package includes:
[`ACT_filter.m`](./ACT_filter.m)  the denoising filter based on ACT
\
[`demo_denoising.m`](./demo_denoising.m) a demo code for denoising of additive stationary Gaussian noise
\
[`./DCuT/fdct_wrapping.m`](./DCuT/fdct_wrapping.m) forward discrete curvelet transform
\
[`./DCuT/ifdct_wrapping.m`](./DCuT/ifdct_wrapping.m) inverse/backward discrete curvelet transform
\
[`./Test_Images`](./Test_Images) some images for testing the denoising performance 
\
\
For a detailed information on each function, please check its corresponding syntax, e.g.,
```
  doc ACT_filter
```



## Visualization of some random tests
https://user-images.githubusercontent.com/48449082/170306455-d4a87a82-0011-4127-a266-cc35021e2c76.mp4









## Disclaimer
Copyright (c) 2006-2022    All rights reserved.
This work should be used for nonprofit purposes and non-commercial use only, see [`LICENSE`](./LICENSE).




## Feedback
If you have any comment, suggestion, or question, please do contact
 [Nasser Eslahi](https://orcid.org/0000-0002-1134-9318)
\
 nasser.eslahi@tuni.fi
\
nasser.eslahi@gmail.com