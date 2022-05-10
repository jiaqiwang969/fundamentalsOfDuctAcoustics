# fundamentalsOfDuctAcoustics
- Aim: Reproduce figs in Rienstra's paper "fundementals of duct acoustics"



### Results List
|  figure  | info  |  link |   Derivatives |
| -- | --| --| ---|
| fig05| Complex axial wave number (m=0, w=5) |  Fig5.m | Fig5_deriv01|



## IDEA LIST

- Fig5_deriv01.m: reconstruct single duct mode

<img src="https://github.com/jiaqiwang969/fundamentalsOfDuctAcoustics/blob/main/results/01-modeEx.png" width="640px">


- Fig5_deriv02.m: reconstruct several modes
 <img src="https://cdn.mathpix.com/snip/images/KbsQ8HWF-WHqnBqArVxMUZRoDyS9RvL_w40SMxk7ZtE.original.fullsize.png" width="640px">


- Fig5_deriv03.m: Fig 07 prove?
 <img src="https://cdn.mathpix.com/snip/images/B4uI1LCB_8HkSKWrftj3-hNinE_G6XitbQcaCOqSiS8.original.fullsize.png" width="640px">
 
 - Fig5_deriv04.m: Does green function is possible to plot? with solid wall.
 
   Not really sure, code is here, to be vertifying...
   
 - Fig5_deriv05.m: fig8 prove?

<img src="https://cdn.mathpix.com/snip/images/upy0Qko5EfVRDb8Ab7GavynHwZMhS5T5fywrybI6tvY.original.fullsize.png" width="640px">

 - Fig5_deriv06.m: fig8 extends to 3D, which includes "n" in z-axis, all this mode add together is real world!
 
<img src="https://cdn.mathpix.com/snip/images/E4TzasLgU3K-1MpaUUcGYMYcJ3CwmS8fIwyKE9l9I6k.original.fullsize.png" width="640px">

- Fig5_deriv07.m: Assume the amplitudes of this modes is random at first, so what? Please think about it, when w=20, which level of "m" at least is cut off? 

  Answer: This derives the cut-off boundary of the duct modes related to "m". So can we plot it?
<!-- <img src="https://cdn.mathpix.com/snip/images/F7-GGLUARejYqMu7StndYSgUOavKQp2DAVR90YNnwQw.original.fullsize.png" width="1040px" height="340px"> -->
<img src="https://cdn.mathpix.com/snip/images/LbCY0kiyZ83jh39DAphqoGRzr7_9EV7ehM3us0ATZvQ.original.fullsize.png"  width="640px">


<!--  x<img src="https://cdn.mathpix.com/snip/images/NIQzwduMnV2r7sTSPvBIwi0Oy1pv3y0enY68T35-OtU.original.fullsize.png" width="965px" height="380px">
 -->
<img src="https://cdn.mathpix.com/snip/images/1-8xsj8IBFTSpWCP0CxxxN6jrKk-5sQYy8IiFjQ6OGQ.original.fullsize.png" width="640px">
How to Apply it in experiment? Aim: seprator the effect of propagation and source!

<img src="https://cdn.mathpix.com/snip/images/IAKHqKXxvfGm1AQIPDPOfGbXNu9Wbrd3eRqSyxiAWqg.original.fullsize.png" width="640px">

Step1: recontruct a propagator: assume all modes' amplitudes is same at z=0; we consider the source in the right, propagting to left with infinite duct, thus only care about left running modes.  
