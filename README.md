# 3dshadow

3dshadow is my research project (KMITL Dr.Eng program). 

It's assume that most of fruits shape is symmetrical object. Measuring volume by means of 3D machine vision is can be carried out.

This project is used fast and simple technique to reconstructs 3D shape from 2D single image, makes cost effecting against another method. 

Object height is computed from object shadow casting, by means of every crossing point of inverse ray tracing from edge shadow and upper edge shadow on object is compute from cosine angle. 

With proper camera calibration, undistorted image is used in whole process with homography estimation for transform image unit to 3D world unit
#
### This repository is contain testing script function based on MATLAB:

#### main script is divided in to two method: 
* cross inverse tracing between lower edge of shadow casted and upper edge of shadow on object.
* cross inverse tracing between lower edge of shadow casted and pseudo skeleton by means of object thinning.
