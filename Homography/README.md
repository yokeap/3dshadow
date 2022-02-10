# Homograhy estimation

Homography matrix is estimate from undistored image, so camera calibration must be done properly. Single square patern is used as plane reference, after loading image, starting click on left upper of square pattern then shifting to right upper (clockwise direction).

```bash
├── Homography
    |-- findHomography.m         // estimate homography
    |-- poseEstimation.m         // measure object with homography estimated
    |-- ..
    |-- .. 
    |-- (tooling functions)
```