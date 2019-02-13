imOrig = imread('DSC_0298.JPG');
[im, newOrigin] = undistortImage(imOrig, cameraParams, 'OutputView', 'full');
imshow(im);


% Detect the checkerboard. (distorted image)
[imagePoints, boardSize] = detectCheckerboardPoints(im);

% Generate the world coordinates of the checkerboard corners in the
% pattern-centric coordinate system, with the upper-left corner at (0,0).
squareSize = 33.10; % in millimeters
worldPoints = generateCheckerboardPoints(boardSize, squareSize);

% Detect the checkerboard. (in undistorted image)
[imagePoints, boardSize] = detectCheckerboardPoints(im);

% Adjust the imagePoints so that they are expressed in the coordinate system
% used in the original image, before it was undistorted.  This adjustment
% makes it compatible with the cameraParameters object computed for the original image.
imagePoints = imagePoints + newOrigin; % adds newOrigin to every row of imagePoints

% Compute rotation and translation of the camera.
[R, t] = extrinsics(imagePoints, worldPoints, cameraParams);