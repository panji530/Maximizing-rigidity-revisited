# Maximizing-rigidity-revisited
This code is for the following paper:

Pan Ji, Hongdong Li, Yuchao Dai, and Ian Reid. "Maximizing rigidity" revisited: a convex programming approach for generic 3D shape reconstruction from multiple perspective views. in ICCV 2017.

The code relies on CVX to solve the SDP problem. To test the code, please install cvx from http://cvxr.com/cvx/download/, and add its path to your Matlab. 

Non-rigid surfaces: run the following files

test_Hulk.m

test_Tshirt.m

test_KINECT_Paper.m

Rigid shapes: run the following files

test_Rigid_ModelHouse.m

test_Rigid_Sythetic.m

Articulated motion: run the following files

test_Articulation_Synthetic.m

test_Articulation_Dance.m

