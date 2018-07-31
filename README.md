# Identifying Oblique Objects in a Cartesian Mesh - C++ Code
This version of the code includes functions for oblique rectangular objects but we will consider general oblique objects with arbitrary cross sections.
Example of the Cartesian Mesh is the Yee Cells of the FDTD Algorithm.
Each rectangular object is defined by one corner (x0,y0,z0), the directional unit vectors u_axis[nn][:],v_axis[nn][:],w_axis[nn][:],and the lengths along the the three directions lu[nn],lv[nn],lw[nn] respectively.
The function ReadObjectInfoaAndFindMinMax is used use the information that were passed from GUI about the rectangular oblique objects to find Xmin,Xmax,Ymin,Ymax,Zmin
As special cases the rectangular object may be a 2D (rectangle) if lu or lv=0 or a 1D (straight line) if lu=lv=0.
W-direction is the Longitudinal direction
This function will calculate finally in the Cartesian coordinate Xmin,Xmax,Ymin,Ymax,Zmin,Zmax.
The function findCornersObliqueRectangular is used to find 7 corners of the rectangular object while the first corner xyz0Corner is already given
The function xyzMinMaxRectangular is used to find Xmin,Xmax,Ymin,Ymax,Zmin,Zmax of the oblique rectangular object
The function checkIfPointInObliqueRectangular is used to to check if the point (x,y,z) is inside the oblique rectangular object or not.
nObjects: is the total number of the Oblique Objects in the computational domain
