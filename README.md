# Identifying Oblique Objects in a Cartesian Mesh - C++ Code
# We will describe the main ideas of this code soon.
This version of the code is written for oblique rectangular objects but we will consider general oblique objects with arbitrary cross sections.
Each rectangular object is defined by one corner (x0,y0,z0), the directional unit vectors u_axis[nn][:],v_axis[nn][:],w_axis[nn][:],and the lengths along the the three directions lu[nn],lv[nn],lw[nn] respectively.
The function ReadObjectInfoaAndFindMinMax is used to use the information that were passed from GUI about the rectangular oblique objects to find Xmin,Xmax,Ymin,Ymax,Zmin for this object. It calls xyzMinMaxRectangular for the case of rectangular objects.
The function findCornersObliqueRectangular is used to find 7 corners of the rectangular object while the first corner xyz0Corner is already given
The function xyzMinMaxRectangular is used to find Xmin,Xmax,Ymin,Ymax,Zmin,Zmax of the oblique rectangular object
