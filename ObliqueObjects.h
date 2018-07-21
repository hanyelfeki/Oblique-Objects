//#include "stdafx.h"
//#include "stdafx.h"
#include <iostream>
 using namespace std;
 double minOfArray(double x[7]);
 double maxOfArray(double x[7]);
 int ne_oblique=1; // ne_oblique is the number of oblique objects.


  double xMesh[100],yMesh[100],zMesh[100],dx[100]; // This is only a virtial mesh. In the original program
                                                  // It is calculated by a mesh generation code
                                                 // Also they should be a dynamic array.


class ObliqueObjects
{
	// This is a class of Rectangular Oblique Objects
	// Each object has (x0,y0,z0) as a refernce corner of the object, lu,lv, lw are the lengths along the
	// curvilinear coordinates. u_axis, v_axis, w_axis are the unit vectors along the 3 axes of the curvilinear coordinates.
public:
	double** u_axis; double** v_axis; double** w_axis;
	double** xyz0Corner; double** xyz7Corner; // First corner and last corner
	double* lu; double* lv; double* lw;
	int *** source_range;
	void constructObliqueObjects(int nObjects);
	void destructObliqueObjects(int nObjects);
	void ReadObjectInfoaAndFindMinMax(int nn,double x0,double y0,double z0,double lu,double lv,double lw,int i_dir,double a1_axis[3],double MyV_axis[3],double MyW_axis[3],double& xmin,double& xmax,double& ymin,double& ymax,double& zmin,double& zmax);
	void  findCornersObliqueRectangular(int nn,double x0,double y0,double z0,double lu,double lv,double lw,double xyz1Corner[3],double xyz2Corner[3],double xyz3Corner[3],double xyz4Corner[3],double xyz5Corner[3],double xyz6Corner[3]);
	void  xyzMinMaxRectangular(int nn,double xyz1Corner[3],double xyz2Corner[3],double xyz3Corner[3],double xyz4Corner[3],double xyz5Corner[3],double xyz6Corner[3],double& xmin,double& xmax,double& ymin,double& ymax,double& zmin,double& zmax);
	void checkIfPointInObliqueRectangular(double x,double y,double z,double x0,double y0,double z0,double u_length,double v_length,double w_length,bool& isit_in,int nn);
}
;
 ObliqueObjects Object1;
