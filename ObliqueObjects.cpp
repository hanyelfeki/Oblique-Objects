// Identifying Oblique Objects in a Cartesian Mesh
# include "ObliqueObjects.h"
#include <math.h>

void ObliqueObjects::constructObliqueObjects(int nObjects)
{
// nObjects: is the total number of the Oblique Objects in the computational domain
	int i,j,k;
// Allocate 2D dynamic arrays
	u_axis = new double*[nObjects];v_axis = new double*[nObjects];w_axis = new double*[nObjects];
	xyz0Corner=new double*[nObjects];
	lu=new double[nObjects]; lv=new double[nObjects];lw=new double[nObjects];

	for (i=0; i<nObjects; i++)
	{
	u_axis[i]=new double[3];v_axis[i]=new double[3];w_axis[i]=new double[3];
	xyz0Corner[i]=new double[3];
	}
// Allocate the 3D dynamic Array "Object_range"
	Object_range=new int**[nObjects];
	for (i=0; i<nObjects; i++)
	{
		Object_range[i]=new int*[3];
		for(j=0;j<3;j++)
		{
			Object_range[i][j]=new int[6];
		}
	}
// Initialize the arrays
	for (i=0; i<nObjects; i++)
	{
		for(j=0;j<3;j++)
		{
			for(k=0;k<6;k++)
			{
			Object_range[i][j][k]=0;
			}
		}
	}

	for (i=0; i<nObjects; i++)
	{
		for (j=0;j<3;j++)
		{
		u_axis[i][j]=0.0;v_axis[i][j]=0.0;w_axis[i][j]=0.0;
		xyz0Corner[i][j]=0.0;
		}
	}

	for (i=0; i<nObjects; i++)
	{
		lu[i]=0.0;lv[i]=0.0;lw[i]=0.0;
	}

}
void ObliqueObjects::destructObliqueObjects(int nObjects)
{
	int i; int j;
	for (i=0; i<nObjects; i++){

	delete[] u_axis[i]; delete[] v_axis[i];delete[] w_axis[i]; //  deallocate the 2D dynamic arrays - 1st part

	}

	delete[] u_axis;  delete[] v_axis;delete[] w_axis;  // deallocate the 2D dynamic arrays - 2nd part

		for (i=0; i<nObjects; i++)
	{

		for(j=0;j<3;j++)
		{
			delete[] Object_range[i][j];
		}
		delete[] Object_range[i];
	}
		delete[] Object_range;

}
void ObliqueObjects::ReadObjectInfoaAndFindMinMax(int nn,double x0,double y0,double z0,double luu,double lvv,double lww,int i_dir,double MyU_axis[3],double MyV_axis[3],double MyW_axis[3],double& xmin,double& xmax,double& ymin,double& ymax,double& zmin,double& zmax)
{
	// In this function the solver use the information that were passed from GUI about the rectangular oblique objects to find Xmin,Xmax,Ymin,Ymax,Zmin
	// As special cases the rectangular object may be a 2D (rectangle) if lu or lv=0 or a 1D (straight line) if lu=lv=0.
	// W-direction is the Longitudinal direction
	// This function will calculate finally in the Cartesian coordinate Xmin,Xmax,Ymin,Ymax,Zmin,Zmax
	// to be used in the Cartesian grid for the solver
	int i;
	double xyz1Corner[3], xyz2Corner[3], xyz3Corner[3], xyz4Corner[3], xyz5Corner[3], xyz6Corner[3],xyz7Corner[3];
	for(i=0;i<3; i++)
	{
	u_axis[nn][i]=MyU_axis[i];v_axis[nn][i]=MyV_axis[i];w_axis[nn][i]=MyW_axis[i];
	}
         xyz0Corner[nn][0]=x0; xyz0Corner[nn][1]=y0; xyz0Corner[nn][2]=z0;
         lu[nn]=luu; lv[nn] = lvv; lw[nn] = lww;
 		 findCornersObliqueRectangular(nn,x0,y0,z0,luu,lvv,lww,xyz1Corner,xyz2Corner,xyz3Corner,xyz4Corner,xyz5Corner,xyz6Corner,xyz7Corner);
		 xyzMinMaxRectangular(nn,xyz1Corner,xyz2Corner,xyz3Corner,xyz4Corner,xyz5Corner,xyz6Corner,xyz7Corner,xmin,xmax,ymin,ymax,zmin,zmax);
}
    void  ObliqueObjects::findCornersObliqueRectangular(int nn,double x0,double y0,double z0,double lu,double lv,double lw,double xyz1Corner[3],double xyz2Corner[3],double xyz3Corner[3],double xyz4Corner[3],double xyz5Corner[3],double xyz6Corner[3],double xyz7Corner[3])
	{
// This function finds 7 corners of the rectangular object while the first corner xyz0Corner is already given
    xyz1Corner[0]=x0+lu*u_axis[nn][0];xyz1Corner[1]=y0+lu*u_axis[nn][1];xyz1Corner[2]=z0+lu*u_axis[nn][2];
    xyz2Corner[0]=x0+lv*v_axis[nn][0];xyz2Corner[1]=y0+lv*v_axis[nn][1];xyz2Corner[2]=z0+lv*v_axis[nn][2];
    xyz3Corner[0]=x0+lw*w_axis[nn][0];xyz3Corner[1]=y0+lw*w_axis[nn][1];xyz3Corner[2]=z0+lw*w_axis[nn][2];
    xyz4Corner[0]=x0+lu*u_axis[nn][0]+lv*v_axis[nn][0]; xyz4Corner[1]=y0+lu*u_axis[nn][1]+lv*v_axis[nn][1];
    xyz4Corner[2]=z0+lu*u_axis[nn][2]+lv*v_axis[nn][2];
    xyz5Corner[0]=x0+lv*v_axis[nn][0]+lw*w_axis[nn][0];xyz5Corner[1]=y0+lv*v_axis[nn][1]+lw*w_axis[nn][1];
    xyz5Corner[2]=z0+lv*v_axis[nn][2]+lw*w_axis[nn][2];
    xyz6Corner[0]=x0+lu*u_axis[nn][0]+lw*w_axis[nn][0]; xyz6Corner[1]=y0+lu*u_axis[nn][1]+lw*w_axis[nn][1];
    xyz6Corner[2]=z0+lu*u_axis[nn][2]+lw*w_axis[nn][2];
    xyz7Corner[0]=x0+lu*u_axis[nn][0]+lv*v_axis[nn][0]+lw*w_axis[nn][0];
    xyz7Corner[1]=y0+lu*u_axis[nn][1]+lv*v_axis[nn][1]+lw*w_axis[nn][1];
    xyz7Corner[2]=z0+lu*u_axis[nn][2]+lv*v_axis[nn][2]+lw*w_axis[nn][2];
	}
  void ObliqueObjects::xyzMinMaxRectangular(int nn,double xyz1Corner[3],double xyz2Corner[3],double xyz3Corner[3],double xyz4Corner[3],double xyz5Corner[3],double xyz6Corner[3],double xyz7Corner[3],double& xmin,double& xmax,double& ymin,double& ymax,double& zmin,double& zmax)
    {
// Find Xmin,Xmax,Ymin,Ymax,Zmin,Zmax of the oblique rectangular object.
		xmin=0.0e+0;xmax=0.0e+0;ymin=0.0e+0;ymax=0.0e+0;zmin=0.0e+0;zmax=0.0e+0;
		double uFindMinMax[8];
		uFindMinMax[0]= xyz0Corner[nn][0]; uFindMinMax[1]=xyz1Corner[0]; uFindMinMax[2]=xyz2Corner[0];
        uFindMinMax[3]=xyz3Corner[0]; uFindMinMax[4]=xyz4Corner[0];
		uFindMinMax[5]=xyz5Corner[0];uFindMinMax[6]=xyz6Corner[0],uFindMinMax[7]=xyz7Corner[0];
		xmin=minOfArray(uFindMinMax);
        xmax= maxOfArray(uFindMinMax);
      uFindMinMax[0]= xyz0Corner[nn][1];uFindMinMax[1]=xyz1Corner[1];uFindMinMax[2]=xyz2Corner[1];
	  uFindMinMax[3]=xyz3Corner[1];uFindMinMax[4]=xyz4Corner[1]; uFindMinMax[5]=xyz5Corner[1];
	  uFindMinMax[6]=xyz6Corner[1],uFindMinMax[7]=xyz7Corner[1];
		ymin=minOfArray(uFindMinMax);
        ymax= maxOfArray(uFindMinMax);
      uFindMinMax[0]= xyz0Corner[nn][2];uFindMinMax[1]=xyz1Corner[2];uFindMinMax[2]=xyz2Corner[2];
	  uFindMinMax[3]=xyz3Corner[2];uFindMinMax[4]=xyz4Corner[2]; uFindMinMax[5]=xyz5Corner[2];
	  uFindMinMax[6]=xyz6Corner[2],uFindMinMax[7]=xyz7Corner[2];
		zmin=minOfArray(uFindMinMax);
        zmax= maxOfArray(uFindMinMax);
		return;
  }
   void ObliqueObjects::checkIfPointInObliqueRectangular(double x,double y,double z,double x0,double y0,double z0,double u_length,double v_length,double w_length,bool& isit_in,int nn)
   {
//  This function to check if the point (x,y,z) is inside the oblique rectangular object or outside.
      isit_in = false;
//
	  double u,v,w;
      u=(x-x0)*u_axis[nn][0] +(y-y0) * u_axis[nn][1] +(z-z0) * u_axis[nn][2];
      v=(x-x0)*v_axis[nn][0] +(y-y0) * v_axis[nn][1] +(z-z0) * v_axis[nn][2];
      w=(x-x0)*w_axis[nn][0] +(y-y0) * w_axis[nn][1] +(z-z0) * w_axis[nn][2];
//
      if(( u >=0. && u<=u_length) &&( v>=0. && v<=v_length) &&( w>=0. && w<=w_length))
		  isit_in= true;
   }

  double minOfArray(double x[8])
  {
	double xmin=1.0e34;
	  for(int i=0;i<8;i++)
	  {
		  if(x[i]<xmin)
			  xmin=x[i];
	  }
	  return xmin;
  }
  double maxOfArray(double x[8])
  {
	double xmax=-1.e34;
	  for(int i=0;i<8;i++)
	  {
		  if(x[i]>xmax)
			  xmax=x[i];
	  }
	  return xmax;
  }



	void constructSimpleMesh()
	{
		int i; int nx=100;
		// double xMesh[100],yMesh[100],zMesh[100],dx[100]
		for(i=0;i<nx;i++)
		{
			dx[i]=0.1;
		}
		for(i=0;i<nx;i++)
		{
			xMesh[i]=(double)i*dx[i];yMesh[i]=(double)i*dx[i];zMesh[i]=(double)i*dx[i];

		}
	}

			//

	   void testPoints()
	   {
		   // This subroutine is used to test points from the mesh

		   int nn,i,j,k;

		   double xi,yj,zk;
		   double x0,y0,z0,lu,lv,lw;
		   bool isit_in=false;

		   for (nn=0;nn<ne_oblique;nn++)
		   {


			   for (i=Object1.Object_range[nn][2][0];i<Object1.Object_range[nn][2][1]+1;i++)
			   {
				   for (j=Object1.Object_range[nn][2][2];j<Object1.Object_range[nn][2][3]+1;j++)
				   {
					   for (k=Object1.Object_range[nn][2][4];k<Object1.Object_range[nn][2][5]+1;k++)
					   {
						   xi=xMesh[i]+dx[i]/2.;
						   yj=yMesh[j];
                           zk=zMesh[k];
					   //
                       x0=Object1.xyz0Corner[nn][0];
                       y0=Object1.xyz0Corner[nn][1];
                       z0=Object1.xyz0Corner[nn][2];
                       lu=Object1.lu[nn];
                       lv=Object1.lv[nn];
                       lw=Object1.lw[nn];

   //
                       Object1.checkIfPointInObliqueRectangular(xi,yj,zk,x0,y0,z0,lu,lv,lw,isit_in,nn);
					   if(isit_in == true)
					   {
                        cout << "xi= " << xi <<endl;cout << "yj= " << yj <<endl; cout << "zk= " << zk <<endl;
					   }
				   }
			   }
			   }

			   }

	   }

int main()
{
	// Example

	constructSimpleMesh();
	Object1.constructObliqueObjects(ne_oblique);
	double x0=5.0, y0=5., z0=4., lu=1.,lv=1.,lw=2.;
	double a1_axis[3],a2_axis[3],a3_axis[3];
	a1_axis[0]=0.;a1_axis[1]=0.;a1_axis[2]=1.0;
	a2_axis[0]=1.0;a2_axis[1]=0.0;a2_axis[2]=0.0;
	a3_axis[0]=0.0;a3_axis[1]=1.0;a3_axis[2]=0.0;
	int idir=3, nn=0;
	double xmin=0.,xmax=0.,ymin=0.,ymax=0.,zmin=0.,zmax=0.;
	Object1.ReadObjectInfoaAndFindMinMax(nn,x0,y0,z0,lu,lv,lw,idir,a1_axis,a2_axis,a3_axis,xmin,xmax,ymin,ymax,zmin,zmax);
	cout<<"xmin,xmax,ymin,ymax,zmin,zmax ="<<endl;
	cout << "Xmin= " << xmin <<endl;
	cout << "Xmax= " << xmax <<endl;
	cout << "Ymin= " << ymin <<endl;
	cout << "Ymax= " << ymax <<endl;
    cout << "Zmin= " << zmin <<endl;
	cout << "Zmax= " << zmax <<endl;
	Object1.Object_range[nn][2][0]=(int)(xmin/dx[0]);Object1.Object_range[nn][2][1]=(int)(xmax/dx[0]);
	Object1.Object_range[nn][2][2]=(int)(ymin/dx[0]);Object1.Object_range[nn][2][3]=(int)(ymax/dx[0]);
	Object1.Object_range[nn][2][4]=(int)(zmin/dx[0]);Object1.Object_range[nn][2][5]=(int)(zmax/dx[0]);

	cout << "dx[0]=" <<dx[0]<<endl;
	cout << "I_Xmin= " << Object1.Object_range[nn][2][0] <<endl;
	cout << "I_Xmax= " << Object1.Object_range[nn][2][1] <<endl;
	cout << "J_Ymin= " << Object1.Object_range[nn][2][2] <<endl;
	cout << "J_Ymax= " << Object1.Object_range[nn][2][3] <<endl;
    cout << "K_Zmin= " << Object1.Object_range[nn][2][4] <<endl;
	cout << "K_Zmax= " << Object1.Object_range[nn][2][5] <<endl;

		testPoints();

	Object1.destructObliqueObjects(ne_oblique);
	return 0;
}
