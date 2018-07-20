// Oblique Objects in general. It can used for field sources
# include "ObliqueObjects.h"
#include <math.h>

void ObliqueObjects::constructObliqueObjects(int nObjects)
{
// nObjects the total number of the Oblique Objects in the computational domain
	int i,j,k;
// Allocate the 2D dynamic arrays
	u_axis = new double*[nObjects];v_axis = new double*[nObjects];w_axis = new double*[nObjects];
	xyz0Corner=new double*[nObjects];xyz7Corner=new double*[nObjects];
	lu=new double[nObjects]; lv=new double[nObjects];lw=new double[nObjects];

	for (i=0; i<nObjects; i++)
	{
	u_axis[i]=new double[3];v_axis[i]=new double[3];w_axis[i]=new double[3];
	xyz0Corner[i]=new double[3];xyz7Corner[i]=new double[3];
	}
// Allocate the 3D dynamic Array "source_range"
	source_range=new int**[nObjects];
	for (i=0; i<nObjects; i++)
	{
		source_range[i]=new int*[3];
		for(j=0;j<3;j++)
		{
			source_range[i][j]=new int[6];
		}
	}
// Initialize the arrays
	for (i=0; i<nObjects; i++)
	{
		for(j=0;j<3;j++)
		{
			for(k=0;k<6;k++)
			{
			source_range[i][j][k]=0;
			}
		}
	}

	for (i=0; i<nObjects; i++)
	{
		for (j=0;j<3;j++)
		{
		u_axis[i][j]=0.0;v_axis[i][j]=0.0;w_axis[i][j]=0.0;
		xyz0Corner[i][j]=0.0;xyz7Corner[i][j]=0.0;
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
			delete[] source_range[i][j];
		}
		delete[] source_range[i];
	}
		delete[] source_range;

}
void ObliqueObjects::ReadObjectInfoaAndFindMinMax(int nn,double x0,double y0,double z0,double luu,double lvv,double lww,int i_dir,double a1_axis_Esource[3],double a2_axis_Esource[3],double a3_axis_Esource[3],double& xmin,double& xmax,double& ymin,double& ymax,double& zmin,double& zmax)
{
	// In this function the solver use the information that were passed from GUI about objects
	// As a special case we have a source or load can be between 2 points, 2 parallel wires, or 2 parallel plates.
	// W-direction is the direction of the source
	// This function will calaculate finally in the Caretizian coordinate Xmin,Xmax,Ymin,Ymax,Zmin,Zmax
	// to be used in the Cartezian grid for the solver
	int i;
	double xyz1Corner[3], xyz2Corner[3], xyz3Corner[3], xyz4Corner[3], xyz5Corner[3], xyz6Corner[3];
	for(i=0;i<3; i++)
	{
	u_axis[nn][i]=a1_axis_Esource[i];v_axis[nn][i]=a2_axis_Esource[i];w_axis[nn][i]=a3_axis_Esource[i];
	}
         xyz0Corner[nn][0]=x0; xyz0Corner[nn][1]=y0; xyz0Corner[nn][2]=z0;
         lu[nn]=luu; lv[nn] = lvv; lw[nn] = lww;
 		 findCornersObliqueRectangular(nn,x0,y0,z0,luu,lvv,lww,xyz1Corner,xyz2Corner,xyz3Corner,xyz4Corner,xyz5Corner,xyz6Corner);
		 xyzMinMaxRectangular(nn,xyz1Corner,xyz2Corner,xyz3Corner,xyz4Corner,xyz5Corner,xyz6Corner,xmin,xmax,ymin,ymax,zmin,zmax);
}
    void  ObliqueObjects::findCornersObliqueRectangular(int nn,double x0,double y0,double z0,double lu,double lv,double lw,double xyz1Corner[3],double xyz2Corner[3],double xyz3Corner[3],double xyz4Corner[3],double xyz5Corner[3],double xyz6Corner[3])
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
    xyz7Corner[nn][0]=x0+lu*u_axis[nn][0]+lv*v_axis[nn][0]+lw*w_axis[nn][0];
    xyz7Corner[nn][1]=y0+lu*u_axis[nn][1]+lv*v_axis[nn][1]+lw*w_axis[nn][1];
    xyz7Corner[nn][2]=z0+lu*u_axis[nn][2]+lv*v_axis[nn][2]+lw*w_axis[nn][2];
	}
  void ObliqueObjects::xyzMinMaxRectangular(int nn,double xyz1Corner[3],double xyz2Corner[3],double xyz3Corner[3],double xyz4Corner[3],double xyz5Corner[3],double xyz6Corner[3],double& xmin,double& xmax,double& ymin,double& ymax,double& zmin,double& zmax)
    {
// Find Xmin,Xmax,Ymin,Ymax,Zmin,Zmax of the oblique rectangular object.
		xmin=0.0e+0;xmax=0.0e+0;ymin=0.0e+0;ymax=0.0e+0;zmin=0.0e+0;zmax=0.0e+0;
		double uFindMinMax[8];
		uFindMinMax[0]= xyz0Corner[nn][0]; uFindMinMax[1]=xyz1Corner[0]; uFindMinMax[2]=xyz2Corner[0];
        uFindMinMax[3]=xyz3Corner[0]; uFindMinMax[4]=xyz4Corner[0];
		uFindMinMax[5]=xyz5Corner[0];uFindMinMax[6]=xyz6Corner[0],uFindMinMax[7]=xyz7Corner[nn][0];
		xmin=minOfArray(uFindMinMax);
        xmax= maxOfArray(uFindMinMax);
      uFindMinMax[0]= xyz0Corner[nn][1];uFindMinMax[1]=xyz1Corner[1];uFindMinMax[2]=xyz2Corner[1];
	  uFindMinMax[3]=xyz3Corner[1];uFindMinMax[4]=xyz4Corner[1]; uFindMinMax[5]=xyz5Corner[1];
	  uFindMinMax[6]=xyz6Corner[1],uFindMinMax[7]=xyz7Corner[nn][1];
		ymin=minOfArray(uFindMinMax);
        ymax= maxOfArray(uFindMinMax);
      uFindMinMax[0]= xyz0Corner[nn][2];uFindMinMax[1]=xyz1Corner[2];uFindMinMax[2]=xyz2Corner[2];
	  uFindMinMax[3]=xyz3Corner[2];uFindMinMax[4]=xyz4Corner[2]; uFindMinMax[5]=xyz5Corner[2];
	  uFindMinMax[6]=xyz6Corner[2],uFindMinMax[7]=xyz7Corner[nn][2];
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

void Source::constructSource(int nObjects)
{
	int i; int j;
	// Allocate index_source
	index_source=new int*[nObjects+2];
	for(i=0;i<nObjects+2;i++)
	{
		index_source[i]=new int[4];
	}
	// Initialize index_source
	for(i=0;i<nObjects+2;i++)
	{
		for(j=0;j<4;j++)   //  j=0,1,2 for x,y,z-directions, j=3 for oblique direction
		{
			index_source[i][j]=0;  // this index is used for sources and is not used for loads
		}
	}
}
void Source::destruct_Source(int nObjects)
{
	int i;
	for(i=0;i<nObjects+2;i++)
	{
		delete[]index_source[i];
	}
	delete[]index_source;
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
			xMesh[i]=(double)i*dx[i];yMesh[i]=i*dx[i];zMesh[i]=(double)i*dx[i];

		}
	}

			//

	   void testPoints()
	   {
		   // This subroutine is used for the field source which can be  E-field sources, H-Field source, voltage lumped source, or Waveguide source

		   int nn,i,j,k;

		   double xi,yj,zk;
		   double x0,y0,z0,lu,lv,lw;
		   bool isit_in=false;
		   // E-field source
		   //  E-field source along X-direction

		   for (nn=0;nn<ne_oblique;nn++)
		   {
	                        //  as it is not the complete code


			   for (i=Object1.source_range[nn][2][0];i<Object1.source_range[nn][2][1]+1;i++)
			   {
				   for (j=Object1.source_range[nn][2][2];j<Object1.source_range[nn][2][3]+1;j++)
				   {
					   for (k=Object1.source_range[nn][2][4];k<Object1.source_range[nn][2][5]+1;k++)
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

			   // The above was for updating the ex component only. The original program includes also updating ey, and ez components.
			   // In this algorithm ex, ey, and ez are not located at the same point in the cell.
			   }

	   }

int main()
{
	// Example
	              // it is usually calculated from the duration that the user defines through the UI.
	constructSimpleMesh();
	Object1.constructObliqueObjects(ne_oblique);
	double x0=5.0, y0=5., z0=4., lu=1.,lv=1.,lw=2.;
	double a1_axis[3],a2_axis[3],a3_axis[3];
	a1_axis[0]=0.;a1_axis[1]=0.;a1_axis[2]=1.0;
	a2_axis[0]=1.0;a2_axis[1]=0.0;a2_axis[2]=0.0;
	a3_axis[0]=0.0;a3_axis[1]=1.0;a3_axis[2]=0.0;
	int idir=3, nn=0;
	double n1_x=0.,n2_x=0.,n3_y=0.,n4_y=0.,n5_z=0.,n6_z=0.;
	Object1.ReadObjectInfoaAndFindMinMax(nn,x0,y0,z0,lu,lv,lw,idir,a1_axis,a2_axis,a3_axis,n1_x,n2_x,n3_y,n4_y,n5_z,n6_z);
	cout<<"n1_x,n2_x,n3_y,n4_y,n5_z,n6_z ="<<endl;
	cout << "Xmin= " << n1_x <<endl;
	cout << "Xmax= " << n2_x <<endl;
	cout << "Ymin= " << n3_y <<endl;
	cout << "Ymax= " << n4_y <<endl;
    cout << "Zmin= " << n5_z <<endl;
	cout << "Zmax= " << n6_z <<endl;
	Object1.source_range[nn][2][0]=(int)(n1_x/dx[0]);Object1.source_range[nn][2][1]=(int)(n2_x/dx[0]);
	Object1.source_range[nn][2][2]=(int)(n3_y/dx[0]);Object1.source_range[nn][2][3]=(int)(n4_y/dx[0]);
	Object1.source_range[nn][2][4]=(int)(n5_z/dx[0]);Object1.source_range[nn][2][5]=(int)(n6_z/dx[0]);

	cout << "dx[0]=" <<dx[0]<<endl;
	cout << "I_Xmin= " << Object1.source_range[nn][2][0] <<endl;
	cout << "I_Xmax= " << Object1.source_range[nn][2][1] <<endl;
	cout << "J_Ymin= " << Object1.source_range[nn][2][2] <<endl;
	cout << "J_Ymax= " << Object1.source_range[nn][2][3] <<endl;
    cout << "K_Zmin= " << Object1.source_range[nn][2][4] <<endl;
	cout << "K_Zmax= " << Object1.source_range[nn][2][5] <<endl;

		testPoints();

	Object1.destructObliqueObjects(ne_oblique);
	return 0;
}
