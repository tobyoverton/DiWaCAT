//
//  RoFSuP.h
//  RoFSuP
//
//  Created by MacBook on 8/28/17.
//  Copyright © 2017 Casto. All rights reserved.
//
//  Edited by T.Overton December 2019

#ifndef RoFSuP_h
#define _USE_MATH_DEFINES
#include <ctime>
#include <math.h>
#include <iostream>
#include <algorithm>
#include <chrono>
#include <complex>
#include <vector>
#include <fstream>
#include <numeric>
#include <functional>

#include <H5Cpp.h>
#include "CircularGreenFunction.h"
#include "PlanarGreensFunctions.h"

/*---------
    Code is set up in sections:
            *Calculating the beam distributions
            *The numbers needed for eigenvalues
            *Individual point wake potential (Green's function)
            *Overall function for calculating the field
 --------*/

using std::vector;
using std::ifstream;

//Find the closest point to a macroparticle from the grid
//Output the array position of the nearest point
int closest3D(double xWanted, double yWanted, double zWanted, vector<double> vx, vector<double> vy, vector<double> vz)
{
    //Return the vector position of the nearest element of a vector to a position
    int j(0);
    double distanceToPoint = pow((xWanted-vx[0]),2) + pow((yWanted-vy[0]),2) + pow((zWanted-vz[0]),2);
    for (int i=1; i<(int) vx.size(); i++)
    {
        double intDist = pow((xWanted-vx[i]),2) + pow((yWanted-vy[i]),2) + pow((zWanted-vz[i]),2);
        if(intDist<distanceToPoint){
            distanceToPoint=intDist;
            j=i;
        }
    }
    return j;
}

//Check if positive and give +/-1 for each
int PositiveNumber(double Number){
    if(Number >= 0) return 1;
    else return -1;
}

//------------------------------------------------------------------------//
//-----------            Trilinear Interpolation          ----------------//
//------------------------------------------------------------------------//


vector<double> TriLinearFormal(vector<double> Fx, vector<double> Fy, vector<double> Ez, vector<double> xB, vector<double> yB, vector<double> zB, double x0, double y0, double z0)
{
    //Find
    double FxInterpolation(0), FyInterpolation(0), EzInterpolation(0);
    vector<int> CubeArrayPositions = {0,0,0,0,0,0,0,0};
    vector<double> InverseDistances = {0,0,0,0,0,0,0,0};
    
    //Use the closest position to get the cube
    //Use the repeating numbers in the grid to get the distance between grid positions
    vector<double> xValues(xB), yValues(yB), zValues(zB);
    std::sort(xValues.begin(),xValues.end());
    auto lastx = std::unique(xValues.begin(),xValues.end());
    xValues.erase(lastx,xValues.end());
    std::sort(yValues.begin(),yValues.end());
    auto lasty = std::unique(yValues.begin(),yValues.end());
    yValues.erase(lasty,yValues.end());
    std::sort(zValues.begin(),zValues.end());
    auto lastz = std::unique(zValues.begin(),zValues.end());
    zValues.erase(lastz,zValues.end());
    double xGrid(abs(xValues[1] - xValues[0])), yGrid(abs(yValues[1] - yValues[0])), zGrid(abs(zValues[1] - zValues[0]));
    
    //Find the closest position
    CubeArrayPositions[0] = closest3D(x0, y0, z0, xB, yB, zB);
    vector<double> c000 = {Fx[CubeArrayPositions[0]],Fy[CubeArrayPositions[0]],Ez[CubeArrayPositions[0]]};
    double ClosestDistance = pow(pow((x0-xB[CubeArrayPositions[0]]),2) + pow((y0-yB[CubeArrayPositions[0]]),2) + pow((z0-zB[CubeArrayPositions[0]]),2),0.5);
    double xD = abs((x0-xB[CubeArrayPositions[0]])/xGrid);
    double yD = abs((y0-yB[CubeArrayPositions[0]])/yGrid);
    double zD = abs((z0-zB[CubeArrayPositions[0]])/zGrid);
    //Just double check the macroparticle isn't exactly at grid position (they're be a 1/0 factor if there is)
    if(ClosestDistance != 0){
        InverseDistances[0]= 1/ClosestDistance;
        double xClosest = xB[CubeArrayPositions[0]];
        double yClosest = yB[CubeArrayPositions[0]];
        double zClosest = zB[CubeArrayPositions[0]];
    
        //Use whether the closest point is ahead or behind the macroparticle and know whether need to add or subtract to get the result
        double xDifference(PositiveNumber(x0-xClosest)*xGrid);
        double yDifference(PositiveNumber(y0-yClosest)*yGrid);
        double zDifference(PositiveNumber(z0-zClosest)*zGrid);
        CubeArrayPositions[1] = closest3D(xClosest+xDifference, yClosest, zClosest, xB, yB, zB);
        vector<double> c100 = {Fx[CubeArrayPositions[1]],Fy[CubeArrayPositions[1]],Ez[CubeArrayPositions[1]]};
        CubeArrayPositions[2] = closest3D(xClosest, yClosest+yDifference, zClosest, xB, yB, zB);
        vector<double> c010 = {Fx[CubeArrayPositions[2]],Fy[CubeArrayPositions[2]],Ez[CubeArrayPositions[2]]};
        CubeArrayPositions[3] = closest3D(xClosest, yClosest, zClosest+zDifference, xB, yB, zB);
        vector<double> c001 = {Fx[CubeArrayPositions[3]],Fy[CubeArrayPositions[3]],Ez[CubeArrayPositions[3]]};
        CubeArrayPositions[4] = closest3D(xClosest+xDifference, yClosest+yDifference, zClosest, xB, yB, zB);
        vector<double> c110 = {Fx[CubeArrayPositions[4]],Fy[CubeArrayPositions[4]],Ez[CubeArrayPositions[4]]};
        CubeArrayPositions[5] = closest3D(xClosest+xDifference, yClosest, zClosest+zDifference, xB, yB, zB);
        vector<double> c101 = {Fx[CubeArrayPositions[5]],Fy[CubeArrayPositions[5]],Ez[CubeArrayPositions[5]]};
        CubeArrayPositions[6] = closest3D(xClosest, yClosest+yDifference, zClosest+zDifference, xB, yB, zB);
        vector<double> c011 = {Fx[CubeArrayPositions[6]],Fy[CubeArrayPositions[6]],Ez[CubeArrayPositions[6]]};
        CubeArrayPositions[7] = closest3D(xClosest+xDifference, yClosest+yDifference, zClosest+zDifference, xB, yB, zB);
        vector<double> c111 = {Fx[CubeArrayPositions[7]],Fy[CubeArrayPositions[7]],Ez[CubeArrayPositions[7]]};
        
        vector<double> c00 = {(1-xD)*c000[0] + xD*c100[0],(1-xD)*c000[1] + xD*c100[1],(1-xD)*c000[2] + xD*c100[2]};
        vector<double> c01 = {(1-xD)*c001[0] + xD*c101[0],(1-xD)*c001[1] + xD*c101[1],(1-xD)*c001[2] + xD*c101[2]};
        vector<double> c10 = {(1-xD)*c010[0] + xD*c110[0],(1-xD)*c010[1] + xD*c110[1],(1-xD)*c010[2] + xD*c110[2]};
        vector<double> c11 = {(1-xD)*c011[0] + xD*c111[0],(1-xD)*c011[1] + xD*c111[1],(1-xD)*c011[2] + xD*c111[2]};
        
        vector<double> c0 = {(1-yD)*c00[0] +yD*c10[0],(1-yD)*c00[1] +yD*c10[1],(1-yD)*c00[2] +yD*c10[2]};
        vector<double> c1 = {(1-yD)*c01[0] +yD*c11[0],(1-yD)*c01[1] +yD*c11[1],(1-yD)*c01[2] +yD*c11[2]};
        
        vector<double> c = {(1-zD)*c0[0] +zD*c1[0],(1-zD)*c0[1] +zD*c1[1],(1-zD)*c0[2] +zD*c1[2]};
        FxInterpolation=c[0];
        FyInterpolation=c[1];
        EzInterpolation=c[2];
        //std::cout<<std::endl;
    }
    else{
        FxInterpolation=Fx[CubeArrayPositions[0]];
        FyInterpolation=Fy[CubeArrayPositions[0]];
        EzInterpolation=Ez[CubeArrayPositions[0]];
    }
    vector<double> InterpolationOutput = {FxInterpolation,FyInterpolation,EzInterpolation};
    return InterpolationOutput;
}

//-------------------------------------------------------------------//
//-------------     HDF5 File Functions         ---------------------//
//-------------------------------------------------------------------//

//Read in HDF5 files and output the vectors to use
//Beam only
void HDF5Reader(std::string filename,vector<double> &x, vector<double> &y, vector<double> &z, vector<double> &px, vector<double> &py, vector<double> &pz, vector<double> &t, vector<double> &q, double y0, double x0){
    using namespace H5;
    //Have the filename (input) and where in the file our data is
    const H5std_string    FILE_NAME(filename);
    const H5std_string    DATASET_NAME("beam/beam");
    //We know how wide the data is (6D + t + charge)
  
    
    H5File fp( FILE_NAME, H5F_ACC_RDONLY );
    DataSet dset = fp.openDataSet(DATASET_NAME);
    
    // get the dataspace
    DataSpace dspace = dset.getSpace();
    // get the dataset type class
    H5T_class_t type_class = dset.getTypeClass();
    
    // get the size of the dataset
    hsize_t rank;
    hsize_t dims[2];
    rank = dspace.getSimpleExtentDims(dims, NULL); // rank = 1

    // Define the memory dataspace
    hsize_t dimsm[2];
    dimsm[0] = dims[0];
    dimsm[1] = dims[1];
    DataSpace memspace ((int) rank,dimsm);
    
    // create a vector the same size as the dataset
    
    vector<float> data;
    data.resize(dims[0]);
     
    int nrow = (int) dims[0];
    int ncol = (int) dims[1];
    
    float *data_out=new float[nrow*ncol];
    
    for (int i=0;i<nrow*ncol;i++)
    {
        data_out[i] = 0;
    }
    
    dset.read(data_out, H5::PredType::NATIVE_FLOAT, memspace,dspace);
    
    fp.close();
    
    x.resize(nrow);
    y.resize(nrow);
    z.resize(nrow);
    px.resize(nrow);
    py.resize(nrow);
    pz.resize(nrow);
    t.resize(nrow);
    q.resize(nrow);
    
    for(int i=0; i<nrow; i++){
        x[i] = (x0 + (data_out[8*i]*100));
        y[i] = (y0 + (data_out[8*i + 1]*100));
        z[i] = data_out[8*i + 2] * 100;
        px[i] = data_out[8*i + 3];
        py[i] = data_out[8*i + 4];
        pz[i] = data_out[8*i + 5];
        t[i] = data_out[8*i + 6];
        q[i] = data_out[8*i + 7];
    }
}

//Meshed beam reading
void HDF5ReaderMacro(std::string filename, vector<double> &x, vector<double> &y, vector<double> &z, vector<double> &px, vector<double> &py, vector<double> &pz, vector<double> &t, vector<double> &q, double y0, double x0, double MaxZ) {
	using namespace H5;
	//Have the filename (input) and where in the file our data is
	const H5std_string    FILE_NAME(filename);
	const H5std_string    DATASET_NAME("MeshValues/macrosAtPoint");
	const H5std_string    DATASET_NAME_XVALS("MeshValues/xPoints");
	const H5std_string    DATASET_NAME_YVALS("MeshValues/yPoints");
	const H5std_string    DATASET_NAME_ZVALS("MeshValues/zPoints");
	const H5std_string    DATASET_NAME_CHARGE("MeshParameters/chargePerMacro");
	//We know how wide the data is (6D + t + charge)


	H5File fp(FILE_NAME, H5F_ACC_RDONLY);
	DataSet dset = fp.openDataSet(DATASET_NAME);
	DataSet dset_xvals = fp.openDataSet(DATASET_NAME_XVALS);
	DataSet dset_yvals = fp.openDataSet(DATASET_NAME_YVALS);
	DataSet dset_zvals = fp.openDataSet(DATASET_NAME_ZVALS);
	DataSet dset_q = fp.openDataSet(DATASET_NAME_CHARGE);

	// get the dataspace
	DataSpace dspace = dset.getSpace();
	DataSpace dspace_xvals = dset_xvals.getSpace();
	DataSpace dspace_yvals = dset_yvals.getSpace();
	DataSpace dspace_zvals = dset_zvals.getSpace();
	DataSpace dspace_q = dset_q.getSpace();

	// get the dataset type class
	H5T_class_t type_class = dset.getTypeClass();

	// get the size of the dataset
	hsize_t rank;
	hsize_t rankVALS;
	hsize_t rankQ;
	hsize_t dims[3];
	hsize_t dimsVALS[1];
	hsize_t dimsQ[1];
	rank = dspace.getSimpleExtentDims(dims, NULL); // rank = 1
	//rank is same for x,y,z values so no need to do this 3 times
	rankVALS = dspace_xvals.getSimpleExtentDims(dimsVALS, NULL);
	rankQ = dspace_q.getSimpleExtentDims(dimsQ, NULL);

	// Define the memory dataspace
	hsize_t dimsm[3];
	dimsm[0] = dims[0];
	dimsm[1] = dims[1];
	dimsm[2] = dims[2];
	hsize_t dimsmXVALS[1];
	dimsmXVALS[0] = dims[0];
	hsize_t dimsmYVALS[1];
	dimsmYVALS[0] = dims[1];
	hsize_t dimsmZVALS[1];
	dimsmZVALS[0] = dims[2];
	hsize_t dimsmQ[1];
	dimsmQ[0] = dimsQ[0];

	DataSpace memspace((int)rank, dimsm);
	DataSpace memspace_xval((int)rankVALS, dimsmXVALS);
	DataSpace memspace_yval((int)rankVALS, dimsmYVALS);
	DataSpace memspace_zval((int)rankVALS, dimsmZVALS);
	DataSpace memspace_q((int)rankQ, dimsmQ);

	// create a vector the same size as the dataset

	vector<float> data;
	data.resize(dims[0]);

	int nrow = (int)dims[0];
	int ncol = (int)dims[1];
	int nz = (int)dims[2];

	float *data_out = new float[nrow*ncol*nz];
	float *xvals = new float[nrow];
	float *yvals = new float[ncol];
	float *zvals = new float[nz];
	float *qval = new float[1];

	for (int i = 0; i < nrow*ncol*nz; i++)
	{
		data_out[i] = 0;
	}
	for (int i = 0; i < nrow; i++)
	{
		xvals[i] = 0;
	}
	for (int i = 0; i < ncol; i++)
	{
		yvals[i] = 0;
	}
	for (int i = 0; i < nz; i++)
	{
		zvals[i] = 0;
	}
	qval[0] = 0;

	dset.read(data_out, H5::PredType::NATIVE_FLOAT, memspace, dspace);
	dset_xvals.read(xvals, H5::PredType::NATIVE_FLOAT, memspace_xval, dspace_xvals);
	dset_yvals.read(yvals, H5::PredType::NATIVE_FLOAT, memspace_yval, dspace_yvals);
	dset_zvals.read(zvals, H5::PredType::NATIVE_FLOAT, memspace_zval, dspace_zvals);
	dset_q.read(qval, H5::PredType::NATIVE_FLOAT, memspace_q, dspace_q);

	fp.close();

	float *zvalsFull;
	int nzBeam(nz);
	if (MaxZ > zvals[nz - 1]) {
		double zDiv = zvals[1] - zvals[0];
		int NumberAdditions = (int)((MaxZ - zvals[nz - 1]) / zDiv);
		zvalsFull = new float[nz + NumberAdditions];
		for (int i = 0; i < nz; i++) {
			zvalsFull[i] = zvals[i];
		}
		for (int i = 0; i < NumberAdditions; i++) {
			zvalsFull[nz] = zvalsFull[nz - 1] + zDiv;
			nz++;
		}
	}
	else {
		zvalsFull = zvals;
	}

	x.resize(nrow*ncol*nz);
	y.resize(nrow*ncol*nz);
	z.resize(nrow*ncol*nz);
	px.resize(nrow*ncol*nz);
	py.resize(nrow*ncol*nz);
	pz.resize(nrow*ncol*nz);
	t.resize(nrow*ncol*nz);
	q.resize(nrow*ncol*nz);

	int V(0);
	int Extras(0);
	for (int i = 0; i < nrow; i++)
	{
		for (int j = 0; j < ncol; j++)
		{
			for (int k = 0; k < nz; k++)
			{
				z[V + Extras] = zvalsFull[k] * 100;
				x[V + Extras] = x0 + (xvals[j] * 100);
				y[V + Extras] = y0 + (yvals[i] * 100);
				q[V + Extras] = data_out[V] * qval[0];
				if (k >= nzBeam) {
					q[V + Extras] = 0;
					Extras++;
				}
				else V++;
			}
		}
	}
}


//Write to a HDF5 file
//Want to copy the initial file format with dielectric parameters
//Also want the forces - we'll put it in another set
//The mode composition will be the last dataset included
void HDF5Write(std::string filename, double x0, double y0, vector<double> xMacro, vector<double> yMacro, vector<double> zMacro, vector<double> px, vector<double> py, vector<double> pz, vector<double> t, vector<double> q, vector<double> xf, vector<double> yf, vector<double> zf, vector<double> Fx, vector<double> Fy, vector<double> Ez, double a, double delta, double w, double Epsilon, double StructureLength, int xModes, int yModes, vector<vector<double>> zeroesES, vector<vector<double>> zeroesEA, vector<vector<double>> zeroesHS, vector<vector<double>> zeroesHA, char Orientation){
    using namespace H5;
    //Have the filename (input) and where in the file our data is
    const H5std_string    FILE_NAME(filename);
    
    H5File file(FILE_NAME,H5F_ACC_TRUNC);
    //Group grp1(file.createGroup("Parameters"));
    //Create the group for each set of data to live in - sensibly named
    Group grp2(file.createGroup("beam"));
    Group grp3(file.createGroup("DielectricForces"));
    Group grp4(file.createGroup("DielectricParameters"));
    Group grp5(file.createGroup("Mode Amplitudes"));
    
    
    //Create the dataset for the macroparticle 6D distribution (+charge)
    double *data = new double[xMacro.size()*8];
    for(int i=0; i<xMacro.size()*8; i+=8) {
        data[i] = (xMacro[i/8]-x0)*0.01;
        data[i+1] = (yMacro[i/8]-y0)*0.01;
        data[i+2] = zMacro[i/8]*0.01;
        data[i+3] = px[i/8];
        data[i+4] = py[i/8];
        data[i+5] = pz[i/8];
        data[i+6] = t[i/8];
        //If the charge is outside the dielectric, set charge to zero
        if(Orientation == 'h' && (abs(yMacro[i/8])>a || abs(xMacro[i/8])>w)) data[i+7] = 0;
        else if(Orientation == 'v' && (abs(xMacro[i/8])>a || abs(yMacro[i/8])>w)) data[i+7] = 0;
        else data[i+7] = q[i/8];
    }
    
    hsize_t dimsf[2] = {xMacro.size(),8};
    DataSpace dataspace(2, dimsf);
    DataSet dataset = grp2.createDataSet("beam", PredType::NATIVE_DOUBLE, dataspace);
    dataset.write(data, PredType::NATIVE_DOUBLE);
    
    char const **DataLabel = new char const *[8];
    for (int i=0; i<8; i++) DataLabel[i]=new char[10];
    DataLabel[0]= "x";
    DataLabel[1]= "y";
    DataLabel[2]= "z";
    DataLabel[3]= "cpx";
    DataLabel[4]= "cpy";
    DataLabel[5]= "cpz";
    DataLabel[6]= "t";
    DataLabel[7]= "q";
    
    char const **UnitLabel = new char const *[8];
    for (int i=0; i<8; i++) UnitLabel[i]=new char[10];
    UnitLabel[0]= "m";
    UnitLabel[1]= "m";
    UnitLabel[2]= "m";
    UnitLabel[3]= "eV";
    UnitLabel[4]= "eV";
    UnitLabel[5]= "eV";
    UnitLabel[6]= "s";
    UnitLabel[7]= "s";
    
    H5::StrType datatype(H5::PredType::C_S1, H5T_VARIABLE);
    hsize_t dimsf3[2] = {8,1};
    DataSpace dataspace3(2, dimsf3);
    DataSet dataset3 = grp2.createDataSet("columns", datatype, dataspace3);
    dataset3.write(DataLabel, datatype);
    DataSet dataset4 = grp2.createDataSet("units", datatype, dataspace3);
    dataset4.write(UnitLabel, datatype);
    
    //Write force values at the mesh points (position in m and force in eV)
    double *ForceData = new double[Fx.size()*6];
    for (int i=0; i<Fx.size()*6; i+=6){
        ForceData[i] = (xf[i/6]-x0) * 0.01;
        ForceData[i + 1] = (yf[i/6]-y0) * 0.01;
        ForceData[i + 2] = zf[i/6] * 0.01;
        ForceData[i + 3] = Fx[i/6];
        ForceData[i + 4] = Fy[i/6];
        ForceData[i + 5] = Ez[i/6];
    }
    
    hsize_t dimsf2[2] = {Fx.size(),6};
    DataSpace dataspace2(2, dimsf2);
    DataSet dataset2 = grp3.createDataSet("ForceField", PredType::NATIVE_DOUBLE, dataspace2);
    dataset2.write(ForceData, PredType::NATIVE_DOUBLE);
    
    char const **ForceLabel = new char const *[6];
    for (int i=0; i<6; i++) ForceLabel[i]=new char[10];
    ForceLabel[0]= "x [m]";
    ForceLabel[1]= "y [m]";
    ForceLabel[2]= "z [m]";
    ForceLabel[3]= "Fx [eV]";
    ForceLabel[4]= "Fy [eV]";
    ForceLabel[5]= "Fz [eV]";
    
    hsize_t dimsf5[2] = {6,1};
    DataSpace dataspace5(2, dimsf5);
    DataSet dataset5 = grp3.createDataSet("columns", datatype, dataspace5);
    dataset5.write(ForceLabel, datatype);
    
    double *Parameters = new double[9];
    Parameters[1] = delta * 1e4;
	if (Orientation == 'v') {
		Parameters[2] = a * 1e4;
		Parameters[0] = w * 1e4;
	}
	else{
		Parameters[0] = a * 1e4;
		Parameters[2] = w * 1e4;
	}
	Parameters[3] = StructureLength;
    Parameters[4] = x0 * 0.01;
    Parameters[5] = y0 * 0.01;
    Parameters[6] = Epsilon;
    Parameters[7] = xModes;
    Parameters[8] = yModes;
    hsize_t dimsf6[2] = {9,1};
    DataSpace dataspace6(2, dimsf6);
    DataSet dataset6 = grp4.createDataSet("Parameters", PredType::NATIVE_DOUBLE, dataspace6);
    dataset6.write(Parameters, PredType::NATIVE_DOUBLE);
    
    char const **ParameterLabel = new char const *[9];
    for (int i=0; i<9; i++) ParameterLabel[i]=new char[25];
    ParameterLabel[1]= "delta [micron]";
	if (Orientation == 'v') {
		ParameterLabel[2] = "a [micron]";
		ParameterLabel[0] = "w [micron]";
	}
	else {
		ParameterLabel[0] = "a [micron]";
		ParameterLabel[2] = "w [micron]";
	}
    ParameterLabel[3]= "Structure Length [m]";
    ParameterLabel[4]= "x0 [m]";
    ParameterLabel[5]= "y0 [m]";
    ParameterLabel[6]= "Dielectric Permitivity";
    ParameterLabel[7]= "nX Modes";
    ParameterLabel[8]= "nY Modes";
    
    hsize_t dimsf7[2] = {9,1};
    DataSpace dataspace7(2, dimsf7);
    DataSet dataset7 = grp4.createDataSet("columns", datatype, dataspace7);
    dataset7.write(ParameterLabel, datatype);
    
    
    double *ModeValuesES = new double[zeroesES.size()*zeroesES[0].size()];
    for (int i=0; i<zeroesES.size(); i++){
        for (int j=0; j<zeroesES[i].size(); j++){
            ModeValuesES[i*zeroesES[0].size() + j] = zeroesES[i][j];
        }
    }
    double *ModeValuesEA = new double[zeroesEA.size()*zeroesEA[0].size()];
    for (int i=0; i<zeroesEA.size(); i++){
        for (int j=0; j<zeroesEA[i].size(); j++){
            ModeValuesEA[i*zeroesEA[0].size() + j] = zeroesEA[i][j];
        }
    }
    double *ModeValuesHS = new double[zeroesHS.size()*zeroesHS[0].size()];
    for (int i=0; i<zeroesHS.size(); i++){
        for (int j=0; j<zeroesHS[i].size(); j++){
            ModeValuesHS[i*zeroesHS[0].size() + j] = zeroesHS[i][j];
        }
    }
    double *ModeValuesHA = new double[zeroesHA.size()*zeroesHA[0].size()];
    for (int i=0; i<zeroesHA.size(); i++){
        for (int j=0; j<zeroesHA[i].size(); j++){
            ModeValuesHA[i*zeroesHA[0].size() + j] = zeroesHA[i][j];
        }
    }
    
    hsize_t dimsf8[2] = {zeroesES.size(),zeroesES[0].size()};
    DataSpace dataspace8(2,dimsf8);
    hsize_t dimsf9[2] = {zeroesEA.size(),zeroesEA[0].size()};
    DataSpace dataspace9(2,dimsf9);
    hsize_t dimsf10[2] = {zeroesHS.size(),zeroesHS[0].size()};
    DataSpace dataspace10(2,dimsf10);
    hsize_t dimsf11[2] = {zeroesHA.size(),zeroesHA[0].size()};
    DataSpace dataspace11(2,dimsf11);
    DataSet dataset8 = grp5.createDataSet("LSE Symmetric", PredType::NATIVE_DOUBLE, dataspace8);
    DataSet dataset9 = grp5.createDataSet("LSE Asymmetric", PredType::NATIVE_DOUBLE, dataspace9);
    DataSet dataset10 = grp5.createDataSet("LSM Symmetric", PredType::NATIVE_DOUBLE, dataspace10);
    DataSet dataset11 = grp5.createDataSet("LSM Asymmetric", PredType::NATIVE_DOUBLE, dataspace11);
    dataset8.write(ModeValuesES, PredType::NATIVE_DOUBLE);
    dataset9.write(ModeValuesEA, PredType::NATIVE_DOUBLE);
    dataset10.write(ModeValuesHS, PredType::NATIVE_DOUBLE);
    dataset11.write(ModeValuesHA, PredType::NATIVE_DOUBLE);
}
//Circular Case - Writing to a file
void HDF5Write_Circ(std::string filename, double x0, double y0, vector<double> xMacro, vector<double> yMacro, vector<double> zMacro, vector<double> px, vector<double> py, vector<double> pz, vector<double> t, vector<double> q, vector<double> xf, vector<double> yf, vector<double> zf, vector<double> Fx, vector<double> Fy, vector<double> Ez, double a, double delta, double Epsilon, double StructureLength, int rModes, int thetaModes, vector<vector<double>> ModeAmp, vector<vector<double>> WaveVec) {
	using namespace H5;
	//Have the filename (input) and where in the file our data is
	const H5std_string    FILE_NAME(filename);

	H5File file(FILE_NAME, H5F_ACC_TRUNC);
	//Group grp1(file.createGroup("Parameters"));
	//Create the group for each set of data to live in - sensibly named
	Group grp2(file.createGroup("beam"));
	Group grp3(file.createGroup("DielectricForces"));
	Group grp4(file.createGroup("DielectricParameters"));
	Group grp5(file.createGroup("Mode Amplitudes"));


	//Create the dataset for the macroparticle 6D distribution (+charge)
	double *data = new double[xMacro.size() * 8];
	for (int i = 0; i<xMacro.size() * 8; i += 8) {
		data[i] = (xMacro[i / 8] - x0)*0.01;
		data[i + 1] = (yMacro[i / 8] - y0)*0.01;
		data[i + 2] = zMacro[i / 8] * 0.01;
		data[i + 3] = px[i / 8];
		data[i + 4] = py[i / 8];
		data[i + 5] = pz[i / 8];
		data[i + 6] = t[i / 8];
		//If the charge is outside the dielectric, set charge to zero
		if ((abs(yMacro[i / 8])>a || abs(xMacro[i / 8])>a)) data[i + 7] = 0;
		else data[i + 7] = q[i / 8];
	}

	hsize_t dimsf[2] = { xMacro.size(),8 };
	DataSpace dataspace(2, dimsf);
	DataSet dataset = grp2.createDataSet("beam", PredType::NATIVE_DOUBLE, dataspace);
	dataset.write(data, PredType::NATIVE_DOUBLE);

	char const **DataLabel = new char const *[8];
	for (int i = 0; i<8; i++) DataLabel[i] = new char[10];
	DataLabel[0] = "x";
	DataLabel[1] = "y";
	DataLabel[2] = "z";
	DataLabel[3] = "cpx";
	DataLabel[4] = "cpy";
	DataLabel[5] = "cpz";
	DataLabel[6] = "t";
	DataLabel[7] = "q";

	char const **UnitLabel = new char const *[8];
	for (int i = 0; i<8; i++) UnitLabel[i] = new char[10];
	UnitLabel[0] = "m";
	UnitLabel[1] = "m";
	UnitLabel[2] = "m";
	UnitLabel[3] = "eV";
	UnitLabel[4] = "eV";
	UnitLabel[5] = "eV";
	UnitLabel[6] = "s";
	UnitLabel[7] = "s";

	H5::StrType datatype(H5::PredType::C_S1, H5T_VARIABLE);
	hsize_t dimsf3[2] = { 8,1 };
	DataSpace dataspace3(2, dimsf3);
	DataSet dataset3 = grp2.createDataSet("columns", datatype, dataspace3);
	dataset3.write(DataLabel, datatype);
	DataSet dataset4 = grp2.createDataSet("units", datatype, dataspace3);
	dataset4.write(UnitLabel, datatype);

	//Write force values at the mesh points (position in m and force in eV)
	double *ForceData = new double[Fx.size() * 6];
	for (int i = 0; i<Fx.size() * 6; i += 6) {
		ForceData[i] = (xf[i / 6] - x0) * 0.01;
		ForceData[i + 1] = (yf[i / 6] - y0) * 0.01;
		ForceData[i + 2] = zf[i / 6] * 0.01;
		ForceData[i + 3] = Fx[i / 6];
		ForceData[i + 4] = Fy[i / 6];
		ForceData[i + 5] = Ez[i / 6];
	}

	hsize_t dimsf2[2] = { Fx.size(),6 };
	DataSpace dataspace2(2, dimsf2);
	DataSet dataset2 = grp3.createDataSet("ForceField", PredType::NATIVE_DOUBLE, dataspace2);
	dataset2.write(ForceData, PredType::NATIVE_DOUBLE);

	char const **ForceLabel = new char const *[6];
	for (int i = 0; i<6; i++) ForceLabel[i] = new char[10];
	ForceLabel[0] = "x [m]";
	ForceLabel[1] = "y [m]";
	ForceLabel[2] = "z [m]";
	ForceLabel[3] = "Fx [eV]";
	ForceLabel[4] = "Fy [eV]";
	ForceLabel[5] = "Fz [eV]";

	hsize_t dimsf5[2] = { 6,1 };
	DataSpace dataspace5(2, dimsf5);
	DataSet dataset5 = grp3.createDataSet("columns", datatype, dataspace5);
	dataset5.write(ForceLabel, datatype);

	double *Parameters = new double[9];
	Parameters[0] = a * 1e4;
	Parameters[1] = delta * 1e4;
	Parameters[2] = 0;
	Parameters[3] = StructureLength;
	Parameters[4] = x0 * 0.01;
	Parameters[5] = y0 * 0.01;
	Parameters[6] = Epsilon;
	Parameters[7] = rModes;
	Parameters[8] = thetaModes;
	hsize_t dimsf6[2] = { 9,1 };
	DataSpace dataspace6(2, dimsf6);
	DataSet dataset6 = grp4.createDataSet("Parameters", PredType::NATIVE_DOUBLE, dataspace6);
	dataset6.write(Parameters, PredType::NATIVE_DOUBLE);

	char const **ParameterLabel = new char const *[9];
	for (int i = 0; i<9; i++) ParameterLabel[i] = new char[25];
	ParameterLabel[0] = "a [micron]";
	ParameterLabel[1] = "delta [micron]";
	ParameterLabel[2] = "w [micron] (N/A)";
	ParameterLabel[3] = "Structure Length [m]";
	ParameterLabel[4] = "x0 [m]";
	ParameterLabel[5] = "y0 [m]";
	ParameterLabel[6] = "Dielectric Permitivity";
	ParameterLabel[7] = "Radial Modes";
	ParameterLabel[8] = "Azimuthal Modes";

	hsize_t dimsf7[2] = { 9,1 };
	DataSpace dataspace7(2, dimsf7);
	DataSet dataset7 = grp4.createDataSet("columns", datatype, dataspace7);
	dataset7.write(ParameterLabel, datatype);


	double *ModeAmplitudes = new double[ModeAmp.size()*ModeAmp[0].size()];
	for (int i = 0; i<ModeAmp.size(); i++) {
		for (int j = 0; j<ModeAmp[i].size(); j++) {
			ModeAmplitudes[i*ModeAmp[0].size() + j] = ModeAmp[i][j];
		}
	}
	double *WaveVectors = new double[WaveVec.size()*WaveVec[0].size()];
	for (int i = 0; i<WaveVec.size(); i++) {
		for (int j = 0; j<WaveVec[i].size(); j++) {
			WaveVectors[i*WaveVec[0].size() + j] = WaveVec[i][j];
		}
	}

	hsize_t dimsf8[2] = { ModeAmp.size(),ModeAmp[0].size() };
	DataSpace dataspace8(2, dimsf8);
	hsize_t dimsf9[2] = { WaveVec.size(),WaveVec[0].size() };
	DataSpace dataspace9(2, dimsf9);
	DataSet dataset8 = grp5.createDataSet("Wave Amplitude", PredType::NATIVE_DOUBLE, dataspace8);
	DataSet dataset9 = grp5.createDataSet("Wave Vector", PredType::NATIVE_DOUBLE, dataspace9);
	dataset8.write(ModeAmplitudes, PredType::NATIVE_DOUBLE);
	dataset9.write(WaveVectors, PredType::NATIVE_DOUBLE);
}


void HDF5Write_PlanarGreen1D(std::string filename, vector<double> zValues, vector<double> Wx, vector<double> Wy, vector<double> Wz, double x0, double y0, double a, double delta, double w, double Epsilon, int xModes, int yModes) {
	using namespace H5;
	//Have the filename (input) and where in the file our data is
	const H5std_string    FILE_NAME(filename);

	H5File file(FILE_NAME, H5F_ACC_TRUNC);
	//Group grp1(file.createGroup("Parameters"));
	//Create the group for each set of data to live in - sensibly named
	Group grp3(file.createGroup("Wake Potential"));
	Group grp4(file.createGroup("Parameters"));
	
	//Write force values at the mesh points (position in m and force in eV)
	double *ForceData = new double[Wx.size() * 4];
	//1e8 is to convert from cgs to SI
	for (int i = 0; i<Wx.size() * 4; i += 4) {
		ForceData[i] = (zValues[i / 4]) * 0.01;
		ForceData[i + 1] = Wx[i / 4] * pow(2.99792458,-1)*1e2;
		ForceData[i + 2] = Wy[i / 4] * pow(2.99792458,-1)*1e2;
		ForceData[i + 3] = Wz[i / 4] * pow(2.99792458,-1)*1e2;
	}

	hsize_t dimsf2[2] = { Wx.size(),4 };
	DataSpace dataspace2(2, dimsf2);
	DataSet dataset2 = grp3.createDataSet("WakePotentials", PredType::NATIVE_DOUBLE, dataspace2);
	dataset2.write(ForceData, PredType::NATIVE_DOUBLE);

	char const **ForceLabel = new char const *[4];
	for (int i = 0; i<4; i++) ForceLabel[i] = new char[10];
	ForceLabel[0] = "z - z0 [m]";
	ForceLabel[1] = "Wx [V/m/C]";
	ForceLabel[2] = "Wy [V/m/C]";
	ForceLabel[3] = "Wz [V/m/C]";

	H5::StrType datatype(H5::PredType::C_S1, H5T_VARIABLE);
	hsize_t dimsf5[2] = { 4,1 };
	DataSpace dataspace5(2, dimsf5);
	DataSet dataset5 = grp3.createDataSet("columns", datatype, dataspace5);
	dataset5.write(ForceLabel, datatype);

	double *Parameters = new double[8];
	Parameters[0] = a * 1e4;
	Parameters[1] = delta * 1e4;
	Parameters[2] = w * 0.01;
	Parameters[3] = x0 * 0.01;
	Parameters[4] = y0 * 0.01;
	Parameters[5] = Epsilon;
	Parameters[6] = xModes;
	Parameters[7] = yModes;
	hsize_t dimsf6[2] = { 8,1 };
	DataSpace dataspace6(2, dimsf6);
	DataSet dataset6 = grp4.createDataSet("Parameters", PredType::NATIVE_DOUBLE, dataspace6);
	dataset6.write(Parameters, PredType::NATIVE_DOUBLE);

	char const **ParameterLabel = new char const *[8];
	for (int i = 0; i<8; i++) ParameterLabel[i] = new char[25];
	ParameterLabel[0] = "a [micron]";
	ParameterLabel[1] = "delta [micron]";
	ParameterLabel[2] = "w [m]";
	ParameterLabel[3] = "x0 [m]";
	ParameterLabel[4] = "y0 [m]";
	ParameterLabel[5] = "Dielectric Permitivity";
	ParameterLabel[6] = "nX Modes";
	ParameterLabel[7] = "nY Modes";

	hsize_t dimsf7[2] = { 8,1 };
	DataSpace dataspace7(2, dimsf7);
	DataSet dataset7 = grp4.createDataSet("columns", datatype, dataspace7);
	dataset7.write(ParameterLabel, datatype);
}


// ----------------------------------------------------------------- //
// ---------------- Kick the Beam Using the Field ------------------ //
// ----------------------------------------------------------------- //

void InterpolationKicker(vector<double>& xB,vector<double>& yB,vector<double>& zB,vector<double>& pxB,vector<double>& pyB,vector<double>& pzB, vector<double> qB, vector<double> x,vector<double> y,vector<double> z,vector<double> Fz,vector<double> Fy,vector<double> Fx, double l){
    //Here we have calculated Fi at (x,y,z) and we want to use interpolation to calculate the force at (xB,yB,zB)
    //Use that to kick the momentum of the beam after passing through structure of length l
    //Also kick the beam based on the original momenta
    int nParticles = (int) xB.size();
    std::cout<<"Interpolating Macroparticle Kicks:"<<std::endl;
    auto InterpolationStartTime = std::chrono::high_resolution_clock::now();
    
    for(int i=0; i<nParticles; i++)
    {
        /*
        int closestPosition = closest3D(xB[i],yB[i],zB[i],x,y,z);
        double FyClosest = Fy[closestPosition];
        double FxClosest = Fx[closestPosition];
        double FzClosest = Fz[closestPosition];
        */
        vector<double> Forces = TriLinearFormal(Fx, Fy, Fz, x, y, z, xB[i], yB[i], zB[i]);
        double FxClosest = Forces[0];
        double FyClosest = Forces[1];
        double FzClosest = Forces[2];
        
        //To kick the beam from initial momentum use p = mv to get v_i and then x += v * l
        xB[i] += l * (pxB[i] / pzB[i]);
        yB[i] += l * (pyB[i] / pzB[i]);
        
        pyB[i] += FyClosest*l;
        pxB[i] += FxClosest*l;
        pzB[i] += FzClosest*l;
        if (i % 5000 == 0) std::cout<<(double) 100*i/nParticles<<"%"<<std::endl;
    }
    auto InterpolationEndTime = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> InterpolationTime = InterpolationEndTime-InterpolationStartTime;
    std::cout<< "Time for Macroparticle Kicks:" <<InterpolationTime.count() <<"s"<<std::endl;
    
}

// ------------------------------------------------------------------------------------- //
// ---------------------------- Complete Field Calculator ------------------------------ //
//-------------------------------------------------------------------------------------- //


void DiWaCATField_Planar(vector <std::string> InputList, std::string FileIn, std::string FileOut)
{

	/*----------------------------------------------------------*/
	/*----------------Set up the simulation space---------------*/
	/*----------------------------------------------------------*/
	//Whether the beam is from "HDF5" file of defined with an analytical "Function"
	std::cout << "Entered DiWaCAT" << std::endl;

	//The file path to the meshed HDF5 file (n.b. mesh points must be in a folder MeshValues)
	std::string HDF5BeamIn = FileIn;

	//File path for the kicked beam (positions changed using initial momenta)
	std::string HDF5BeamOut = FileOut;

	char PlateOrientation = InputList[0][0]; //Plate orientation: either 'h' or 'v'
	double MaxZ = std::stold(InputList[14]); //Maximum value of longitudinal coordinate in units of sigmat
									  //Timing info for root finder
	auto start = std::chrono::high_resolution_clock::now();

	//---------------------- READ IN THE BEAM AND STRUCTURE PARAMETERS --------------//
	
	double xi = 0;    //Position of force calculation in x if xdiv = 1
	double yi = 0;   //Position of force calculation in y if ydiv = 1
	double zi = 0;   //Position of force calculation in zeta if zetadiv = 1
	double x0 = std::stold(InputList[1]);    //Beam position x
	double y0 = std::stold(InputList[2]);    //Beam position y
	double z0 = 0;
	double Ep = std::stold(InputList[3]);    //Relative Permittivity of dielectric
	double Mu = std::stold(InputList[4]);   // Relative Permeability of dielectric
	double b = std::stold(InputList[5]);   //Distance from cavity center to dielectric
	double c = std::stold(InputList[5]) + std::stold(InputList[6]);   //Distance from cavity center to conducting wall (-c < y < c)
	double w = std::stold(InputList[7]);   //Width of Cavity (0 < x < w)
	int sN = std::stoi(InputList[8]);     //Number of frequcy modes in x
	int sI = std::stoi(InputList[9]);     //Number of frequency modes in y for each mode in x
	double acc = std::stold(InputList[10]); //Percision of root finder
	double ModeAccuracy = std::stold(InputList[11]) / 100;
	double StructureLength = std::stold(InputList[12]);
	int NSteps = std::stold(InputList[13]);
	if (NSteps == 0) {
		std::cout << "Number of steps set to zero. No kicks will be applied." << std::endl;
	}
	//CHANGE TO FALSE IF USING USER INPUT MODE NUMBERS
	bool ConvergenceCalculate = true;
	if (InputList[15] == "f") {
		ConvergenceCalculate = false;
	}

	std::string fileIn = HDF5BeamIn;
	std::string fileOut = HDF5BeamOut;

	//Need to add half the width to the initial position
	x0 += w / 2;
	if (PlateOrientation == 'v') {
		//Swap initial positions for vertical setup - y0 is always towards the dielectric
		y0 = x0;
		x0 = std::stold(InputList[2]);
	}


	/*--------------Read in the macro file----------------*/

	//Set up the vectors
	vector<double> x;
	vector<double> y;
	vector<double> z;
	vector<double> px;
	vector<double> py;
	vector<double> pz;
	vector<double> t;
	vector<double> q;

	//Read from the HDF5 file and write the data to the vectors
	HDF5ReaderMacro(fileIn, x, y, z, px, py, pz, t, q, y0, x0, MaxZ);
	//If Max(Z) < MaxZ --> need to add mesh points, i.e. find the repeated x and y and add more z points

	std::cout << "Read HDF5 File" << std::endl;


	/*--------------Read in the full beam file----------------*/

	//Set up the vectors
	vector<double> xB;
	vector<double> yB;
	vector<double> zB;
	vector<double> pxB;
	vector<double> pyB;
	vector<double> pzB;
	vector<double> tB;
	vector<double> qB;

	//Read from the HDF5 file and write the data to the vectors
	HDF5Reader(fileIn, xB, yB, zB, pxB, pyB, pzB, tB, qB, y0, x0);

	//Calculate the average momentum to calculate beta
	double E(0);
	for (int i = 0; i<xB.size(); i++) {
		E += pzB[i] / (xB.size() * 1e6);
	}
	std::cout << "Average momentum: " << E << " MeV/c" << std::endl;
	//Calculate a z mesh interval (we'll use it later to make sure modes are right)
	double zMeshInterval(0);
	int zLoop(0);
	while (zMeshInterval <= 0) {
		zMeshInterval = zB[zLoop] - zB[0];
		zLoop++;
	}
	double yMax = yB[0];
	double xMax = xB[0];
	for (int i = 0; i<xB.size(); i++) {
		if (yB[i]>yMax) yMax = yB[i];
		if (xB[i]>xMax) xMax = xB[i];
	}

	double B = sqrt(1 - pow(E / 0.510999, -2)); // Use the energy to calculate beta
	double B2 = B*B;

	double zetamin = *std::min_element(z.begin(), z.end());
	double zetamax = *std::max_element(z.begin(), z.end());

	//------------------ SET UP 3D FORCE MATRICES  -------------------//

	int NParticle = (int)x.size();

	//Create and size vectors to hold forces.
	//Need to make all values 0 too
	vector <double> Fz(NParticle);
	vector <double> Fx(NParticle);
	vector <double> Fy(NParticle);
	std::fill(Fz.begin(), Fz.end(), 0);
	std::fill(Fx.begin(), Fx.end(), 0);
	std::fill(Fy.begin(), Fy.end(), 0);

	//Set up smaller subset of macroparticles (9*9*30)

	//--------------- CALCULATE THE EIGENVALUES OF FREQUENCY MODES ------------//
	//Create and resize vectors to hold eigenvalues of allowed frequency modes
	vector<vector<double> > zeroesES;
	vector<vector<double> > zeroesEA;
	vector<vector<double> > zeroesHS;
	vector<vector<double> > zeroesHA;

	//Calculate the modes -- if last mode != 1% of the total, add more modes
	bool yModesAccurate = false;
	bool xModesAccurate = false;
	if (ConvergenceCalculate == true) {
		while (yModesAccurate == false || xModesAccurate == false)
		{
			EigenvaluesCalculator(sI, sN, zeroesES, zeroesEA, zeroesHS, zeroesHA, Ep, Mu, B2, b, c, w, acc);

			//Calculate a beam centre Green's function
			//If the function with one fewer mode is less than 1% lower we have enough modes
			double dFzFullModes(0), dFyFullModes(0), dFxFullModes(0);
			double dFzFewerModesX(0), dFyFewerModesX(0), dFxFewerModesX(0);
			double dFzFewerModesY(0), dFyFewerModesY(0), dFxFewerModesY(0);
			if (PlateOrientation == 'h') {
				//If beam is beyond the structure edge, set the checking point as 5micron from edge
				if (yMax > b) yMax = b - 0.0005;
				CalcWakeElement(dFzFullModes, dFxFullModes, dFyFullModes, xMax, x0, yMax, y0, zMeshInterval, 0, zeroesES, zeroesEA, zeroesHS, zeroesHA, sN, sI, Ep, Mu, B2, b, c, w);
				if (isnan(dFzFullModes)) {
					EigenvaluesCalculator(sI, sN, zeroesES, zeroesEA, zeroesHS, zeroesHA, Ep, Mu, B2, b, c, w, acc);
					CalcWakeElement(dFzFullModes, dFxFullModes, dFyFullModes, xMax, x0, yMax, y0, zMeshInterval, 0, zeroesES, zeroesEA, zeroesHS, zeroesHA, sN, sI, Ep, Mu, B2, b, c, w);
				}
				CalcWakeElement(dFzFewerModesX, dFxFewerModesX, dFyFewerModesX, xMax, x0, yMax, y0, zMeshInterval, 0, zeroesES, zeroesEA, zeroesHS, zeroesHA, sN - 5, sI, Ep, Mu, B2, b, c, w);
				CalcWakeElement(dFzFewerModesY, dFxFewerModesY, dFyFewerModesY, xMax, x0, yMax, y0, zMeshInterval, 0, zeroesES, zeroesEA, zeroesHS, zeroesHA, sN, sI - 5, Ep, Mu, B2, b, c, w);
			}
			if (PlateOrientation == 'v') {
				//If beam is beyond the structure edge, set the checking point as 5micron from edge
				if (xMax > b) xMax = b - 0.0005;
				CalcWakeElement(dFzFullModes, dFxFullModes, dFyFullModes, yMax, y0, xMax, x0, zMeshInterval, 0, zeroesES, zeroesEA, zeroesHS, zeroesHA, sN, sI, Ep, Mu, B2, b, c, w);
				if (isnan(dFzFullModes)) {
					EigenvaluesCalculator(sI, sN, zeroesES, zeroesEA, zeroesHS, zeroesHA, Ep, Mu, B2, b, c, w, acc);
					CalcWakeElement(dFzFullModes, dFxFullModes, dFyFullModes, yMax, y0, xMax, x0, zMeshInterval, 0, zeroesES, zeroesEA, zeroesHS, zeroesHA, sN, sI, Ep, Mu, B2, b, c, w);
				}
				CalcWakeElement(dFzFewerModesX, dFxFewerModesX, dFyFewerModesX, yMax, y0, xMax, x0, zMeshInterval, 0, zeroesES, zeroesEA, zeroesHS, zeroesHA, sN - 5, sI, Ep, Mu, B2, b, c, w);
				CalcWakeElement(dFzFewerModesY, dFxFewerModesY, dFyFewerModesY, yMax, y0, xMax, x0, zMeshInterval, 0, zeroesES, zeroesEA, zeroesHS, zeroesHA, sN, sI - 5, Ep, Mu, B2, b, c, w);
			}

			if (((dFzFullModes - dFzFewerModesX) / dFzFullModes < ModeAccuracy) && ((dFyFullModes - dFyFewerModesX) / dFyFullModes < ModeAccuracy) && ((dFxFullModes - dFxFewerModesX) / dFxFullModes < ModeAccuracy)) {
				xModesAccurate = true;
			}
			else sN += 2;

			if (((dFzFullModes - dFzFewerModesY) / dFzFullModes < ModeAccuracy) && ((dFyFullModes - dFyFewerModesY) / dFyFullModes < ModeAccuracy) && ((dFxFullModes - dFxFewerModesY) / dFxFullModes < ModeAccuracy)) {
				yModesAccurate = true;
			}
			else {
				sI += 2;
			}
			/*
			if(sI>sN){
			sN++;
			sI=20;
			}
			*/
		}
	}
	EigenvaluesCalculator(sI, sN, zeroesES, zeroesEA, zeroesHS, zeroesHA, Ep, Mu, B2, b, c, w, acc);
	std::cout << "Modes used: \n Nx = " << sN << "\n Ny = " << sI << std::endl;
	//---------------------    PRE-RUN ADMIN    -----------------------//
	//Timing info for root finder
	auto finish = std::chrono::high_resolution_clock::now();

	std::chrono::duration<double> elapsed = finish - start;
	std::cout << "Elapsed time for Root Finder: " << elapsed.count() << " s\n";

	//Timing info for force calculator
	auto start2 = std::chrono::high_resolution_clock::now();


	/*----------------------------------------------------------*/
	/*-----------------Run code on meshed macros----------------*/
	/*----------------------------------------------------------*/


	//Run field calculator at the macroparticle positions
	//---------------         CODE RUNNING       -------------------//

	//Calculate forces at each desired point in mesh
	if (PlateOrientation == 'h') {
		std::cout << "Horizontal Plate" << std::endl;
		TotalForceMeshHDF5(Fz, Fx, Fy, NParticle, xi, x0, yi, y0, zi, zetamin, zetamax, z0, x, y, z, zeroesES, zeroesEA, zeroesHS, zeroesHA, sN, sI, Ep, Mu, B2, b, c, w, q);
	}

	if (PlateOrientation == 'v') {
		std::cout << "Vertical Plate" << std::endl;
		TotalForceMeshHDF5VerticalPlate(Fz, Fx, Fy, NParticle, xi, x0, yi, y0, zi, zetamin, zetamax, z0, x, y, z, zeroesES, zeroesEA, zeroesHS, zeroesHA, sN, sI, Ep, Mu, B2, b, c, w, q);
	}


	/*----------------------------------------------------------*/
	/*----------Use force on mesh on full beam sample-----------*/
	/*----------------------------------------------------------*/
	//Apply the kick to the big sample of particles
	//Need to apply the kick to the positions from the original beam not the smaller mesh
	//This is the most simple kick possible - i.e. the one time interval


	/* We need to:
	* read in the beam positions from the HDF5 file (got a function to do that
	* get the approximate force at each position using the force calculated from the mesh
	* i.e interpolation
	* change the phase space of the particles using this force
	* write the beam, along with it's new phase space info to a HDF5

	*/
	if (StructureLength>0 && NSteps>0) {
		for (int i = 0; i<NSteps; i++) {
			std::cout << "Kick Step: " << i + 1 << " of " << NSteps << std::endl;
			InterpolationKicker(xB, yB, zB, pxB, pyB, pzB, qB, x, y, z, Fz, Fy, Fx, StructureLength / NSteps);
		}
	}
	else std::cout << "DLW Length = 0. No Kicks applied" << std::endl;


	/*----------------------------------------------------------*/
	/*----------------------Output results----------------------*/
	/*----------------------------------------------------------*/
	//Read the simulated beam back to a HDF5 file
	//Add another path to have the field strengths in so we can plot that if wanted

	HDF5Write(fileOut, x0, y0, xB, yB, zB, pxB, pyB, pzB, tB, qB, x, y, z, Fx, Fy, Fz, b, c - b, w, Ep, StructureLength, sN, sI, zeroesES, zeroesEA, zeroesHS, zeroesHA, PlateOrientation);

	//Timing info for force calculator
	auto finish2 = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed2 = finish2 - start2;
	std::cout << "Elapsed time for Simulation: " << elapsed2.count() << " s\n";
}


void DiWaCATField_Circ(vector <std::string> InputList, std::string FileIn, std::string FileOut) {
	/*----------------------------------------------------------*/
	/*----------------Set up the simulation space---------------*/
	/*----------------------------------------------------------*/
	//Whether the beam is from "HDF5" file of defined with an analytical "Function"
	std::cout << "Entered DiWaCAT" << std::endl;

	//The file path to the meshed HDF5 file (n.b. mesh points must be in a folder MeshValues)

	std::string HDF5BeamIn = FileIn;

	
	//File path for the kicked beam (positions changed using initial momenta)
	std::string HDF5BeamOut = FileOut;

	double MaxZ = std::stold(InputList[12]); //Maximum value of longitudinal coordinate in units of sigmat
											 //Timing info for root finder
	auto start = std::chrono::high_resolution_clock::now();

	//---------------------- READ IN THE BEAM AND STRUCTURE PARAMETERS --------------//

	double xi = 0;    //Position of force calculation in x if xdiv = 1
	double yi = 0;   //Position of force calculation in y if ydiv = 1
	double zi = 0;   //Position of force calculation in zeta if zetadiv = 1
	double x0 = std::stold(InputList[0]);    //Beam position x
	double y0 = std::stold(InputList[1]);    //Beam position y
	double z0 = 0;
	double epsilon = std::stold(InputList[2]);    //Relative Permittivity of dielectric
	double mu = std::stold(InputList[3]);   // Relative Permeability of dielectric
	double b = std::stold(InputList[4]);   //Distance from cavity center to dielectric
	double c = (std::stold(InputList[4]) + std::stold(InputList[5]));   //Distance from cavity center to conducting wall (-c < y < c)
	int nR = std::stoi(InputList[6]);     //Number of frequcy modes in r
	int nT = std::stoi(InputList[7]);     //Number of frequency modes in theta
	double acc = std::stold(InputList[8]); //Percision of root finder
	double ModeAccuracy = std::stold(InputList[9]) / 100;
	double StructureLength = std::stold(InputList[10]);
	int NSteps = std::stold(InputList[11]);
	if (NSteps == 0) {
		std::cout << "Number of steps set to zero. No kicks will be applied." << std::endl;
	}
	//CHANGE TO FALSE IF USING USER INPUT MODE NUMBERS
	bool ConvergenceCalculate = true;
	if (InputList[13] == "f") {
		ConvergenceCalculate = false;
	}

	std::string fileIn = HDF5BeamIn;
	std::string fileOut = HDF5BeamOut;

	/*--------------Read in the macro file----------------*/

	//Set up the vectors
	vector<double> x;
	vector<double> y;
	vector<double> z;
	vector<double> px;
	vector<double> py;
	vector<double> pz;
	vector<double> t;
	vector<double> q;

	//Read from the HDF5 file and write the data to the vectors
	HDF5ReaderMacro(fileIn, x, y, z, px, py, pz, t, q, y0, x0, MaxZ);
	//Convert the transverse positions to radial
	vector<double> r(x.size());
	vector<double> theta(x.size());
	CartMeshToCyl(r, theta, x, y);

	std::cout << "Read HDF5 File" << std::endl;

	/*--------------Read in the full beam file----------------*/

	//Set up the vectors
	vector<double> xB;
	vector<double> yB;
	vector<double> zB;
	vector<double> pxB;
	vector<double> pyB;
	vector<double> pzB;
	vector<double> tB;
	vector<double> qB;

	//Read from the HDF5 file and write the data to the vectors
	HDF5Reader(fileIn, xB, yB, zB, pxB, pyB, pzB, tB, qB, y0, x0);

	//Calculate the average momentum to calculate beta
	double E(0);
	for (int i = 0; i<xB.size(); i++) {
		E += pzB[i] / (xB.size() * 1e6);
	}
	std::cout << "Average momentum: " << E << " MeV/c" << std::endl;

	double rMax(0);
	for (int i = 0; i<r.size(); i++) {
		if (r[i]>rMax) rMax = r[i];
	}
	if (rMax > b) {
		rMax = 0.9*b;
	}

	double zetamin = *std::min_element(z.begin(), z.end());
	double zetamax = *std::max_element(z.begin(), z.end());

	//------------------ SET UP 3D FORCE MATRICES  -------------------//

	int NParticle = (int)x.size();

	//Create and size vectors to hold forces.
	//Need to make all values 0 too
	vector <double> Fz(NParticle);
	vector <double> Fx(NParticle);
	vector <double> Fy(NParticle);
	vector <double> Fr(NParticle);
	vector <double> Ftheta(NParticle);
	std::fill(Fz.begin(), Fz.end(), 0);
	std::fill(Fx.begin(), Fx.end(), 0);
	std::fill(Fy.begin(), Fy.end(), 0);
	std::fill(Fr.begin(), Fr.end(), 0);
	std::fill(Ftheta.begin(), Ftheta.end(), 0);

	//----------------- Calculate Mode Amplitude --------------------//
	vector<vector<double>> WaveVectors;
	vector<vector<double>> ModeAmplitude;
	//Find the modes
	//Find the number of modes needed for convergence if set to true
	std::cout << "Calculating Mode Frequencies and Amplitudes" << std::endl;
	if (ConvergenceCalculate == true) {
		ModeConvergence(ModeAmplitude, WaveVectors, rMax*0.01, b*0.01, c*0.01, mu, epsilon, acc, ModeAccuracy);
	}
	else FindModes(WaveVectors, ModeAmplitude, b*0.01, c*0.01, mu, epsilon, nT, nR, acc,true);

	
	double fx(0), fy(0), ez(0);

	//---------------------    PRE-RUN ADMIN    -----------------------//
	//Timing info for root finder
	auto finish = std::chrono::high_resolution_clock::now();

	std::chrono::duration<double> elapsed = finish - start;
	std::cout << "Elapsed time for Root Finder: " << elapsed.count() << " s\n";


	//------------------- Calculate Forces ---------------------- //
	TotalForceMeshCircular(Fz, Fr, Ftheta, r, theta, z, q, ModeAmplitude, WaveVectors, b*0.01, c*0.01, mu, epsilon);
	//--------------Convert Cylindrical Forces to Cartesian ------------//
	CylForceToCart(Fx, Fy, Fr, Ftheta, theta);

	if (StructureLength>0 && NSteps>0) {
		for (int i = 0; i<NSteps; i++) {
			std::cout << "Kick Step: " << i + 1 << " of " << NSteps << std::endl;
			InterpolationKicker(xB, yB, zB, pxB, pyB, pzB, qB, x, y, z, Fz, Fy, Fx, StructureLength / NSteps);
		}
	}
	else std::cout << "DLW Length = 0. No Kicks applied" << std::endl;


	/*----------------------------------------------------------*/
	/*----------------------Output results----------------------*/
	/*----------------------------------------------------------*/
	//Read the simulated beam back to a HDF5 file
	//Add another path to have the field strengths in so we can plot that if wanted

	HDF5Write_Circ(fileOut, x0, y0, xB, yB, zB, pxB, pyB, pzB, tB, qB, x, y, z, Fx, Fy, Fz, b, c - b, epsilon, StructureLength, nR, nT, ModeAmplitude, WaveVectors);

	//Timing info for force calculator
	auto finish2 = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed2 = finish2 - start;
	std::cout << "Elapsed time for Simulation: " << elapsed2.count() << " s\n";

}

void DiWaCATField_1DGreen_Planar(vector <std::string> InputList,std::string FileOut) {
	//File path for the kicked beam (positions changed using initial momenta)
	std::string HDF5BeamOut = FileOut;

	double MaxZ = std::stold(InputList[12]); //Maximum value of longitudinal coordinate in units of sigmat
											 //Timing info for root finder

	auto start = std::chrono::high_resolution_clock::now();

	//---------------------- READ IN THE BEAM AND STRUCTURE PARAMETERS --------------//
	double x0 = std::stold(InputList[0]);    //Beam position x
	double y0 = std::stold(InputList[1]);    //Beam position y
	double z0 = 0;
	double Ep = std::stold(InputList[2]);    //Relative Permittivity of dielectric
	double Mu = std::stold(InputList[3]);   // Relative Permeability of dielectric
	double b = std::stold(InputList[4]);   //Distance from cavity center to dielectric
	double c = std::stold(InputList[4]) + std::stold(InputList[5]);   //Distance from cavity center to conducting wall (-c < y < c)
	double w = std::stold(InputList[6]);   //Width of Cavity (0 < x < w)
	int sN = std::stoi(InputList[7]);     //Number of frequcy modes in x
	int sI = std::stoi(InputList[8]);     //Number of frequency modes in y for each mode in x
	double acc = std::stold(InputList[9]); //Percision of root finder
	double ModeAccuracy = std::stold(InputList[10]) / 100;
	int NPoints = std::stold(InputList[11]);
	//CHANGE TO FALSE IF USING USER INPUT MODE NUMBERS
	bool ConvergenceCalculate = true;
	if (InputList[13] == "f") {
		ConvergenceCalculate = false;
	}
	std::string fileOut = HDF5BeamOut;
	
	//Need to add half the width to the initial position
	x0 += w / 2;

	//Set up the vector of Z positions
	vector<double> zValues(NPoints);
	vector<double> Wx(NPoints), Wy(NPoints), Wz(NPoints);
	double zInterval = MaxZ / double(NPoints);
	
	std::generate(zValues.begin(), zValues.end(), [n = 0, &zInterval]() mutable { return n++ * zInterval; });
	
	//Set momentum as constant - high enough to ensure relativistic
	double E(250);
	double B = sqrt(1 - pow(E / 0.510999, -2)); // Use the energy to calculate beta
	double B2 = B*B;

	//--------------- CALCULATE THE EIGENVALUES OF FREQUENCY MODES ------------//
	//Create and resize vectors to hold eigenvalues of allowed frequency modes
	vector<vector<double> > zeroesES;
	vector<vector<double> > zeroesEA;
	vector<vector<double> > zeroesHS;
	vector<vector<double> > zeroesHA;

	//Calculate the modes -- if last mode != 1% of the total, add more modes
	bool yModesAccurate = false;
	bool xModesAccurate = false;
	if (ConvergenceCalculate == true) {
		while (yModesAccurate == false || xModesAccurate == false)
		{
			EigenvaluesCalculator(sI, sN, zeroesES, zeroesEA, zeroesHS, zeroesHA, Ep, Mu, B2, b, c, w, acc);
			//Calculate a beam centre Green's function
			//If the function with one fewer mode is less than 1% lower we have enough modes
			double dFzFullModes(0), dFyFullModes(0), dFxFullModes(0);
			double dFzFewerModesX(0), dFyFewerModesX(0), dFxFewerModesX(0);
			double dFzFewerModesY(0), dFyFewerModesY(0), dFxFewerModesY(0);
			//If beam is beyond the structure edge, set the checking point as 5micron from edge
			CalcWakeElement(dFzFullModes, dFxFullModes, dFyFullModes, x0, x0, y0, y0, zInterval, 0, zeroesES, zeroesEA, zeroesHS, zeroesHA, sN, sI, Ep, Mu, B2, b, c, w);
			if (isnan(dFzFullModes)) {
				EigenvaluesCalculator(sI, sN, zeroesES, zeroesEA, zeroesHS, zeroesHA, Ep, Mu, B2, b, c, w, acc);
				CalcWakeElement(dFzFullModes, dFxFullModes, dFyFullModes, x0, x0, y0, y0, zInterval, 0, zeroesES, zeroesEA, zeroesHS, zeroesHA, sN, sI, Ep, Mu, B2, b, c, w);
			}
			CalcWakeElement(dFzFewerModesX, dFxFewerModesX, dFyFewerModesX, x0, x0, y0, y0, zInterval, 0, zeroesES, zeroesEA, zeroesHS, zeroesHA, sN - 5, sI, Ep, Mu, B2, b, c, w);
			CalcWakeElement(dFzFewerModesY, dFxFewerModesY, dFyFewerModesY, x0, x0, y0, y0, zInterval, 0, zeroesES, zeroesEA, zeroesHS, zeroesHA, sN, sI - 5, Ep, Mu, B2, b, c, w);
			if (((dFzFullModes - dFzFewerModesX) / dFzFullModes < ModeAccuracy) && ((dFyFullModes - dFyFewerModesX) / dFyFullModes < ModeAccuracy) && ((dFxFullModes - dFxFewerModesX) / dFxFullModes < ModeAccuracy)) {
				xModesAccurate = true;
			}
			else sN += 2;
			if (((dFzFullModes - dFzFewerModesY) / dFzFullModes < ModeAccuracy) && ((dFyFullModes - dFyFewerModesY) / dFyFullModes < ModeAccuracy) && ((dFxFullModes - dFxFewerModesY) / dFxFullModes < ModeAccuracy)) {
				yModesAccurate = true;
			}
			else {
				sI += 2;
			}
		}
	}
	EigenvaluesCalculator(sI, sN, zeroesES, zeroesEA, zeroesHS, zeroesHA, Ep, Mu, B2, b, c, w, acc);
	std::cout << "Modes used: \n Nx = " << sN << "\n Ny = " << sI << std::endl;
	//---------------------    PRE-RUN ADMIN    -----------------------//
	//Timing info for root finder
	auto finish = std::chrono::high_resolution_clock::now();

	std::chrono::duration<double> elapsed = finish - start;
	std::cout << "Elapsed time for Root Finder: " << elapsed.count() << " s\n";

	//Timing info for force calculator
	auto start2 = std::chrono::high_resolution_clock::now();
	//Conversion from cgs to SI
	for (int i = 0; i < NPoints; i++) {
		CalcWakeElement(Wz[i], Wx[i], Wy[i], x0, x0, y0, y0, zValues[i], 0, zeroesES, zeroesEA, zeroesHS, zeroesHA, sN, sI, Ep, Mu, B2, b, c, w);
	}

	/*----------------------------------------------------------*/
	/*----------------------Output results----------------------*/
	/*----------------------------------------------------------*/
	HDF5Write_PlanarGreen1D(fileOut, zValues,Wx,Wy,Wz,x0, y0, b, c - b, w, Ep, sN, sI);

	//Timing info for force calculator
	auto finish2 = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed2 = finish2 - start2;
	std::cout << "Elapsed time for Simulation: " << elapsed2.count() << " s\n";


}

#define RoFSuP_h

#endif /* RoFSuP_h */

