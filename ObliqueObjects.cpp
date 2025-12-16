// Oblique_Objects_Modern.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <algorithm>

using namespace std;

class ObliqueObjects {
public:
    // --- Constructor & Destructor ---
    explicit ObliqueObjects(int nObjects);
    ~ObliqueObjects() = default;

    // --- Public Methods ---
    void ReadObjectInfoAndFindMinMax(
        int nn, double x0, double y0, double z0,
        double luu, double lvv, double lww,
        int i_dir,
        const array<double, 3>& MyU_axis,
        const array<double, 3>& MyV_axis,
        const array<double, 3>& MyW_axis,
        double& xmin, double& xmax,
        double& ymin, double& ymax,
        double& zmin, double& zmax);

    void findCornersObliqueRectangular(
        int nn, double x0, double y0, double z0,
        double lu, double lv, double lw,
        array<double, 3>& xyz1,
        array<double, 3>& xyz2,
        array<double, 3>& xyz3,
        array<double, 3>& xyz4,
        array<double, 3>& xyz5,
        array<double, 3>& xyz6,
        array<double, 3>& xyz7);

    void xyzMinMaxRectangular(
        int nn,
        const array<double, 3>& xyz1,
        const array<double, 3>& xyz2,
        const array<double, 3>& xyz3,
        const array<double, 3>& xyz4,
        const array<double, 3>& xyz5,
        const array<double, 3>& xyz6,
        const array<double, 3>& xyz7,
        double& xmin, double& xmax,
        double& ymin, double& ymax,
        double& zmin, double& zmax);

    void checkIfPointInObliqueRectangular(
        double x, double y, double z,
        double x0, double y0, double z0,
        double u_length, double v_length, double w_length,
        bool& isit_in, int nn);

    // --- Data ---
    vector<array<double, 3>> u_axis, v_axis, w_axis;
    vector<array<double, 3>> xyz0Corner;
    vector<double> lu, lv, lw;
    vector<vector<vector<int>>> Object_range;

private:
    // Utility
    static double minOfArray(const array<double, 8>& x);
    static double maxOfArray(const array<double, 8>& x);
};

// ---------------------------------------------------------------
// Constructor
// ---------------------------------------------------------------
ObliqueObjects::ObliqueObjects(int nObjects)
{
    u_axis.assign(nObjects, { 0.0, 0.0, 0.0 });
    v_axis.assign(nObjects, { 0.0, 0.0, 0.0 });
    w_axis.assign(nObjects, { 0.0, 0.0, 0.0 });
    xyz0Corner.assign(nObjects, { 0.0, 0.0, 0.0 });

    lu.assign(nObjects, 0.0);
    lv.assign(nObjects, 0.0);
    lw.assign(nObjects, 0.0);

    Object_range.assign(nObjects, vector<vector<int>>(3, vector<int>(6, 0)));
}

// ---------------------------------------------------------------
// Read Object Info and Find Min/Max
// ---------------------------------------------------------------
void ObliqueObjects::ReadObjectInfoAndFindMinMax(
    int nn, double x0, double y0, double z0,
    double luu, double lvv, double lww,
    int /*i_dir*/,
    const array<double, 3>& MyU_axis,
    const array<double, 3>& MyV_axis,
    const array<double, 3>& MyW_axis,
    double& xmin, double& xmax,
    double& ymin, double& ymax,
    double& zmin, double& zmax)
{
    u_axis[nn] = MyU_axis;
    v_axis[nn] = MyV_axis;
    w_axis[nn] = MyW_axis;
    xyz0Corner[nn] = { x0, y0, z0 };
    lu[nn] = luu;
    lv[nn] = lvv;
    lw[nn] = lww;

    array<double, 3> xyz1, xyz2, xyz3, xyz4, xyz5, xyz6, xyz7;
    findCornersObliqueRectangular(nn, x0, y0, z0, luu, lvv, lww,
        xyz1, xyz2, xyz3, xyz4, xyz5, xyz6, xyz7);
    xyzMinMaxRectangular(nn, xyz1, xyz2, xyz3, xyz4, xyz5, xyz6, xyz7,
        xmin, xmax, ymin, ymax, zmin, zmax);
}

// ---------------------------------------------------------------
// Compute 7 Corners
// ---------------------------------------------------------------
void ObliqueObjects::findCornersObliqueRectangular(
    int nn, double x0, double y0, double z0,
    double lu, double lv, double lw,
    array<double, 3>& xyz1, array<double, 3>& xyz2,
    array<double, 3>& xyz3, array<double, 3>& xyz4,
    array<double, 3>& xyz5, array<double, 3>& xyz6,
    array<double, 3>& xyz7)
{
    const auto& U = u_axis[nn];
    const auto& V = v_axis[nn];
    const auto& W = w_axis[nn];

    xyz1 = { x0 + lu * U[0], y0 + lu * U[1], z0 + lu * U[2] };
    xyz2 = { x0 + lv * V[0], y0 + lv * V[1], z0 + lv * V[2] };
    xyz3 = { x0 + lw * W[0], y0 + lw * W[1], z0 + lw * W[2] };
    xyz4 = { x0 + lu * U[0] + lv * V[0],
            y0 + lu * U[1] + lv * V[1],
            z0 + lu * U[2] + lv * V[2] };
    xyz5 = { x0 + lv * V[0] + lw * W[0],
            y0 + lv * V[1] + lw * W[1],
            z0 + lv * V[2] + lw * W[2] };
    xyz6 = { x0 + lu * U[0] + lw * W[0],
            y0 + lu * U[1] + lw * W[1],
            z0 + lu * U[2] + lw * W[2] };
    xyz7 = { x0 + lu * U[0] + lv * V[0] + lw * W[0],
            y0 + lu * U[1] + lv * V[1] + lw * W[1],
            z0 + lu * U[2] + lv * V[2] + lw * W[2] };
}

// ---------------------------------------------------------------
// Find XYZ Min/Max
// ---------------------------------------------------------------
void ObliqueObjects::xyzMinMaxRectangular(
    int nn,
    const array<double, 3>& xyz1,
    const array<double, 3>& xyz2,
    const array<double, 3>& xyz3,
    const array<double, 3>& xyz4,
    const array<double, 3>& xyz5,
    const array<double, 3>& xyz6,
    const array<double, 3>& xyz7,
    double& xmin, double& xmax,
    double& ymin, double& ymax,
    double& zmin, double& zmax)
{
    array<double, 8> vals;

    // X
    vals = { xyz0Corner[nn][0], xyz1[0], xyz2[0], xyz3[0], xyz4[0], xyz5[0], xyz6[0], xyz7[0] };
    xmin = minOfArray(vals);
    xmax = maxOfArray(vals);

    // Y
    vals = { xyz0Corner[nn][1], xyz1[1], xyz2[1], xyz3[1], xyz4[1], xyz5[1], xyz6[1], xyz7[1] };
    ymin = minOfArray(vals);
    ymax = maxOfArray(vals);

    // Z
    vals = { xyz0Corner[nn][2], xyz1[2], xyz2[2], xyz3[2], xyz4[2], xyz5[2], xyz6[2], xyz7[2] };
    zmin = minOfArray(vals);
    zmax = maxOfArray(vals);
}

// ---------------------------------------------------------------
// Check if a point is inside an oblique rectangular object
// ---------------------------------------------------------------
void ObliqueObjects::checkIfPointInObliqueRectangular(
    double x, double y, double z,
    double x0, double y0, double z0,
    double u_length, double v_length, double w_length,
    bool& isit_in, int nn)
{
    isit_in = false;

    const auto& U = u_axis[nn];
    const auto& V = v_axis[nn];
    const auto& W = w_axis[nn];

    double u = (x - x0) * U[0] + (y - y0) * U[1] + (z - z0) * U[2];
    double v = (x - x0) * V[0] + (y - y0) * V[1] + (z - z0) * V[2];
    double w = (x - x0) * W[0] + (y - y0) * W[1] + (z - z0) * W[2];

    if ((u >= 0. && u <= u_length) &&
        (v >= 0. && v <= v_length) &&
        (w >= 0. && w <= w_length))
        isit_in = true;
}

// ---------------------------------------------------------------
// Utility helpers
// ---------------------------------------------------------------
double ObliqueObjects::minOfArray(const array<double, 8>& x) {
    return *min_element(x.begin(), x.end());
}
double ObliqueObjects::maxOfArray(const array<double, 8>& x) {
    return *max_element(x.begin(), x.end());
}

int main() {
    int ne_oblique = 1;
    ObliqueObjects obj(ne_oblique);

    double x0 = 5.0, y0 = 5.0, z0 = 4.0;
    double lu = 1.0, lv = 1.0, lw = 2.0;

    array<double, 3> a1_axis = { 0.0, 0.0, 1.0 };
    array<double, 3> a2_axis = { 1.0, 0.0, 0.0 };
    array<double, 3> a3_axis = { 0.0, 1.0, 0.0 };

    double xmin, xmax, ymin, ymax, zmin, zmax;
    obj.ReadObjectInfoAndFindMinMax(0, x0, y0, z0, lu, lv, lw, 3,
        a1_axis, a2_axis, a3_axis,
        xmin, xmax, ymin, ymax, zmin, zmax);

    cout << "Xmin=" << xmin << " Xmax=" << xmax << endl;
    cout << "Ymin=" << ymin << " Ymax=" << ymax << endl;
    cout << "Zmin=" << zmin << " Zmax=" << zmax << endl;

    bool inside;
    obj.checkIfPointInObliqueRectangular(5.5, 5.5, 4.5, x0, y0, z0, lu, lv, lw, inside, 0);
    cout << "Point inside? " << boolalpha << inside << endl;

    return 0;
}

