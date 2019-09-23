#ifndef _COMPADRE_MISC_HPP_
#define _COMPADRE_MISC_HPP_

namespace Compadre {

struct XYZ {

    KOKKOS_INLINE_FUNCTION
    XYZ() : x(0), y(0), z(0) {}

    KOKKOS_INLINE_FUNCTION
    XYZ(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {}

    double x;
    double y;
    double z;

    KOKKOS_INLINE_FUNCTION
    double& operator [](const int i) {
        switch (i) {
            case 0:
                return x;
            case 1:
                return y;
            default:
                return z;
        }
    }

    KOKKOS_INLINE_FUNCTION
    XYZ operator *(double scalar) {
        XYZ result;
        result.x = scalar*x;
        result.y = scalar*y;
        result.z = scalar*z;
        return result;
    }
}; // XYZ

KOKKOS_INLINE_FUNCTION
void getRHSDims(DenseSolverType dense_solver_type, BoundaryType boundary_type, const int M, const int N, int* dims) {
    if (_dense_solver_type != DenseSolverType::LU) {
        dims[0] = max_num_rows;
        dims[1] = max_num_rows;
    } else {
        dims[0] = this_num_cols;
        dims[1] = max_num_rows;
    }
}

}; // Compadre

#endif
