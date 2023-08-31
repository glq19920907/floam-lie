#pragma once

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <pcl/common/transforms.h>

#include <sophus/se3.hpp>

// transform
namespace utils_tool_transform
{
    // double ( lie <--> Isometry )
    // double xi[6]; // t_x,t_y,t_z,phi_x,phi_y,phi_z
    inline void ArrayToIsometry3d(const double * xi, Eigen::Isometry3d & trans)
    {
        Eigen::Map<const Eigen::Matrix<double, 6, 1>> lie(xi);
        Sophus::SE3d SE3_T = Sophus::SE3d::exp(lie);
        trans = SE3_T.matrix();
    }

    inline void Isometry3dToArray(const Eigen::Isometry3d trans, double * xi)
    {
        Eigen::Map<Eigen::Matrix<double, 6, 1>> lie(xi);
        Sophus::SE3d SE3_T(trans.matrix());
        lie = Sophus::SE3d::log(SE3_T);
    }
}

namespace utils_tool_math
{
    
}

namespace utils_tool_visualization
{
    template<typename T>
    inline void printArray(T arr[], int size) {
        for (int i = 0; i < size; i++) {
            std::cout << arr[i] << " ";
        }
        std::cout << std::endl;
    }
}
