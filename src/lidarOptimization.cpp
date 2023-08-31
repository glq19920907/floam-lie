// Author of FLOAM: Wang Han 
// Email wh200720041@gmail.com
// Homepage https://wanghan.pro

#include "lidarOptimization.h"

EdgeAnalyticCostFunction::EdgeAnalyticCostFunction(Eigen::Vector3d curr_point_, Eigen::Vector3d last_point_a_, Eigen::Vector3d last_point_b_)
        : curr_point(curr_point_), last_point_a(last_point_a_), last_point_b(last_point_b_){

}

bool EdgeAnalyticCostFunction::Evaluate(double const *const *parameters, double *residuals, double **jacobians) const
{ 
    Eigen::Isometry3d trans = Eigen::Isometry3d::Identity();
    utils_tool_transform::ArrayToIsometry3d(parameters[0], trans);
    Eigen::Vector3d lp = trans * curr_point;

    Eigen::Vector3d nu = (lp - last_point_a).cross(lp - last_point_b);
    Eigen::Vector3d de = last_point_a - last_point_b;
    double de_norm = de.norm();
    residuals[0] = nu.norm()/de_norm;
    
    if(jacobians != NULL)
    {
        if(jacobians[0] != NULL)
        {
            Eigen::Matrix3d skew_lp = skew(lp);
            Eigen::Matrix<double, 3, 6> dp_by_se3;
            dp_by_se3.block<3,3>(0,3) = -skew_lp;
            (dp_by_se3.block<3,3>(0, 0)).setIdentity();
            Eigen::Map<Eigen::Matrix<double, 1, 6, Eigen::RowMajor> > J_se3(jacobians[0]);
            J_se3.setZero();
            Eigen::Matrix3d skew_de = skew(de);
            J_se3.block<1,6>(0,0) = - nu.transpose() / nu.norm() * skew_de * dp_by_se3/de_norm;
      
        }
    }  

    return true;
 
}   


SurfNormAnalyticCostFunction::SurfNormAnalyticCostFunction(Eigen::Vector3d curr_point_, Eigen::Vector3d plane_unit_norm_, double negative_OA_dot_norm_) 
                                                        : curr_point(curr_point_), plane_unit_norm(plane_unit_norm_), negative_OA_dot_norm(negative_OA_dot_norm_){

}

bool SurfNormAnalyticCostFunction::Evaluate(double const *const *parameters, double *residuals, double **jacobians) const
{
    Eigen::Isometry3d trans = Eigen::Isometry3d::Identity();
    utils_tool_transform::ArrayToIsometry3d(parameters[0], trans);
    Eigen::Vector3d point_w = trans * curr_point;

    residuals[0] = plane_unit_norm.dot(point_w) + negative_OA_dot_norm;

    if(jacobians != NULL)
    {
        if(jacobians[0] != NULL)
        {
            Eigen::Matrix3d skew_point_w = skew(point_w);
            Eigen::Matrix<double, 3, 6> dp_by_se3;
            dp_by_se3.block<3,3>(0,3) = -skew_point_w;
            (dp_by_se3.block<3,3>(0, 0)).setIdentity();
            Eigen::Map<Eigen::Matrix<double, 1, 6, Eigen::RowMajor> > J_se3(jacobians[0]);
            J_se3.setZero();
            J_se3.block<1,6>(0,0) = plane_unit_norm.transpose() * dp_by_se3;
   
        }
    }
    return true;

}   


bool PoseSE3Parameterization::Plus(const double *x, const double *delta, double *x_plus_delta) const
{
    Eigen::Map<const Eigen::Matrix<double, 6, 1>> lie(x);
    Eigen::Map<const Eigen::Matrix<double, 6, 1>> delta_lie(delta);

    Sophus::SE3d T = Sophus::SE3d::exp(lie);
    Sophus::SE3d delta_T = Sophus::SE3d::exp(delta_lie);
    Eigen::Matrix<double, 6, 1> x_plus_delta_lie = (delta_T * T).log();

    for(int i = 0; i < 6; ++i) x_plus_delta[i] = x_plus_delta_lie(i, 0);

    return true;
}

bool PoseSE3Parameterization::ComputeJacobian(const double *x, double *jacobian) const
{
    ceres::MatrixRef(jacobian, 6, 6) = ceres::Matrix::Identity(6, 6);
    return true;
}

Eigen::Matrix<double,3,3> skew(Eigen::Matrix<double,3,1>& mat_in){
    Eigen::Matrix<double,3,3> skew_mat;
    skew_mat.setZero();
    skew_mat(0,1) = -mat_in(2);
    skew_mat(0,2) =  mat_in(1);
    skew_mat(1,2) = -mat_in(0);
    skew_mat(1,0) =  mat_in(2);
    skew_mat(2,0) = -mat_in(1);
    skew_mat(2,1) =  mat_in(0);
    return skew_mat;
}