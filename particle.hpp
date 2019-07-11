#pragma once

#include "Eigen/Dense"
#include "Eigen/SVD"
#include "GL/gl.h"
#include "GL/glut.h"

class Particle
{
public:
    double mass;
    double volume;
    
    Eigen::Vector2d position;
    Eigen::Vector2d velocity;

    Eigen::Matrix2d B = Eigen::Matrix2d::Zero();
    Eigen::Matrix2d FE = Eigen::Matrix2d::Identity();
    Eigen::Matrix2d FP = Eigen::Matrix2d::Identity();

    int gx, gy;
    double diffx, diffy;
    Eigen::Vector3d wx, wy;
    Eigen::Vector3d dwx, dwy;

    // sand
    double q = 0, phi = 30, alpha = sqrt(0.6667) / 2.5;

    void renderPosition();
};