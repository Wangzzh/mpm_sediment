#pragma once

#include "Eigen/Dense"
#include "Eigen/SVD"
#include "GL/gl.h"
#include "GL/glut.h"

class Particle
{
public:
    double mass;
    Eigen::Vector2d position;
    Eigen::Vector2d velocity;

    Eigen::Matrix2d B = Eigen::Matrix2d::Zero();

    int gx, gy;
    double diffx, diffy;
    Eigen::Vector3d wx, wy;

    void renderPosition();
};