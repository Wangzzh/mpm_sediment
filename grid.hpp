#pragma once

#include "Eigen/Dense"
#include "Eigen/SVD"
#include "GL/gl.h"
#include "GL/glut.h"

class Grid
{
public:
    double mass;
    Eigen::Vector2d position;
    Eigen::Vector2d momentum;

    void renderPosition();
};