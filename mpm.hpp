#pragma once

#include <cstdlib>
#include <cmath>
#include <ctime>
#include <random>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>

#include "Eigen/Dense"
#include "Eigen/SVD"
#include "GL/gl.h"
#include "GL/glut.h"

#include "particle.hpp"
#include "grid.hpp"
#include "material.hpp"

class MPM
{
public:
    // Time step
    double dt;
    int fluidStepPerDT;
    int sandStepPerDT;

    // MPM class for sediment simulation
    std::vector<Grid> grids;
    std::vector<Particle> pFluid;
    std::vector<Particle> pSand;

    void initFromConfig(std::string);
    void render();
    
};