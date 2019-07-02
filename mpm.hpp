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
    int sedimentStepPerDT;
    double dtFluid, dtSediment;

    // Grid
    int nGrid;
    double sGrid;

    // Parameters
    Eigen::Vector2d gravity;

    // MPM class for sediment simulation
    std::vector<std::vector<Grid*>> grids;
    std::vector<Particle*> pFluid;
    std::vector<Particle*> pSediment;
    
    // Materials
    Material fluidMaterial;
    Material sedimentMaterial;

    void initFromConfig(std::string);
    void render();

    void step();
    void fluidStep();
    void sedimentStep();

private:
    ~MPM();

};