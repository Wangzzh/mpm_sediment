#include "mpm.hpp"

void MPM::step() {
    for (int i = 0; i < fluidStepPerDT; i++) {
        fluidStep();
    }
    for (int i = 0; i < sedimentStepPerDT; i++) {
        sedimentStep();
    }
}

void MPM::fluidStep() {
    for (int i = 0; i < nGrid; i++) {
        for (int j = 0; j < nGrid; j++) {
            grids[i][j]->mass = 0;
            grids[i][j]->fMomentum << 0,0;
        }
    }

    // Particle to grid
    for (auto &p : pFluid) {
        p->gx = floor(p->position(0) / sGrid);
        p->gy = floor(p->position(1) / sGrid);
        p->diffx = p->position(0) / sGrid - (double)p->gx - 0.5;
        p->diffy = p->position(1) / sGrid - (double)p->gy - 0.5;

        p->wx << 0.5*(0.5-p->diffx)*(0.5-p->diffx),
                 0.75-p->diffx*p->diffx,
                 0.5*(0.5+p->diffx)*(0.5+p->diffx);
        p->wy << 0.5*(0.5-p->diffy)*(0.5-p->diffy),
                 0.75-p->diffy*p->diffy,
                 0.5*(0.5+p->diffy)*(0.5+p->diffy);

        for (int dx = -1; dx <= 1; dx++) {
            for (int dy = -1; dy <= 1; dy ++) {
                if (p->gx+dx>=0 && p->gx+dx<nGrid && p->gy+dy>=0 && p->gy+dy<nGrid) {
                    grids[p->gx+dx][p->gy+dy]->mass += p->wx(dx+1) * p->wy(dy+1) * p->mass;
                    grids[p->gx+dx][p->gy+dy]->fMomentum += p->wx(dx+1) * p->wy(dy+1) * p->mass * p->velocity;
                }
            }
        }
    }

    // forces
    for (int i = 1; i + 1 < nGrid; i++) {
        for (int j = 1; j + 1 < nGrid; j++) {
            grids[i][j]->fMomentum += gravity * grids[i][j]->mass * dtFluid;
        }
        if (grids[i][0]->fMomentum(1) < 0) {grids[i][0]->fMomentum(1) = 0;}
    }

    // projection
    for (int it = 0; it < 10; it ++) {
        for (int i = 1; i + 1 < nGrid; i++) {
            for (int j = 1; j + 1 < nGrid; j++) {
                double inflow = grids[i+1][j]->fMomentum(0) + grids[i][j+1]->fMomentum(1) - grids[i][j-1]->fMomentum(1) - grids[i-1][j]->fMomentum(0);
                if (i-1>0) grids[i-1][j]->fMomentum(0) += 0.25 * inflow;
                if (i+1<nGrid-1) grids[i+1][j]->fMomentum(0) -= 0.25 * inflow;
                if (j-1>0) grids[i][j-1]->fMomentum(1) += 0.25 * inflow;
                if (j-1<nGrid-1) grids[i][j+1]->fMomentum(1) -= 0.25 * inflow;
            }
        }
    }

    // Grid to Particle
    for (auto &p : pFluid) {
        p->velocity << 0,0;
        for (int dx = -1; dx <= 1; dx++) {
            for (int dy = -1; dy <= 1; dy ++) {
                if (p->gx+dx>=0 && p->gx+dx<nGrid && p->gy+dy>=0 && p->gy+dy<nGrid) {
                    p->velocity += grids[p->gx+dx][p->gy+dy]->fMomentum * p->wx(dx+1) * p->wy(dy+1) / grids[p->gx+dx][p->gy+dy]->mass;
                }
            }
        }
        p->position += dtFluid * p->velocity;
    }
    
}

void MPM::sedimentStep() {
    for (int i = 0; i < nGrid; i++) {
        for (int j = 0; j < nGrid; j++) {
            grids[i][j]->mass = 0;
            grids[i][j]->sMomentum << 0,0;
        }
    }

    // Particle to grid
    for (auto &p : pSediment) {
        p->gx = floor(p->position(0) / sGrid);
        p->gy = floor(p->position(1) / sGrid);
        p->diffx = p->position(0) / sGrid - (double)p->gx - 0.5;
        p->diffy = p->position(1) / sGrid - (double)p->gy - 0.5;

        p->wx << 0.5*(0.5-p->diffx)*(0.5-p->diffx),
                 0.75-p->diffx*p->diffx,
                 0.5*(0.5+p->diffx)*(0.5+p->diffx);
        p->wy << 0.5*(0.5-p->diffy)*(0.5-p->diffy),
                 0.75-p->diffy*p->diffy,
                 0.5*(0.5+p->diffy)*(0.5+p->diffy);

        for (int dx = -1; dx <= 1; dx++) {
            for (int dy = -1; dy <= 1; dy ++) {
                if (p->gx+dx>=0 && p->gx+dx<nGrid && p->gy+dy>=0 && p->gy+dy<nGrid) {
                    grids[p->gx+dx][p->gy+dy]->mass += p->wx(dx+1) * p->wy(dy+1) * p->mass;
                    grids[p->gx+dx][p->gy+dy]->sMomentum += p->wx(dx+1) * p->wy(dy+1) * p->mass * 
                        (p->velocity + 4. * sGrid * sGrid * p->B * (grids[p->gx+dx][p->gy+dy]->position - p->position));
                }
            }
        }
    }

    // forces
    for (int i = 0; i < nGrid; i++) {
        for (int j = 0; j < nGrid; j++) {
            grids[i][j]->sMomentum += gravity * grids[i][j]->mass * dtSediment;
        }
    }

    // Grid to Particle
    for (auto &p : pSediment) {
        p->velocity << 0,0;
        for (int dx = -1; dx <= 1; dx++) {
            for (int dy = -1; dy <= 1; dy ++) {
                if (p->gx+dx>=0 && p->gx+dx<nGrid && p->gy+dy>=0 && p->gy+dy<nGrid) {
                    p->velocity += grids[p->gx+dx][p->gy+dy]->sMomentum * p->wx(dx+1) * p->wy(dy+1) / grids[p->gx+dx][p->gy+dy]->mass;
                    p->B += grids[p->gx+dx][p->gy+dy]->sMomentum * p->wx(dx+1) * p->wy(dy+1) * (grids[p->gx+dx][p->gy+dy]->position - p->position).transpose(); 
                }
            }
        }
        p->position += dtSediment * p->velocity;
    }
}

void MPM::render() {
    // Particles
    glColor3f(0.5, 0.5, 1);
    for (auto &p : pFluid) {
        p->renderPosition();
    }
    glColor3f(1, 0.8, 0.5);
    for (auto &p : pSediment) {
        p->renderPosition();
    }

    // Grids
    glColor3f(0.3, 0.3, 0.3);
    for (int i = 0; i < nGrid; i++) {
        for (int j = 0; j < nGrid; j++) {
            grids[i][j]->renderPosition();
        }
    }
}

void MPM::initFromConfig(std::string filename) {
    std::ifstream f;
    f.open(filename, std::fstream::in);
    if (!f) {
        std::cout << "Unable to load file: " << filename << std::endl;
        throw "File error";
    }

    enum {NORMAL, MATERIAL} state = NORMAL;
    Material* materialPtr;
    std::string line;

    int i, j; double d, x, y, z;
    while(std::getline(f, line)) {
        std::istringstream iss(line);
        std::string word;
        iss >> word;

        if (word == "") {}
        else if (word.size() >= 1 && word[0] == '#') {}
        else if (word == "startFluidMaterial") {state = MATERIAL; materialPtr = &fluidMaterial;}
        else if (word == "startSedimentMaterial") {state = MATERIAL; materialPtr = &sedimentMaterial;}

        else if (state == NORMAL) {
            if (word == "dt") {iss >> d; dt = d;} 
            else if (word == "fluidStepPerDT") {iss >> i; fluidStepPerDT = i; dtFluid = dt / (double)fluidStepPerDT;}
            else if (word == "sedimentStepPerDT") {iss >> i; sedimentStepPerDT = i; dtSediment = dt / (double)sedimentStepPerDT;}
            else if (word == "gravity") {iss >> x >> y; gravity << x, y;}
            else if (word == "nGrid") {
                iss >> i; nGrid = i; sGrid = 1. / (double)i;
                grids = std::vector<std::vector<Grid*>>(i, std::vector<Grid*>(i));
                for (i = 0; i < nGrid; i++) {
                    for (j = 0; j < nGrid; j++) {
                        grids[i][j] = new Grid();
                        grids[i][j]->position << (((double)i + 0.5) / nGrid), (((double)j + 0.5) / nGrid);
                        grids[i][j]->fMomentum << 0, 0;
                        grids[i][j]->sMomentum << 0, 0;
                    }
                }
            }
            else if (word == "addFluid" || word == "addSediment") {
                Particle* p = new Particle();
                iss >> x;
                p->mass = x;
                iss >> x >> y;
                p->position << x, y;
                if (iss >> x >> y) {
                    p->velocity << x, y;
                } else {
                    p->velocity << 0, 0;
                } 
                if (word == "addFluid") {pFluid.push_back(p);}
                else {pSediment.push_back(p);}
            }
            else if (word == "addFluidBlock" || word == "addSedimentBlock") {
                iss >> x >> y;
                Eigen::Vector2d center;
                center << x, y;
                iss >> x >> y;
                Eigen::Vector2d size;
                size << x, y;

                iss >> x;
                double interval = x;

                for (x = center[0] - 0.5*size[0]; x <= center[0] + 0.5*size[0]; x += interval) {
                    for (y = center[1] - 0.5*size[1]; y <= center[1] + 0.5*size[1]; y += interval) {
                        Particle* p = new Particle();
                        p->mass = 1;
                        p->position << x, y;
                        p->velocity << 0, 0;
                        if (word == "addFluidBlock") {
                            p->mass = interval * interval * fluidMaterial.density;
                            pFluid.push_back(p);
                        }
                        else {
                            p->mass = interval * interval * sedimentMaterial.density;
                            pSediment.push_back(p);
                        }
                    }
                }
            }
            else {std::cout << "Unknown param: " << word << std::endl;}
        }
        
        else if (state == MATERIAL) {
            if (word == "endMaterial") {state = NORMAL;}
            else if (word == "density") {iss >> d; materialPtr->density = d;}
            else {std::cout << "Unknown param: " << word << std::endl;}
        }

        else {std::cout << "Unknown param: " << word << std::endl;}
    }

    f.close();
}

MPM::~MPM() {
    for (auto &p : pFluid) {delete p;}
    for (auto &p : pSediment) {delete p;}
    for (auto &gv : grids) {for (auto &g : gv) {delete g;}}
}