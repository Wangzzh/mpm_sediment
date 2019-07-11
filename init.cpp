#include "mpm.hpp"

void MPM::initFromConfig(std::string filename) {
    srand(0);

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

                Eigen::Vector2d unitX; unitX << 1, 0;
                Eigen::Vector2d unitY; unitY << 0, 1;
                if (iss >> d) {
                    unitX << cos(d), sin(d);
                    unitY << -sin(d), cos(d);
                }

                for (x = -0.5*size[0]; x <= 0.5*size[0]; x += interval) {
                    for (y = -0.5*size[1]; y <= 0.5*size[1]; y += interval) {
                        Particle* p = new Particle();
                        p->mass = 1;
                        p->volume = interval * interval;
                        p->position = center + x * unitX + y * unitY;
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
            else if (word == "addFluidBlockRand" || word == "addSedimentBlockRand") {
                iss >> x >> y;
                Eigen::Vector2d center;
                center << x, y;
                iss >> x >> y;
                Eigen::Vector2d size;
                size << x, y;

                iss >> x;
                double interval = x;

                Eigen::Vector2d unitX; unitX << 1, 0;
                Eigen::Vector2d unitY; unitY << 0, 1;
                if (iss >> d) {
                    unitX << cos(d), sin(d);
                    unitY << -sin(d), cos(d);
                }

                for (x = -0.5*size[0]; x <= 0.5*size[0]; x += interval) {
                    for (y = -0.5*size[1]; y <= 0.5*size[1]; y += interval) {
                        Particle* p = new Particle();
                        p->mass = 1;
                        p->volume = interval * interval;
                        p->position = ((rand() % 10000 / 10000. - 0.5) * size(0)) * unitX + 
                                      ((rand() % 10000 / 10000. - 0.5) * size(1)) * unitY + center;
                        p->velocity << 0, 0;
                        if (word == "addFluidBlockRand") {
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
            else if (word == "E") {iss >> d; materialPtr->E = d; materialPtr->calc();}
            else if (word == "nu") {iss >> d; materialPtr->nu = d; materialPtr->calc();}
            else {std::cout << "Unknown param: " << word << std::endl;}
        }

        else {std::cout << "Unknown param: " << word << std::endl;}
    }

    f.close();
}