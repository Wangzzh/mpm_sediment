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
            grids[i][j]->pressure = 0;
        }
        if (grids[i][0]->fMomentum(1) < 0) {grids[i][0]->fMomentum(1) = 0;}
    }

    // projection
    for (int it = 0; it < 30; it ++) {
        for (int i = 1; i + 1 < nGrid; i++) {
            for (int j = 1; j + 1 < nGrid; j++) {
                double inflow = grids[i+1][j]->fMomentum(0) / (grids[i+1][j] -> mass + 0.0001)+ 
                    grids[i][j+1]->fMomentum(1) / (grids[i][j+1]->mass + 0.0001) - 
                    grids[i][j-1]->fMomentum(1) / (grids[i][j-1]->mass + 0.0001) - 
                    grids[i-1][j]->fMomentum(0) / (grids[i-1][j]->mass + 0.0001);

                grids[i][j] -> pressure = 0.25 * (
                    (grids[i+1][j]->pressure + grids[i-1][j]->pressure + grids[i][j+1]->pressure + grids[i][j-1]->pressure)
                    - 0.5 / dtFluid * inflow  * grids[i][j]->mass
                );
            }
        }
        // set boundaries
        for (int i = 1; i + 1 < nGrid; i++) {
            grids[0][i] -> pressure = grids[1][i] -> pressure;
            grids[nGrid-1][i]->pressure = grids[nGrid-2][i]->pressure;
            grids[i][0] -> pressure = grids[i][1]->pressure;
            grids[i][nGrid-1]->pressure = grids[i][nGrid-2]->pressure;
        }
    }
    for (int i = 1; i + 1 < nGrid; i++) {
        for (int j = 1; j + 1 < nGrid; j++) {
            grids[i][j] -> fMomentum(0) += 0.5 * dtFluid * (grids[i-1][j]->pressure - grids[i+1][j]->pressure) * grids[i][j]->mass;
            grids[i][j] -> fMomentum(1) += 0.5 * dtFluid * (grids[i][j-1]->pressure - grids[i][j+1]->pressure) * grids[i][j]->mass;
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

        p->wx << 0.5*(1.5-(1+p->diffx))*(1.5-(1+p->diffx)),
                 0.75-p->diffx*p->diffx,
                 0.5*(1.5-(1-p->diffx))*(1.5-(1-p->diffx));
        p->wy << 0.5*(1.5-(1+p->diffy))*(1.5-(1+p->diffy)),
                 0.75-p->diffy*p->diffy,
                 0.5*(1.5-(1-p->diffy))*(1.5-(1-p->diffy));

        p->dwx << -1.5+(1+p->diffx),
                  -2. * p->diffx,
                  1.5-(1-p->diffx);
        p->dwy << -1.5+(1+p->diffy),
                  -2. * p->diffy,
                  1.5-(1-p->diffy);

        for (int dx = -1; dx <= 1; dx++) {
            for (int dy = -1; dy <= 1; dy ++) {
                if (p->gx+dx>=0 && p->gx+dx<nGrid && p->gy+dy>=0 && p->gy+dy<nGrid) {
                    grids[p->gx+dx][p->gy+dy]->mass += p->wx(dx+1) * p->wy(dy+1) * p->mass;
                    grids[p->gx+dx][p->gy+dy]->sMomentum += p->wx(dx+1) * p->wy(dy+1) * p->mass * 
                        (p->velocity + 4. * p->B * (grids[p->gx+dx][p->gy+dy]->position - p->position));
                }
            }
        }
        
        // stress force
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(p->FE, Eigen::ComputeThinU | Eigen::ComputeThinV);
        Eigen::Matrix2d U = svd.matrixU();
        Eigen::Matrix2d V = svd.matrixV();
        Eigen::Matrix2d Sig = U.transpose() * p->FE * V;
        Eigen::Matrix2d LogSig; LogSig << log(Sig(0, 0)), 0, 0, log(Sig(1, 1));
        // std::cout << p->FE << std::endl  << LogSig << std::endl << std::endl;

        // Snow
        Eigen::Matrix2d RE = U * V.transpose();
        double JE = p->FE.determinant();
        Eigen::Matrix2d PF = 2 * sedimentMaterial.mu * (p->FE - RE) + sedimentMaterial.lambda * (JE - 1) * JE * p->FE.transpose().inverse();
        // Sand
        // Eigen::Matrix2d PF = U * (2 * sedimentMaterial.mu * Sig.inverse() * LogSig 
        //     + sedimentMaterial.lambda * LogSig.trace() * Sig.inverse()) * V.transpose();
        
        for (int dx = -1; dx <= 1; dx++) {
            for (int dy = -1; dy <= 1; dy ++) {
                if (p->gx+dx>=0 && p->gx+dx<nGrid && p->gy+dy>=0 && p->gy+dy<nGrid) {
                    Eigen::Vector2d dWeight; dWeight << p->dwx(dx+1) * p->wy(dy+1), p->wx(dx+1) * p->dwy(dy+1);
                    grids[p->gx+dx][p->gy+dy] -> sMomentum -= dtSediment * 
                        p->volume * PF * p->FE.transpose() * dWeight;
                }
            }
        }
    }

    // forces
    double friction_mu = 0.5;
    for (int i = 0; i < nGrid; i++) {
        for (int j = 0; j < nGrid; j++) {
            grids[i][j]->sMomentum += gravity * grids[i][j]->mass * dtSediment;
            
            // collision
            if (grids[i][j]->position(1) < 0.1) {
                if (abs(grids[i][j] -> sMomentum(0)) < friction_mu * abs(grids[i][j] -> sMomentum(1))) {
                    grids[i][j] -> sMomentum(1) = 0;
                    grids[i][j] -> sMomentum(0) = 0;
                } else {
                    grids[i][j] -> sMomentum(0) = grids[i][j] -> sMomentum(0) / abs(grids[i][j] -> sMomentum(0)) * 
                        (abs(grids[i][j] -> sMomentum(0)) - friction_mu * abs(grids[i][j] -> sMomentum(1)));
                    grids[i][j] -> sMomentum(1) = 0;
                }
            }
        }
    }

    // Grid to Particle
    for (auto &p : pSediment) {
        p->velocity << 0,0;
        Eigen::Matrix2d deltaF = Eigen::Matrix2d::Zero();
        p->B = Eigen::Matrix2d::Zero();

        for (int dx = -1; dx <= 1; dx++) {
            for (int dy = -1; dy <= 1; dy ++) {
                if (p->gx+dx>=0 && p->gx+dx<nGrid && p->gy+dy>=0 && p->gy+dy<nGrid) {
                    Eigen::Vector2d dWeight; dWeight << p->dwx(dx+1) * p->wy(dy+1), p->wx(dx+1) * p->dwy(dy+1);
                    p->velocity += grids[p->gx+dx][p->gy+dy]->sMomentum * p->wx(dx+1) * p->wy(dy+1) / grids[p->gx+dx][p->gy+dy]->mass;
                    p->B += grids[p->gx+dx][p->gy+dy]->sMomentum  / grids[p->gx+dx][p->gy+dy]->mass * 
                        p->wx(dx+1) * p->wy(dy+1) * (grids[p->gx+dx][p->gy+dy]->position - p->position).transpose(); 

                    deltaF += grids[p->gx+dx][p->gy+dy]->sMomentum / grids[p->gx+dx][p->gy+dy]->mass * dWeight.transpose();
                }
            }
        }
        p->position += dtSediment * p->velocity;
        // if (p->position(1) < 0.1) {
        //     p->position(1) = 0.1;
        //     if (p->velocity(1) < 0) {
        //         p->velocity(1) = 0;   
        //     }
        // }

        Eigen::Matrix2d Fnew = (Eigen::Matrix2d::Identity() + dtSediment * deltaF) * p->FE;
        
        // Project to yield space
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(Fnew, Eigen::ComputeThinU | Eigen::ComputeThinV);
        Eigen::Matrix2d U = svd.matrixU();
        Eigen::Matrix2d V = svd.matrixV();
        Eigen::Matrix2d Sig = U.transpose() * Fnew * V;
        Eigen::Matrix2d e; e << log(Sig(0, 0)), 0, 0, log(Sig(1, 1));

        Eigen::Matrix2d ehat = e - e.trace() / 2. * Eigen::Matrix2d::Identity();
        double normF_ehat = sqrt((ehat * ehat).trace());
        double dgamma = normF_ehat + 
            (sedimentMaterial.lambda + sedimentMaterial.mu) / sedimentMaterial.mu * 
            e.trace() * p->alpha;
        if (dgamma < 0) {
            // p->color = Eigen::Vector3d(0, 1, 0);
        } 
        else if (normF_ehat == 0 || e.trace()>0) {
            Sig = Eigen::Matrix2d::Identity();
            p->q += sqrt((e * e).trace());
            // p->color = Eigen::Vector3d(1, 0, 0);
        } 
        else {
            Eigen::Matrix2d H = e - dgamma * ehat / normF_ehat;
            Sig << exp(H(0, 0)), 0, 0, exp(H(1, 1));
            p->q += dgamma;
            // p->color = Eigen::Vector3d(0, 0, 1);
        }
        p->phi = sedimentMaterial.h0 + (sedimentMaterial.h1 * p->q - sedimentMaterial.h3) * 
            exp(-sedimentMaterial.h2 * p->q);
        p->phi = 30;
        p->alpha = sqrt(2./3.) * 2 * sin(p->phi*3.1416/180) / (3 - sin(p->phi*3.1416/180));

        p->FP = (V * Sig.inverse() * U.transpose() * Fnew) * p->FP;
        p->FE = U * Sig * V.transpose();

        // p->FE = Fnew;
    }
}

void MPM::render() {
    // Particles
    for (auto &p : pFluid) {
        p->renderPosition();
    }
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

MPM::~MPM() {
    for (auto &p : pFluid) {delete p;}
    for (auto &p : pSediment) {delete p;}
    for (auto &gv : grids) {for (auto &g : gv) {delete g;}}
}