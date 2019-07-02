#include "mpm.hpp"

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
            else if (word == "fluidStepPerDT") {iss >> i; fluidStepPerDT = i;}
            else if (word == "sedimentStepPerDT") {iss >> i; sedimentStepPerDT = i;}
            else if (word == "nGrid") {
                iss >> i; nGrid = i; 
                grids = std::vector<std::vector<Grid*>>(i, std::vector<Grid*>(i));
                for (i = 0; i < nGrid; i++) {
                    for (j = 0; j < nGrid; j++) {
                        grids[i][j] = new Grid();
                        grids[i][j]->position << (((double)i + 0.5) / nGrid), (((double)j + 0.5) / nGrid);
                        grids[i][j]->momentum << 0, 0;
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
                    p->momentum << x, y;
                } else {
                    p->momentum << 0, 0;
                } 
                if (word == "addFluid") {pFluid.push_back(p);}
                else {pSediment.push_back(p);}
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