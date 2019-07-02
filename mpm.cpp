#include "mpm.hpp"

void MPM::render() {
    
}

void MPM::initFromConfig(std::string filename) {
    std::ifstream f;
    f.open(filename, std::fstream::in);
    if (!f) {
        std::cout << "Unable to load file: " << filename << std::endl;
        throw "File error";
    }

    std::string line;

    int i; double d;
    while(std::getline(f, line)) {
        std::istringstream iss(line);
        std::string word;
        iss >> word;

        if (word.size() >= 1 && word[0] == '#') {}
        else if (word == "dt") {iss >> d; dt = d;} 
        else if (word == "fluidStepPerDT") {iss >> i; fluidStepPerDT = i;}
        else if (word == "sandStepPerDT") {iss >> i; fluidStepPerDT = i;}
        else {std::cout << "Unknown param: " << word << std::endl;}
    }



    f.close();
}