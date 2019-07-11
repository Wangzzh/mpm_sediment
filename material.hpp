#pragma once

struct Material
{
    double density;
    
    double E = 10000;
    double nu = 0.2;
    double mu;
    double lambda;

    double h0 = 35, h1 = 9, h2 = 0.2, h3 = 10;

    void calc();
};