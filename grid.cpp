#include "grid.hpp"

void Grid::renderPosition() {
    glRectd(position(0) - 0.002, position(1) - 0.002, position(0) + 0.002, position(1) + 0.002);
    
    double velocityFactor = 100;
    glBegin(GL_LINES);
    glVertex3f(position[0], position[1], 0.0);
    glVertex3f(position[0] + fMomentum[0] * velocityFactor, 
        position[1] + fMomentum[1] * velocityFactor, 0.);
    glEnd();
}