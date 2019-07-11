#include "particle.hpp"

void Particle::renderPosition() {
    glRectd(position(0) - 0.002, position(1) - 0.002, position(0) + 0.002, position(1) + 0.002);
    
    // double velocityFactor = 2;
    // glBegin(GL_LINES);
    // glVertex3f(position[0], position[1], 0.0);
    // glVertex3f(position[0] + velocity[0] * mass * velocityFactor, 
    //     position[1] + velocity[1] * mass * velocityFactor, 0.);
    // glEnd();
}