#include "particle.hpp"

void Particle::renderPosition() {
    glRectd(position(0) - 0.002, position(1) - 0.002, position(0) + 0.002, position(1) + 0.002);
}