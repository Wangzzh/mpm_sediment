#include "grid.hpp"

void Grid::renderPosition() {
    glRectd(position(0) - 0.002, position(1) - 0.002, position(0) + 0.002, position(1) + 0.002);
}