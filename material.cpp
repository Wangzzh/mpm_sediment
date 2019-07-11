#include "material.hpp"

void Material::calc() {
    mu = E / 2. / (1 + nu);
    lambda = E * nu / (1 + nu) / (1 - 2 * nu);
}