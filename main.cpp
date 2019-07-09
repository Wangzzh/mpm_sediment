#include <fstream>
#include <cstdlib>
#include <iostream>

#include "GL/gl.h"
#include "GL/glut.h"

#include "particle.hpp"
#include "grid.hpp"
#include "mpm.hpp"
#include "material.hpp"

bool simulate = false;

MPM* mpm;

void display() 
{
    glClearColor(0, 0, 0, 0);
    glClear(GL_COLOR_BUFFER_BIT);

    mpm->render();

    glFlush();
}

void idle()
{
    if (simulate) {
        mpm->step();
    }
    glutPostRedisplay();
}

void keyboard(unsigned char key, int x, int y)
{
    if (key == ' ') {simulate = !simulate; std::cout << "Simulate: " << simulate << std::endl;}
    if (key == 'g' && !simulate) {mpm->step(); std::cout << "Step" << std::endl;}
}


int main(int argc, char* argv[]) 
{
    std::string configFilename = "config.txt";

    if (argc >= 2) {
        configFilename = std::string(argv[1]);
    }

    mpm = new MPM();
    mpm->initFromConfig(configFilename);
    
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_SINGLE);
    glutInitWindowSize(600, 600);
    glutInitWindowPosition(100, 100);
    glutCreateWindow("MPM Sediment Simulation");
    glutDisplayFunc(display);
    glutKeyboardFunc(keyboard);
    glutIdleFunc(idle);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0, 1, 0, 1, 0, 100);

    glutMainLoop();
    return 0;
}