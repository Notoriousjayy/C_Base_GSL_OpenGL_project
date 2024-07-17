//
// Created by jorda on 7/16/2024.
//

#include "graphics2d.h"
#include "../Linear-Algebra/svd_sparse.h"
#include <GL/gl.h>
#include <GL/glu.h>
#include <cmath>
#include <iostream>

void Vector2::normalize() {
    double sizeSq = x * x + y * y;
    if (sizeSq < 0.0000001) {
        std::cerr << "\nnormalize() sees vector (0,0)!";
        return; // does nothing to zero vectors
    }
    float scaleFactor = 1.0 / static_cast<float>(sqrt(sizeSq));
    x *= scaleFactor;
    y *= scaleFactor;
}

Canvas::Canvas(int width, int height, char* title) {
    // Constructor implementation
}

void Canvas::setWindow(float l, float r, float b, float t) {
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(static_cast<GLdouble>(l), static_cast<GLdouble>(r), static_cast<GLdouble>(b), static_cast<GLdouble>(t));
    if (t == b) return;
    windowAspect = (r - l) / (t - b);
}

void Canvas::setViewport(int l, int r, int b, int t) {
    glViewport(static_cast<GLint>(l), static_cast<GLint>(b), static_cast<GLint>(r - l), static_cast<GLint>(t - b));
}

void Canvas::lineTo(float x, float y) {
    glBegin(GL_LINES);
    glVertex2f(static_cast<GLfloat>(CP.x), static_cast<GLfloat>(CP.y));
    CP.x = x; CP.y = y;
    glVertex2f(static_cast<GLfloat>(CP.x), static_cast<GLfloat>(CP.y));
    glEnd();
    glFlush();
}

void Canvas::forward(float dist, int vis) {
    float RadPerDeg = 0.017453393; // radians per degree
    float x = CP.x + dist * cos(RadPerDeg * dist);
    float y = CP.y + dist * sin(RadPerDeg * dist);
    if (vis) lineTo(x, y);
    else moveTo(x, y);
    CP.x = x; CP.y = y;
}

void Canvas::initCT() {
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

void Canvas::rotate2D(double angle) {
    glMatrixMode(GL_MODELVIEW);
    glRotated(angle, 0.0, 0.0, 1.0);
}

void Canvas::translate2D(double dx, double dy) {
    glMatrixMode(GL_MODELVIEW);
    glTranslated(dx, dy, 0.0);
}

void Canvas::scale2D(double sx, double sy) {
    glMatrixMode(GL_MODELVIEW);
    glScaled(sx, sy, 1.0);
}

void Canvas::pushCT() {
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
}

void Canvas::popCT() {
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
}

void Canvas::ngon(int n, float cx, float cy, float radius) {
    float RadPerDeg = 0.017453393; // radians per degree
    if (n < 3) return; // bad number of sides
    double angle = 0, angleInc = 2 * 3.14159265 / n; // angle increment
    moveTo(cx + radius, cy);
    for (int k = 1; k <= n; k++) {
        angle += angleInc;
        lineTo(radius * cos(angle) + cx, radius * sin(angle) + cy);
    }
}

int RGBpixmap::readBMPFile(string fname) {
    ifstream inf(fname.c_str(), ios::in | ios::binary);
    if (!inf) {
        cout << "can't open file: " << fname << endl;
        return 0;
    }
    char ch1, ch2;
    inf.get(ch1); inf.get(ch2);
    if (ch1 != 'B' || ch2 != 'M') {
        inf.close(); return 0;
    }
    ulong fileSize = getLong(inf);
    ushort reserved1 = getShort(inf), reserved2 = getShort(inf);
    ulong offBits = getLong(inf), headerSize = getLong(inf);
    ulong numCols = getLong(inf), numRows = getLong(inf);
    ushort planes = getShort(in &#8203;:citation[oaicite:0]{index=0}&#8203;
