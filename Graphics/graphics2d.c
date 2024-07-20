//
// Created by jorda on 7/19/2024.
//

#include "graphics2d.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// Normalize the vector
void normalize(Vector2 *v) {
    double sizeSq = v->x * v->x + v->y * v->y;
    if (sizeSq < 0.0000001) {
        fprintf(stderr, "\nnormalize() sees vector (0,0)!");
        return;
    }
    float scaleFactor = 1.0 / (float)sqrt(sizeSq);
    v->x *= scaleFactor;
    v->y *= scaleFactor;
}

// Initialize the Canvas
void Canvas_init(Canvas *c, int width, int height, char* title) {
    // Constructor implementation (not needed for this example)
}

// Set the window for the Canvas
void Canvas_setWindow(Canvas *c, float l, float r, float b, float t) {
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D((GLdouble)l, (GLdouble)r, (GLdouble)b, (GLdouble)t);
    if (t == b) return;
    c->windowAspect = (r - l) / (t - b);
}

// Set the viewport for the Canvas
void Canvas_setViewport(Canvas *c, int l, int r, int b, int t) {
    glViewport((GLint)l, (GLint)b, (GLint)(r - l), (GLint)(t - b));
}

// Draw a line to a point in the Canvas
void Canvas_lineTo(Canvas *c, float x, float y) {
    glBegin(GL_LINES);
    glVertex2f((GLfloat)c->CP.x, (GLfloat)c->CP.y);
    c->CP.x = x;
    c->CP.y = y;
    glVertex2f((GLfloat)c->CP.x, (GLfloat)c->CP.y);
    glEnd();
    glFlush();
}

// Move forward in the Canvas
void Canvas_forward(Canvas *c, float dist, int vis) {
    float RadPerDeg = 0.017453393f;
    float x = c->CP.x + dist * cos(RadPerDeg * dist);
    float y = c->CP.y + dist * sin(RadPerDeg * dist);
    if (vis) Canvas_lineTo(c, x, y);
    else {
        c->CP.x = x;
        c->CP.y = y;
    }
}

// Initialize the transformation matrix
void Canvas_initCT() {
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

// Rotate the Canvas
void Canvas_rotate2D(double angle) {
    glMatrixMode(GL_MODELVIEW);
    glRotated(angle, 0.0, 0.0, 1.0);
}

// Translate the Canvas
void Canvas_translate2D(double dx, double dy) {
    glMatrixMode(GL_MODELVIEW);
    glTranslated(dx, dy, 0.0);
}

// Scale the Canvas
void Canvas_scale2D(double sx, double sy) {
    glMatrixMode(GL_MODELVIEW);
    glScaled(sx, sy, 1.0);
}

// Push the current transformation matrix
void Canvas_pushCT() {
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
}

// Pop the current transformation matrix
void Canvas_popCT() {
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
}

// Draw an n-gon on the Canvas
void Canvas_ngon(Canvas *c, int n, float cx, float cy, float radius) {
    float RadPerDeg = 0.017453393f;
    if (n < 3) return;
    double angle = 0, angleInc = 2 * 3.14159265 / n;
    c->CP.x = cx + radius;
    c->CP.y = cy;
    for (int k = 1; k <= n; k++) {
        angle += angleInc;
        Canvas_lineTo(c, radius * cos(angle) + cx, radius * sin(angle) + cy);
    }
}

// Read a BMP file
int RGBpixmap_readBMPFile(const char *fname) {
    FILE *inf = fopen(fname, "rb");
    if (!inf) {
        printf("can't open file: %s\n", fname);
        return 0;
    }
    char ch1, ch2;
    fread(&ch1, sizeof(char), 1, inf);
    fread(&ch2, sizeof(char), 1, inf);
    if (ch1 != 'B' || ch2 != 'M') {
        fclose(inf);
        return 0;
    }
    unsigned long fileSize;
    fread(&fileSize, sizeof(unsigned long), 1, inf);
    unsigned short reserved1, reserved2;
    fread(&reserved1, sizeof(unsigned short), 1, inf);
    fread(&reserved2, sizeof(unsigned short), 1, inf);
    unsigned long offBits, headerSize, numCols, numRows;
    fread(&offBits, sizeof(unsigned long), 1, inf);
    fread(&headerSize, sizeof(unsigned long), 1, inf);
    fread(&numCols, sizeof(unsigned long), 1, inf);
    fread(&numRows, sizeof(unsigned long), 1, inf);
    unsigned short planes;
    fread(&planes, sizeof(unsigned short), 1, inf);
    // Additional reading code goes here...

    fclose(inf);
    return 1;
}
