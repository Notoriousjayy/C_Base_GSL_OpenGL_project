//
// Created by jorda on 7/19/2024.
//

#ifndef GRAPHICS2D_H
#define GRAPHICS2D_H

#include <GL/gl.h>
#include <GL/glu.h>

// Define the Vector2 structure
typedef struct {
    float x;
    float y;
} Vector2;

// Define the Canvas structure
typedef struct {
    Vector2 CP;
    float windowAspect;
} Canvas;

void normalize(Vector2 *v);
void Canvas_init(Canvas *c, int width, int height, char* title);
void Canvas_setWindow(Canvas *c, float l, float r, float b, float t);
void Canvas_setViewport(Canvas *c, int l, int r, int b, int t);
void Canvas_lineTo(Canvas *c, float x, float y);
void Canvas_forward(Canvas *c, float dist, int vis);
void Canvas_initCT();
void Canvas_rotate2D(double angle);
void Canvas_translate2D(double dx, double dy);
void Canvas_scale2D(double sx, double sy);
void Canvas_pushCT();
void Canvas_popCT();
void Canvas_ngon(Canvas *c, int n, float cx, float cy, float radius);
int RGBpixmap_readBMPFile(const char *fname);

#endif // GRAPHICS2D_H
