//
// Created by jorda on 7/19/2024.
//

#ifndef _GRAPHICS2D_H_
#define _GRAPHICS2D_H_

#include <stdio.h>
#include <stdlib.h>

typedef struct {
    int x, y;
} IntPoint;

typedef struct {
    float x, y;
} Point2;

typedef struct {
    int left, top, right, bott;
} IntRect;

typedef struct {
    float x, y;
} Vector2;

typedef struct {
    unsigned char r, g, b;
} mRGB;

typedef struct {
    mRGB* pixel; // array of pixels
    int nRows, nCols; // dimensions of the pix map
} RGBpixmap;

typedef struct {
    int num;
    Point2 pt[80]; // may need larger arrays in some circumstances
} PolyLine;

typedef struct {
    Point2 CP; // current position in world
    float windowAspect;
} Canvas;

// Function prototypes for drawing functions
void drawLine(double x0, double y0, double x1, double y1);
void drawCircle(double x, double y, double r);
void drawFilledCircle(double x, double y, double r);
void drawRectangle(double x, double y, double width, double height);
void drawFilledRectangle(double x, double y, double width, double height);

// IntPoint functions
void setIntPoint(IntPoint* p, int dx, int dy);
void setIntPointFromPoint(IntPoint* p, IntPoint* src);

// Point2 functions
void setPoint2(Point2* p, float dx, float dy);
void setPoint2FromPoint(Point2* p, Point2* src);

// IntRect functions
void setIntRect(IntRect* r, int l, int t, int ri, int b);
void setIntRectFromRect(IntRect* r, IntRect* src);

// Vector2 functions
void setVector2(Vector2* v, float dx, float dy);
void setVector2FromVector(Vector2* v, Vector2* src);
void setVector2Diff(Vector2* v, Point2* a, Point2* b);
void normalizeVector2(Vector2* v);
float dotVector2(Vector2* v, Vector2* b);
void perpVector2(Vector2* v);
float perpDotVector2(Vector2* v, Vector2* b);

// Canvas functions
void initCanvas(Canvas* canvas, int width, int height, char* title);
void setCanvasWindow(Canvas* canvas, float l, float r, float b, float t);
void setCanvasViewport(Canvas* canvas, int l, int r, int b, int t);
float getCanvasWindowAspect(Canvas* canvas);
void canvasLineTo(Canvas* canvas, float x, float y);
void canvasMoveTo(Canvas* canvas, float x, float y);
void canvasForward(Canvas* canvas, float dist, int vis);
void canvasInitCT(Canvas* canvas);
void canvasRotate2D(Canvas* canvas, double angle);
void canvasTranslate2D(Canvas* canvas, double dx, double dy);
void canvasScale2D(Canvas* canvas, double sx, double sy);
void canvasPushCT(Canvas* canvas);
void canvasPopCT(Canvas* canvas);
void canvasNgon(Canvas* canvas, int n, float cx, float cy, float radius);

// RGBpixmap functions
int readBMPFile(RGBpixmap* pixmap, char* fname);
void freeRGBpixmap(RGBpixmap* pixmap);
void copyRGBpixmap(RGBpixmap* pixmap, IntPoint from, IntPoint to, int width, int height);
void drawRGBpixmap(RGBpixmap* pixmap);
int readRGBpixmap(RGBpixmap* pixmap, int x, int y, int width, int height);
int readRGBpixmapFromRect(RGBpixmap* pixmap, IntRect r);
void setRGBpixmapPixel(RGBpixmap* pixmap, int x, int y, mRGB color);
mRGB getRGBpixmapPixel(RGBpixmap* pixmap, int x, int y);

// Utility functions
unsigned short getShort(FILE* fp);
unsigned long getLong(FILE* fp);

#endif // _GRAPHICS2D_H_

