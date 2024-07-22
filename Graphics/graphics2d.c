//
// Created by jorda on 7/19/2024.
//

#include "graphics2d.h"
#include <GL/gl.h>
#include <GL/glu.h>
#include <math.h>
#include <stdio.h>


// Drawing functions
void drawLine(double x0, double y0, double x1, double y1) {
    glBegin(GL_LINES);
    glVertex2d(x0, y0);
    glVertex2d(x1, y1);
    glEnd();
}

void drawCircle(double x, double y, double r) {
    glBegin(GL_LINE_LOOP);
    for (int i = 0; i < 360; i++) {
        double theta = i * 3.1415926 / 180;
        double cx = r * cos(theta);
        double cy = r * sin(theta);
        glVertex2d(x + cx, y + cy);
    }
    glEnd();
}

void drawFilledCircle(double x, double y, double r) {
    glBegin(GL_TRIANGLE_FAN);
    glVertex2d(x, y);
    for (int i = 0; i <= 360; i++) {
        double theta = i * 3.1415926 / 180;
        double cx = r * cos(theta);
        double cy = r * sin(theta);
        glVertex2d(x + cx, y + cy);
    }
    glEnd();
}

void drawRectangle(double x, double y, double width, double height) {
    glBegin(GL_LINE_LOOP);
    glVertex2d(x, y);
    glVertex2d(x + width, y);
    glVertex2d(x + width, y + height);
    glVertex2d(x, y + height);
    glEnd();
}

void drawFilledRectangle(double x, double y, double width, double height) {
    glBegin(GL_QUADS);
    glVertex2d(x, y);
    glVertex2d(x + width, y);
    glVertex2d(x + width, y + height);
    glVertex2d(x, y + height);
    glEnd();
}

// IntPoint functions
void setIntPoint(IntPoint* p, int dx, int dy) {
    p->x = dx;
    p->y = dy;
}

void setIntPointFromPoint(IntPoint* p, IntPoint* src) {
    p->x = src->x;
    p->y = src->y;
}

// Point2 functions
void setPoint2(Point2* p, float dx, float dy) {
    p->x = dx;
    p->y = dy;
}

void setPoint2FromPoint(Point2* p, Point2* src) {
    p->x = src->x;
    p->y = src->y;
}

// IntRect functions
void setIntRect(IntRect* r, int l, int t, int ri, int b) {
    r->left = l;
    r->top = t;
    r->right = ri;
    r->bott = b;
}

void setIntRectFromRect(IntRect* r, IntRect* src) {
    r->left = src->left;
    r->top = src->top;
    r->right = src->right;
    r->bott = src->bott;
}

// Vector2 functions
void setVector2(Vector2* v, float dx, float dy) {
    v->x = dx;
    v->y = dy;
}

void setVector2FromVector(Vector2* v, Vector2* src) {
    v->x = src->x;
    v->y = src->y;
}

void setVector2Diff(Vector2* v, Point2* a, Point2* b) {
    v->x = a->x - b->x;
    v->y = a->y - b->y;
}

void normalizeVector2(Vector2* v) {
    double sizeSq = v->x * v->x + v->y * v->y;
    if (sizeSq < 0.0000001) {
        fprintf(stderr, "normalize() sees vector (0,0)!\n");
        return; // does nothing to zero vectors
    }
    float scaleFactor = 1.0 / sqrt(sizeSq);
    v->x *= scaleFactor;
    v->y *= scaleFactor;
}

float dotVector2(Vector2* v, Vector2* b) {
    return v->x * b->x + v->y * b->y;
}

void perpVector2(Vector2* v) {
    float tmp = v->x;
    v->x = -v->y;
    v->y = tmp;
}

float perpDotVector2(Vector2* v, Vector2* b) {
    return v->x * b->x - v->y * b->y;
}

// Canvas functions
void initCanvas(Canvas* canvas, int width, int height, char* title) {
    // Initialization code for canvas
}

void setCanvasWindow(Canvas* canvas, float l, float r, float b, float t) {
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(l, r, b, t);
    if (t == b) return;
    canvas->windowAspect = (r - l) / (t - b);
}

void setCanvasViewport(Canvas* canvas, int l, int r, int b, int t) {
    glViewport(l, b, r - l, t - b);
}

float getCanvasWindowAspect(Canvas* canvas) {
    return canvas->windowAspect;
}

void canvasLineTo(Canvas* canvas, float x, float y) {
    glBegin(GL_LINES);
    glVertex2f(canvas->CP.x, canvas->CP.y);
    canvas->CP.x = x;
    canvas->CP.y = y;
    glVertex2f(canvas->CP.x, canvas->CP.y);
    glEnd();
    glFlush();
}

void canvasMoveTo(Canvas* canvas, float x, float y) {
    canvas->CP.x = x;
    canvas->CP.y = y;
}

void canvasForward(Canvas* canvas, float dist, int vis) {
    float RadPerDeg = 0.017453393; // radians per degree
    float x = canvas->CP.x + dist * cos(RadPerDeg * dist);
    float y = canvas->CP.y + dist * sin(RadPerDeg * dist);
    if (vis) canvasLineTo(canvas, x, y);
    else canvasMoveTo(canvas, x, y);
    canvas->CP.x = x;
    canvas->CP.y = y;
}

void canvasInitCT(Canvas* canvas) {
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

void canvasRotate2D(Canvas* canvas, double angle) {
    glMatrixMode(GL_MODELVIEW);
    glRotated(angle, 0.0, 0.0, 1.0);
}

void canvasTranslate2D(Canvas* canvas, double dx, double dy) {
    glMatrixMode(GL_MODELVIEW);
    glTranslated(dx, dy, 0.0);
}

void canvasScale2D(Canvas* canvas, double sx, double sy) {
    glMatrixMode(GL_MODELVIEW);
    glScaled(sx, sy, 1.0);
}

void canvasPushCT(Canvas* canvas) {
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
}

void canvasPopCT(Canvas* canvas) {
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
}

void canvasNgon(Canvas* canvas, int n, float cx, float cy, float radius) {
    float RadPerDeg = 0.017453393; // radians per degree
    if (n < 3) return; // bad number of sides
    double angle = 0, angleInc = 2 * 3.14159265 / n; // angle increment
    canvasMoveTo(canvas, cx + radius, cy);
    for (int k = 1; k <= n; k++) {
        angle += angleInc;
        canvasLineTo(canvas, radius * cos(angle) + cx, radius * sin(angle) + cy);
    }
}

// RGBpixmap functions
int readBMPFile(RGBpixmap* pixmap, char* fname) {
    FILE* inf = fopen(fname, "rb");
    if (!inf) {
        printf("can't open file: %s\n", fname);
        return 0;
    }
    char ch1, ch2;
    fread(&ch1, 1, 1, inf);
    fread(&ch2, 1, 1, inf);
    if (ch1 != 'B' || ch2 != 'M') {
        fclose(inf);
        return 0;
    }
    unsigned long fileSize = getLong(inf);
    unsigned short reserved1 = getShort(inf);
    unsigned short reserved2 = getShort(inf);
    unsigned long offBits = getLong(inf);
    unsigned long headerSize = getLong(inf);
    unsigned long numCols = getLong(inf);
    unsigned long numRows = getLong(inf);
    unsigned short planes = getShort(inf);
    unsigned short bitsPerPixel = getShort(inf);
    unsigned long compression = getLong(inf);
    unsigned long imageSize = getLong(inf);
    unsigned long xPels = getLong(inf);
    unsigned long yPels = getLong(inf);
    unsigned long numLUTentries = getLong(inf);
    unsigned long impColors = getLong(inf);
    if (bitsPerPixel != 24) {
        printf("not a 24-bit pixel image, or is compressed!\n");
        fclose(inf);
        return 0;
    }
    pixmap->nRows = numRows;
    pixmap->nCols = numCols;
    pixmap->pixel = (mRGB*)malloc(pixmap->nRows * pixmap->nCols * sizeof(mRGB));
    if (!pixmap->pixel) return 0;
    long count = 0;
    char dum;
    for (int row = 0; row < pixmap->nRows; row++) {
        for (int col = 0; col < pixmap->nCols; col++) {
            unsigned char r, g, b;
            fread(&b, 1, 1, inf);
            fread(&g, 1, 1, inf);
            fread(&r, 1, 1, inf);
            pixmap->pixel[count].r = r;
            pixmap->pixel[count].g = g;
            pixmap->pixel[count].b = b;
            count++;
        }
        for (int k = 0; k < (4 - (pixmap->nCols * 3) % 4) % 4; k++) {
            fread(&dum, 1, 1, inf);
        }
    }
    fclose(inf);
    return 1;
}

void freeRGBpixmap(RGBpixmap* pixmap) {
    free(pixmap->pixel);
    pixmap->nRows = 0;
    pixmap->nCols = 0;
}

void copyRGBpixmap(RGBpixmap* pixmap, IntPoint from, IntPoint to, int width, int height) {
    if (pixmap->nRows == 0 || pixmap->nCols == 0) return;
    glCopyPixels(from.x, from.y, width, height, GL_COLOR);
}

void drawRGBpixmap(RGBpixmap* pixmap) {
    if (pixmap->nRows == 0 || pixmap->nCols == 0) return;
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glDrawPixels(pixmap->nCols, pixmap->nRows, GL_RGB, GL_UNSIGNED_BYTE, pixmap->pixel);
}

int readRGBpixmap(RGBpixmap* pixmap, int x, int y, int width, int height) {
    pixmap->nRows = height;
    pixmap->nCols = width;
    pixmap->pixel = (mRGB*)malloc(pixmap->nRows * pixmap->nCols * sizeof(mRGB));
    if (!pixmap->pixel) return -1;
    glPixelStorei(GL_PACK_ALIGNMENT, 1);
    glReadPixels(x, y, pixmap->nCols, pixmap->nRows, GL_RGB, GL_UNSIGNED_BYTE, pixmap->pixel);
    return 0;
}

int readRGBpixmapFromRect(RGBpixmap* pixmap, IntRect r) {
    return readRGBpixmap(pixmap, r.left, r.top, r.right - r.left, r.top - r.bott);
}

void setRGBpixmapPixel(RGBpixmap* pixmap, int x, int y, mRGB color) {
    if (x >= 0 && x < pixmap->nCols && y >= 0 && y < pixmap->nRows) {
        pixmap->pixel[pixmap->nCols * y + x] = color;
    }
}

mRGB getRGBpixmapPixel(RGBpixmap* pixmap, int x, int y) {
    mRGB bad = {255, 255, 255};
    if (x < 0 || x >= pixmap->nCols || y < 0 || y >= pixmap->nRows) return bad;
    return pixmap->pixel[pixmap->nCols * y + x];
}

unsigned short getShort(FILE* fp) {
    char ic;
    unsigned short ip;
    fread(&ic, 1, 1, fp); ip = ic;
    fread(&ic, 1, 1, fp); ip |= ((unsigned short)ic << 8);
    return ip;
}

unsigned long getLong(FILE* fp) {
    unsigned long ip = 0;
    unsigned char uc = 0;
    for (int i = 0; i < 4; i++) {
        fread(&uc, 1, 1, fp);
        ip |= ((unsigned long)uc << (i * 8));
    }
    return ip;
}
