//
// Created by jorda on 7/16/2024.
//

#ifndef _GRAPHICS2D
#define _GRAPHICS2D

#include <string>
#include <iostream>
#include <fstream>

using namespace std;

class IntPoint {
        public:
        int x, y;
        void set(int dx, int dy) { x = dx; y = dy; }
        void set(IntPoint& p) { x = p.x; y = p.y; }
        IntPoint(int xx, int yy) { x = xx; y = yy; }
        IntPoint() { x = y = 0; }
};

class Point2 {
        public:
        float x, y;
        void set(float dx, float dy) { x = dx; y = dy; }
        void set(Point2& p) { x = p.x; y = p.y; }
        Point2(float xx, float yy) { x = xx; y = yy; }
        Point2() { x = y = 0; }
};

class IntRect {
        public:
        int left, top, right, bott;
        IntRect() { left = top = right = bott = 0; }
        IntRect(int l, int t, int r, int b) { left = l; top = t; right = r; bott = b; }
        void set(int l, int t, int r, int b) { left = l; top = t; right = r; bott = b; }
        void set(IntRect& r) { left = r.left; top = r.top; right = r.right; bott = r.bott; }
};

class Vector2 {
        public:
        float x, y;
        void set(float dx, float dy) { x = dx; y = dy; }
        void set(Vector2& v) { x = v.x; y = v.y; }
        void setDiff(Point2& a, Point2& b) { x = a.x - b.x; y = a.y - b.y; }
        void normalize(); // adjust this vector to unit length
        float dot(Vector2 b) { return x * b.x + y * b.y; }
        void perp() { float tmp = x; x = -y; y = tmp; }
        float perpDot(Vector2& v) { return x * v.x - y * v.y; }

        Vector2(float xx, float yy) { x = xx; y = yy; }
        Vector2(Vector2& v) { x = v.x; y = v.y; }
        Vector2() { x = y = 0; }
};

class Canvas {
        private:
        Point2 CP; // current position in world
        float windowAspect;

        public:
        Canvas(int width, int height, char* title);
        void setWindow(float l, float r, float b, float t);
        void setViewport(int l, int r, int b, int t);
        float getWindowAspect() { return windowAspect; }
        void lineTo(float x, float y);
        void moveTo(float x, float y) { CP.x = x; CP.y = y; }
        void forward(float dist, int vis);
        void initCT(); // initialize the CT (model view matrix)
        void rotate2D(double angle);
        void translate2D(double dx, double dy);
        void scale2D(double sx, double sy);
        void pushCT();
        void popCT();
        void ngon(int n, float cx, float cy, float radius);
};

class mRGB {
        public:
        unsigned char r, g, b;
        mRGB(unsigned char rr = 0, unsigned char gg = 0, unsigned char bb = 0) : r(rr), g(gg), b(bb) {}
        void set(unsigned char rr, unsigned char gg, unsigned char bb) { r = rr; g = gg; b = bb; }
};

class RGBpixmap {
        private:
        mRGB* pixel; // array of pixels
        public:
        int nRows, nCols; // dimensions of the pix map
        RGBpixmap() { nRows = nCols = 0; pixel = 0; }
        RGBpixmap(int rows, int cols) {
            nRows = rows;
            nCols = cols;
            pixel = new mRGB[rows * cols];
        }
        ~RGBpixmap() { delete[] pixel; }

        int readBMPFile(string fname);
        void freeIt() { delete[] pixel; nRows = nCols = 0; }
        void copy(IntPoint from, IntPoint to, int width, int height);
        void draw();
        int read(int x, int y, int width, int height);
        int read(IntRect r);
        void setPixel(int x, int y, mRGB color);
        mRGB getPixel(int x, int y);
};

class PolyLine {
        public:
        int num;
        Point2 pt[80]; // may need larger arrays in some circumstances
        PolyLine() { num = 0; }
};

#endif // _GRAPHICS2D
