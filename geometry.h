// geometry.h

#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <math.h>
#include <float.h>
#include <stdbool.h>

#define CMP(x, y) \
   (fabsf((x) - (y)) <= FLT_EPSILON * \
   fmaxf(1.0f, fmaxf(fabsf(x), fabsf(y))))

#define RectangleCircle(rectangle, circle) CircleRectangle(circle, rectangle)

#define OrientedRectangleCircle(rectangle, circle) CircleOrientedRectangle(circle, rectangle)

#define RAD2DEG(x) ((x) * 57.295754f)

#define DEG2RAD(x) ((x) * 0.0174533f)

#define OVERLAP(aMin, aMax, bMin, bMax) \
   ((bMin <= aMax) && (aMin <= bMax))

#define OrientedRectangleRectangle(oriented, regular) RectangleOrientedRectangle(regular, oriented)

typedef struct {
    float x, y, z;
} Point, Vector;

typedef struct {
    Point P0, P1;
} Line, Segment;

typedef struct {
    Segment position;
    float radius;
} Circle;

typedef struct {
    Point origin;
    Vector size;
} Rectangle2D;

typedef struct {
    Point position;
    Vector halfExtents;
    float rotation;
} OrientedRectangle;

typedef struct {
    float min;
    float max;
} Interval2D;

typedef struct {
    Point V0;
    Vector n;
} Plane;

#define dot(u,v)   ((u).x * (v).x + (u).y * (v).y + (u).z * (v).z)
#define norm(v)    sqrt(dot(v,v))   // norm = length of vector
#define d(u,v)     norm(sub(u,v))   // distance = norm of difference

int isLeft(Point P0, Point P1, Point P2);
int orientation2D_Triangle(Point V0, Point V1, Point V2);
float area2D_Triangle(Point V0, Point V1, Point V2);
int orientation2D_Polygon(int n, Point* V);
float area2D_Polygon(int n, Point* V);
float area3D_Polygon(int n, Point* V, Point N);

Vector sub(Vector u, Vector v);
Vector add(Vector u, Vector v);
Vector scalar_multiply(float s, Vector v);

int closest2D_Point_to_Line(Point P[], int n, Line L);
float dist_Point_to_Line(Point P, Line L);
float dist_Point_to_Segment(Point P, Segment S);

int cn_PnPoly(Point P, Point* V, int n);
int wn_PnPoly(Point P, Point* V, int n);

Vector subtract(Point a, Point b);
Point add_point(Point a, Vector b);
Vector scale(Vector v, float scalar);
float dist_Point_to_Plane(Point P, Plane PL, Point *B);
void initRectangle2D(Rectangle2D* rect);
void initRectangle2DWithValues(Rectangle2D* rect, const Point * o, const Vector * s);

void initOrientedRectangle(OrientedRectangle* rect);
void initOrientedRectangleWithParams(OrientedRectangle* rect, const Point * position, const Vector * halfExtents, float rotation);

float lengthSq(const Line * line);
void initLine2D(Line * line, const Point* start, const Point* end);

bool PointOnLine(const Point* point, const Line* line);
bool PointInCircle(const Point* point, const Circle* circle);
bool PointInRectangle(const Point* point, const Rectangle2D* rectangle);
bool PointInOrientedRectangle(const Point* point, const OrientedRectangle* rectangle);

bool CircleCircle(const Circle* c1, const Circle* c2);
bool CircleRectangle(const Circle* circle, const Rectangle2D* rectangle);

Point GetMin(const Rectangle2D* rect);
Point GetMax(const Rectangle2D* rect);

bool CircleOrientedRectangle(const Circle* circle, const OrientedRectangle* rect);

bool RectangleRectangle(const Rectangle2D* rect1, const Rectangle2D* rect2);

Interval2D GetInterval(const OrientedRectangle* rect, const Vector* axis);
bool OverlapOnAxis(const Rectangle2D* rect1, const OrientedRectangle* rect2, const Vector* axis);
bool RectangleOrientedRectangle(const Rectangle2D* rect1, const OrientedRectangle* rect2);

bool OrientedRectangleOrientedRectangle(const OrientedRectangle* r1, const OrientedRectangle* r2);

#endif // GEOMETRY_H
