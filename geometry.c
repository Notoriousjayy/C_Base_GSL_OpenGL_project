// geometry.c

#include "geometry.h"
#include <stdio.h>

int isLeft(Point P0, Point P1, Point P2) {
    return ((P1.x - P0.x) * (P2.y - P0.y) - (P2.x - P0.x) * (P1.y - P0.y));
}

int orientation2D_Triangle(Point V0, Point V1, Point V2) {
    return isLeft(V0, V1, V2);
}

float area2D_Triangle(Point V0, Point V1, Point V2) {
    return (float)isLeft(V0, V1, V2) / 2.0;
}

int orientation2D_Polygon(int n, Point* V) {
    int rmin = 0;
    int xmin = V[0].x;
    int ymin = V[0].y;

    for (int i = 1; i < n; i++) {
        if (V[i].y > ymin)
            continue;
        if (V[i].y == ymin) { // just as low
            if (V[i].x < xmin) // and to left
                continue;
        }
        rmin = i; // a new rightmost lowest vertex
        xmin = V[i].x;
        ymin = V[i].y;
    }

    if (rmin == 0)
        return isLeft(V[n-1], V[0], V[1]);
    else
        return isLeft(V[rmin-1], V[rmin], V[rmin+1]);
}

float area2D_Polygon(int n, Point* V) {
    float area = 0;
    int i, j, k;

    if (n < 3) return 0;

    for (i = 1, j = 2, k = 0; i < n; i++, j++, k++)
        area += V[i].x * (V[j].y - V[k].y);

    area += V[n].x * (V[1].y - V[n-1].y);
    return area / 2.0;
}

float area3D_Polygon(int n, Point* V, Point N) {
    float area = 0;
    float an, ax, ay, az;
    int coord;
    int i, j, k;

    if (n < 3) return 0;

    ax = (N.x > 0.0 ? N.x : -N.x);
    ay = (N.y > 0.0 ? N.y : -N.y);
    az = (N.z > 0.0 ? N.z : -N.z);

    coord = 3;
    if (ax > ay) {
        if (ax > az)
            coord = 1;
    } else if (ay > az)
        coord = 2;

    switch (coord) {
        case 1:
            for (i = 1, j = 2, k = 0; i < n; i++, j++, k++)
                area += V[i].y * (V[j].z - V[k].z);
            break;
        case 2:
            for (i = 1, j = 2, k = 0; i < n; i++, j++, k++)
                area += V[i].z * (V[j].x - V[k].x);
            break;
        case 3:
            for (i = 1, j = 2, k = 0; i < n; i++, j++, k++)
                area += V[i].x * (V[j].y - V[k].y);
            break;
    }

    switch (coord) {
        case 1:
            area += V[n].y * (V[1].z - V[n-1].z);
            break;
        case 2:
            area += V[n].z * (V[1].x - V[n-1].x);
            break;
        case 3:
            area += V[n].x * (V[1].y - V[n-1].y);
            break;
    }

    an = sqrt(ax*ax + ay*ay + az*az);
    switch (coord) {
        case 1:
            area *= (an / (2 * N.x));
            break;
        case 2:
            area *= (an / (2 * N.y));
            break;
        case 3:
            area *= (an / (2 * N.z));
            break;
    }
    return area;
}

Vector sub(Vector u, Vector v) {
    Vector result = { u.x - v.x, u.y - v.y, u.z - v.z };
    return result;
}

Vector add(Vector u, Vector v) {
    Vector result = { u.x + v.x, u.y + v.y, u.z + v.z };
    return result;
}

Vector scalar_multiply(float s, Vector v) {
    Vector result = { s * v.x, s * v.y, s * v.z };
    return result;
}

int closest2D_Point_to_Line(Point P[], int n, Line L) {
    float a = L.P0.y - L.P1.y;
    float b = L.P1.x - L.P0.x;
    float c = L.P0.x * L.P1.y - L.P1.x * L.P0.y;

    int mi = 0;
    float min = a * P[0].x + b * P[0].y + c;
    if (min < 0) min = -min;

    for (int i = 1; i < n; i++) {
        float dist = a * P[i].x + b * P[i].y + c;
        if (dist < 0) dist = -dist;
        if (dist < min) {
            mi = i;
            min = dist;
        }
    }
    return mi;
}

float dist_Point_to_Line(Point P, Line L) {
    Vector v = sub(L.P1, L.P0);
    Vector w = sub(P, L.P0);

    double c1 = dot(w, v);
    double c2 = dot(v, v);
    double b = c1 / c2;

    Point Pb = add(L.P0, scalar_multiply(b, v));
    return d(P, Pb);
}

float dist_Point_to_Segment(Point P, Segment S) {
    Vector v = sub(S.P1, S.P0);
    Vector w = sub(P, S.P0);

    double c1 = dot(w, v);
    if (c1 <= 0) {
        return d(P, S.P0);
    }

    double c2 = dot(v, v);
    if (c2 <= c1) {
        return d(P, S.P1);
    }

    double b = c1 / c2;
    Point Pb = add(S.P0, scalar_multiply(b, v));
    return d(P, Pb);
}

int cn_PnPoly(Point P, Point* V, int n) {
    int cn = 0;    // the crossing number counter

    // loop through all edges of the polygon
    for (int i = 0; i < n; i++) {    // edge from V[i] to V[i+1]
        if (((V[i].y <= P.y) && (V[i+1].y > P.y))     // an upward crossing
            || ((V[i].y > P.y) && (V[i+1].y <= P.y))) {  // a downward crossing
            // compute the actual edge-ray intersect x-coordinate
            float vt = (float)(P.y - V[i].y) / (V[i+1].y - V[i].y);
            if (P.x < V[i].x + vt * (V[i+1].x - V[i].x)) // P.x < intersect
                ++cn;    // a valid crossing of y=P.y right of P.x
        }
    }
    return (cn & 1);    // 0 if even (out), and 1 if odd (in)
}

int wn_PnPoly(Point P, Point* V, int n) {
    int wn = 0;    // the winding number counter

    // loop through all edges of the polygon
    for (int i = 0; i < n; i++) {    // edge from V[i] to V[i+1]
        if (V[i].y <= P.y) {         // start y <= P.y
            if (V[i+1].y > P.y)      // an upward crossing
                if (isLeft(V[i], V[i+1], P) > 0) // P left of edge
                    ++wn;             // have a valid up intersect
        }
        else {                       // start y > P.y (no test needed)
            if (V[i+1].y <= P.y)     // a downward crossing
                if (isLeft(V[i], V[i+1], P) < 0) // P right of edge
                    --wn;             // have a valid down intersect
        }
    }
    return wn;
}

// Function to subtract two points or vectors
Vector subtract(Point a, Point b) {
    Vector result;
    result.x = a.x - b.x;
    result.y = a.y - b.y;
    result.z = a.z - b.z;
    return result;
}

// Function to add a point and a vector
Point add_point(Point a, Vector b) {
    Point result;
    result.x = a.x + b.x;
    result.y = a.y + b.y;
    result.z = a.z + b.z;
    return result;
}

// Function to scale a vector by a scalar
Vector scale(Vector v, float scalar) {
    Vector result;
    result.x = v.x * scalar;
    result.y = v.y * scalar;
    result.z = v.z * scalar;
    return result;
}

// dist_Point_to_Plane(): get the distance from a point to a plane
//     Input:  P  = a 3D point
//             PL = a plane with point V0 and normal n
//     Output: *B = base point on PL of perpendicular from P
//     Return: the distance from P to the plane PL
float dist_Point_to_Plane(Point P, Plane PL, Point *B) {
    float sb, sn, sd;
    Vector diff = subtract(P, PL.V0);

    sn = -dot(PL.n, diff);
    sd = dot(PL.n, PL.n);
    sb = sn / sd;

    Vector scaled_n = scale(PL.n, sb);
    *B = add_point(P, scaled_n);
    return d(P, *B);
}

// Initialize a Rectangle2D with default values
void initRectangle2D(Rectangle2D* rect) {
    rect->origin.x = 0;
    rect->origin.y = 0;
    rect->size.x = 1;
    rect->size.y = 1;
}

// Initialize a Rectangle2D with specified values
void initRectangle2DWithValiues(Rectangle2D* rect, const Point * o, const Vector * s) {
    rect->origin.x = o->x;
    rect->origin.y = o->y;
    rect->size.x = s->x;
    rect->size.y = s->y;
}

void initOrientedRectangle(OrientedRectangle* rect) {
    rect->position.x = 0.0f;
    rect->position.y = 0.0f;
    rect->halfExtents.x = 1.0f;
    rect->halfExtents.y = 1.0f;
    rect->rotation = 0.0f;
}

void initOrientedRectangleWithParams(OrientedRectangle* rect, const Point * position, const Vector * halfExtents, float rotation) {
    rect->position.x = position->x;
    rect->position.y = position->y;
    rect->halfExtents.x = halfExtents->x;
    rect->halfExtents.y = halfExtents->y;
    rect->rotation = rotation;
}

float lengthSq(const Line* line) {
    float dx = line->P1.x - line->P0.x;
    float dy = line->P1.y - line->P0.y;
    return dx * dx + dy * dy;
}

void initLine(Line* line, const Point* start, const Point* end) {
    line->P0.x = start->x;
    line->P0.y = start->y;
    line->P1.x = end->x;
    line->P1.y = end->y;
}

bool PointOnLine(const Point* point, const Line* line) {
    float dx = line->P1.x - line->P0.x;
    float dy = line->P1.y - line->P0.y;
    float dx1 = point->x - line->P0.x;
    float dy1 = point->y - line->P0.y;
    float cross = dx * dy1 - dy * dx1;
    return fabs(cross) < 1e-6;
}

bool PointInCircle(const Point* point, const Circle* circle) {
    Line line;
    initLine(&line, point, &circle->position);
    if (lengthSq(&line) < circle->radius * circle->radius) {
        return true;
    }
    return false;
}

bool PointInRectangle(const Point* point, const Rectangle2D* rectangle) {
    if (point->x < rectangle->origin.x ||
        point->x > rectangle->origin.x + rectangle->size.x ||
        point->y < rectangle->origin.y ||
        point->y > rectangle->origin.y + rectangle->size.y) {
        return false;
    }
    return true;
}

bool PointInOrientedRectangle(const Point* point, const OrientedRectangle* rectangle) {
    // Translate point back to origin
    float cos_theta = cos(-rectangle->rotation);
    float sin_theta = sin(-rectangle->rotation);

    float translated_x = point->x - rectangle->position.x;
    float translated_y = point->y - rectangle->position.y;

    // Rotate point
    float rotated_x = translated_x * cos_theta - translated_y * sin_theta;
    float rotated_y = translated_x * sin_theta + translated_y * cos_theta;

    // Check if point is within half extents
    if (fabs(rotated_x) > rectangle->halfExtents.x || fabs(rotated_y) > rectangle->halfExtents.y) {
        return false;
    }
    return true;
}

bool CircleCircle(const Circle *c1, const Circle *c2) {
    Line line;
    initLine(&line, &c1->position.P0, &c2->position.P0);
    float radiiSum = c1->radius + c2->radius;
    return lengthSq(&line) <= radiiSum * radiiSum;
}

Point GetMin(const Rectangle2D* rect) {
    Point min;
    min.x = rect->origin.x;
    min.y = rect->origin.y;
    return min;
}

Point GetMax(const Rectangle2D* rect) {
    Point max;
    max.x = rect->origin.x + rect->size.x;
    max.y = rect->origin.y + rect->size.y;
    return max;
}

bool CircleRectangle(const Circle* circle, const Rectangle2D* rect) {
    Point min = GetMin(rect);
    Point max = GetMax(rect);
    Point closestPoint = circle->position.P0;
    if (closestPoint.x < min.x) {
        closestPoint.x = min.x;
    } else if (closestPoint.x > max.x) {
        closestPoint.x = max.x;
    }
    if (closestPoint.y < min.y) {
        closestPoint.y = min.y;
    } else if (closestPoint.y > max.y) {
        closestPoint.y = max.y;
    }

    Line line;
    initLine(&line, &circle->position.P0, &closestPoint);
    return lengthSq(&line) <= circle->radius * circle->radius;
}

void initRectangle2DWithValues(Rectangle2D* rect, const Point* origin, const Vector* size) {
    rect->origin = *origin;
    rect->size = *size;
}

Vector vector_sub(Vector a, Vector b) {
    Vector result = {a.x - b.x, a.y - b.y};
    return result;
}
Vector vector_add(Vector a, Vector b) {
    Vector result = {a.x + b.x, a.y + b.y};
    return result;
}

Vector vector_multiply(Vector v, float s) {
    Vector result = {v.x * s, v.y * s};
    return result;
}

void mat2_multiply(float* result, const float* mat, const float* vec) {
    result[0] = mat[0] * vec[0] + mat[1] * vec[1];
    result[1] = mat[2] * vec[0] + mat[3] * vec[1];
}

float dot2(Vector a, Vector b) {
    return a.x * b.x + a.y * b.y;
}


bool CircleOrientedRectangle(const Circle* circle, const OrientedRectangle* rect) {
    Vector r = vector_sub(*(Vector *)&circle->position.P0, *(Vector*)&rect->position);
    float theta = -DEG2RAD(rect->rotation);
    float zRotation2x2[4] = {
            cosf(theta), sinf(theta),
            -sinf(theta), cosf(theta)
    };

    Vector rotatedR;
    mat2_multiply((float*)&rotatedR, zRotation2x2, (float*)&r);
    Vector lCirclePos = vector_add(rotatedR, *(Vector *)&rect->halfExtents);
    Circle lCircle = {{{lCirclePos.x, lCirclePos.y, 0.0f}, {lCirclePos.x, lCirclePos.y, 0.0f}}, circle->radius};
    Rectangle2D lRect = {{0.0f, 0.0f, 0.0f}, {rect->halfExtents.x * 2.0f, rect->halfExtents.y * 2.0f, 0.0f}};

    return CircleRectangle(&lCircle, &lRect);
}

bool RectangleRectangle(const Rectangle2D* rect1, const Rectangle2D* rect2) {
    Vector aMin = GetMin(rect1);
    Vector aMax = GetMax(rect1);
    Vector bMin = GetMin(rect2);
    Vector bMax = GetMax(rect2);
    bool overX = OVERLAP(aMin.x, aMax.x, bMin.x, bMax.x);
    bool overY = OVERLAP(aMin.y, aMax.y, bMin.y, bMax.y);
    return overX && overY;
}

Interval2D GetInterval(const OrientedRectangle* rect, const Vector* axis) {
    Rectangle2D r;
    Point origin = {rect->position.x - rect->halfExtents.x, rect->position.y - rect->halfExtents.y, 0};
    Vector size = {rect->halfExtents.x * 2.0f, rect->halfExtents.y * 2.0f, 0};
    initRectangle2DWithValues(&r, &origin, &size);

    Vector min = GetMin(&r);
    Vector max = GetMax(&r);
    Vector verts[] = { min, max, {min.x, max.y}, {max.x, min.y} };

    float t = DEG2RAD(rect->rotation);
    float zRot[] = { cosf(t), sinf(t), -sinf(t), cosf(t) };

    for (int i = 0; i < 4; ++i) {
        Vector v = vector_sub(verts[i], *(Vector*)&rect->position);
        mat2_multiply((float*)&v, zRot, (float*)&v);
        verts[i] = vector_add(v, *(Vector*)&rect->position);
    }

    Interval2D res;
    res.min = res.max = dot2(*axis, verts[0]);
    for (int i = 1; i < 4; ++i) {
        float proj = dot2(*axis, verts[i]);
        if (proj < res.min) res.min = proj;
        if (proj > res.max) res.max = proj;
    }
    return res;
}

bool OverlapOnAxis(const Rectangle2D* rect1, const OrientedRectangle* rect2, const Vector* axis) {
    Interval2D a = {dot2(*axis, GetMin(rect1)), dot2(*axis, GetMax(rect1))};
    Interval2D b = GetInterval(rect2, axis);
    return OVERLAP(a.min, a.max, b.min, b.max);
}

bool RectangleOrientedRectangle(const Rectangle2D* rect1, const OrientedRectangle* rect2) {
    Vector axisToTest[] = { {1, 0}, {0, 1}, {0, 0}, {0, 0} };

    float t = DEG2RAD(rect2->rotation);
    float zRot[] = { cosf(t), sinf(t), -sinf(t), cosf(t) };

    Vector axis = {rect2->halfExtents.x, 0};
    Vector normalizedAxis;
    mat2_multiply((float*)&normalizedAxis, zRot, (float*)&axis);
    axisToTest[2] = normalizedAxis;

    axis = (Vector){0, rect2->halfExtents.y};
    mat2_multiply((float*)&normalizedAxis, zRot, (float*)&axis);
    axisToTest[3] = normalizedAxis;

    for (int i = 0; i < 4; ++i) {
        if (!OverlapOnAxis(rect1, rect2, &axisToTest[i])) {
            return false; // No collision has taken place
        }
    }
    return true; // We have a collision
}

bool OrientedRectangleOrientedRectangle(const OrientedRectangle* r1, const OrientedRectangle* r2) {
    Rectangle2D local1;
    Vector size = {r1->halfExtents.x * 2.0f, r1->halfExtents.y * 2.0f, 0.0f};
    Point origin = {0.0f, 0.0f, 0.0f};
    initRectangle2DWithValues(&local1, &origin, &size);

    Vector r = vector_sub(*(Vector*)&r2->position, *(Vector*)&r1->position);

    OrientedRectangle local2;
    initOrientedRectangleWithParams(&local2, &r2->position, &r2->halfExtents, r2->rotation);
    local2.rotation = r2->rotation - r1->rotation;

    float t = -DEG2RAD(r1->rotation);
    float z[] = {
            cosf(t), sinf(t),
            -sinf(t), cosf(t)
    };

    Vector rotated_r;
    mat2_multiply((float*)&rotated_r, z, (float*)&r);
    *(Vector*)&local2.position = vector_add(rotated_r, *(Vector*)&r1->halfExtents);

    return RectangleOrientedRectangle(&local1, &local2);
}