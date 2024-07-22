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
    initLine(&line, point, &circle->position.P0);
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

Vector Vector_add(const Vector* l, const Vector* r) {
    Vector result = { l->x + r->x, l->y + r->y };
    return result;
}

Vector Vector_sub(const Vector* l, const Vector* r) {
    Vector result = { l->x - r->x, l->y - r->y };
    return result;
}


Vector Vector_mul(const Vector* l, const Vector* r) {
    Vector result = { l->x * r->x, l->y * r->y };
    return result;
}


Vector Vector_mul_scalar(const Vector* l, float r) {
    Vector result = { l->x * r, l->y * r };
    return result;
}

bool Vector_eq(const Vector* l, const Vector* r) {
    return CMP(l->x, r->x) && CMP(l->y, r->y);
}


bool Vector_neq(const Vector* l, const Vector* r) {
    return !Vector_eq(l, r);
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

float dot_Vector(const Vector* l, const Vector* r) {
    return l->x * r->x + l->y * r->y + l->z * r->z;
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

float magnitude_Vector(const Vector* v) {
    return sqrtf(dot_Vector(v, v));
}

float magnitude_sq_Vector(const Vector* v) {
    return dot_Vector(v, v);
}

float distance_Vector(const Vector* p1, const Vector* p2) {
    Vector t = Vector_sub(p1, p2);
    return magnitude_Vector(&t);
}

void normalize_Vector(Vector* v) {
    float mag = magnitude_Vector(v);
    *v = Vector_mul_scalar(v, 1.0f / mag);
}

Vector normalized_Vector(const Vector* v) {
    Vector result = Vector_mul_scalar(v, 1.0f / magnitude_Vector(v));
    return result;
}

Vector cross_Vector(const Vector* l, const Vector* r) {
    Vector result;
    result.x = l->y * r->z - l->z * r->y;
    result.y = l->z * r->x - l->x * r->z;
    result.z = l->x * r->y - l->y * r->x;
    return result;
}

float angle_Vector(const Vector* l, const Vector* r) {
    float m = sqrtf(magnitude_sq_Vector(l) * magnitude_sq_Vector(r));
    return acos(dot_Vector(l, r) / m);
}


Vector project_Vector(const Vector* length, const Vector* direction) {
    float dot = dot_Vector(length, direction);
    float magSq = magnitude_sq_Vector(direction);
    return Vector_mul_scalar(direction, dot / magSq);
}

Vector perpendicular_Vector(const Vector* len, const Vector* dir) {
    Vector proj = project_Vector(len, dir);
    return Vector_sub(len, &proj);
}

Vector reflection_Vector(const Vector* vec, const Vector* normal) {
    float d = dot_Vector(vec, normal);
    Vector temp = Vector_mul_scalar(normal, d * 2.0f);
    Vector result = Vector_sub(vec, &temp);
    return result;
}

// Transpose functions
void Transpose(const float *srcMat, float *dstMat, int srcRows, int srcCols) {
    for (int i = 0; i < srcRows * srcCols; i++) {
        int row = i / srcCols;
        int col = i % srcCols;
        dstMat[i] = srcMat[srcCols * col + row];
    }
}

mat2 Transpose2(const mat2* matrix) {
    mat2 result;
    Transpose(matrix->asArray, result.asArray, 2, 2);
    return result;
}

mat3 Transpose3(const mat3* matrix) {
    mat3 result;
    Transpose(matrix->asArray, result.asArray, 3, 3);
    return result;
}

mat4 Transpose4(const mat4* matrix) {
    mat4 result;
    Transpose(matrix->asArray, result.asArray, 4, 4);
    return result;
}

// Scalar multiplication functions
mat2 mat2_scalar_mult(const mat2* matrix, float scalar) {
    mat2 result;
    for (int i = 0; i < 4; ++i) {
        result.asArray[i] = matrix->asArray[i] * scalar;
    }
    return result;
}

mat3 mat3_scalar_mult(const mat3* matrix, float scalar) {
    mat3 result;
    for (int i = 0; i < 9; ++i) {
        result.asArray[i] = matrix->asArray[i] * scalar;
    }
    return result;
}

mat4 mat4_scalar_mult(const mat4* matrix, float scalar) {
    mat4 result;
    for (int i = 0; i < 16; ++i) {
        result.asArray[i] = matrix->asArray[i] * scalar;
    }
    return result;
}

// Matrix multiplication functions
bool Multiply(float* out, const float* matA, int aRows, int aCols, const float* matB, int bRows, int bCols) {
    if (aCols != bRows) {
        return false;
    }
    for (int i = 0; i < aRows; ++i) {
        for (int j = 0; j < bCols; ++j) {
            out[bCols * i + j] = 0.0f;
            for (int k = 0; k < bRows; ++k) {
                int a = aCols * i + k;
                int b = bCols * k + j;
                out[bCols * i + j] += matA[a] * matB[b];
            }
        }
    }
    return true;
}

mat2 mat2_mult(const mat2* matA, const mat2* matB) {
    mat2 res;
    Multiply(res.asArray, matA->asArray, 2, 2, matB->asArray, 2, 2);
    return res;
}

mat3 mat3_mult(const mat3* matA, const mat3* matB) {
    mat3 res;
    Multiply(res.asArray, matA->asArray, 3, 3, matB->asArray, 3, 3);
    return res;
}

mat4 mat4_mult(const mat4* matA, const mat4* matB) {
    mat4 res;
    Multiply(res.asArray, matA->asArray, 4, 4, matB->asArray, 4, 4);
    return res;
}

// Determinant functions
float Determinant_mat2(const mat2* matrix) {
    return matrix->_11 * matrix->_22 - matrix->_12 * matrix->_21;
}

mat2 Cut_mat3(const mat3* mat, int row, int col) {
    mat2 result;
    int index = 0;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            if (i == row || j == col) {
                continue;
            }
            int target = index++;
            int source = 3 * i + j;
            result.asArray[target] = mat->asArray[source];
        }
    }
    return result;
}

mat2 Minor_mat2(const mat2* mat) {
    return (mat2){
            mat->_22, -mat->_21,
            -mat->_12, mat->_11
    };
}

mat3 Minor_mat3(const mat3* mat) {
    mat3 result;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            mat2 cut = Cut_mat3(mat, i, j);  // Store the result in a variable
            result.asArray[3 * i + j] = Determinant_mat2(&cut);
        }
    }
    return result;
}

void Cofactor(float* out, const float* minor, int rows, int cols) {
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            int t = cols * i + j;
            float sign = powf(-1.0f, i + j);
            out[t] = minor[t] * sign;
        }
    }
}

mat2 Cofactor_mat2(const mat2* mat) {
    mat2 result;
    Cofactor(result.asArray, Minor_mat2(mat).asArray, 2, 2);
    return result;
}

mat3 Cofactor_mat3(const mat3* mat) {
    mat3 result;
    Cofactor(result.asArray, Minor_mat3(mat).asArray, 3, 3);
    return result;
}

float Determinant_mat3(const mat3* mat) {
    float result = 0.0f;
    mat3 cofactor = Cofactor_mat3(mat);
    for (int j = 0; j < 3; ++j) {
        int index = 3 * 0 + j;
        result += mat->asArray[index] * cofactor.asArray[index];
    }
    return result;
}

mat3 Cut_mat4(const mat4* mat, int row, int col) {
    mat3 result;
    int index = 0;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            if (i == row || j == col) {
                continue;
            }
            int target = index++;
            int source = 4 * i + j;
            result.asArray[target] = mat->asArray[source];
        }
    }
    return result;
}

mat4 Minor_mat4(const mat4* mat) {
    mat4 result;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            mat3 cut = Cut_mat4(mat, i, j);  // Store the result in a variable
            result.asArray[4 * i + j] = Determinant_mat3(&cut);
        }
    }
    return result;
}

mat4 Cofactor_mat4(const mat4* mat) {
    mat4 result;
    Cofactor(result.asArray, Minor_mat4(mat).asArray, 4, 4);
    return result;
}

float Determinant_mat4(const mat4* mat) {
    float result = 0.0f;
    mat4 cofactor = Cofactor_mat4(mat);
    for (int j = 0; j < 4; ++j) {
        result += mat->asArray[4 * 0 + j] * cofactor.asArray[4 * 0 + j];
    }
    return result;
}

mat2 Adjugate_mat2(const mat2* mat) {
    mat2 cofactor = Cofactor_mat2(mat);
    return Transpose2(&cofactor);
}

mat3 Adjugate_mat3(const mat3* mat) {
    mat3 cofactor = Cofactor_mat3(mat);
    return Transpose3(&cofactor);
}

mat4 Adjugate_mat4(const mat4* mat) {
    mat4 cofactor = Cofactor_mat4(mat);
    return Transpose4(&cofactor);
}

mat2 Inverse_mat2(const mat2* mat) {
    float det = Determinant_mat2(mat);
    if (CMP(det, 0.0f)) {
        return (mat2){0};
    }
    mat2 adjugate = Adjugate_mat2(mat);
    return mat2_scalar_mult(&adjugate, 1.0f / det);
}

mat3 Inverse_mat3(const mat3* mat) {
    float det = Determinant_mat3(mat);
    if (CMP(det, 0.0f)) {
        return (mat3){0};
    }
    mat3 adjugate = Adjugate_mat3(mat);
    return mat3_scalar_mult(&adjugate, 1.0f / det);
}

mat4 Inverse_mat4(const mat4* mat) {
    float det = Determinant_mat4(mat);
    if (CMP(det, 0.0f)) {
        return (mat4){0};
    }
    mat4 adjugate = Adjugate_mat4(mat);
    return mat4_scalar_mult(&adjugate, 1.0f / det);
}

mat4 Translation_xyz(float x, float y, float z) {
    return (mat4){
            1.0f, 0.0f, 0.0f, 0.0f,
            0.0f, 1.0f, 0.0f, 0.0f,
            0.0f, 0.0f, 1.0f, 0.0f,
            x, y, z, 1.0f
    };
}

mat4 Translation_vec(const Vector* pos) {
    return (mat4){
            1.0f, 0.0f, 0.0f, 0.0f,
            0.0f, 1.0f, 0.0f, 0.0f,
            0.0f, 0.0f, 1.0f, 0.0f,
            pos->x, pos->y, pos->z, 1.0f
    };
}

Vector GetTranslation(const mat4* mat) {
    return (Vector){ mat->_41, mat->_42, mat->_43 };
}

mat4 Scale_xyz(float x, float y, float z) {
    return (mat4){
            x, 0.0f, 0.0f, 0.0f,
            0.0f, y, 0.0f, 0.0f,
            0.0f, 0.0f, z, 0.0f,
            0.0f, 0.0f, 0.0f, 1.0f
    };
}

mat4 Scale_vec(const Vector* vec) {
    return (mat4){
            vec->x, 0.0f, 0.0f, 0.0f,
            0.0f, vec->y, 0.0f, 0.0f,
            0.0f, 0.0f, vec->z, 0.0f,
            0.0f, 0.0f, 0.0f, 1.0f
    };
}

Vector GetScale(const mat4* mat) {
    return (Vector){ mat->_11, mat->_22, mat->_33 };
}

mat4 Rotation(float pitch, float yaw, float roll) {
    mat4 xRot = XRotation(pitch);
    mat4 yRot = YRotation(yaw);
    mat4 zRot = ZRotation(roll);
    mat4 xyRot = mat4_mult(&xRot, &yRot);
    return mat4_mult(&zRot, &xyRot);
}

mat3 Rotation3x3(float pitch, float yaw, float roll) {
    mat3 xRot = XRotation3x3(pitch);
    mat3 yRot = YRotation3x3(yaw);
    mat3 zRot = ZRotation3x3(yaw);
    mat3 xyRot = mat3_mult(&xRot, &yRot);
    return mat3_mult(&zRot, &xyRot);
}

mat4 ZRotation(float angle) {
    angle = DEG2RAD(angle);
    return (mat4){
            cosf(angle), sinf(angle), 0.0f, 0.0f,
            -sinf(angle), cosf(angle), 0.0f, 0.0f,
            0.0f, 0.0f, 1.0f, 0.0f,
            0.0f, 0.0f, 0.0f, 1.0f
    };
}

mat3 ZRotation3x3(float angle) {
    angle = DEG2RAD(angle);
    return (mat3){
            cosf(angle), sinf(angle), 0.0f,
            -sinf(angle), cosf(angle), 0.0f,
            0.0f, 0.0f, 1.0f
    };
}

mat4 XRotation(float angle) {
    angle = DEG2RAD(angle);
    return (mat4){
            1.0f, 0.0f, 0.0f, 0.0f,
            0.0f, cosf(angle), sinf(angle), 0.0f,
            0.0f, -sinf(angle), cosf(angle), 0.0f,
            0.0f, 0.0f, 0.0f, 1.0f
    };
}

mat3 XRotation3x3(float angle) {
    angle = DEG2RAD(angle);
    return (mat3){
            1.0f, 0.0f, 0.0f,
            0.0f, cosf(angle), sinf(angle),
            0.0f, -sinf(angle), cosf(angle)
    };
}

mat4 YRotation(float angle) {
    angle = DEG2RAD(angle);
    return (mat4){
            cosf(angle), 0.0f, -sinf(angle), 0.0f,
            0.0f, 1.0f, 0.0f, 0.0f,
            sinf(angle), 0.0f, cosf(angle), 0.0f,
            0.0f, 0.0f, 0.0f, 1.0f
    };
}

mat3 YRotation3x3(float angle) {
    angle = DEG2RAD(angle);
    return (mat3){
            cosf(angle), 0.0f, -sinf(angle),
            0.0f, 1.0f, 0.0f,
            sinf(angle), 0.0f, cosf(angle)
    };
}

mat4 AxisAngle(const Vector* axis, float angle) {
    angle = DEG2RAD(angle);
    float c = cosf(angle);
    float s = sinf(angle);
    float t = 1.0f - cosf(angle);

    float x = axis->x;
    float y = axis->y;
    float z = axis->z;
    if (!CMP((x * x + y * y + z * z), 1.0f)) {
        float inv_len = 1.0f / sqrtf(x * x + y * y + z * z);
        x *= inv_len;
        y *= inv_len;
        z *= inv_len;
    }

    return (mat4){
            t * (x * x) + c, t * x * y + s * z, t * x * z - s * y, 0.0f,
            t * x * y - s * z, t * (y * y) + c, t * y * z + s * x, 0.0f,
            t * x * z + s * y, t * y * z - s * x, t * (z * z) + c, 0.0f,
            0.0f, 0.0f, 0.0f, 1.0f
    };
}

mat3 AxisAngle3x3(const Vector* axis, float angle) {
    angle = DEG2RAD(angle);
    float c = cosf(angle);
    float s = sinf(angle);
    float t = 1.0f - cosf(angle);

    float x = axis->x;
    float y = axis->y;
    float z = axis->z;
    if (!CMP((x * x + y * y + z * z), 1.0f)) {
        float inv_len = 1.0f / sqrtf(x * x + y * y + z * z);
        x *= inv_len;
        y *= inv_len;
        z *= inv_len;
    }

    return (mat3){
            t * (x * x) + c, t * x * y + s * z, t * x * z - s * y,
            t * x * y - s * z, t * (y * y) + c, t * y * z + s * x,
            t * x * z + s * y, t * y * z - s * x, t * (z * z) + c
    };
}

Vector MultiplyPoint(const Vector* vec, const mat4* mat) {
    return (Vector){
            vec->x * mat->_11 + vec->y * mat->_21 + vec->z * mat->_31 + 1.0f * mat->_41,
            vec->x * mat->_12 + vec->y * mat->_22 + vec->z * mat->_32 + 1.0f * mat->_42,
            vec->x * mat->_13 + vec->y * mat->_23 + vec->z * mat->_33 + 1.0f * mat->_43
    };
}

Vector MultiplyVector(const Vector* vec, const mat4* mat) {
    return (Vector){
            vec->x * mat->_11 + vec->y * mat->_21 + vec->z * mat->_31 + 0.0f * mat->_41,
            vec->x * mat->_12 + vec->y * mat->_22 + vec->z * mat->_32 + 0.0f * mat->_42,
            vec->x * mat->_13 + vec->y * mat->_23 + vec->z * mat->_33 + 0.0f * mat->_43
    };
}

Vector MultiplyVector3x3(const Vector* vec, const mat3* mat) {
    return (Vector){
            vec->x * mat->_11 + vec->y * mat->_21 + vec->z * mat->_31,
            vec->x * mat->_12 + vec->y * mat->_22 + vec->z * mat->_32,
            vec->x * mat->_13 + vec->y * mat->_23 + vec->z * mat->_33
    };
}

mat4 Transform(const Vector* scale, const Vector* eulerRotation, const Vector* translate) {
    mat4 scaled = Scale_vec(scale);
    mat4 rotated = Rotation(eulerRotation->x, eulerRotation->y, eulerRotation->z);
    mat4 translated = Translation_vec(translate);
    mat4 temp = mat4_mult(&scaled, &rotated);
    return mat4_mult(&temp, &translated);
}

mat4 TransformAxisAngle(const Vector* scale, const Vector* rotationAxis, float rotationAngle, const Vector* translate) {
    mat4 scaled = Scale_vec(scale);
    mat4 rotated = AxisAngle(rotationAxis, rotationAngle);
    mat4 translated = Translation_vec(translate);
    mat4 temp = mat4_mult(&scaled, &rotated);
    return mat4_mult(&temp, &translated);
}

mat4 LookAt(const Vector* position, const Vector* target, const Vector* up) {
    Vector forwardVec = Vector_sub(target, position);
    Vector forward = normalized_Vector(&forwardVec);

    Vector rightVec = cross_Vector(up, &forward);
    Vector right = normalized_Vector(&rightVec);

    Vector newUp = cross_Vector(&forward, &right);

    return (mat4){
            right.x, newUp.x, forward.x, 0.0f,
            right.y, newUp.y, forward.y, 0.0f,
            right.z, newUp.z, forward.z, 0.0f,
            -dot_Vector(&right, position),
            -dot_Vector(&newUp, position),
            -dot_Vector(&forward, position), 1.0f
    };
}


mat4 Projection(float fov, float aspect, float zNear, float zFar) {
    float tanHalfFov = tanf(DEG2RAD((fov * 0.5f)));
    float fovY = 1.0f / tanHalfFov;
    float fovX = fovY / aspect;

    mat4 result = {0};
    result._11 = fovX;
    result._22 = fovY;
    result._33 = zFar / (zFar - zNear);
    result._34 = 1.0f;
    result._43 = -zNear * result._33;
    return result;
}

mat4 Ortho(float left, float right, float bottom, float top, float zNear, float zFar) {
    float _11 = 2.0f / (right - left);
    float _22 = 2.0f / (top - bottom);
    float _33 = 1.0f / (zFar - zNear);
    float _41 = (left + right) / (left - right);
    float _42 = (top + bottom) / (bottom - top);
    float _43 = zNear / (zNear - zFar);

    return (mat4){
            _11, 0.0f, 0.0f, 0.0f,
            0.0f, _22, 0.0f, 0.0f,
            0.0f, 0.0f, _33, 0.0f,
            _41, _42, _43, 1.0f
    };
}