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

typedef struct mat2 {
    union {
        struct {
            float _11, _12,
                    _21, _22;
        };
        float asArray[4];
    };
} mat2;

typedef struct mat3 {
    union {
        struct {
            float _11, _12, _13,
                    _21, _22, _23,
                    _31, _32, _33;
        };
        float asArray[9];
    };
} mat3;

typedef struct mat4 {
    union {
        struct {
            float _11, _12, _13, _14,
                    _21, _22, _23, _24,
                    _31, _32, _33, _34,
                    _41, _42, _43, _44;
        };
        float asArray[16];
    };
} mat4;

Vector Vector_add(const Vector* l, const Vector* r);
Vector Vector3_add(const Vector* l, const Vector* r);
Vector Vector_sub(const Vector* l, const Vector* r);
Vector Vector3_sub(const Vector* l, const Vector* r);
Vector Vector_mul(const Vector* l, const Vector* r);
Vector Vector3_mul(const Vector* l, const Vector* r);
Vector Vector_mul_scalar(const Vector* l, float r);
Vector Vector3_mul_scalar(const Vector* l, float r);
bool Vector_eq(const Vector* l, const Vector* r);
bool Vector3_eq(const Vector* l, const Vector* r);
bool Vector_neq(const Vector* l, const Vector* r);
bool Vector3_neq(const Vector* l, const Vector* r);

float dot_Vector(const Vector* l, const Vector* r);
float magnitude_Vector(const Vector* v);
float magnitude_Vector(const Vector* v);
float magnitude_sq_Vector(const Vector* v);
float magnitude_sq_Vector(const Vector* v);
float distance_Vector(const Vector* p1, const Vector* p2);
void normalize_Vector(Vector* v);
void normalize_Vector(Vector* v);
Vector normalized_Vector(const Vector* v);
Vector normalized_Vector(const Vector* v);
Vector cross_Vector(const Vector* l, const Vector* r);
float angle_Vector(const Vector* l, const Vector* r);
float angle_Vector(const Vector* l, const Vector* r);
Vector project_Vector(const Vector* length, const Vector* direction);
Vector project_Vector(const Vector* length, const Vector* direction);
Vector perpendicular_Vector(const Vector* len, const Vector* dir);
Vector perpendicular_Vector(const Vector* len, const Vector* dir);
Vector reflection_Vector(const Vector* vec, const Vector* normal);
Vector reflection_Vector(const Vector* vec, const Vector* normal);

#define dot(u,v)   ((u).x * (v).x + (u).y * (v).y + (u).z * (v).z)
#define norm(v)    sqrt(dot(v,v))   // norm = length of Vectortor
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


void Transpose(const float *srcMat, float *dstMat, int srcRows, int srcCols);
mat2 Transpose2(const mat2* matrix);
mat3 Transpose3(const mat3* matrix);
mat4 Transpose4(const mat4* matrix);

mat2 mat2_scalar_mult(const mat2* matrix, float scalar);
mat3 mat3_scalar_mult(const mat3* matrix, float scalar);
mat4 mat4_scalar_mult(const mat4* matrix, float scalar);

bool Multiply(float* out, const float* matA, int aRows, int aCols, const float* matB, int bRows, int bCols);
mat2 mat2_mult(const mat2* matA, const mat2* matB);
mat3 mat3_mult(const mat3* matA, const mat3* matB);
mat4 mat4_mult(const mat4* matA, const mat4* matB);

float Determinant_mat2(const mat2* matrix);
mat2 Minor_mat2(const mat2* mat);
mat3 Minor_mat3(const mat3* mat);
mat4 Minor_mat4(const mat4* mat);

mat2 Cofactor_mat2(const mat2* mat);
mat3 Cofactor_mat3(const mat3* mat);
mat4 Cofactor_mat4(const mat4* mat);

mat2 Adjugate_mat2(const mat2* mat);
mat3 Adjugate_mat3(const mat3* mat);
mat4 Adjugate_mat4(const mat4* mat);

mat2 Inverse_mat2(const mat2* mat);
mat3 Inverse_mat3(const mat3* mat);
mat4 Inverse_mat4(const mat4* mat);

mat2 mat2_scalar_mult(const mat2* matrix, float scalar);
mat3 mat3_scalar_mult(const mat3* matrix, float scalar);
mat4 mat4_scalar_mult(const mat4* matrix, float scalar);

mat2 mat2_mult(const mat2* matA, const mat2* matB);
mat3 mat3_mult(const mat3* matA, const mat3* matB);
mat4 mat4_mult(const mat4* matA, const mat4* matB);

// Transformation functions
mat4 Translation_xyz(float x, float y, float z);
mat4 Translation_vec(const Vector* pos);
Vector GetTranslation(const mat4* mat);

mat4 Scale_xyz(float x, float y, float z);
mat4 Scale_vec(const Vector* vec);
Vector GetScale(const mat4* mat);

mat4 Rotation(float pitch, float yaw, float roll);
mat3 Rotation3x3(float pitch, float yaw, float roll);

mat4 ZRotation(float angle);
mat3 ZRotation3x3(float angle);
mat4 XRotation(float angle);
mat3 XRotation3x3(float angle);
mat4 YRotation(float angle);
mat3 YRotation3x3(float angle);

mat4 AxisAngle(const Vector* axis, float angle);
mat3 AxisAngle3x3(const Vector* axis, float angle);

Vector MultiplyPoint(const Vector* vec, const mat4* mat);
Vector MultiplyVector(const Vector* vec, const mat4* mat);
Vector MultiplyVector3x3(const Vector* vec, const mat3* mat);

mat4 Transform(const Vector* scale, const Vector* eulerRotation, const Vector* translate);
mat4 TransformAxisAngle(const Vector* scale, const Vector* rotationAxis, float rotationAngle, const Vector* translate);

mat4 LookAt(const Vector* position, const Vector* target, const Vector* up);
mat4 Projection(float fov, float aspect, float zNear, float zFar);
mat4 Ortho(float left, float right, float bottom, float top, float zNear, float zFar);

#endif // GEOMETRY_H
