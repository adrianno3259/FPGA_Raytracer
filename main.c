#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>

#define db_color(A) fprintf(stderr, "%s = (%f, %f, %f)\n", #A, A.r, A.g, A.b)
#define PC(A) printf("%s = (%f, %f, %f)\n", #A, A.r, A.g, A.b)
#define PV(A) printf("%s = (%f, %f, %f)\n", #A, A.x, A.y, A.z)

#define true 1
#define false 0
typedef int Bool;

typedef struct {
    double r, g, b;
} Color;

double max3(double a, double b, double c){
    return (a>b ? (a>c? a : c) : (b>c? b : c));
}

Color colorNew(double r, double g, double b){
    Color res;
    res.r = r, res.g = g, res.b = b;
    return res;
}

Color colorClamp(Color c){
    Color res;
    double max = max3(c.r, c.g, c.b);
    if(max > 1.0)   res = colorNew(c.r/max, c.g/max, c.b/max);
    else            res = c;
    return res;
}

Color colorAdd(Color c1, Color c2){
    return (colorNew(c1.r + c2.r, c1.g + c2.g, c1.b + c2.b));
}

Color colorSubtract(Color c1, Color c2){
    return (colorNew(c1.r - c2.r, c1.g - c2.g, c1.b - c2.b));
}

Color colorMultiply(Color c1, Color c2){
    return (colorNew(c1.r * c2.r, c1.g * c2.g, c1.b * c2.b));
}

Color colorScalarMultiply(double k, Color c){
    return colorNew(c.r*k, c.g*k, c.b*k);
}

Color colorDivide(Color c1, Color c2){
    Color res;
    if(c2.r == 0) res.r = 0.0; else res.r = c1.r/c2.r;
    if(c2.g == 0) res.g = 0.0; else res.g = c1.g/c2.g;
    if(c2.b == 0) res.b = 0.0; else res.b = c1.b/c2.b;
    return (res);
}

Color colorPower(Color c, double p){
    return (colorNew(pow(c.r, p), pow(c.g, p), pow(c.b, p)));
}

void test_color(){
    Color c1 = colorNew(0.4, 0.4, 4.0), c2 = colorNew(1.0, 0.0, 0.0), c3 = colorNew(2, 3, 1);
    db_color(c1); db_color(c2);
    db_color(colorAdd(c1, c2));
    db_color(colorMultiply(c1, c2));
    db_color(colorClamp(c3));
    Color div1 = colorNew(0.5, 2, 1), div2 = colorNew(1, 2, 4);
    double k = 2, p = 2;
    db_color(colorScalarMultiply( k,div1));
    db_color(colorDivide(div2, div1));
    db_color(colorPower(div2, p));
}

/************************************ VECTOR ************************************/

typedef struct{
    float x, y, z;
} Vec3d;

Vec3d vecNew(double x, double y, double z){
    Vec3d v;
    v.x = x, v.y = y, v.z = z;
    return v;
}

Vec3d vecAdd(Vec3d a, Vec3d b){
    return vecNew(a.x + b.x, a.y + b.y, a.z + b.z);
}

Vec3d vecSubtract(Vec3d a, Vec3d b){
    return vecNew(a.x - b.x, a.y - b.y, a.z - b.z);
}

Vec3d vecScalarMultiply(double k, Vec3d v){
    return vecNew(k*v.x, k*v.y, k*v.z);
}

Vec3d vecScalarDivide(double k, Vec3d v){
    return vecNew(v.x/k, v.y/k, v.z/k);
}

double vecSquaredLength(Vec3d v){
    return (v.x*v.x + v.y*v.y + v.z*v.z);
}

double vecLength(Vec3d v){
    return sqrt(vecSquaredLength(v));
}

double vecDot(Vec3d v1, Vec3d v2){
    return (v1.x*v2.x + v1.y*v2.y + v1.z*v2.z);
}

Vec3d vecCross(Vec3d u, Vec3d v){
    return vecNew(u.y*v.z - u.z*v.y, u.z*v.x - u.x*v.z, u.x*v.y - u.y*v.x);
}

Vec3d vecInvert(Vec3d v){
    return vecNew(-v.x, -v.y, -v.z);
}

Vec3d vecNormalize(Vec3d v){
    return vecScalarDivide(vecLength(v), v);
}

/********************** RAY *******************************/

typedef struct{
    Vec3d origin, direction; //origin e direction
} Ray;

Ray rayNew(Vec3d o, Vec3d d){
    Ray r;
    r.origin = o, r.direction = d;
    return r;
}

Vec3d rayPoint(Ray r, double t){
    return vecAdd(r.origin, vecScalarMultiply(t, r.direction));
}

/*************************** IMAGE ***********************************/

typedef struct{
    int w, h;
    Color** pixels;
} ImagePPM;

ImagePPM imageNew(int w, int h){
    ImagePPM res; int i;
    res.w = w, res.h = h;
    res.pixels = (Color**) malloc(h*sizeof(Color*));
    for(i = 0; i < h; i++)
        res.pixels[i] = (Color*) malloc(w*sizeof(Color));
    return res;
}

void imageSetPixel(ImagePPM im, int i, int j, double r, double g, double b){
    im.pixels[i][j].r = r;
    im.pixels[i][j].g = g;
    im.pixels[i][j].b = b;
}

void imageSetPixelColor(ImagePPM im, int i, int j, Color c){
    im.pixels[i][j] = c;
}

void imageSave(ImagePPM im){
    int i, j;
    FILE* file = fopen("image.ppm", "wb");
    fprintf(file, "P6\n%d %d\n255\n", im.w, im.h);
    for(i = im.h-1; i >= 0; i--)
    for(j = im.w-1; j >= 0; j--){
        char col[3];
        col[0] = 255*im.pixels[i][j].r;
        col[1] = 255*im.pixels[i][j].g;
        col[2] = 255*im.pixels[i][j].b;
        fwrite(col, sizeof(char), 3, file);
    }
    fclose(file);
}

/********************* CAMERA **************************/

typedef struct{
    Vec3d eyePoint, lookPoint, upVector;
    Vec3d u, v, w;
    double distance, exposureTime, pSize;
    int hres, vres;
} Camera;


Camera cameraNew(int hres, int vres, Vec3d eye, Vec3d lkp, Vec3d upv, double dist, double size){
    Camera c;
    c.eyePoint = eye, c.lookPoint = lkp,
    c.upVector = upv, c.distance = dist,
    c.pSize = size, c.hres = hres,
    c.vres = vres;
    c.w = vecNormalize(vecSubtract(eye, lkp));
    c.u = vecNormalize(vecCross(upv, c.w));
    c.v = vecCross(c.w, c.u);
    return c;
}

Ray cameraGetRay(Camera cam, int r, int c){
    Ray res;
    res.origin = cam.eyePoint;

    double x, y;
    x = cam.pSize*(c - cam.hres/2),
    y = cam.pSize*(r - cam.vres/2);

    Vec3d xu = vecScalarMultiply(x, cam.u);
    Vec3d yv = vecScalarMultiply(y, cam.v);
    Vec3d dw = vecScalarMultiply(cam.distance, cam.w);
    res.direction = vecAdd(xu, vecSubtract(yv, dw)); // xw + yv - dw
    res.direction = vecNormalize(res.direction);
    return res;
}

// INTERSECTION INFO
typedef struct{
    Bool hit, entering;
    double t;
    Vec3d hitPoint;
    Vec3d normal;
    Ray r;
    Color color;            // resíduo
    unsigned long int id;   // para acessar o material
} Intersect;

/*************************** BRDF and MATERIALS ***************************/

#define invPI 0.318309886183

typedef struct{
    float kd;
    Color cd;
} Lambertian;

typedef struct{
    float kr;
    Color cr;
} PerfectSpecular;

typedef struct{
    float ks, exp;
    Color cs;
} GlossySpecular;

typedef enum { LambertianBRDF, PerfectSpecularBRDF, GlossySpecularBRDF } BRDFType;

typedef struct{
    BRDFType brdfType;
    union{
        Lambertian      l;
        PerfectSpecular p;
        GlossySpecular  g;
    } brdfData;
} BRDF;

#define PDouble(A) printf("%s = %f\n", #A, A);
#define P_BRDF(A) printf("%s ---- \n", #A); PDouble(A.brdfData.l.kd); PC(A.brdfData.l.cd);

Color brdfF(BRDF brdf, Intersect it, Vec3d wi, Vec3d wo){
    if(brdf.brdfType == LambertianBRDF){
        //P_BRDF(brdf)
        Color res =(colorScalarMultiply(brdf.brdfData.l.kd, brdf.brdfData.l.cd));
        //PC(res);
        return res;
    }
    else
        return colorNew(1,0,0);
}

Color brdfRho(BRDF brdf, Intersect it, Vec3d wo){
    if(brdf.brdfType == LambertianBRDF){
        //printf("%s = %f\n", "brdf.kd", brdf.brdfData.l.kd);
        //PC(brdf.brdfData.l.cd);
        return ( colorScalarMultiply(brdf.brdfData.l.kd,brdf.brdfData.l.cd));
    }
    else
        return colorNew(1,0,0);
}

/******************** MATERIALS ***********************/

typedef struct{
    BRDF ambientBRDF;
    BRDF diffuseBDRF;
} Matte;

#define PMatte(M) printf("%s----\n", #M); P_BRDF(M.materialData.m.ambientBDRF); P_BRDF(M.materialData.m.diffuseBRDF)

typedef enum {MatteMaterial} MaterialType;

typedef struct {
    MaterialType type;
    union{
        Matte m;
    } materialData;
} Material;

Material materialNew(MaterialType type, Color cd, Color cs, Color cr){
    Material res;
    if(type == MatteMaterial){
        BRDF ambient, diffuse;
        ambient.brdfType = diffuse.brdfType = LambertianBRDF;
        ambient.brdfData.l.kd = 0.25, diffuse.brdfData.l.kd = 0.65;
        ambient.brdfData.l.cd = diffuse.brdfData.l.cd = cd;
        res.materialData.m.ambientBRDF = ambient, res.materialData.m.diffuseBDRF = diffuse;
        P_BRDF(res.materialData.m.ambientBRDF); P_BRDF(res.materialData.m.diffuseBDRF);
        //PMatte(res);
        return res;
    }
}

/******************* OBJECTS **************************/

int NEXT_ID = 0;

#define PM(M) printf("matrix %s ----\n", #M); for(i = 0; i < 4; i++){ for(j = 0; j < 4; j++) printf("%lf ", M.data[i][j]); printf("\n"); }

typedef struct{
    double data[4][4];
} TransformationMatrix;

TransformationMatrix matrixNewIdentity(){
    TransformationMatrix res;
    int i, j;
    for(i = 0; i < 4; i++) for(j = 0; j < 4; j++)
        if(i == j) res.data[i][j] = 1; else res.data[i][j] = 0;
    return res;
}

Vec3d transformPoint(TransformationMatrix m, Vec3d p){
    return vecNew(m.data[0][0]*p.x + m.data[0][1]*p.y + m.data[0][2]*p.z + m.data[0][3],
                  m.data[1][0]*p.x + m.data[1][1]*p.y + m.data[1][2]*p.z + m.data[1][3],
                  m.data[2][0]*p.x + m.data[2][1]*p.y + m.data[2][2]*p.z + m.data[2][3]);
}

Vec3d transformDirection(TransformationMatrix m, Vec3d p){
    return vecNew(m.data[0][0]*p.x + m.data[0][1]*p.y + m.data[0][2]*p.z,
                  m.data[1][0]*p.x + m.data[1][1]*p.y + m.data[1][2]*p.z,
                  m.data[2][0]*p.x + m.data[2][1]*p.y + m.data[2][2]*p.z);
}

Vec3d transformNormal(TransformationMatrix m, Vec3d p){
    return vecNew(m.data[0][0]*p.x + m.data[1][0]*p.y + m.data[2][0]*p.z,
                  m.data[0][1]*p.x + m.data[1][1]*p.y + m.data[2][1]*p.z,
                  m.data[0][2]*p.x + m.data[1][2]*p.y + m.data[2][2]*p.z);
}

TransformationMatrix matrixMultiply(TransformationMatrix m1, TransformationMatrix m2){
    TransformationMatrix res;
    int i, j, k;

    for(i = 0; i < 4; i++)
    for(j = 0; j < 4; j++)
    for(k = 0; k < 4; k++)
        res.data[i][j] += m1.data[i][k] * m2.data[k][j];
    return res;
}

#define PI 3.141592653589793

#define TRANSFORMATION_ROTATION_X 0
#define TRANSFORMATION_ROTATION_Y 1
#define TRANSFORMATION_ROTATION_Z 2

TransformationMatrix matrixGenerateInverseSingleRotation(int axis, double deg){
    TransformationMatrix res = matrixNewIdentity();
    double degrees = (deg*PI)/180;
    if(degrees != 0){
        if(axis == TRANSFORMATION_ROTATION_X)
            res.data[1][1] = cos(-degrees), res.data[1][2] = -sin(-degrees),
            res.data[2][1] = sin(-degrees), res.data[2][2] = cos(-degrees);
        else if(axis == TRANSFORMATION_ROTATION_Y)
            res.data[0][0] = cos(-degrees), res.data[0][2] = sin(-degrees),
            res.data[2][0] = -sin(-degrees), res.data[2][2] = cos(-degrees);
        else if(axis == TRANSFORMATION_ROTATION_Z)
            res.data[0][0] = cos(-degrees), res.data[0][1] = -sin(-degrees),
            res.data[1][0] = sin(-degrees), res.data[1][1] = cos(-degrees);
    }
    return res;
}

TransformationMatrix matrixGenerateInverseTranslation(double tx, double ty, double tz){
    TransformationMatrix res = matrixNewIdentity();
    res.data[0][3] = -tx, res.data[1][3] = -ty, res.data[2][3] = -tz;
    return res;
}

TransformationMatrix matrixGenerateInverseGlobalScale(double sx, double sy, double sz){
    TransformationMatrix res = matrixNewIdentity();
    res.data[0][0] = 1/sx, res.data[1][1] = 1/sy, res.data[2][2] = 1/sz;
    return res;
}

// Prototype
/*
TransformationMatrix matrixInverseLocalScale(Vec3d origin, double sx, double sy, double sz){
    TransformationMatrix centralize, scaleChange, returnToPosition, res;
    centralize =  matrixGenerateInverseTranslation(origin.x, origin.y, origin.z);
    scaleChange = matrixGenerateInverseGlobalScale(sx, sy, sz);
    returnToPosition = matrixGenerateInverseTranslation(-origin.x, -origin.y, -origin.z);
    res = matrixMultiply(returnToPosition, scaleChange);
    res = matrixMultiply(res, )
}*/

TransformationMatrix matrixUberTransformation(double tx, double ty, double tz,
                                              double rx, double ry, double rz,
                                              double sx, double sy, double sz)
{
    TransformationMatrix res, rotx, roty, rotz, trans, scale;
    TransformationMatrix tmp[5];
    /*
    scale = matrixGenerateInverseGlobalScale(sx, sy, sz);
    rotx = matrixGenerateInverseSingleRotation(TRANSFORMATION_ROTATION_X, rx);
    roty = matrixGenerateInverseSingleRotation(TRANSFORMATION_ROTATION_Y, ry);
    rotz = matrixGenerateInverseSingleRotation(TRANSFORMATION_ROTATION_Z, rz);
    trans = matrixGenerateInverseTranslation(tx, ty, tz);
    */

    tmp[0] = matrixGenerateInverseGlobalScale(sx, sy, sz);
    tmp[1] = matrixGenerateInverseSingleRotation(TRANSFORMATION_ROTATION_X, rx);
    tmp[2] = matrixGenerateInverseSingleRotation(TRANSFORMATION_ROTATION_Y, ry);
    tmp[3] = matrixGenerateInverseSingleRotation(TRANSFORMATION_ROTATION_Z, rz);
    tmp[4] = matrixGenerateInverseTranslation(tx, ty, tz);
    res = tmp[0];

    int i, j;
    PM(tmp[i])
    for(i = 1; i < 5; i++)
        {PM(tmp[i]); res = matrixMultiply(res, tmp[i]);}
    return res;
}

typedef struct{
    unsigned long int id;
    double radius;
    Vec3d center;
    Material material;
    TransformationMatrix transformation;
    Bool transform;
} Sphere;

Sphere sphereNew(Vec3d center, double radius, Color diffuseColor, Bool transform){
    Sphere r;
    r.id = NEXT_ID;
    r.center = center, r.radius = radius;
    r.material = materialNew(MatteMaterial, diffuseColor, colorNew(0,0,0),colorNew(0,0,0) );
    r.transformation = matrixNewIdentity();
    r.transform = transform;
    return r;
}

#define K_EPSILON 0.00001

Intersect sphereHit(Sphere s, Ray ray){
    Intersect i;
    double t;

    if(s.transform){
        ray.origin = transformPoint(s.transformation, ray.origin);
        ray.direction = transformDirection(s.transformation, ray.direction);
    }

    Vec3d temp = vecSubtract(ray.origin, s.center);
    double a = vecDot(ray.direction, ray.direction);
    double b = 2.0 * vecDot(temp, ray.direction);
    double c = vecDot(temp, temp) - s.radius*s.radius;
    double disc = b * b - (4.0 * a * c);

    i.hit = false;
    i.entering = false;
    i.t = 0.0;
    i.hitPoint = vecNew(0,0,0);
    i.normal = vecNew(0,0,0);
    i.r = ray;
    if(disc < 0.0) {i.hit = false; return i;}
    else {
        double e = sqrt(disc);
        double denom = 2.0*a;
        t = ((-b) -e)/denom;
        if (t > K_EPSILON ) {
            i.t= t;
            i.entering = true;
            i.hit = true;
            i.r = ray;
            i.hitPoint = rayPoint(ray, t);
            Vec3d tmp = vecNew(i.hitPoint.x - s.center.x,
                               i.hitPoint.y - s.center.y,
                               i.hitPoint.z - s.center.z);
            if(s.transform){
                i.normal = transformNormal(s.transformation, vecNormalize(tmp));
                i.normal = vecNormalize(i.normal);
            }else i.normal = vecNormalize(tmp);
            return i;
        }
        t = (-b + e)/denom;
        if (t > K_EPSILON) {
            i.t= t;
            i.entering = false;
            i.hit = true;
            i.r = ray;
            i.hitPoint = rayPoint(ray, t);
            Vec3d tmp = vecNew(i.hitPoint.x - s.center.x,
                               i.hitPoint.y - s.center.y,
                               i.hitPoint.z - s.center.z);
            if(s.transform){
                i.normal = transformNormal(s.transformation, vecNormalize(tmp));
                i.normal = vecNormalize(i.normal);
            }else i.normal = vecNormalize(tmp);
            return i;
        }
    }
    i.hit = false;
    return i;
}

typedef struct{
    Vec3d v1, v2, v3, normal;
    //Material material;
    //Bool transform;
} Triangle;

Triangle triangleNew(Vec3d v1, Vec3d v2, Vec3d v3){
    Triangle t;
    t.v1 = v1, t.v2 = v2, t.v3 = v3;
    t.normal = vecCross(vecSubtract(v2, v1), vecSubtract(v3, v1));
    t.normal = vecNormalize(t.normal);
    //t.material = materialNew(MatteMaterial, c, colorNew(0,0,0), colorNew(0,0,0));
    return t;
}

Intersect triangleIntersect(Triangle t, Ray ray){

    //printf("teste___\n"); fflush(stdout);
    Intersect it;
    it.hit = false;
    it.entering = false;
    it.t = 0.0;
    it.hitPoint = vecNew(0,0,0);
    it.normal = vecNew(0,0,0);
    it.r = ray;

    double a = t.v1.x - t.v2.x, b = t.v1.x - t.v3.x, c = ray.direction.x, d = t.v1.x - ray.origin.x;
    double e = t.v1.y - t.v2.y, f = t.v1.y - t.v3.y, g = ray.direction.y, h = t.v1.y - ray.origin.y;
    double i = t.v1.z - t.v2.z, j = t.v1.z - t.v3.z, k = ray.direction.z, l = t.v1.z - ray.origin.z;

    double m = f*k - g*j, n = h*k - g*l, p = f*l - h*j;
    double q = g*i - e*k, s = e*j - f*i;

    double inv_denom = 1.0/(a*m + b*q + c*s);

    double e1 = d*m - b*n - c*p;
    double beta = e1*inv_denom;

    if(beta < 0.0)      return it;

    double r = e*l - h*i;
    double e2 = a*n + d*q + c*r;
    double gamma = e2 * inv_denom;

    if(gamma < 0.0)         return it;
    if(beta + gamma > 1.0)  return it;

    double e3 = a*p - b*r + d*s;
    double tmin = e3 * inv_denom;

    if(tmin < K_EPSILON)    return it;

    //printf("teste___\n"); fflush(stdout);

    it.t = tmin;
    it.normal = t.normal;
    it.hit = true;
    it.entering = true;
    it.hitPoint = rayPoint(ray, tmin);

    return it;
}

#define MAX_TRIANGLES 50000

typedef struct {
    Triangle triangles[MAX_TRIANGLES];
    int numTriangles;
    Material material;
    TransformationMatrix transformation;
    Bool transform;
} Mesh;

Mesh meshNew(Color c){
    Mesh res;
    res.numTriangles = 0;
    res.material = materialNew(MatteMaterial, c, c, c);
    P_BRDF(res.material.materialData.m.ambientBRDF); P_BRDF(res.material.materialData.m.diffuseBDRF);
    return res;
}

typedef struct{
    double x, y, z;
} vertex;

#define MAX_VERTICES 30000
#define PT(T) PV(T.v1); PV(T.v2); PV(T.v3)

Mesh meshImporter(char filename[], Color c){
    Mesh res;
    res.numTriangles = 0;
    res.material = materialNew(MatteMaterial, c,c,c);


    FILE* file = fopen(filename, "r");
    char line[200];
    int v_count = 0, vt_count = 0, vn_count = 0, f_count = 0;
    Vec3d v[MAX_VERTICES];//, vn[MAX_VERTICES];
    //Vec2d vt[MAX_VERTICES];

    //printf("teste 1---------------\n"); fflush(stdout);


    while(!feof(file)){
        if(fgets(line, 200, file) == NULL) break;
        //printf("teste---- %s -----------\n", line); fflush(stdout);
        if(line[0] == 'v'){
            if(line[1] == ' '){ //vertex
                sscanf(line+1, "%f %f %f\n", &v[v_count].x, &v[v_count].y, &v[v_count].z);
                //PV(v[v_count]);
                //if(v_count % 5000 == 0) PV3(v[v_count]);
                v_count++;
            } else if(line[1] == 'n'){
                //sscanf(line+2, "%f %f %f\n", &vn[vn_count].x, &vn[vn_count].y, &vn[vn_count].z);
                //if(vn_count % 5000 == 0) PV3(vn[vn_count]);
                vn_count++;
            }else if(line[1] == 't'){
                //sscanf(line+2, "%lf %lf\n", &vt[vt_count].x, &vt[vt_count].y);
                //if(vt_count % 5000 == 0) PV2(vt[vt_count]);
                vt_count++;
            }

        }else if(line[0] == 'f'){
            //printf("%s\n",line);
            if(res.numTriangles<MAX_TRIANGLES){
                int a, b, c, d, e, f, g, h, i;
                //sscanf(line+1, "%d/%d/%d %d/%d/%d %d/%d/%d\n", &a,&g,&h,&b,&e,&f,&c,&h,&i);
                sscanf(line+1, "%d %d %d\n", &a,&b,&c);
                //printf("%s -----\na, d, g = %d, %d, %d\n", line+1, a, d, g);
                //PV(v[a-1]); PV(v[b-1]); PV(v[c-1]);
                res.triangles[res.numTriangles] = triangleNew(v[a-1], v[b-1], v[c-1]);
                //PT(res.triangles[res.numTriangles]);
                res.numTriangles++;
            }
            f_count++;
            //printf("%d/%d/%d nt: %d\n", a,d,g,res.numTriangles);
        }
    }

    int i;
    //for(i = 0; i < v_count; i++) PV(v[i]);

    printf("Counts:\nv: %d, vt: %d, vn: %d\n", v_count, vt_count, vn_count);
    printf("triangle number: %d  triangles in file %d\n", res.numTriangles, f_count);
    fclose(file);
    return res;
}

Intersect meshHit(Mesh mesh, Ray ray){

    if(mesh.transform){
        ray.origin = transformPoint(mesh.transformation, ray.origin);
        ray.direction = transformDirection(mesh.transformation, ray.direction);
    }

    //printf("teste_hit_enter_\n"); fflush(stdout);
    Intersect it, itMin;
    itMin.hit = false;
    double tmin = 10000000.;
    int i;

    for(i = 0; i < mesh.numTriangles; i++){
        //printf("teste_triangle_intersect_\n"); fflush(stdout);
        it = triangleIntersect(mesh.triangles[i], ray);
        //printf("teste_triangle_end_\n"); fflush(stdout);
        if(it.hit && it.t < tmin)
            itMin = it, tmin = it.t;
    }
    //printf("teste_hit_out_\n"); fflush(stdout);
    itMin.normal = transformNormal(mesh.transformation, itMin.normal);
    return itMin;
}



#define N_OBJS 1
Mesh scene_m[N_OBJS];
Triangle scene_t[N_OBJS];
Sphere scene[N_OBJS];

/*********************** LIGHT ***************************/

typedef enum {AmbientLight, PointLight, DirectionalLight} LightType;

typedef struct{
    LightType type;     // type of light
    Vec3d posDir;       // (x, y, z) Dir <- directional light, pos <- point light
    Color color;        // (r, g, b)
    double intensity;   // [0, 1]
    int shadows;        // [0-1]    <- point and directional lights
} Light;

Light lightNew(LightType type, Vec3d pos_dir, Color c, double intensity, int shadows){
    Light l;
    l.type = type, l.posDir = pos_dir, l.color = c;
    l.intensity = intensity, l.shadows = shadows;
    return l;
}

Vec3d lightGetDirection(Light l, Vec3d hitPoint){
    if(l.type == PointLight){
        return vecNormalize(vecSubtract(l.posDir, hitPoint));
    }else
        return vecNew(0,0,0);
}

Color lightGetL(Light l){
    return colorScalarMultiply(l.intensity, l.color);
}

#define N_LIGHTS 1
Light lights[N_LIGHTS];
Light ambient;

/*************************** TRACER ****************************/

Color background;

Color matteShade(Material m, Intersect i){
    Vec3d wo = vecInvert(i.r.direction);
    Color L = colorMultiply(brdfRho(m.materialData.m.ambientBRDF, i, wo), lightGetL(ambient));
    int j;
    for(j = 0; j < N_LIGHTS; j++){
        Vec3d wi = lightGetDirection(lights[j], i.hitPoint);
        double ndotwi = vecDot(i.normal, wi);
        if(ndotwi>0.0) L = colorAdd(L, colorScalarMultiply(ndotwi, colorAdd(brdfF(m.materialData.m.diffuseBDRF, i, wi, wo), lightGetL(lights[j]))));
    }
    return colorClamp(L);
}

Intersect rayHitObjects(Ray r){
    Intersect it, res;
    double tmin = 1000000.0;
    int i, id = 0;
    res.hit = false;
    for(i = 0; i < N_OBJS; i++){
        it = meshHit(scene_m[i], r);
        if(it.hit && it.t < tmin){
            res = it;
            tmin = res.t;
            id = i;
        }
    }
    res.id = id;
    return res;
}

Color rayTraceSimple(Ray r){
    double tmin = INT_MAX; int i;
    Intersect it; Color res = colorNew(0.1, 0.1, 0.1);
    it = rayHitObjects(r);
    if(it.hit){ res = matteShade(scene_m[it.id].material, it);}// PC(res);}
    //else  res = background;
    //PC(background);
    //PC(res);
    return res;
}

int main(){
    const int IMAGE_VRES = 500, IMAGE_HRES = 500;// NUM_OBJS = 10, NUM_LIGHTS = 10;
    int i, j;
    ImagePPM im = imageNew(IMAGE_HRES, IMAGE_VRES);
    Camera cam = cameraNew(IMAGE_HRES, IMAGE_VRES, vecNew(60, 0, 60),
                           vecNew(0,0,0), vecNew(0,0,1), 200, 1.0);

    background = colorNew(0.1, 0.1, 0.1);
    scene[0] = sphereNew(vecNew(0, 0, 0), 10, colorNew(1,0,1), true);
    scene[0].transformation = matrixGenerateInverseGlobalScale(1, 2, 3);
    scene[1] = sphereNew(vecNew(0, 40, -40), 20, colorNew(1,0,0), false);

    scene_m[0] = meshImporter("albertosaurus2.obj", colorNew(1, 0, 0));

    scene_m[0].transform = false;
    scene_m[0].transformation = matrixGenerateInverseSingleRotation(TRANSFORMATION_ROTATION_X, 90);

    printf("triangulos cubo:\n");
    for(i = 0; i < scene_m[0].numTriangles; i++){
        //PT(scene_m[0].triangles[i]);
    }


/*
    scene_m[0] = meshNew(colorNew(1,0,0));
    scene_m[0].triangles[0] = triangleNew(vecNew(-10,-10,0), vecNew(0, 10, 0), vecNew(0,0,10));
    scene_m[0].triangles[1] = triangleNew(vecNew(-10,-30,30), vecNew(0, 10, 30), vecNew(0,0,40));
    scene_m[0].numTriangles = 2;*/



    lights[0] = lightNew(PointLight, vecNew(100,40,0), colorNew(1,1,1), .3, 0);
    ambient = lightNew(AmbientLight, vecNew(0,0,0), colorNew(1,1,1), 0.1, 0);


    for(i = 0; i < im.h; i++)
        for(j = 0; j < im.w; j++){
            Ray ray = cameraGetRay(cam, i, j);
            //printf("teste___\n"); fflush(stdout);
            Color col = rayTraceSimple(ray);
            //printf("teste___\n"); fflush(stdout);
            imageSetPixelColor(im, i, j, col);
        }
    imageSave(im);

    return 0;
}

