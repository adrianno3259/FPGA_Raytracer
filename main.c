#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
//#include "helper.c"
//#include "image.h"

#define db_color(A) fprintf(stderr, "%s = (%f, %f, %f)\n", #A, A.r, A.g, A.b)

typedef struct {
    double r, g, b;
} Color;

double max3(double a, double b, double c){
    return (a>b ? (a>c? a : c) : (b>c? b : c));
}

Color colorNew(double r, double g, double b)
{
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
    double x, y, z;
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

void imageSave(ImagePPM im)
{
    FILE* file = fopen("image.ppm", "wb");

    fprintf(file, "P6\n%d %d\n255\n", im.w, im.h);
    //file<<'P'<<'6'<<'\n'<<im.w<<' '<<im.h<<'\n'<<'2'<<'5'<<'5'<<'\n';
    for(int i = im.h-1; i >= 0; i--)
    for(int j = im.w-1; j >= 0; j--)
    {
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
    //c.u = vecInvert()
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



#define true 1
#define false 0
typedef int Bool;

// INTERSECTION INFO
typedef struct{
    Bool hit, entering;
    double t;
    Vec3d hitPoint;
    Vec3d normal;
    Ray r;
    Color color;
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

Color brdfF(BRDF brdf, Intersect it, Vec3d wi, Vec3d wo){
    if(brdf.brdfType == LambertianBRDF)
        return (colorScalarMultiply(brdf.brdfData.l.kd*invPI, brdf.brdfData.l.cd));
    else
        return colorNew(1,0,0);
}

Color brdfRho(BRDF brdf, Intersect it, Vec3d wo){
    if(brdf.brdfType == LambertianBRDF)
        return ( colorScalarMultiply(brdf.brdfData.l.kd,brdf.brdfData.l.cd));
    else
        return colorNew(1,0,0);
}

/******************** MATERIALS ***********************/

typedef struct{
    BRDF ambientBRDF;
    BRDF diffuseBDRF;
} Matte;

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
        return res;
    }
}



/******************* OBJECTS **************************/

typedef struct{
    double radius;
    Vec3d center;
    Material material;
} Sphere;

Sphere sphereNew(Vec3d center, double radius, Color c){
    Sphere r;
    r.center = center, r.radius = radius, r.color = c;
    return r;
}

#define K_EPSILON 0.00001

Intersect sphereHit(Sphere s, Ray ray){
    Intersect i;
    double t;
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
    //i.obj = NULL;
    i.r = ray;
    //printVar(a); printVar(b); printVar(c); printVar(disc);
    if(disc < 0.0) {i.hit = false; return i;}
    else {
        double e = sqrt(disc);
        double denom = 2.0*a;
        t = ((-b) -e)/denom;
        if (t > K_EPSILON ) {
            i.t= t;
            i.entering = true;
            i.hit = true;
            //i.obj = s;
            i.color = s.color;
            //i.m = this->m;
            i.r = ray;
            i.hitPoint = rayPoint(ray, t);
            Vec3d tmp = vecNew(i.hitPoint.x - s.center.x,
                               i.hitPoint.y - s.center.y,
                               i.hitPoint.z - s.center.z);
            i.normal = vecNormalize(tmp);
            return i;
        }
        t = (-b + e)/denom;
        if (t > K_EPSILON) {
            i.t= t;
            i.entering = false;
            i.hit = true;
            //i.obj = s;
            i.color = s.color;
            //i.m = this->m;
            i.r = ray;
            i.hitPoint = rayPoint(ray, t);
            Vec3d tmp = vecNew(i.hitPoint.x - s.center.x,
                               i.hitPoint.y - s.center.y,
                               i.hitPoint.z - s.center.z);
            i.normal = vecNormalize(tmp);
            return i;
        }
    }
    i.hit = false;
    return i;
}

#define N_OBJS 10
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

#define N_LIGHTS 10
Light lights[N_LIGHTS];
Light ambient;

/*************************** TRACER ****************************/

Color matteShade(Material m, Intersect i){
    Vec3d wo = vecInvert(i.r.direction);
    Color L = colorMultiply(brdfRho(m.materialData.m.ambientBRDF, i, wo), lightGetL(ambient));
    int j;
    for(j = 0; j < N_LIGHTS; j++){
        Vec3d wi = lightGetDirection(lights[j], i.hitPoint);
        float ndotwi = vecDot(i.normal, wi);
        if(ndotwi>0.0) colorAdd(L, colorScalarMultiply(ndotwi, colorAdd(brdfF(m.materialData.m.diffuseBDRF, i, wi, wo), lightGetL(lights[j]))));
    }
    return L;
}


Intersect rayHitObjects(Ray r){
    Intersect it, res;
    double tmin = 1000000.0;
    int i;
    res.hit = false;
    for(i = 0; i < N_OBJS; i++){
        it = sphereHit(scene[i], r);
        if(it.hit && it.t < tmin){
            res = it;
            tmin = res.t;
        }
    }

    return res;
}

Color rayTraceSimple(Ray r){
    double tmin = INT_MAX; int i;
    Intersect it; Color res = colorNew(0, 0, 0); //Background color

    it = rayHitObjects(r);

    if(it.hit){
        return matteShade(it.)
    }

    /*for(i = 0; i < numSpheres; i++){
        it = sphereHit(scene[i], r);
        if(it.hit){
            if(it.t < tmin){
                tmin = it.t;
                //imin = i;
                res = it.color;
            }
        }
    }*/

    return res;
}


int main(){

    //test_color();

    const int IMAGE_VRES = 1000, IMAGE_HRES = 1000, NUM_OBJS = 10, NUM_LIGHTS = 10;

    ImagePPM im = imageNew(IMAGE_HRES, IMAGE_VRES);
    Camera cam = cameraNew(IMAGE_HRES, IMAGE_VRES, vecNew(500, 0, 0),
                           vecNew(0,0,0), vecNew(0,0,1), 200, 1.0);


    scene[0] = sphereNew(vecNew(0, 0, 0), 20, colorNew(1,0,1));
    scene[1] = sphereNew(vecNew(100, 100, 100), 20, colorNew(1,0,0));


    lights[1] = lightNew(PointLight, vecNew(10,0,20), colorNew(1,1,1), 1.0, 0);

    int i, j;
    for(i = 0; i < im.h; i++)
        for(j = 0; j < im.w; j++){

            Ray ray = cameraGetRay(cam, i, j);

            Color col = rayTraceSimple(ray, scene, 2);

            imageSetPixelColor(im, i, j, col);

            /*int rad = (int)sqrt((pow(i-im.h/2, 2) + pow(j-im.w/2, 2)));
            if( rad <= 500 && rad >400)
                imageSetPixel(im, i, j, 0, 0, 0);
            else
              imageSetPixel(im, i, j, 1, 1, 1);*/
        }
    imageSave(im);
    return 0;

}

