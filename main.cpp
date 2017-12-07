#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <cmath>


using namespace std;

#define printVec(A) std::cout << #A << " = (" << (A).x <<", "<<(A).y<<", "<<(A).z<<")"<<std::endl

#define PV(A) printf("%s = (%f, %f, %f)\n", #A, A.x, A.y, A.z)

class Vec3d
{
public:
    float x, y, z;
    Vec3d(): x(0.0), y(0.0), z(0.0) {} // KAIOSHIN
    Vec3d(const float a) : x(a), y(a), z(a) {} // KAIOSHIN
    Vec3d(const float a, const float b, const float c) : x(a), y(b), z(c) {} // KAIOSHIN
    Vec3d(const Vec3d& v) : x(v.x), y(v.y), z(v.z) {} // KAIOSHIN
    ~Vec3d (void){}; // HAKAISHIN
    Vec3d& operator= (const Vec3d& rhs); //Sobrecarrega os operadores de atribuição
    Vec3d operator- (void) const; // Operador menos unário
    Vec3d operator* (const float a) const; //multiplicação por escalar
    Vec3d operator/ (const float a) const; // divisão por escalar
    Vec3d operator+ (const Vec3d& v) const; // soma de vetores
    Vec3d& operator+= (const Vec3d& v); // soma e atribuição de vetores
    Vec3d operator- (const Vec3d& v) const; // subtração de vetores
    float operator* (const Vec3d& b) const; // produto escalar - dot product de dois vetores
    Vec3d operator^ (const Vec3d& v) const; // produto vetorial - cross product de dois vetores
    float length(void); // retorna a magnitude do vetor
    float len_squared(void); // retorna o quadrado da magnitude do vetor
    void normalize(void); // transforma o vetor em vetor unitário
    Vec3d& hat(void); // retorna o vetor normalizado e transforma o vetor em unitário
};

// ----- multiplicação com o escalar do lado esquerdo -------
Vec3d operator* (const double a, const Vec3d& v);
Vec3d operator* (const int a, const Vec3d& v);


inline Vec3d Vec3d::operator- (void) const { return Vec3d(-x, -y, -z); }
inline Vec3d Vec3d::operator* (const float a) const { return Vec3d(a*x, a*y, a*z); }
inline Vec3d Vec3d::operator/ (const float a) const { return Vec3d(x/a, y/a, z/a); }
inline Vec3d Vec3d::operator+ (const Vec3d& c) const { return Vec3d(x + c.x, y + c.y, z + c.z); }
inline Vec3d& Vec3d::operator+= (const Vec3d& c) { x += c.x; y += c.y; z += c.z; return (*this); }
inline Vec3d Vec3d::operator- (const Vec3d& c) const { return Vec3d(x - c.x, y - c.y, z - c.z); }
inline float Vec3d::operator* (const Vec3d& c) const { return ( (x * c.x)+ (y * c.y) + (z * c.z) ) ; }
inline Vec3d Vec3d::operator^ (const Vec3d& v) const { return Vec3d(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x); }

inline Vec3d operator* (const double a, const Vec3d& v) { return (Vec3d(a * v.x, a * v.y, a * v.z)); }
inline Vec3d operator* (const int a, const Vec3d& v) { return (Vec3d(a * v.x, a * v.y, a * v.z)); }
Vec3d& Vec3d::operator= (const Vec3d& rhs){
    if (this == &rhs)
        return (*this);
    x = rhs.x; y = rhs.y; z = rhs.z;
    return (*this);
}

float Vec3d::length(void) { return sqrt((x*x) + (y*y) + (z*z)); }
float Vec3d::len_squared(void) { return (x*x) + (y*y) + (z*z); }

void Vec3d::normalize(void){
	double length = sqrt(x * x + y * y + z * z);
	x /= length; y /= length; z /= length;
}

Vec3d& Vec3d::hat(void){
	double length = sqrt(x * x + y * y + z * z);
	x /= length; y /= length; z /= length;
	return (*this);
}

/******************************** RAY *****************************************************/

#define printRay(R) std::cout<< #R << " ----------" <<std::endl; printVec(R.origin); printVec(R.direction);

class Ray
{
public:
    Vec3d origin, direction;
    Ray(const Vec3d& o, const Vec3d& d) : origin(o), direction(d) {}
    inline Vec3d rayPoint(float t) const { return (t * direction + origin); }
};



/**************************************** CAMERA *******************************************/


class Camera{
public:

    Camera(const int hres, const int vres, const Vec3d& eye, const Vec3d& lkp,
           const Vec3d& upv, const double dist, const double s){
        setup(hres, vres, eye, lkp, upv, dist, s);

    }


    Ray getRay(int r, int c){
        Ray res = Ray(Vec3d(), Vec3d());
        res.origin = eyePoint;
        double x, y;
        x = pSize*(c - hres/2),
        y = pSize*(r - vres/2);

        Vec3d xu = x*u;// vecScalarMultiply(x, cam.u);
        Vec3d yv = y*v;//vecScalarMultiply(y, cam.v);
        Vec3d dw = distance*w;//vecScalarMultiply(cam.distance, cam.w);
        res.direction = xu + yv - dw; // xw + yv - dw
        res.direction.normalize();
        return res;
    }

private:

    void setup(const int hres, const int vres, const Vec3d& eye, const Vec3d& lkp,
           const Vec3d& upv, const double dist, const double s){
        eyePoint = eye, lookPoint = lkp,
        upVector = upv, distance = dist,
        pSize = s, this->hres = hres, this->vres = vres;
        w = eye - lkp;
        w.normalize();
        u = upv ^ w;
        u.normalize();
        v = w ^ u;

    }

    Vec3d eyePoint, lookPoint, upVector;
    Vec3d u, v, w;
    double distance, exposureTime, pSize;
    int hres, vres;

};





/**************************************** COLOR ******************************************/


#define printCol(A) std::cout << #A << " = (" << (A).r <<", "<<(A).g<<", "<<(A).b<<")"<<std::endl

#define BLACK Color(0.0)
#define WHITE Color(1.0)
#define RED Color(1.0, 0.0, 0.0)
#define GREEN Color(0.0, 1.0, 0.0)
#define BLUE Color(0.0, 0.0, 1.0)
#define YELLOW Color(1.0, 1.0, 0.0)
#define MAGENTA Color(0.0, 1.0, 1.0)
#define ORANGE Color(1.0, 0.5, 0.0)


class Color
{
public:
    float r, g, b;
    Color(): r(0.0), g(0.0), b(0.0) {}                                          // KAIOSHIN
    Color(const float a) : r(a), g(a), b(a) {}                                  // KAIOSHIN
    Color(const float a, const float b, const float c) : r(a), g(b), b(c) {}    // KAIOSHIN
    Color(const Color& v) : r(v.r), g(v.g), b(v.b) {}                           // KAIOSHIN
    ~Color (void) {}                                                            // HAKAISHIN
    Color& operator= (const Color& rhs);    //Sobrecarrega os operadores de atribuição
    Color operator- (void) const;           // Operador menos unário
    Color operator* (const float a) const;  //multiplicação por escalar
    Color operator/ (const float a) const;  // divisão por escalar
    Color operator+ (const Color& v) const; // soma de cores
    Color& operator+= (const Color& v);     // soma e atribuição de cores
    Color operator- (const Color& v) const; // subtração de cores
    Color operator* (const Color& b) const; // produto de cores
    Color operator/ (const Color& b) const; // divisão de cores
    Color clamp();
};



inline Color Color::operator- (void) const { return Color(-r, -g, -b); }
inline Color Color::operator* (const float a) const { return Color(a*r, a*g, a*b); }
inline Color Color::operator/ (const float a) const { return Color(r/a, g/a, b/a); }
inline Color Color::operator+ (const Color& c) const { return Color(r + c.r, g + c.g, b + c.b); }
inline Color& Color::operator+= (const Color& c) { r += c.r; g += c.g; b += c.b; return (*this); }
inline Color Color::operator- (const Color& c) const { return Color(r - c.r, g - c.g, b - c.b); }
inline Color Color::operator* (const Color& c) const { return Color(r * c.r, g * c.g, b * c.b); }
inline Color Color::operator/ (const Color& c) const { return Color(r / c.r, g / c.g, b / c.b); }

// ----- multiplicação com o escalar do lado esquerdo -------
inline Color operator* (const double a, const Color& v) { return (Color(a * v.r, a * v.g, a * v.b)); }
inline Color operator* (const int a, const Color& v) { return (Color(a * v.r, a * v.g, a * v.b)); }

Color& Color::operator= (const Color& rhs){
    if (this == &rhs) return (*this);
    r = rhs.r; g = rhs.g; b = rhs.b;
    return (*this);
}

Color Color::clamp()
{
    float maxVal = std::max(r, std::max(g, b));
    if(maxVal > 1.0)
        r /= maxVal, g /= maxVal, b /= maxVal;
    return Color(r, g, b);
}



/************************************* MATERIAL ****************************************/

/*
    Classe generalizada para um material em RayTracer sem BRDFs

    kd -> coeficiente difuso
    ks -> coeficiente especular
    exp-> expoente especular
    kr -> coeficiente reflexivo
    cd -> cor difusa
    cs -> cor da reflexão especular
    cr -> cor da reflexão perfeita

*/

class Material{

public:
    Color cd, cs, cr;
    float kd, ks, exp, kr;
    Material(const Color& cd, const float kd, const Color& cs, const float ks, const float exp, const Color& cr, const float kr) :
        cd(cd) , cs(cs), cr(cr), kd(kd), ks(ks), exp(exp), kr(kr){}
};

/************************************** LIGHT ***************************************/


class Light
{
public:
    Vec3d position;
    float strength;
    Color color;
    Light(const Vec3d& p, const float s, const Color& c) : position(p), strength(s), color(c){};
    Color getL() const;
    Vec3d getDirection(const Vec3d& pt) const;
};


/**************************************  BBox ********************************************/


class BBox{
public:
    Vec3d minp, maxp;
    BBox(Vec3d v0, Vec3d v1) : minp(v0), maxp(v1){}
};




/************************************** TRIANGLE ****************************************/

class Triangle{
public:
    Vec3d p1, p2, p3, normal;
    int materialID;
    Triangle(const Vec3d& v1, const Vec3d& v2, const Vec3d& v3, const int material) : p1(v1), p2(v2), p3(v3), materialID(material) {
        normal = (v2 - v1) ^ (v3 - v1);
        normal.normalize();
    }

    Vec3d getMaxPoint(){
        Vec3d p(-1e10);

        // maior coordenada de cada vértice

        if(p1.x > p.x) p.x = p1.x;
        if(p1.y > p.y) p.y = p1.y;
        if(p1.z > p.z) p.z = p1.z;

        if(p2.x > p.x) p.x = p2.x;
        if(p2.y > p.y) p.y = p2.y;
        if(p2.z > p.z) p.z = p2.z;

        if(p3.x > p.x) p.x = p3.x;
        if(p3.y > p.y) p.y = p3.y;
        if(p3.z > p.z) p.z = p3.z;

        return p;

    }

    Vec3d getMinPoint(){
        Vec3d p(1e10);

        // menor coordenada de cada pértice

        if(p1.x < p.x) p.x = p1.x;
        if(p1.y < p.y) p.y = p1.y;
        if(p1.z < p.z) p.z = p1.z;

        if(p2.x < p.x) p.x = p2.x;
        if(p2.y < p.y) p.y = p2.y;
        if(p2.z < p.z) p.z = p2.z;

        if(p3.x < p.x) p.x = p3.x;
        if(p3.y < p.y) p.y = p3.y;
        if(p3.z < p.z) p.z = p3.z;

        return p;

    }

};



class Grid{
public:
    Grid() : boundingBox(BBox(Vec3d(), Vec3d())){
        //boundingBox = BBox(Vec3d(), Vec3d());
    }

    BBox getBoundingBox(){}
    void setup(const vector<Triangle>& tris){}


private:
    vector<Triangle> orderedCells; // Voxels com repetição de triângulos para serem enviados diretamente para a FPGA
    BBox boundingBox;
    int nx, ny, nz;
    Vec3d minCoordinates();
    Vec3d maxCoordinates();
};




#define MAX_TRIANGLES 50000

vector<Triangle> meshes;
vector<Material> materials;
vector<Light> lights;
vector<Ray> rays;

// matID -> id do material

#define MAX_VERTICES 30000
#define PT(T) PV(T.v1); PV(T.v2); PV(T.v3)

vector<Triangle> meshImporter(char filename[], int mID){

    FILE* file = fopen(filename, "r");
    char line[200];
    int v_count = 0, f_count = 0, numTriangles = 0;
    float x, y, z;
    vector<Vec3d> v;
    vector<Triangle> tris;

    while(!feof(file)){
        if(fgets(line, 200, file) == NULL) break;

        if(line[0] == 'v'){
            if(line[1] == ' '){ //vertex
                sscanf(line+1, "%f %f %f\n", &x, &y, &z);
                v.push_back(Vec3d(x, y, z));
                v_count++;
            }
        }else if(line[0] == 'f'){

            if(v.size() < MAX_TRIANGLES){
                int a, b, c, d, e, f, g, h, i;
                //sscanf(line+1, "%d/%d/%d %d/%d/%d %d/%d/%d\n", &a,&g,&h,&b,&e,&f,&c,&h,&i);
                sscanf(line+1, "%d %d %d\n", &a,&b,&c);
                tris.push_back( Triangle(v[a-1], v[b-1], v[c-1], mID) );
                numTriangles++;
            }
            f_count++;
        }
    }

    printf("Counts:\nv: %d\n", v_count);
    printf("triangle number: %d  triangles in file %d\n", numTriangles, f_count);
    fclose(file);
    return tris;
}




int main(int argc, char**argv){

    char* filename;
    if(argc == 2)
        filename = argv[1];
    else
        filename = "teddy.obj";

    Material red_matte(Color(1,0,0),0.65, Color(),0,0,Color(),0);

    //Meshes = meshImporter(filename, )

    Triangle t(Vec3d(-5,0,0), Vec3d(0,0,10), Vec3d(7, 5, 2), 0);

    PV(t.getMinPoint()); PV(t.getMaxPoint());


    return 0;
}
