#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <cmath>


#define MAX_TRIANGLES 50000
#define MAX_VERTICES 30000

#define TRIANGLE_ATTRIBUTE_NUMBER 13
#define LIGHT_ATTRIBUTE_NUMBER 7
#define MATERIAL_ATTRIBUTE_NUMBER 4
#define RAY_ATTRIBUTE_NUMBER 6

#define PT(T) PV(T.p1); PV(T.p2); PV(T.p3)
#define PV(A) printf("%s = (%f, %f, %f)\n", #A, A.x, A.y, A.z)



using namespace std;


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

inline float clamp(float x, float minimum, float maximum){
    return (x < minimum ? minimum : (x > maximum ? maximum : x));
}

#define pVar(V) cout<<#V<<" "<<V<<endl

class Grid{
public:
    Grid(const vector<Triangle>& src) :totalTriangles(0), boundingBox(BBox(Vec3d(), Vec3d())), source(src){ }
    BBox getBoundingBox(){}
    void setup(){
        Vec3d p0 = minCoordinates(), p1 = maxCoordinates();
        boundingBox.maxp = p1, boundingBox.minp = p0;


        // obter número de voxels
        int nObjects = source.size();
        float wx = p1.x - p0.x, wy = p1.y - p0.y, wz = p1.z - p0.z;
        float multiplier = 1.0;
        float s = pow(wx*wy*wz/nObjects, 0.3333333);
        nx = multiplier*wx/s+1;
        ny = multiplier*wy/s+1;
        nz = multiplier*wz/s+1;

        pVar(nx); pVar(ny); pVar(nz);
        cout<<endl<<endl<<endl;

        // contagem de elementos por voxel
        int nCells = nx*ny*nz;

        //orderedCells = vector<vector<Triangle> >(nCells, vector<Triangle>());
        for(int i = 0; i < nCells; i++) orderedCells.push_back(vector<Triangle>());

        for(int i = 0; i < nCells; i++) counts.push_back(0);
        pVar(orderedCells.size());
        pVar(counts.size());

        BBox tmp_box = BBox(Vec3d(), Vec3d());
        int index;

        for(int j = 0; j < nObjects; j++){
            tmp_box = BBox(source[j].getMinPoint(), source[j].getMaxPoint());

            int ixmin = clamp((tmp_box.minp.x - p0.x) * nx / (p1.x - p0.x), 0, nx - 1),
                iymin = clamp((tmp_box.minp.y - p0.y) * ny / (p1.y - p0.y), 0, ny - 1),
                izmin = clamp((tmp_box.minp.z - p0.z) * nz / (p1.z - p0.z), 0, nz - 1),
                ixmax = clamp((tmp_box.maxp.x - p0.x) * nx / (p1.x - p0.x), 0, nx - 1),
                iymax = clamp((tmp_box.maxp.y - p0.y) * ny / (p1.y - p0.y), 0, ny - 1),
                izmax = clamp((tmp_box.maxp.z - p0.z) * nz / (p1.z - p0.z), 0, nz - 1);

            pVar(ixmin); pVar(iymin); pVar(izmin);
            pVar(ixmax); pVar(iymax); pVar(izmax);
            cout<<endl<<endl;

            for(int iz = izmin; iz <= izmax; iz++)
            for(int iy = iymin; iy <= iymax; iy++)
            for(int ix = ixmin; ix <= ixmax; ix++){
                index = ix + nx*iy + nx*ny*iz;
                pVar(index);
                PT(source[j]);
                orderedCells[index].push_back(source[j]);
                PT(orderedCells[index][counts[j]]);
                counts[index]++;
                totalTriangles++;
            }

        }
        pVar(orderedCells.size());
        pVar(counts.size());

    }

    vector<vector<Triangle> > orderedCells; // Voxels com repetição de triângulos para serem enviados diretamente para a FPGA
    vector<int> counts;
    int totalTriangles;
private:
    vector<Triangle> source;
    BBox boundingBox;
    int nx, ny, nz;

    inline Vec3d mergeMinPoint(Vec3d v1, Vec3d v2){
        return Vec3d(min(v1.x, v2.x), min(v1.y, v2.y), min(v1.z, v2.z));
    }

    inline Vec3d mergeMaxPoint(Vec3d v1, Vec3d v2){
        return Vec3d(max(v1.x, v2.x), max(v1.y, v2.y), max(v1.z, v2.z));
    }

    Vec3d minCoordinates(){
        Vec3d tmp(1e10);
        for(int i = 0; i < source.size(); i++)
            tmp = mergeMinPoint(source[i].getMinPoint(), tmp);
        return tmp;
    }

    Vec3d maxCoordinates(){
        Vec3d tmp(-1e10);
        for(int i = 0; i < source.size(); i++)
            tmp = mergeMaxPoint(source[i].getMaxPoint(), tmp);
        return tmp;
    }
};


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

vector<Triangle> meshes;
vector<Material> materials;
vector<Light> lights;
vector<Ray> rays;

float *material_data, *light_data, *ray_data, *voxel_data;
vector<float> triangle_data;

const int IMAGE_VRES = 600, IMAGE_HRES = 600;

int main(int argc, char**argv){

    char* filename;
    if(argc == 2)
        filename = argv[1];
    else
        filename = "teddy.obj";


    Camera cam = Camera(IMAGE_HRES, IMAGE_VRES, Vec3d(2, 3, 1),
                       Vec3d(0,1,0), Vec3d(0,0,1), 200, 1.0);

    Material red_matte(Color(1,0,0),0.65, Color(),0,0,Color(),0);
    materials.push_back(red_matte);

    Triangle t(Vec3d(-5,0,0), Vec3d(0,0,10), Vec3d(7, 5, 0), 0);
    Triangle t2(Vec3d(6,4,6), Vec3d(6.5,4.5,9), Vec3d(5.5,3,9.5), 0);

    meshes.push_back(t);
    meshes.push_back(t2);

    PV(t.getMinPoint()); PV(t.getMaxPoint());

    Grid grid(meshes);
    grid.setup();

    /*
    cout << "teste"<<endl;
    cout<<endl<<endl<<endl;
    //pVar(grid.counts.size());
    //PT(grid.orderedCells[1][0]);
    for(int i = 0; i < grid.counts.size(); i++)
    for(int j = 0; j < grid.counts[i]; j++)
    {
        cout<<"------------------------------"<<endl;
        PT(grid.orderedCells[i][j]);
    }*/

    for(int r = 0; r < cam.vres; r++)
    for(int c = 0; c < cam.hres; c++)
        rays.push_back(cam.getRay(r, c));

    /*
        pVar(grid.totalTriangles);
        for(int j = 0; j < grid.totalTriangles*TRIANGLE_ATTRIBUTE_NUMBER; j++){cout<<j<<"---"<<triangle_data[j]<<endl;}
        for(int j = 0; j < grid.counts.size(); j++){cout<<"voxel "<<j<<" --> "<<endl<<endl; pVar(grid.counts[j]); }
    */

    //int gridNumValues = TRIANGLE_ATTRIBUTE_NUMBER*grid.totalTriangles,

    int rayNumValues = RAY_ATTRIBUTE_NUMBER*rays.size(),
        materialNumValues = MATERIAL_ATTRIBUTE_NUMBER*materials.size(),
        lightNumValues = LIGHT_ATTRIBUTE_NUMBER*lights.size();

    ray_data        = (float*) malloc(sizeof(float)*rayNumValues);
    material_data   = (float*) malloc(sizeof(float)*materialNumValues);
    light_data      = (float*) malloc(sizeof(float)*lightNumValues);
    voxel_data      = (float*) malloc(sizeof(float)*grid.counts.size());

    for(int i = 0; i < lights.size(); i++){
        light_data[i*LIGHT_ATTRIBUTE_NUMBER] = lights[i].position.x;
        light_data[i*LIGHT_ATTRIBUTE_NUMBER+1] = lights[i].position.y;
        light_data[i*LIGHT_ATTRIBUTE_NUMBER+2] = lights[i].position.z;

        light_data[i*LIGHT_ATTRIBUTE_NUMBER+3] = lights[i].intensity;

        light_data[i*LIGHT_ATTRIBUTE_NUMBER+4] = lights[i].color.r;
        light_data[i*LIGHT_ATTRIBUTE_NUMBER+5] = lights[i].color.g
        light_data[i*LIGHT_ATTRIBUTE_NUMBER+6] = lights[i].color.b;
    }


    for(int i = 0; i < materials.size(); i++){
        // k diffuse
        material_data[i*MATERIAL_ATTRIBUTE_NUMBER]      = materials[i].kd;
        material_data[i*MATERIAL_ATTRIBUTE_NUMBER+1]    = materials[i].cd.r;
        material_data[i*MATERIAL_ATTRIBUTE_NUMBER+2]    = materials[i].cd.g;
        material_data[i*MATERIAL_ATTRIBUTE_NUMBER+3]    = materials[i].cd.b;
    }

    for(int i = 0; i < rays.size(); i++){
        ray_data[i*RAY_ATTRIBUTE_NUMBER]    = rays[i].origin.x;
        ray_data[i*RAY_ATTRIBUTE_NUMBER+1]  = rays[i].origin.y;
        ray_data[i*RAY_ATTRIBUTE_NUMBER+2]  = rays[i].origin.z;

        ray_data[i*RAY_ATTRIBUTE_NUMBER+3]  = rays[i].direction.x;
        ray_data[i*RAY_ATTRIBUTE_NUMBER+4]  = rays[i].direction.y;
        ray_data[i*RAY_ATTRIBUTE_NUMBER+5]  = rays[i].direction.z;
    }

    for(int i = 0; i < grid.counts.size(); i++){
        // soma cumulativa de elementos em cada voxel
        if(i == 0)  voxel_data[i] = grid.counts[i];
        else        voxel_data[i] = voxel_data[i-1] + grid.counts[i];
        // copiando os triangulos para
        for(int j = 0; j < grid.counts[i]; j++){
            cout<<endl;
            pVar(i);
            PT(grid.orderedCells[i][j]);
            PV(grid.orderedCells[i][j].normal);
            cout<<"materialID = "<<grid.orderedCells[i][j].materialID<<endl<<endl;

            triangle_data.push_back(grid.orderedCells[i][j].p1.x);
            triangle_data.push_back(grid.orderedCells[i][j].p1.y);
            triangle_data.push_back(grid.orderedCells[i][j].p1.z);

            triangle_data.push_back(grid.orderedCells[i][j].p2.x);
            triangle_data.push_back(grid.orderedCells[i][j].p2.y);
            triangle_data.push_back(grid.orderedCells[i][j].p2.z);

            triangle_data.push_back(grid.orderedCells[i][j].p3.x);
            triangle_data.push_back(grid.orderedCells[i][j].p3.y);
            triangle_data.push_back(grid.orderedCells[i][j].p3.z);

            triangle_data.push_back(grid.orderedCells[i][j].normal.x);
            triangle_data.push_back(grid.orderedCells[i][j].normal.y);
            triangle_data.push_back(grid.orderedCells[i][j].normal.z);

            triangle_data.push_back(grid.orderedCells[i][j].materialID);

        }
    }




    return 0;
}
