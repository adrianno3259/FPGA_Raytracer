#include "tracerHLS.h"

#define TRIANGLE_ATTRIBUTE_NUMBER 13
#define LIGHT_ATTRIBUTE_NUMBER 7
#define MATERIAL_ATTRIBUTE_NUMBER 4
#define RAY_ATTRIBUTE_NUMBER 6

float *m_data, *l_data, *r_data, *v_data, *t_data;
int num_m, num_l, num_r, num_v, num_t;

int hit, entering;
double t;
float hitPoint[3];
float normal[3];
int irayID;
int materialID;

//interseção  temporária para hitObjects()
int localHit, localEntering;
double localT;
float localHitPoint[3];
float localNormal[3];
int localRayID;
int localMaterialID;


float colorOut[1024*1024*3];
int nColor;


float max3(float a, float b, float c){
    return (a>b ? (a>c? a : c) : (b>c? b : c));
}


void setlocalHitPoint(int rayID, float t){
    localHitPoint[0] = t*r_data[rayID+3] + r_data[rayID];
    localHitPoint[1] = t*r_data[rayID+4] + r_data[rayID+1];
    localHitPoint[2] = t*r_data[rayID+5] + r_data[rayID+2];
}

#define tri_v1x(T) t_data[triBaseAddr]
#define tri_v1y(T) t_data[triBaseAddr+1]
#define tri_v1z(T) t_data[triBaseAddr+2]

#define tri_v2x(T) t_data[triBaseAddr+3]
#define tri_v2y(T) t_data[triBaseAddr+4]
#define tri_v2z(T) t_data[triBaseAddr+5]

#define tri_v3x(T) t_data[triBaseAddr+6]
#define tri_v3y(T) t_data[triBaseAddr+7]
#define tri_v3z(T) t_data[triBaseAddr+8]

#define tri_normalx(T) t_data[triBaseAddr+9]
#define tri_normaly(T) t_data[triBaseAddr+10]
#define tri_normalz(T) t_data[triBaseAddr+11]

#define tri_material(T) t_data[triBaseAddr+12]

void triangleIntersect(int tID, int rayID){

    localHit = 1;
    localEntering = 1;
    localT = 0.0;
    localHitPoint[0] = localHitPoint[1]  = localHitPoint[2] = 0;
    localNormal[0] = localNormal[1] = localNormal[2] = 0;
    localRayID = rayID;
    /*
    double a = t.v1.x - t.v2.x, b = t.v1.x - t.v3.x, c = ray.direction.x, d = t.v1.x - ray.origin.x;
    double e = t.v1.y - t.v2.y, f = t.v1.y - t.v3.y, g = ray.direction.y, h = t.v1.y - ray.origin.y;
    double i = t.v1.z - t.v2.z, j = t.v1.z - t.v3.z, k = ray.direction.z, l = t.v1.z - ray.origin.z;*/

    int triBaseAddr = tID*TRIANGLE_ATTRIBUTE_NUMBER;

    float  a = tri_v1x(tID) - tri_v2x(tID),
            b = tri_v1x(tID) - tri_v3x(tID),
            c = r_data[rayID+3],
            d = tri_v1x(tID) - r_data[rayID];

    float  e = tri_v1y(tID) - tri_v2y(tID),
            f = tri_v1y(tID) - tri_v3y(tID),
            g = r_data[rayID+4],
            h = tri_v1y(tID) - r_data[rayID+1];

    float  i = tri_v1z(tID) - tri_v2z(tID),
            j = tri_v1z(tID) - tri_v3z(tID),
            k = r_data[rayID+5],
            l = tri_v1z(tID) - r_data[rayID+2];

    float m = f*k - g*j, n = h*k - g*l, p = f*l - h*j;
    float q = g*i - e*k, s = e*j - f*i;
    float inv_denom = 1.0/(a*m + b*q + c*s);
    float e1 = d*m - b*n - c*p;
    float beta = e1*inv_denom;
    if(beta > 0.0)//      return;
    {
        float r = e*l - h*i;
        float e2 = a*n + d*q + c*r;
        float gamma = e2 * inv_denom;
        if(gamma > 0.0){
            if(beta + gamma < 1.0){
                float e3 = a*p - b*r + d*s;
                float tmin = e3 * inv_denom;
                if(tmin > 0.00001){
                    localT = tmin;
                    setlocalHitPoint(rayID, t);
                    normal[0] = tri_normalx(tID),
                    normal[1] = tri_normaly(tID),
                    normal[2] = tri_normalz(tID);
                    localHit = 1;
                    localEntering = 1;
                }
            }
        }
    }
}

void colorClamp()
{
    float maxVal = (r, std::max(g, b));
    if(maxVal > 1.0)
        r /= maxVal, g /= maxVal, b /= maxVal;
    return Color(r, g, b);
}

Vec3d lightGetDirection(Light l, Vec3d hitPoint){
    if(l.type == PointLight){
        //l.pos
        return vecNormalize(vecSubtract(l.posDir, hitPoint));
    }else
        return vecNew(0,0,0);
}

Color lightGetL(Light l){
    return colorScalarMultiply(l.intensity, l.color);
}

void matteShade(int matID){
    float wo[] = {r_data[rayID+3], r_data[rayID+4], r_data[rayID+5]}; // vecInvert(i.r.direction);
    float L[3];
    //Color L = colorMultiply(brdfRho(m.materialData.m.ambientBRDF, i, wo), lightGetL(ambient));
    int j;
    for(j = 0; j < num_l; j++){
        float wi[3];
        wi[0] = l_data[j], wi[1] = wi[2];// lightGetDirection(lights[j], i.hitPoint);
        double ndotwi = vecDot(i.normal, wi);
        if(ndotwi>0.0) L = colorAdd(L, colorScalarMultiply(ndotwi, colorAdd(brdfF(m.materialData.m.diffuseBDRF, i, wi, wo), lightGetL(lights[j]))));
    }
    return colorClamp(L);
}

void rayHitObjects(int rayID){
    Intersect it, res;
    double tmin = 1000000.0;
    int i, id = 0;
    localHit = false;

    for(i = 0; i < num_t; i++){
        triangleIntersect(rayID);
        if(localHit && localT < tmin){
            setlocalHitPoint(rayID, tmin);
            t = localT
            tmin = localT;

            hitPoint[0] = localHitPoint[0],
            hitPoint[1] = localHitPoint[1],
            hitPoint[2] = localHitPoint[2];

            normal[0] = localNormal[0],
            normal[1] = localNormal[1],
            normal[2] = localNormal[2];

            hit = true;
            entering = true;
            id = i;
        }
    }
}
/*
void rayHitObjectsVoxel(int rayID, int voxelID){

    double tmin = 1000000.0;
    int i, id = 0;
    localHit = false;

    for(i = 0; i < v_data[voxelID]; i++){
        int triangleID = (voxelID > 0) ? v_data[voxelID-1] + 1 + i : i;

        triangleIntersect(triangleID, rayID);

        if(localHit && localT < tmin){
            setlocalHitPoint(rayID, tmin);
            t = localT
            tmin = localT;

            hitPoint[0] = localHitPoint[0],
            hitPoint[1] = localHitPoint[1],
            hitPoint[2] = localHitPoint[2];

            normal[0] = localNormal[0],
            normal[1] = localNormal[1],
            normal[2] = localNormal[2];

            hit = true;
            entering = true;
            id = i;
        }
    }
}

double max3(double a, double b, double c){
    return (a > b ? (a > c ? a : c) : (b > c ? b : c));
}

int boundingBoxHit(float minx, float miny, float minz,
                    float maxx, float maxy, float maxz, int rayID){

    double ox = r_data[rayID], oy = r_data[rayID+1], oz = r_data[rayID+2],
           dx = r_data[rayID+3], dy = r_data[rayID+4], dz = r_data[rayID+5];

    double tx_min, ty_min, tz_min;
    double tx_max, ty_max, tz_min;

    double a = 1.0/dx;
    if(a >= 0){
        tx_min = (minx - ox)*a;
        tx_max = (maxx - ox)*a;
    } else {
        tx_max = (minx - ox)*a;
        tx_min = (maxx - ox)*a;
    }

    double b = 1.0/dy;
    if(b >= 0){
        ty_min = (miny - oy)*b;
        ty_max = (maxy - oy)*b;
    } else {
        ty_max = (miny - oy)*b;
        ty_min = (maxy - oy)*b;
    }

    double c = 1.0/dz;
    if(c >= 0){
        tz_min = (minz - oz)*c;
        tz_max = (maxz - oz)*c;
    } else {
        tz_max = (minz - oz)*c;
        tz_min = (maxz - oz)*c;
    }

    double t0, t1;

    t0 = max3(tx_min, ty_min, tz_min);
    t1 = min3(tx_max, ty_max, tz_max);

    return (t0 < t1 && t1 > K_EPSILON);
}

void hitGrid(int rayID){




}

// Salva no array global de cores o resultado do tracing
void rayTraceSimple(int rayID){
    double tmin = 1e15; int i;
    rayHitObjects(rayID);
    if(hit) matteShade(materialID);
    else colorOut[rayID*3] = colorOut[rayID*3+1] = colorOut[rayID*3+2] = 0.1;
}

void tracer(int nVoxels, float* voxel_data, int nTriangles, float* triangle_data,
            int nMaterials, float* material_data, int nLights, float* light_data, int nRays, float* ray_data,
            float gridMinPoint_x, float gridMinPoint_y, float gridMinPoint_z,
            float gridMaxPoint_x, float gridMaxPoint_y, float gridMaxPoint_z)
{
    v_data = voxel_data, t_data = triangle_data, m_data = material_data, l_data = light_data;
    num_v = nVoxels, num_t = nTriangles, num_m = nMaterials, num_l = nLights;

    int i;
    for(i = 0; i < nRays; i++){
        rayTraceSimple(i);
    }

}
