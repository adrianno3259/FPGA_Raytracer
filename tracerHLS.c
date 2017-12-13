#include "tracerHLS.h"
#include <stdio.h>

#define TRIANGLE_ATTRIBUTE_NUMBER 13
#define LIGHT_ATTRIBUTE_NUMBER 7
#define MATERIAL_ATTRIBUTE_NUMBER 4
#define RAY_ATTRIBUTE_NUMBER 6
// data
float *m_data, *l_data, *r_data, *t_data;
int num_m, num_l, num_r, num_v, num_t, *v_data;

//bounding box da grid
float gridMin_x, gridMin_y, gridMin_z, gridMax_x, gridMax_y, gridMax_z;
int nx, ny, nz;

// interseção final
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


float* colorOut;
//int nColor;


void setlocalHitPoint(int rayID, float t){
    //printf("t = %f\n\n\n",t);
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
            c = r_data[rayID*RAY_ATTRIBUTE_NUMBER*RAY_ATTRIBUTE_NUMBER+3],
            d = tri_v1x(tID) - r_data[rayID*RAY_ATTRIBUTE_NUMBER];

    float  e = tri_v1y(tID) - tri_v2y(tID),
            f = tri_v1y(tID) - tri_v3y(tID),
            g = r_data[rayID*RAY_ATTRIBUTE_NUMBER*RAY_ATTRIBUTE_NUMBER+4],
            h = tri_v1y(tID) - r_data[rayID*RAY_ATTRIBUTE_NUMBER+1];

    float  i = tri_v1z(tID) - tri_v2z(tID),
            j = tri_v1z(tID) - tri_v3z(tID),
            k = r_data[rayID*RAY_ATTRIBUTE_NUMBER+5],
            l = tri_v1z(tID) - r_data[rayID*RAY_ATTRIBUTE_NUMBER+2];

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
                    printf("---- %f\n", tri_normaly(tID));
                    localT = tmin;
                    setlocalHitPoint(rayID, localT);
                    localNormal[0] = tri_normalx(tID),
                    localNormal[1] = tri_normaly(tID),
                    localNormal[2] = tri_normalz(tID);
                    localHit = 1;
                    localEntering = 1;
                    localMaterialID = tri_material(tID)
                }
            }
        }
    }
}


#define pInt(A) printf("%s = %d\n", #A, A); fflush(stdout);
#define pFloat(A) printf("%s = %f\n", #A, A); fflush(stdout);

int rayHitObjectsVoxel(int rayID, int voxelID){

    double tmin = 1000000.0;
    int i, id = 0;
    localHit = 0;
    hit = 0;
    // ultimo triangulo do voxel
    pInt(voxelID)
    int end = (voxelID > 0) ? (v_data[voxelID] - v_data[voxelID-1]): v_data[voxelID];
    pInt(end)

    for(i = 0; i < end; i++){
        int triangleID = (voxelID > 0) ? v_data[voxelID-1] + i : i;
        //printf("-------------t---------\n"); fflush(stdout);
        pInt(triangleID);
        triangleIntersect(triangleID, rayID);
        //printf("-------------t---------\n"); fflush(stdout);
        if(localHit && localT < tmin){
            //setlocalHitPoint(rayID, tmin);
            t = localT;
            tmin = localT;

            hitPoint[0] = localHitPoint[0],
            hitPoint[1] = localHitPoint[1],
            hitPoint[2] = localHitPoint[2];

            normal[0] = localNormal[0],
            normal[1] = localNormal[1],
            normal[2] = localNormal[2];

            materialID = localMaterialID;
            hit = 1;
            entering = 1;
            id = i;

        }
    }
    return hit;
}

float max3(float a, float b, float c){
    return (a > b ? (a > c ? a : c) : (b > c ? b : c));
}

float min3(float a, float b, float c){
    return (a < b ? (a < c ? a : c) : (b < c ? b : c));
}

int boundingBoxHit(float minx, float miny, float minz,
                    float maxx, float maxy, float maxz, int rayID){

    double  ox = r_data[ rayID*RAY_ATTRIBUTE_NUMBER ],
            oy = r_data[rayID*RAY_ATTRIBUTE_NUMBER+1],
            oz = r_data[rayID*RAY_ATTRIBUTE_NUMBER+2],
            dx = r_data[rayID*RAY_ATTRIBUTE_NUMBER+3],
            dy = r_data[rayID*RAY_ATTRIBUTE_NUMBER+4],
            dz = r_data[rayID*RAY_ATTRIBUTE_NUMBER+5];

    double tx_min, ty_min, tz_min;
    double tx_max, ty_max, tz_max;

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
    localT = t0;
    setlocalHitPoint(rayID, t0);
    return (t0 < t1 && t1 > 0.00001);
}

int rayInsideBox(int rayID){
    return  (r_data[ rayID*RAY_ATTRIBUTE_NUMBER ] > gridMin_x) &&
            (r_data[ rayID*RAY_ATTRIBUTE_NUMBER ] < gridMax_x) &&
            (r_data[rayID*RAY_ATTRIBUTE_NUMBER+1] > gridMin_y) &&
            (r_data[rayID*RAY_ATTRIBUTE_NUMBER+1] < gridMax_y) &&
            (r_data[rayID*RAY_ATTRIBUTE_NUMBER+2] > gridMin_z) &&
            (r_data[rayID*RAY_ATTRIBUTE_NUMBER+2] < gridMax_z);
}

float clampFloat(float x, float min, float max){
    return (x < min ? min : (x > max ? max : x));
}

void hitGrid(int rayID){

    int ix, iy, iz;
    float   ox = r_data[ rayID*RAY_ATTRIBUTE_NUMBER ],
            oy = r_data[rayID*RAY_ATTRIBUTE_NUMBER+1],
            oz = r_data[rayID*RAY_ATTRIBUTE_NUMBER+2],
            dx = r_data[rayID*RAY_ATTRIBUTE_NUMBER+3],
            dy = r_data[rayID*RAY_ATTRIBUTE_NUMBER+4],
            dz = r_data[rayID*RAY_ATTRIBUTE_NUMBER+5];



    if(rayInsideBox(rayID)){
        ix = clampFloat((ox - gridMin_x) * nx / (gridMax_x - gridMin_x), 0, nx - 1);
        iy = clampFloat((oy - gridMin_y) * ny / (gridMax_y - gridMin_y), 0, ny - 1);
        iz = clampFloat((oz - gridMin_z) * nz / (gridMax_z - gridMin_z), 0, nz - 1);
    } else {
        ix = clampFloat((localHitPoint[0] - gridMin_x) * nx / (gridMax_x - gridMin_x), 0, nx - 1);
        iy = clampFloat((localHitPoint[1] - gridMin_y) * ny / (gridMax_y - gridMin_y), 0, ny - 1);
        iz = clampFloat((localHitPoint[2] - gridMin_z) * nz / (gridMax_z - gridMin_z), 0, nz - 1);
    }

    float tx_min = gridMin_x, ty_min = gridMin_y, tz_min = gridMin_z;
    float tx_max = gridMax_x, ty_max = gridMax_y, tz_max = gridMax_z;

    float   dtx = (tx_max - tx_min)/nx,
            dty = (ty_max - ty_min)/ny,
            dtz = (tz_max - tz_min)/nz;

    float tx_next, ty_next, tz_next;
    int ix_step, iy_step, iz_step;
    int ix_stop, iy_stop, iz_stop;

    tx_next = tx_min + (ix + 1) * dtx;
    ix_step = 1;
    ix_stop = nx;

    ty_next = ty_min + (iy + 1) * dty;
    iy_step = 1;
    iy_stop = ny;

    tz_next = tz_min + (iz + 1) * dtz;
    iz_step = 1;
    iz_stop = nz;


    while(1){
        int voxelID = ix + nx*iy + nx*ny*iz;
        int size = (voxelID > 0) ? (v_data[voxelID] - v_data[voxelID-1]): v_data[voxelID];

        if(tx_next < ty_next && tx_next < tz_next){

            if(size > 0){
                rayHitObjectsVoxel(rayID, voxelID);
                if(hit && (t < tx_next)){
                    return 1;
                }
            }
            tx_next += dtx;
            ix += ix_step;

            if(ix == ix_stop) return 0;
        } else {

            if(ty_next < tz_next){
                if(size > 0){
                    rayHitObjectsVoxel(rayID, voxelID);
                    if(hit && (t < ty_next)){
                        return 1;
                    }
                }
                ty_next += dty;
                iy += iy_step
                if(iy == iy_stop) return 0;
            } else {
                if(size > 0){
                    rayHitObjectsVoxel(rayID, voxelID);
                    if(hit && (t < tz_next)){
                        return 1;
                    }
                }
                tz_next += dtz;
                iz += iz_step
                if(iz == iz_stop) return 0;
            }
        }

    }

}
/*
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
*/

#define setVec(V,X,Y,Z) V[0] = X, V[1] = Y, V[2] = Z
#define lenVec(V) sqrt(V[0]*V[0] + V[1]*V[1] + V[2]*V[2])
#define divideVec(V,S) V[0] /= S, V[1] /= S, V[2] /= S
#define dotVec(U,V) U[0]*V[0] + U[1]*V[1] + U[2]*V[2]

#define l_pos(LID) LIGHT_ATTRIBUTE_NUMBER*LID

#define setColor(C,R,G,B) colorOut[C*3] = R, colorOut[C*3 +1] = G, colorOut[C*3+2] = B

void matteShade(int matID){
    float wo[] = {-r_data[rayID+3], -r_data[rayID+4], -r_data[rayID+5]}; // vecInvert(i.r.direction);
    float L[] = {0, 0, 0} ;
    //Color L = colorMultiply(brdfRho(m.materialData.m.ambientBRDF, i, wo), lightGetL(ambient));
    int j;
    for(j = 0; j < num_l; j++){
        float wi[3];
        // light get direction to hit point
        int baseLight = l_pos(j);
        setVec(wi,  l_data[ baseLight ] - hitPoint[0],
                    l_data[baseLight+1] - hitPoint[1],
                    l_data[baseLight+2] - hitPoint[2]);

        float len = lenVec(wi);
        divideVec(wi, len); // normalize direction


        float ndotwi = dotVec(normal, wi); //vecDot(i.normal, wi);
        if(ndotwi > 0.0) L = colorAdd(L, colorScalarMultiply(ndotwi, colorAdd(brdfF(m.materialData.m.diffuseBDRF, i, wi, wo), lightGetL(lights[j]))));
    }
    return colorClamp(L);
}


// Salva no array global de cores o resultado do tracing
void rayTraceSimple(int rayID){
    double tmin = 1e15; int i;
    hitGrid(rayID);
    if(hit) matteShade(materialID);
    else colorOut[rayID*3] = colorOut[rayID*3+1] = colorOut[rayID*3+2] = 0.1;
}



void  tracer(int nVoxels, float* voxel_data, int nTriangles, float* triangle_data,
            int nMaterials, float* material_data, int nLights, float* light_data, int nRays, float* ray_data,
            float gridMinPoint_x, float gridMinPoint_y, float gridMinPoint_z,
            float gridMaxPoint_x, float gridMaxPoint_y, float gridMaxPoint_z,
            int n_x, int n_y, int n_z,
            float* colors) // return
{
    v_data = voxel_data, t_data = triangle_data, m_data = material_data, l_data = light_data;
    num_v = nVoxels, num_t = nTriangles, num_m = nMaterials, num_l = nLights;
    gridMax_x = gridMaxPoint_x, gridMax_y = gridMaxPoint_y, gridMax_z = gridMaxPoint_z;
    gridMin_x = gridMinPoint_x, gridMin_y = gridMinPoint_y, gridMin_z = gridMinPoint_z;
    nx = n_x, ny = n_y, nz = n_z;
    colorOut = colors;

    int i;
    for(i = 0; i < nRays; i++){
        rayTraceSimple(i);
    }

}



int main(){

    // dados externos:  v_data, r_data, t_data, m_data, l_data
    //                  num_v, num_r, num_t, num_m, num_l
    //                  grid_min_max_xyz, nx, ny, nz

    float triangleData[TRIANGLE_ATTRIBUTE_NUMBER*2];
    float rayData[RAY_ATTRIBUTE_NUMBER];
    int voxelData[2];

    voxelData[0] = 0;
    voxelData[1] = 2;
    /*

    Ray testRay = rayNew(vecNew(0, -1, 0), vecNew(0, 1, 0));
    Triangle testTriangle = triangleNew(vecNew(1,2,-1), vecNew(-1, 2, -1), vecNew(0, 2, 1));
    Intersect inter = triangleIntersect(testTriangle, testRay);
    PV(inter.hitPoint); printf("%s %f", "inter.t", inter.t);
    */

    triangleData[0] = 1, triangleData[1] = 2, triangleData[2] = -1; // point 1
    triangleData[3] = -1, triangleData[4] = 2, triangleData[5] = -1; // point 2
    triangleData[6] = 0, triangleData[7] = 2, triangleData[8] = 1; // point 3
    triangleData[9] = 0, triangleData[10] = -1, triangleData[11] = 0; // normal
    triangleData[13] = 0;

    triangleData[TRIANGLE_ATTRIBUTE_NUMBER+0] = 1,
    triangleData[TRIANGLE_ATTRIBUTE_NUMBER+1] = 3,
    triangleData[TRIANGLE_ATTRIBUTE_NUMBER+2] = -1; // point 1

    triangleData[TRIANGLE_ATTRIBUTE_NUMBER+3] = -1,
    triangleData[TRIANGLE_ATTRIBUTE_NUMBER+4] = 3,
    triangleData[TRIANGLE_ATTRIBUTE_NUMBER+5] = -1; // point 2

    triangleData[TRIANGLE_ATTRIBUTE_NUMBER+6] = 0,
    triangleData[TRIANGLE_ATTRIBUTE_NUMBER+7] = 3,
    triangleData[TRIANGLE_ATTRIBUTE_NUMBER+8] = 1; // point 3

    triangleData[TRIANGLE_ATTRIBUTE_NUMBER+9] = 0,
    triangleData[TRIANGLE_ATTRIBUTE_NUMBER+10] = -1,
    triangleData[TRIANGLE_ATTRIBUTE_NUMBER+11] = 0; // normal
    triangleData[TRIANGLE_ATTRIBUTE_NUMBER+13] = 0;

    rayData[0] = 0, rayData[1] = -1, rayData[2] = 0;
    rayData[3] = 0, rayData[4] = 1, rayData[5] = 0;

    r_data = rayData, t_data = triangleData, v_data = voxelData;
    num_r = 1, num_t = 2, num_v = 2;

    //triangleIntersect(0, 0);
    rayHitObjectsVoxel(0, 1);

    pFloat(t) pFloat(hitPoint[0]) pFloat(hitPoint[1]) pFloat(hitPoint[2])
    gridMin_x = -10, gridMin_y = -10, gridMin_z = -10;
    gridMax_x = 10, gridMax_y = 10, gridMax_z = 10;
    pInt(rayInsideBox(0));
    printf("t = %f\nhitPoint = %f, %f, %f\n", localT, localHitPoint[0], localHitPoint[1], localHitPoint[2]);
    printf("Normal = %f, %f, %f\n", localNormal[0], localNormal[1], localNormal[2]);

    //pInt(boundingBoxHit(-1, 0, -1, 1, 1, 1, 0))

    return 0;
}
