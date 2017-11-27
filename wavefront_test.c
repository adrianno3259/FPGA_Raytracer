#include<stdio.h>
#include<stdlib.h>


#define MAX_TRIANGLES 10000
#define MAX_VERTICES 30000

typedef struct {
    double x, y, z;
} Vec3d;

typedef struct{
    double x, y;
} Vec2d;

#define PV3(A) printf("%s = (%lf, %lf, %lf)\n", #A, A.x, A.y, A.z)
#define PV2(A) printf("%s = (%lf, %lf)\n", #A, A.x, A.y)

int main(){
    FILE* file = fopen("albertosaurus.obj", "r");
    char line[200];
    int v_count = 0, vt_count = 0, vn_count = 0;
    Vec3d v[MAX_VERTICES], vn[MAX_VERTICES];
    Vec2d vt[MAX_VERTICES];

    printf("teste 1---------------\n"); fflush(stdout);

    while(!feof(file)){
        if(fgets(line, 100, file) == NULL) break;
        printf("teste---- %s -----------\n", line); fflush(stdout);
        if(line[0] == 'v'){
            if(line[1] == ' '){ //vertex
                sscanf(line+1, "%lf %lf %lf\n", &v[v_count].x, &v[v_count].y, &v[v_count].z);
                if(v_count % 5000 == 0) PV3(v[v_count]);
                v_count++;
            } else if(line[1] == 'n'){
                sscanf(line+2, "%lf %lf %lf\n", &vn[vn_count].x, &vn[vn_count].y, &vn[vn_count].z);
                if(vn_count % 5000 == 0) PV3(vn[vn_count]);
                vn_count++;
            }else if(line[1] == 't'){
                sscanf(line+2, "%lf %lf\n", &vt[vt_count].x, &vt[vt_count].y);
                if(vt_count % 5000 == 0) PV2(vt[vt_count]);
                vt_count++;
            }else if(line[1] == 'f'){
                double a, b, c, d, e, f, g, h, i;
                sscanf(line+1, "%lf/%lf/%lf lf/%lf/%lf lf/%lf/%lf\n", &a,&b,&c,&d,&e,&f,&g,&h,&i);
                if(vt_count % 5000 == 0) printf("%lf/%lf/%lf lf/%lf/%lf lf/%lf/%lf\n", a,b,c,d,e,f,g,h,i);
            }
        }
    }

    printf("Counts:\nv: %d, vt: %d, vn: %d\n", v_count, vt_count, vn_count);
    fclose(file);


}
