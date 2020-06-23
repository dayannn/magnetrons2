#include <iostream>
#include <math.h>
#include <vector>
#include "omp.h"
#include <ctime>
#include <time.h>

using namespace std;

//радиус барабана
double Rb = 10;

// радиус объекта
double Ro = 1;

// высота объекта
double Ho = 10;

// угловая скорость барабана
double Wb = 0.175;                    // рад/с ~10 град в сек

// угловая скорость объекта
double Wo = 0.35;                  // рад/с ~10 град в сек

// расстояние до мишени
double L = 20;

// ширина мишени
double Wm = 8;

// высота мишени
double Hm = 15;

// количество объектов
double obj_num = 8;

//время полного оборота барабана
double turn_time = 6.283/ abs(Wb);

//# количество измерений за 1 оборот
int time_n = 40;

// # количество оборотов барабана
int turns_num = 20;

double time_delta = turn_time / time_n;

double calc_num = time_n * turns_num;

struct point2D {
    double x;
    double y;
};

struct point3D {
    double x;
    double y;
    double z;
};

point2D starting_point_b = {0, Rb};
point2D starting_point_o = {0, Ro};




// angle in radians
point2D rotate_point_0 (point2D point, double angle){
    point2D res;
    double Xnew = point.x * cos(angle) - point.y * sin(angle);
    double Ynew = point.x * sin(angle) + point.y * cos(angle);

    res.x = Xnew;
    res.y = Ynew;
    return res;
}

double get_current_angle (double w, double t){
    return fmod(w*t, 2*M_PI);
}

point2D get_object_center (double Wb, point2D point_b, double t){
    double alpha_b = get_current_angle(Wb, t);
    point2D vec_b = rotate_point_0(point_b, alpha_b);
    return vec_b;
}


// point_b  начальное опложение центра объекта, point_o - начальное положение точки на объекте относительно его центра
point2D get_point_coords (double Wb, double Wo, point2D point_b, point2D point_o, double t){
    point2D vec_b = get_object_center(Wb, point_b, t);
    double alpha_o = get_current_angle(Wo, t);
    point2D vec_o = rotate_point_0(point_o, alpha_o);
    point2D res = {vec_b.x + vec_o.x, vec_b.y + vec_o.y};
    return res;
}

point2D get_normal_vec(point2D center, point2D point){
    point2D res = {point.x - center.x, point.y - center.y};
    return res;
}

point3D get_normal_vec(point3D center, point3D point){
    point3D res = {point.x - center.x, point.y - center.y, point.z - center.z};
    return res;
}

double get_radians(double degrees){
    return  ( degrees * M_PI ) / 180;
}

point3D from_2d(point2D p) {
    point3D res = {p.x, p.y, 0};
    return res;
}

point3D from_2d(point2D p, double z) {
    point3D res = {p.x, p.y, z};
    return res;
}

double get_cosin(point2D v1, point2D v2) {
    double scalprod = v1.x * v2.x + v1.y * v2.y;
    double modX = sqrt(v1.x * v1.x + v1.y * v1.y);
    double modY = sqrt(v2.x * v2.x + v2.y * v2.y);
    return scalprod/modX/modY;
}

double get_cosin(point3D v1, point3D v2) {
    double scalprod = v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
    double modX = sqrt(v1.x * v1.x + v1.y * v1.y + v1.z * v1.z);
    double modY = sqrt(v2.x * v2.x + v2.y * v2.y + v2.z + v2.z);
    return scalprod/modX/modY;
}

double get_sin(point2D v1, point2D v2) {
    return sqrt(1 - pow(get_cosin(v1,v2), 2));
}

double get_sin(point3D v1, point3D v2) {
    return sqrt(1 - pow(get_cosin(v1,v2), 2));
}

point3D diff (point3D a, point3D b){
    point3D res;
    res.x = a.x - b.x;
    res.y = a.y - b.y;
    res.z = a.z - b.z;
    return res;
}

point3D sum (point3D a, point3D b){
    point3D res;
    res.x = a.x + b.x;
    res.y = a.y + b.y;
    res.z = a.z + b.z;
    return res;
}

point2D sum (point2D a, point2D b){
    point2D res;
    res.x = a.x + b.x;
    res.y = a.y + b.y;
    return res;
}

point3D mult (point3D p, double num){
    point3D res;
    res.x = p.x * num;
    res.y = p.y * num;
    res.z = p.z * num;
    return res;
}

double get_function_value (point3D target_p, point3D object_p, point3D n){


    point3D n1 = {0, 0, -1};   // нормаль к плоскости мишени
    point3D s = diff(target_p, object_p);    // вектор он точки на объекте к точке на мишени
    point3D s1 = diff(object_p, target_p);
    double cos_teta = get_cosin(n, s);

    if (cos_teta < 0)            // точка на объекте не видима с точки на мишени
        return 0;

    double cos_teta1 = abs(get_cosin(n1, s1));

    double sin_teta1 = abs(get_sin(n1, s1));
    double res = (sin_teta1 * cos_teta);
    if (res < 0)
        cout << "Achtung! res = " << res << endl;
    return res; //# * Im0

}

double integrate_one_point(point3D point, point3D n, point3D t1, point3D t2){
    int N = 100;         //or a large number, the discretization step
    double integral = 0.0;

    double x_target = t1.x;
    double dy = (t2.y - t1.y)/N;
    double dz = (t2.z - t1.z)/N;

    for(int yi = 0; yi < N; yi++) {
        double y = t1.y  + yi * dy;     //the y-value


        for(int zi = 0; zi < N; zi++) {
            double z = t1.z + zi*dz;    //the z-value
            point3D curr_point = {x_target, y, z};

            integral += get_function_value(curr_point, point, n);  //function call
        }
    }
    return integral*dy*dz;
}

vector<point3D> get_target_coords (double L, double Wm, double Hm){
    point3D t1 = {L, Wm/2, Hm/2};
    point3D t2 = {L, -Wm/2, -Hm/2};

    vector<point3D> res(2);
    res[0] = t1;
    res[1] = t2;
    return res;
}

// point_o - начальное положение точки на объекте относительно его центра
point2D get_point_coords (point2D point_center, double Wo, double start_angle, double t, double Ro){
    double alpha_o = fmod(get_current_angle(Wo, t) + start_angle, 2*M_PI);
    point2D point_o = {0, Ro};
    point2D vec_o = rotate_point_0(point_o, alpha_o);
    return sum(point_center, vec_o);
}

vector<point3D> get_object_points(point2D point_center, double num, double Wo, double t, double Ro){
    vector<point3D> res;
    double delta = (2*M_PI) / (num);
    for (int i = 0; i < num; i++){
        point2D p = get_point_coords (point_center, Wo, delta*i, t, Ro);
        point3D pp = {p.x, p.y, 0};
        res.push_back(pp);
    }
    return res;
}


void print_res_vect(vector<double> v){
    cout << "vect: ";
    for (int i = 0; i < v.size(); i++){
        cout << v[i] << " ";
    }
    cout << endl;
}

int main()
{

    unsigned int start_time =  clock(); // начальное время
    double cur_time = 0;
    int obj_point_num = 48;      // колчество точек на цилиндре, в которых вычисляется толщина


    double start_omp;
    double end_omp;
    #if defined(_OPENMP)
       start_omp = omp_get_wtime();
    #endif

    vector<double> res(obj_point_num);

    vector<point3D> target_coords = get_target_coords (L, Wm, Hm);

    for (int i = 0; i < calc_num; i++){

        point2D cent = get_object_center(Wb, starting_point_b, cur_time);
        point3D cent_3d = from_2d(cent);
        vector<point3D> pts = get_object_points(cent, obj_point_num, Wo, cur_time, Ro);

        #pragma omp parallel for
        for (int j = 0; j < obj_point_num; j++) {
            point3D p = pts[j];
            point3D normal = get_normal_vec(cent_3d, p);

            res[j] = res[j] + time_delta * integrate_one_point(pts[j], normal, target_coords[0], target_coords[1]);
        }

        cur_time += time_delta;

       // print_res_vect(res);
        cout << "processed " << i << " out of " << calc_num << endl;
    }

    unsigned int end_time = clock(); // конечное время
    unsigned int time = end_time - start_time; // искомое время

    #if defined(_OPENMP)
       end_omp = omp_get_wtime();
    #endif

    cout << "time= " << ((float)time)/CLOCKS_PER_SEC << endl;

    #if defined(_OPENMP)
           cout << "omp time= " << end_omp - start_omp<< endl;
    #endif


    return 0;


    // 6.14  10
    // 15

}
