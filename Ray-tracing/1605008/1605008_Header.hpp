#include <windows.h>
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <bits/stdc++.h>
#include <string.h>
#include <list>
#include <math.h>


#define pi (2*acos(0.0))

using namespace std;

struct point
{
	double x,y,z;
};

class Light
{
public:

    point light_pos;
    double color[3];

    Light(double x , double y , double z)
    {
        light_pos.x = x;
        light_pos.y = y;
        light_pos.z = z;
    }

    void setColor(double red , double green , double blue)
    {
        this->color[0] = red;
        this->color[1] = green;
        this->color[2] = blue;
    }

    void draw()
    {
        glPointSize(5);
        glColor3f(color[0] , color[1] , color[2]);

        glBegin(GL_POINTS);

        glVertex3f(light_pos.x , light_pos.y , light_pos.z);

        glEnd();
    }

    void print()
    {
        cout<<light_pos.x<<" "<<light_pos.y<<" "<<light_pos.z<<" "<<endl;
        int i;
        for(i=0 ; i<3 ; i++)
        {
            cout<<color[i]<<" ";
        }
        cout<<endl;
    }
};
extern vector<Light> lights;

class Object;
extern vector<Object*> objects;

extern int level_of_recursion;


double dotProduct(double vect_A[], double vect_B[])
{
    double product = 0;
    for (int i = 0; i < 3; i++)
        product = product + vect_A[i] * vect_B[i];
    return product;
}

point crossProduct(double vect_A[] , double vect_B[])
{
    point cross_P;

    cross_P.x = vect_A[1] * vect_B[2] - vect_A[2] * vect_B[1];
    cross_P.y = vect_A[2] * vect_B[0] - vect_A[0] * vect_B[2];
    cross_P.z = vect_A[0] * vect_B[1] - vect_A[1] * vect_B[0];

    return cross_P;
}

class Ray
{
public:

    point start;
    point dir;

    Ray(double start_x , double start_y , double start_z , double dir_x , double dir_y , double dir_z)
    {
        this->start.x = start_x;
        this->start.y = start_y;
        this->start.z = start_z;

        this->dir.x = dir_x;
        this->dir.y = dir_y;
        this->dir.z = dir_z;

        double value = sqrt(dir.x*dir.x + dir.y*dir.y + dir.z*dir.z);
        dir.x = dir.x / value;
        dir.y = dir.y / value;
        dir.z = dir.z / value;

        //cout<<start_x<<" "<<start_y<<" "<<start_z<<endl;
    }
};

class Object
{
public:

    point reference_point;
    double height , width , length;
    double color[3];
    double coEfficients[4];
    int shine;

    Object()
    {

    }

    virtual void draw()
    {

    }

    virtual void print()
    {

    }

    void setColor(double red , double green , double blue)
    {
        this->color[0] = red;
        this->color[1] = green;
        this->color[2] = blue;
    }

    void setShine(int shine)
    {
        this->shine = shine;
    }

    void setCoEfficients(double ambient , double diffuse , double specular ,double recursive_reflection)
    {
        this->coEfficients[0] = ambient;
        this->coEfficients[1] = diffuse;
        this->coEfficients[2] = specular;
        this->coEfficients[3] = recursive_reflection;
    }

    virtual double intersect(Ray r, double *colored, int level)
    {
        colored = NULL;
        return -1.0;
    }
};

class Sphere : public Object
{
public:
    double radius;

    Sphere(double center_x , double center_y , double center_z , double radius)
    {
        this->reference_point.x = center_x;
        this->reference_point.y = center_y;
        this->reference_point.z = center_z;
        this->radius = radius;
    }

    void draw()
    {
        double slices = 60;
        double stacks = 20;

        struct point points[100][100];
        int i,j;
        double h,r;
        //generate points
        for(i=0;i<=stacks;i++)
        {
            h=radius*sin(((double)i/(double)stacks)*(pi/2));
            r=radius*cos(((double)i/(double)stacks)*(pi/2));
            for(j=0;j<=slices;j++)
            {
                points[i][j].x=r*cos(((double)j/(double)slices)*2*pi);
                points[i][j].y=r*sin(((double)j/(double)slices)*2*pi);
                points[i][j].z=h;
            }
        }
        //draw quads using generated points
        for(i=0;i<stacks;i++)
        {
            glColor3f(color[0] , color[1] , color[2]);
            for(j=0;j<slices;j++)
            {
                glBegin(GL_QUADS);{
			    //upper hemisphere
                    glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
                    glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
                    glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
                    glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);
                //lower hemisphere
                    glVertex3f(points[i][j].x,points[i][j].y,-points[i][j].z);
                    glVertex3f(points[i][j+1].x,points[i][j+1].y,-points[i][j+1].z);
                    glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,-points[i+1][j+1].z);
                    glVertex3f(points[i+1][j].x,points[i+1][j].y,-points[i+1][j].z);
                }glEnd();
            }
        }
    }

    void print()
    {
        cout<<reference_point.x<<" "<<reference_point.y<<" "<<reference_point.z<<" "<<endl;
        cout<<radius<<endl;
        int i;
        for(i=0 ; i<3 ; i++)
        {
            cout<<color[i]<<" ";
        }
        cout<<endl;
        for(i=0 ; i<4 ; i++)
        {
            cout<<coEfficients[i]<<" ";
        }
        cout<<endl;
        cout<<shine<<endl;

    }

    double intersect(Ray r , double* colored , int level)
    {
        /*colored[0] = this->color[0];
        colored[1] = this->color[1];
        colored[2] = this->color[2];*/

        double t;

        double catch_color[3];

        point oc;
        oc.x = r.start.x - reference_point.x;
        oc.y = r.start.y - reference_point.y;
        oc.z = r.start.z - reference_point.z;

        double a,b,c;
        a = 1;

        double vect_a[] = {oc.x , oc.y , oc.z};
        double vect_b[] = {r.dir.x , r.dir.y , r.dir.z};
        b = 2.0 * dotProduct(vect_a , vect_b);

        //double vect_a[] = {oc.x , oc.y , oc.z};
        c = dotProduct(vect_a , vect_a) - radius*radius;

        double discriminant = b*b - 4*a*c;

        if(discriminant < 0.0)
        {
            return -1.0;
            //t = -1.0;
        }
        else
        {
            double numerator = -b - sqrt(discriminant);
            if (numerator > 0.0)
            {
                //return numerator / (2.0 * a);
                t = numerator / (2.0 * a);
            }

            else
            {
                numerator = -b + sqrt(discriminant);
                if (numerator > 0.0)
                {
                    //return numerator / (2.0 * a);
                    t = numerator / (2.0 * a);
                }
                else
                {
                    return -1.0;
                    //t = -1.0;
                }
            }
        }

        if(level == 0){
            return t;
        }
        else
        {
            point intersection_point;
            intersection_point.x = r.start.x + t*r.dir.x;
            intersection_point.y = r.start.y + t*r.dir.y;
            intersection_point.z = r.start.z + t*r.dir.z;

            catch_color[0] = this->color[0] * coEfficients[0];
            catch_color[1] = this->color[1] * coEfficients[0];
            catch_color[2] = this->color[2] * coEfficients[0];

            point normal_at_intersection_point;
            normal_at_intersection_point.x = intersection_point.x - reference_point.x;
            normal_at_intersection_point.y = intersection_point.y - reference_point.y;
            normal_at_intersection_point.z = intersection_point.z - reference_point.z;
            double value = sqrt(normal_at_intersection_point.x*normal_at_intersection_point.x + normal_at_intersection_point.y*normal_at_intersection_point.y + normal_at_intersection_point.z*normal_at_intersection_point.z);
            normal_at_intersection_point.x = normal_at_intersection_point.x / value;
            normal_at_intersection_point.y = normal_at_intersection_point.y / value;
            normal_at_intersection_point.z = normal_at_intersection_point.z / value;

            double check_normal_direction;
            double vect_check_1[] = {normal_at_intersection_point.x , normal_at_intersection_point.y , normal_at_intersection_point.z};
            double vect_check_2[] = {r.dir.x , r.dir.y , r.dir.z};
            check_normal_direction = dotProduct(vect_check_1 , vect_check_2);
            if(check_normal_direction > 0)
            {
                normal_at_intersection_point.x = -normal_at_intersection_point.x;
                normal_at_intersection_point.y = -normal_at_intersection_point.y;
                normal_at_intersection_point.z = -normal_at_intersection_point.z;
            }

            //cout<<catch_color[0]<<" "<<catch_color[1]<<" "<<catch_color[2]<<" "<<endl;
            for(int i=0 ; i<lights.size() ; i++)
            {
                //cout<<lights[0].light_pos.x<<endl;
                Ray rayl(lights[i].light_pos.x , lights[i].light_pos.y , lights[i].light_pos.z , intersection_point.x-lights[i].light_pos.x , intersection_point.y-lights[i].light_pos.y , intersection_point.z-lights[i].light_pos.z);

                double t_now;
                t_now = (intersection_point.x - rayl.start.x) / rayl.dir.x;
                double* dummy_again = new double[3];
                bool obscured = false;

                for(int k=0 ; k<objects.size() ; k++)
                {
                    double t_l = 0.0;
                    t_l = objects[k]->intersect(rayl , dummy_again , 0);
                    //cout<<t_l<<endl;
                    if(t_l > 0 && floor(t_l) < floor(t_now))
                    {
                        obscured = true;
                        break;
                    }
                }

                delete dummy_again;

                if(obscured)
                {
                    continue;
                }

                double lambart_value;
                double vect_c[] = {normal_at_intersection_point.x , normal_at_intersection_point.y , normal_at_intersection_point.z};
                double vect_d[] = {rayl.dir.x , rayl.dir.y , rayl.dir.z};
                lambart_value = dotProduct(vect_c , vect_d);
                //cout<<lambart_value<<endl;

                /*point rayr;
                rayr.x = 2*lambart_value*normal_at_intersection_point.x - rayl.dir.x;
                rayr.y = 2*lambart_value*normal_at_intersection_point.y - rayl.dir.y;
                rayr.z = 2*lambart_value*normal_at_intersection_point.z - rayl.dir.z;*/
                point rayr;
                rayr.x = rayl.dir.x - 2*lambart_value*normal_at_intersection_point.x;
                rayr.y = rayl.dir.y - 2*lambart_value*normal_at_intersection_point.y;
                rayr.z = rayl.dir.z - 2*lambart_value*normal_at_intersection_point.z;

                double phong_value;
                double vect_e[] = {rayr.x , rayr.y , rayr.z};
                double vect_f[] = {-r.dir.x , -r.dir.y , -r.dir.z};
                phong_value = dotProduct(vect_e , vect_f);

                if(lambart_value < 0)
                {
                    lambart_value = 0;
                }
                if(phong_value < 0)
                {
                    phong_value = 0;
                }

                catch_color[0] = catch_color[0] + lights[i].color[0]*coEfficients[1]*abs(lambart_value)*this->color[0];
                catch_color[1] = catch_color[1] + lights[i].color[1]*coEfficients[1]*abs(lambart_value)*this->color[1];
                catch_color[2] = catch_color[2] + lights[i].color[2]*coEfficients[1]*abs(lambart_value)*this->color[2];

                /*catch_color[0] = catch_color[0] + lights[i].color[0]*coEfficients[1]*max(lambart_value , 0)*this->color[0];
                catch_color[1] = catch_color[1] + lights[i].color[1]*coEfficients[1]*max(lambart_value , 0)*this->color[1];
                catch_color[2] = catch_color[2] + lights[i].color[2]*coEfficients[1]*max(lambart_value , 0)*this->color[2];*/

                //cout<<lambart_value<<endl;

                /*catch_color[0] = catch_color[0] + lights[i].color[0]*coEfficients[2]* pow(abs(phong_value) , shine) *this->color[0];
                catch_color[1] = catch_color[1] + lights[i].color[1]*coEfficients[2]* pow(abs(phong_value) , shine) *this->color[1];
                catch_color[2] = catch_color[2] + lights[i].color[2]*coEfficients[2]* pow(abs(phong_value) , shine) *this->color[2];*/

                catch_color[0] = catch_color[0] + lights[i].color[0]*coEfficients[2]* pow(abs(phong_value) , shine);
                catch_color[1] = catch_color[1] + lights[i].color[1]*coEfficients[2]* pow(abs(phong_value) , shine);
                catch_color[2] = catch_color[2] + lights[i].color[2]*coEfficients[2]* pow(abs(phong_value) , shine);

                /*catch_color[0] = catch_color[0] + lights[i].color[0]*coEfficients[2]* pow(max(phong_value , 0) , shine);
                catch_color[1] = catch_color[1] + lights[i].color[1]*coEfficients[2]* pow(max(phong_value , 0) , shine);
                catch_color[2] = catch_color[2] + lights[i].color[2]*coEfficients[2]* pow(max(phong_value , 0) , shine);*/

            }

            //cout<<catch_color[0]<<" "<<catch_color[1]<<" "<<catch_color[2]<<" "<<endl;

            colored[0] = catch_color[0];
            colored[1] = catch_color[1];
            colored[2] = catch_color[2];

            //return t;
            if(level >= level_of_recursion)
            {
                return t;
            }

            point reflection_direction;
            double lambart_value_like_before;
            double vect_g[] = {normal_at_intersection_point.x , normal_at_intersection_point.y , normal_at_intersection_point.z};
            double vect_h[] = {r.dir.x , r.dir.y , r.dir.z};
            lambart_value_like_before = dotProduct(vect_g , vect_h);
            reflection_direction.x = r.dir.x - 2*lambart_value_like_before*normal_at_intersection_point.x;
            reflection_direction.y = r.dir.y - 2*lambart_value_like_before*normal_at_intersection_point.y;
            reflection_direction.z = r.dir.z - 2*lambart_value_like_before*normal_at_intersection_point.z;

            point advanced;
            advanced.x = intersection_point.x + reflection_direction.x;
            advanced.y = intersection_point.y + reflection_direction.y;
            advanced.z = intersection_point.z + reflection_direction.z;

            Ray ray_reflected(advanced.x , advanced.y , advanced.z , reflection_direction.x , reflection_direction.y , reflection_direction.z);

            int nearest = -1;
            double t_check , t_min = 100000;
            double* color_of_reflection = new double[3];
            double* dummy_yes_I_am = new double[3];

            for(int l=0 ; l<objects.size() ; l++)
            {
                t_check = objects[l]->intersect(ray_reflected , dummy_yes_I_am , 0);
                if(t_check > 0 && t_check < t_min)
                {
                    t_min = t_check;
                    nearest = l;
                    //cout<<"bingo"<<endl;
                }
            }

            if(nearest != -1)
            {
                t_min = objects[nearest]->intersect(ray_reflected , color_of_reflection , level+1);

                colored[0] = colored[0] + color_of_reflection[0]*coEfficients[3];
                colored[1] = colored[1] + color_of_reflection[1]*coEfficients[3];
                colored[2] = colored[2] + color_of_reflection[2]*coEfficients[3];
            }

            delete color_of_reflection;
            delete dummy_yes_I_am;

            return t;

        }
    }
};

class Triangle : public Object
{
public:

    point first_point;
    point second_point;
    point third_point;

    Triangle(double x1 , double y1 , double z1 , double x2 , double y2 , double z2 , double x3 , double y3 , double z3)
    {
        this->first_point.x = x1;
        this->first_point.y = y1;
        this->first_point.z = z1;

        this->second_point.x = x2;
        this->second_point.y = y2;
        this->second_point.z = z2;

        this->third_point.x = x3;
        this->third_point.y = y3;
        this->third_point.z = z3;

        this->reference_point.x = 0;
        this->reference_point.y = 0;
        this->reference_point.z = 0;
    }

    void draw()
    {
        glColor3f(color[0] , color[1] , color[2]);

        glBegin(GL_TRIANGLES);
        {
            glVertex3f(first_point.x , first_point.y , first_point.z);
			glVertex3f(second_point.x , second_point.y , second_point.z);
			glVertex3f(third_point.x , third_point.y , third_point.z);
        }
        glEnd();
    }

    void print()
    {
        cout<<first_point.x<<" "<<first_point.y<<" "<<first_point.z<<" "<<endl;
        cout<<second_point.x<<" "<<second_point.y<<" "<<second_point.z<<" "<<endl;
        cout<<third_point.x<<" "<<third_point.y<<" "<<third_point.z<<" "<<endl;
        int i;
        for(i=0 ; i<3 ; i++)
        {
            cout<<color[i]<<" ";
        }
        cout<<endl;
        for(i=0 ; i<4 ; i++)
        {
            cout<<coEfficients[i]<<" ";
        }
        cout<<endl;
        cout<<shine<<endl;

    }

    double intersect(Ray r , double* colored , int level)
    {
        /*colored[0] = this->color[0];
        colored[1] = this->color[1];
        colored[2] = this->color[2];*/

        double catch_color[3];

        point edge1;
        point edge2;

        edge1.x = second_point.x - first_point.x;
        edge1.y = second_point.y - first_point.y;
        edge1.z = second_point.z - first_point.z;

        edge2.x = third_point.x - first_point.x;
        edge2.y = third_point.y - first_point.y;
        edge2.z = third_point.z - first_point.z;

        point h;
        double vect_a[] = {r.dir.x , r.dir.y , r.dir.z};
        double vect_b[] = {edge2.x , edge2.y , edge2.z};
        h = crossProduct(vect_a , vect_b);

        double a;
        double vect_c[] = {edge1.x , edge1.y , edge1.z};
        double vect_d[] = {h.x , h.y , h.z};
        a = dotProduct(vect_c , vect_d);

        if (a > -0.00001 && a < 0.00001)
        {
            return -1.0;
        }

        double f;
        f = 1/a;

        point s;
        s.x = r.start.x - first_point.x;
        s.y = r.start.y - first_point.y;
        s.z = r.start.z - first_point.z;

        double u;
        double vect_e[] = {s.x , s.y , s.z};
        u = f * dotProduct(vect_e , vect_d);

        if (u < 0.0 || u > 1.0){
            return -1.0;
        }

        point q;
        q = crossProduct(vect_e , vect_c);

        double v;
        double vect_f[] = {q.x , q.y , q.z};
        v = f * dotProduct(vect_a , vect_f);

        if (v < 0.0 || (u + v) > 1.0){
            return -1.0;
        }

        double t;
        t = f * dotProduct(vect_b , vect_f);

        if (t > 0.00001){
            //return t;
        }
        else{
            return -1.0;
        }

        if(level == 0){
            return t;
        }
        else
        {
            point intersection_point;
            intersection_point.x = r.start.x + t*r.dir.x;
            intersection_point.y = r.start.y + t*r.dir.y;
            intersection_point.z = r.start.z + t*r.dir.z;

            catch_color[0] = this->color[0] * coEfficients[0];
            catch_color[1] = this->color[1] * coEfficients[0];
            catch_color[2] = this->color[2] * coEfficients[0];

            point normal_at_intersection_point;
            normal_at_intersection_point = crossProduct(vect_c , vect_b);
            double value = sqrt(normal_at_intersection_point.x*normal_at_intersection_point.x + normal_at_intersection_point.y*normal_at_intersection_point.y + normal_at_intersection_point.z*normal_at_intersection_point.z);
            normal_at_intersection_point.x = normal_at_intersection_point.x / value;
            normal_at_intersection_point.y = normal_at_intersection_point.y / value;
            normal_at_intersection_point.z = normal_at_intersection_point.z / value;

            double check_normal_direction;
            double vect_check_1[] = {normal_at_intersection_point.x , normal_at_intersection_point.y , normal_at_intersection_point.z};
            double vect_check_2[] = {r.dir.x , r.dir.y , r.dir.z};
            check_normal_direction = dotProduct(vect_check_1 , vect_check_2);
            if(check_normal_direction > 0)
            {
                normal_at_intersection_point.x = -normal_at_intersection_point.x;
                normal_at_intersection_point.y = -normal_at_intersection_point.y;
                normal_at_intersection_point.z = -normal_at_intersection_point.z;
            }

            //cout<<catch_color[0]<<" "<<catch_color[1]<<" "<<catch_color[2]<<" "<<endl;
            for(int i=0 ; i<lights.size() ; i++)
            {
                //cout<<lights[0].light_pos.x<<endl;
                Ray rayl(lights[i].light_pos.x , lights[i].light_pos.y , lights[i].light_pos.z , intersection_point.x-lights[i].light_pos.x , intersection_point.y-lights[i].light_pos.y , intersection_point.z-lights[i].light_pos.z);

                double t_now;
                t_now = (intersection_point.x - rayl.start.x) / rayl.dir.x;
                double* dummy_again = new double[3];
                bool obscured = false;

                for(int k=0 ; k<objects.size() ; k++)
                {
                    double t_l = 0.0;
                    t_l = objects[k]->intersect(rayl , dummy_again , 0);
                    //cout<<t_l<<endl;
                    if(t_l > 0 && floor(t_l) < floor(t_now))
                    {
                        obscured = true;
                        break;
                    }
                }

                delete dummy_again;

                if(obscured)
                {
                    continue;
                }

                double lambart_value;
                double vect_u[] = {normal_at_intersection_point.x , normal_at_intersection_point.y , normal_at_intersection_point.z};
                double vect_v[] = {rayl.dir.x , rayl.dir.y , rayl.dir.z};
                lambart_value = dotProduct(vect_u , vect_v);
                //cout<<lambart_value<<endl;

                /*point rayr;
                rayr.x = 2*lambart_value*normal_at_intersection_point.x - rayl.dir.x;
                rayr.y = 2*lambart_value*normal_at_intersection_point.y - rayl.dir.y;
                rayr.z = 2*lambart_value*normal_at_intersection_point.z - rayl.dir.z;*/
                point rayr;
                rayr.x = rayl.dir.x - 2*lambart_value*normal_at_intersection_point.x;
                rayr.y = rayl.dir.y - 2*lambart_value*normal_at_intersection_point.y;
                rayr.z = rayl.dir.z - 2*lambart_value*normal_at_intersection_point.z;

                double phong_value;
                double vect_w[] = {rayr.x , rayr.y , rayr.z};
                double vect_x[] = {-r.dir.x , -r.dir.y , -r.dir.z};
                phong_value = dotProduct(vect_w , vect_x);

                if(lambart_value < 0)
                {
                    lambart_value = 0;
                }
                if(phong_value < 0)
                {
                    phong_value = 0;
                }

                catch_color[0] = catch_color[0] + lights[i].color[0]*coEfficients[1]*abs(lambart_value)*this->color[0];
                catch_color[1] = catch_color[1] + lights[i].color[1]*coEfficients[1]*abs(lambart_value)*this->color[1];
                catch_color[2] = catch_color[2] + lights[i].color[2]*coEfficients[1]*abs(lambart_value)*this->color[2];

                //cout<<lights[i].color[0]<<endl;

                /*catch_color[0] = catch_color[0] + lights[i].color[0]*coEfficients[2]* pow(abs(phong_value) , shine) *this->color[0];
                catch_color[1] = catch_color[1] + lights[i].color[1]*coEfficients[2]* pow(abs(phong_value) , shine) *this->color[1];
                catch_color[2] = catch_color[2] + lights[i].color[2]*coEfficients[2]* pow(abs(phong_value) , shine) *this->color[2];*/

                catch_color[0] = catch_color[0] + lights[i].color[0]*coEfficients[2]* pow(abs(phong_value) , shine);
                catch_color[1] = catch_color[1] + lights[i].color[1]*coEfficients[2]* pow(abs(phong_value) , shine);
                catch_color[2] = catch_color[2] + lights[i].color[2]*coEfficients[2]* pow(abs(phong_value) , shine);

            }

            //cout<<catch_color[0]<<" "<<catch_color[1]<<" "<<catch_color[2]<<" "<<endl;

            colored[0] = catch_color[0];
            colored[1] = catch_color[1];
            colored[2] = catch_color[2];

            //return t;
            if(level >= level_of_recursion)
            {
                return t;
            }

            point reflection_direction;
            double lambart_value_like_before;
            double vect_g[] = {normal_at_intersection_point.x , normal_at_intersection_point.y , normal_at_intersection_point.z};
            double vect_h[] = {r.dir.x , r.dir.y , r.dir.z};
            lambart_value_like_before = dotProduct(vect_g , vect_h);
            reflection_direction.x = r.dir.x - 2*lambart_value_like_before*normal_at_intersection_point.x;
            reflection_direction.y = r.dir.y - 2*lambart_value_like_before*normal_at_intersection_point.y;
            reflection_direction.z = r.dir.z - 2*lambart_value_like_before*normal_at_intersection_point.z;

            point advanced;
            advanced.x = intersection_point.x + reflection_direction.x;
            advanced.y = intersection_point.y + reflection_direction.y;
            advanced.z = intersection_point.z + reflection_direction.z;

            Ray ray_reflected(advanced.x , advanced.y , advanced.z , reflection_direction.x , reflection_direction.y , reflection_direction.z);

            int nearest = -1;
            double t_check , t_min = 100000;
            double* color_of_reflection = new double[3];
            double* dummy_yes_I_am = new double[3];

            for(int l=0 ; l<objects.size() ; l++)
            {
                t_check = objects[l]->intersect(ray_reflected , dummy_yes_I_am , 0);
                if(t_check > 0 && t_check < t_min)
                {
                    t_min = t_check;
                    nearest = l;
                    //cout<<"bingo"<<endl;
                }
            }

            if(nearest != -1)
            {
                t_min = objects[nearest]->intersect(ray_reflected , color_of_reflection , level+1);

                colored[0] = colored[0] + color_of_reflection[0]*coEfficients[3];
                colored[1] = colored[1] + color_of_reflection[1]*coEfficients[3];
                colored[2] = colored[2] + color_of_reflection[2]*coEfficients[3];
            }

            delete color_of_reflection;
            delete dummy_yes_I_am;

            return t;
        }

    }
};

class Floor : public Object
{
public:
    double floor_size;
    double tile_size;

    Floor(double floorWidth , double tileWidth)
    {
        this->reference_point.x = -floorWidth/2;
        this->reference_point.y = -floorWidth/2;
        this->reference_point.z = 0;
        //this->length = tileWidth;
        this->floor_size = floorWidth;
        this->tile_size = tileWidth;

        this->coEfficients[0] = 0.3;
        this->coEfficients[1] = 0.2;
        this->coEfficients[2] = 0.2;
        this->coEfficients[3] = 0.3;

        this->shine = 5;
    }

    void draw()
    {
        bool f = false;

        for(int i=reference_point.x ; i<=floor_size/2 ; i+=tile_size)
        {
            for(int j=reference_point.y ; j<=floor_size/2 ; j+=tile_size)
            {
                if(f)
                {
                    glColor3f(0 , 0 , 0);
                    f = !f;
                }
                else
                {
                    glColor3f(1 , 1 , 1);
                    f = !f;
                }

                glBegin(GL_QUADS);
                {
                    glVertex3f(i , j , 0);
                    glVertex3f(i , j+tile_size , 0);
                    glVertex3f(i+tile_size , j+tile_size , 0);
                    glVertex3f(i+tile_size , j , 0);
                }
                glEnd();
            }
        }
    }

    void print()
    {
        cout<<floor_size<<" "<<tile_size<<" "<<endl;
    }

    double intersect(Ray r, double *colored, int level)
    {
        //cout<<"go"<<endl;
        /*colored[0] = 1;
        colored[1] = 1;
        colored[2] = 1;*/

        double catch_color[3];

        point intersection_point;

        point p0;
        p0.x = 0;
        p0.y = 0;
        p0.z = 0;

        point n;
        n.x = 0;
        n.y = 0;
        n.z = 1;

        //cout<<r.start.x<<" "<<r.start.y<<" "<<r.start.z<<endl;

        double denom;
        double vect_a[] = {n.x , n.y , n.z};
        double vect_b[] = {r.dir.x , r.dir.y , r.dir.z};
        denom = dotProduct(vect_a , vect_b);

        double t = -1.0;

        if(abs(denom) > .000001)
        {
            point p0_start;
            p0_start.x = p0.x - r.start.x;
            p0_start.y = p0.y - r.start.y;
            p0_start.z = p0.z - r.start.z;

            double vect_c[] = {p0_start.x , p0_start.y , p0_start.z};
            t = (dotProduct(vect_c , vect_a)) / denom;

            //point intersection_point;
            intersection_point.x = r.start.x + t*r.dir.x;
            intersection_point.y = r.start.y + t*r.dir.y;
            intersection_point.z = 0;

            if(intersection_point.x < -floor_size/2 || intersection_point.x > floor_size/2 || intersection_point.y < -floor_size/2 || intersection_point.y > floor_size/2)
            {
                return -1.0;
            }
            else
            {
                int row = floor((intersection_point.x - reference_point.x)/tile_size);
                int column = floor((intersection_point.y - reference_point.y)/tile_size);

                if((row + column)%2 != 0)
                {
                    this->color[0] = 0;
                    this->color[1] = 0;
                    this->color[2] = 0;
                }
                else
                {
                    this->color[0] = 1;
                    this->color[1] = 1;
                    this->color[2] = 1;
                }
            }
        }
        else
        {
            return -1;
        }

        //cout<<t<<endl;

        if(level == 0){
            return t;
        }
        else
        {
            /*point intersection_point;
            intersection_point.x = r.start.x + t*r.dir.x;
            intersection_point.y = r.start.y + t*r.dir.y;
            intersection_point.z = r.start.z + t*r.dir.z;*/

            catch_color[0] = this->color[0] * coEfficients[0];
            catch_color[1] = this->color[1] * coEfficients[0];
            catch_color[2] = this->color[2] * coEfficients[0];

            point normal_at_intersection_point;
            normal_at_intersection_point.x = 0;
            normal_at_intersection_point.y = 0;
            normal_at_intersection_point.z = 1;
            /*double value = sqrt(normal_at_intersection_point.x*normal_at_intersection_point.x + normal_at_intersection_point.y*normal_at_intersection_point.y + normal_at_intersection_point.z*normal_at_intersection_point.z);
            normal_at_intersection_point.x = normal_at_intersection_point.x / value;
            normal_at_intersection_point.y = normal_at_intersection_point.y / value;
            normal_at_intersection_point.z = normal_at_intersection_point.z / value;*/

            double check_normal_direction;
            double vect_check_1[] = {normal_at_intersection_point.x , normal_at_intersection_point.y , normal_at_intersection_point.z};
            double vect_check_2[] = {r.dir.x , r.dir.y , r.dir.z};
            check_normal_direction = dotProduct(vect_check_1 , vect_check_2);
            if(check_normal_direction > 0)
            {
                normal_at_intersection_point.x = -normal_at_intersection_point.x;
                normal_at_intersection_point.y = -normal_at_intersection_point.y;
                normal_at_intersection_point.z = -normal_at_intersection_point.z;
            }

            //cout<<catch_color[0]<<" "<<catch_color[1]<<" "<<catch_color[2]<<" "<<endl;
            for(int i=0 ; i<lights.size() ; i++)
            {
                //cout<<lights[0].light_pos.x<<endl;
                Ray rayl(lights[i].light_pos.x , lights[i].light_pos.y , lights[i].light_pos.z , intersection_point.x-lights[i].light_pos.x , intersection_point.y-lights[i].light_pos.y , intersection_point.z-lights[i].light_pos.z);

                double t_now;
                t_now = (intersection_point.x - rayl.start.x) / rayl.dir.x;
                double* dummy_again = new double[3];
                bool obscured = false;

                for(int k=0 ; k<objects.size() ; k++)
                {
                    double t_l = 0.0;
                    t_l = objects[k]->intersect(rayl , dummy_again , 0);
                    //cout<<t_l<<endl;
                    if(t_l > 0 && floor(t_l) < floor(t_now))
                    {
                        obscured = true;
                        break;
                    }
                }

                delete dummy_again;

                if(obscured)
                {
                    continue;
                }

                double lambart_value;
                double vect_u[] = {normal_at_intersection_point.x , normal_at_intersection_point.y , normal_at_intersection_point.z};
                double vect_v[] = {rayl.dir.x , rayl.dir.y , rayl.dir.z};
                lambart_value = dotProduct(vect_u , vect_v);
                //cout<<lambart_value<<endl;

                /*point rayr;
                rayr.x = 2*lambart_value*normal_at_intersection_point.x - rayl.dir.x;
                rayr.y = 2*lambart_value*normal_at_intersection_point.y - rayl.dir.y;
                rayr.z = 2*lambart_value*normal_at_intersection_point.z - rayl.dir.z;*/
                point rayr;
                rayr.x = rayl.dir.x - 2*lambart_value*normal_at_intersection_point.x;
                rayr.y = rayl.dir.y - 2*lambart_value*normal_at_intersection_point.y;
                rayr.z = rayl.dir.z - 2*lambart_value*normal_at_intersection_point.z;

                double phong_value;
                double vect_w[] = {rayr.x , rayr.y , rayr.z};
                double vect_x[] = {-r.dir.x , -r.dir.y , -r.dir.z};
                phong_value = dotProduct(vect_w , vect_x);

                if(lambart_value < 0)
                {
                    lambart_value = 0;
                }
                if(phong_value < 0)
                {
                    phong_value = 0;
                }

                catch_color[0] = catch_color[0] + lights[i].color[0]*coEfficients[1]*abs(lambart_value)*this->color[0];
                catch_color[1] = catch_color[1] + lights[i].color[1]*coEfficients[1]*abs(lambart_value)*this->color[1];
                catch_color[2] = catch_color[2] + lights[i].color[2]*coEfficients[1]*abs(lambart_value)*this->color[2];

                //cout<<this->color[0]<<" "<<lights[i].color[0]*coEfficients[1]*abs(lambart_value)*this->color[0]<<endl;

                /*catch_color[0] = catch_color[0] + lights[i].color[0]*coEfficients[2]* pow(abs(phong_value) , shine) *this->color[0];
                catch_color[1] = catch_color[1] + lights[i].color[1]*coEfficients[2]* pow(abs(phong_value) , shine) *this->color[1];
                catch_color[2] = catch_color[2] + lights[i].color[2]*coEfficients[2]* pow(abs(phong_value) , shine) *this->color[2];*/

                catch_color[0] = catch_color[0] + lights[i].color[0]*coEfficients[2]* pow(abs(phong_value) , shine);
                catch_color[1] = catch_color[1] + lights[i].color[1]*coEfficients[2]* pow(abs(phong_value) , shine);
                catch_color[2] = catch_color[2] + lights[i].color[2]*coEfficients[2]* pow(abs(phong_value) , shine);

            }

            //cout<<catch_color[0]<<" "<<catch_color[1]<<" "<<catch_color[2]<<" "<<endl;

            colored[0] = catch_color[0];
            colored[1] = catch_color[1];
            colored[2] = catch_color[2];

            //return t;
            if(level >= level_of_recursion)
            {
                return t;
            }

            point reflection_direction;
            double lambart_value_like_before;
            double vect_g[] = {normal_at_intersection_point.x , normal_at_intersection_point.y , normal_at_intersection_point.z};
            double vect_h[] = {r.dir.x , r.dir.y , r.dir.z};
            lambart_value_like_before = dotProduct(vect_g , vect_h);
            reflection_direction.x = r.dir.x - 2*lambart_value_like_before*normal_at_intersection_point.x;
            reflection_direction.y = r.dir.y - 2*lambart_value_like_before*normal_at_intersection_point.y;
            reflection_direction.z = r.dir.z - 2*lambart_value_like_before*normal_at_intersection_point.z;

            point advanced;
            advanced.x = intersection_point.x + reflection_direction.x;
            advanced.y = intersection_point.y + reflection_direction.y;
            advanced.z = intersection_point.z + reflection_direction.z;

            Ray ray_reflected(advanced.x , advanced.y , advanced.z , reflection_direction.x , reflection_direction.y , reflection_direction.z);

            int nearest = -1;
            double t_check , t_min = 100000;
            double* color_of_reflection = new double[3];
            double* dummy_yes_I_am = new double[3];

            for(int l=0 ; l<objects.size() ; l++)
            {
                t_check = objects[l]->intersect(ray_reflected , dummy_yes_I_am , 0);
                if(t_check > 0 && t_check < t_min)
                {
                    t_min = t_check;
                    nearest = l;
                    //cout<<"bingo"<<endl;
                }
            }

            if(nearest != -1)
            {
                t_min = objects[nearest]->intersect(ray_reflected , color_of_reflection , level+1);

                colored[0] = colored[0] + color_of_reflection[0]*coEfficients[3];
                colored[1] = colored[1] + color_of_reflection[1]*coEfficients[3];
                colored[2] = colored[2] + color_of_reflection[2]*coEfficients[3];
            }

            delete color_of_reflection;
            delete dummy_yes_I_am;

            return t;
        }
    }
};

class General : public Object
{
public:

    double A,B,C,D,E,F,G,H,I,J;
    //point cube_reference;

    General(double a , double b , double c , double d , double e , double f , double g , double h , double i , double j , double cube_x , double cube_y , double cube_z , double l , double w , double he)
    {
        this->A = a;
        this->B = b;
        this->C = c;
        this->D = d;
        this->E = e;
        this->F = f;
        this->G = g;
        this->H = h;
        this->I = i;
        this->J = j;

        this->reference_point.x = cube_x;
        this->reference_point.y = cube_y;
        this->reference_point.z = cube_z;

        this->length = l;
        this->width = w;
        this->height = he;

        this->reference_point.x = 0;
        this->reference_point.y = 0;
        this->reference_point.z = 0;
    }

    void draw()
    {

    }

    void print()
    {
        cout<<A<<" "<<B<<" "<<C<<" "<<D<<" "<<E<<" "<<F<<" "<<G<<" "<<H<<" "<<I<<" "<<J<<" "<<endl;
        cout<<reference_point.x<<" "<<reference_point.y<<" "<<reference_point.z<<" "<<endl;
        cout<<length<<" "<<width<<" "<<height<<" "<<endl;
        int i;
        for(i=0 ; i<3 ; i++)
        {
            cout<<color[i]<<" ";
        }
        cout<<endl;
        for(i=0 ; i<4 ; i++)
        {
            cout<<coEfficients[i]<<" ";
        }
        cout<<endl;
        cout<<shine<<endl;

    }

    double intersect(Ray r , double *colored , int level)
    {
        /*colored[0] = this->color[0];
        colored[1] = this->color[1];
        colored[2] = this->color[2];*/

        double catch_color[3];

        double rx = r.start.x;
        double ry = r.start.y;
        double rz = r.start.z;

        double dx = r.dir.x;
        double dy = r.dir.y;
        double dz = r.dir.z;

        double a = A*dx*dx + B*dy*dy + C*dz*dz + D*dx*dy + E*dy*dz + F*dz*dx;

        double b = 2*(A*rx*dx + B*ry*dy + C*rz*dz)
                 + D*(rx*dy + ry*dx) + E*(ry*dz + rz*dy) +  F*(rz*dx + rx*dz)
                 + G*dx + H*dy + I*dz;

        double c = A*rx*rx + B*ry*ry + C*rz*rz
                 + D*rx*ry + E*ry*rz + F*rz*rx
                 + G*rx + H*ry + I*rz + J;


        double discriminant = b*b - 4*a*c;


        if(discriminant < 0.0)
        {
            return -1.0;
        }

        double numerator1 = -b - sqrt(discriminant);
        double numerator2 = -b + sqrt(discriminant);

        double t1 = numerator1/(2.0 * a);
        double t2 = numerator2/(2.0 * a);


        double t = -1.0;

        if(t1 < 0 && t2 < 0)
        {
            return -1.0;
        }

        if(t1 < 0)
        {
            //t = t2;
            point intersecting_point2;
            intersecting_point2.x = rx + t2*dx;
            intersecting_point2.y = ry + t2*dy;
            intersecting_point2.z = rz + t2*dz;

            point temp2;
            temp2.x = intersecting_point2.x - reference_point.x;
            temp2.y = intersecting_point2.y - reference_point.y;
            temp2.z = intersecting_point2.z - reference_point.z;

            if(length > 0)
            {
                if(abs(temp2.x) > length)
                    return -1.0;

                else
                {
                    t = t2;
                }
            }

            if(width > 0)
            {
                if(abs(temp2.y) > width)
                    return -1.0;

                else
                {
                    t = t2;
                }
            }

            if(height > 0)
            {
                if(abs(temp2.z) > height)
                    return -1.0;

                else
                {
                    t = t2;
                }
            }
        }

        if(t2 < 0)
        {
            t = t1;
        }

        if(t1 >= 0 && t2 >= 0)
        {
            //point intersecting_point1(rx + t1*dx, ry + t1*dy, rz + t1*dz);
            //point intersecting_point2(rx + t2*dx, ry + t2*dy, rz + t2*dz);
            point intersecting_point1;
            intersecting_point1.x = rx + t1*dx;
            intersecting_point1.y = ry + t1*dy;
            intersecting_point1.z = rz + t1*dz;

            point intersecting_point2;
            intersecting_point2.x = rx + t2*dx;
            intersecting_point2.y = ry + t2*dy;
            intersecting_point2.z = rz + t2*dz;

            //Point temp1 = intersecting_point1.sub(reference_point);
            //Point temp2 = intersecting_point2.sub(reference_point);

            point temp1;
            temp1.x = intersecting_point1.x - reference_point.x;
            temp1.y = intersecting_point1.y - reference_point.y;
            temp1.z = intersecting_point1.z - reference_point.z;

            point temp2;
            temp2.x = intersecting_point2.x - reference_point.x;
            temp2.y = intersecting_point2.y - reference_point.y;
            temp2.z = intersecting_point2.z - reference_point.z;

            if(length > 0)
            {
                if(abs(temp1.x) > length && abs(temp2.x) > length)
                    return -1.0;

                if(abs(temp1.x) > length)
                    t = t2;

                else if(abs(temp2.x) > length)
                    t = t1;

                else
                    t = min(t1, t2);
            }

            if(width > 0)
            {
                if(abs(temp1.y) > width && abs(temp2.y) > width)
                    return -1.0;

                if(abs(temp1.y) > width)
                    t = t2;

                else if(abs(temp2.y) > width)
                    t = t1;

                else
                    t = min(t1, t2);
            }

            if(height > 0)
            {
                if(abs(temp1.z) > height && abs(temp2.z) > height)
                    return -1.0;

                if(abs(temp1.z) > height)
                    t = t2;

                else if(abs(temp2.z) > height)
                    t = t1;

                else
                    t = min(t1, t2);
            }

        }



        if(level == 0){
            return t;
        }
        else
        {
            point intersection_point;
            intersection_point.x = r.start.x + t*r.dir.x;
            intersection_point.y = r.start.y + t*r.dir.y;
            intersection_point.z = r.start.z + t*r.dir.z;

            double x, y, z;

            x  = intersection_point.x;
            y  = intersection_point.y;
            z  = intersection_point.z;

            catch_color[0] = this->color[0] * coEfficients[0];
            catch_color[1] = this->color[1] * coEfficients[0];
            catch_color[2] = this->color[2] * coEfficients[0];

            point normal_at_intersection_point;
            normal_at_intersection_point.x = 2*A*x + D*y + E*z + G;
            normal_at_intersection_point.y = 2*B*y + D*x + F*z + H;
            normal_at_intersection_point.z = 2*C*z + E*x + F*y + I;
            double value = sqrt(normal_at_intersection_point.x*normal_at_intersection_point.x + normal_at_intersection_point.y*normal_at_intersection_point.y + normal_at_intersection_point.z*normal_at_intersection_point.z);
            normal_at_intersection_point.x = normal_at_intersection_point.x / value;
            normal_at_intersection_point.y = normal_at_intersection_point.y / value;
            normal_at_intersection_point.z = normal_at_intersection_point.z / value;

            double check_normal_direction;
            double vect_check_1[] = {normal_at_intersection_point.x , normal_at_intersection_point.y , normal_at_intersection_point.z};
            double vect_check_2[] = {r.dir.x , r.dir.y , r.dir.z};
            check_normal_direction = dotProduct(vect_check_1 , vect_check_2);
            if(check_normal_direction > 0)
            {
                normal_at_intersection_point.x = -normal_at_intersection_point.x;
                normal_at_intersection_point.y = -normal_at_intersection_point.y;
                normal_at_intersection_point.z = -normal_at_intersection_point.z;
            }

            //cout<<catch_color[0]<<" "<<catch_color[1]<<" "<<catch_color[2]<<" "<<endl;
            for(int i=0 ; i<lights.size() ; i++)
            {
                //cout<<lights[0].light_pos.x<<endl;
                Ray rayl(lights[i].light_pos.x , lights[i].light_pos.y , lights[i].light_pos.z , intersection_point.x-lights[i].light_pos.x , intersection_point.y-lights[i].light_pos.y , intersection_point.z-lights[i].light_pos.z);

                double t_now;
                t_now = (intersection_point.x - rayl.start.x) / rayl.dir.x;
                double* dummy_again = new double[3];
                bool obscured = false;

                for(int k=0 ; k<objects.size() ; k++)
                {
                    double t_l = 0.0;
                    t_l = objects[k]->intersect(rayl , dummy_again , 0);
                    //cout<<t_l<<endl;
                    if(t_l > 0 && floor(t_l) < floor(t_now))
                    {
                        obscured = true;
                        break;
                    }
                }

                delete dummy_again;

                if(obscured)
                {
                    continue;
                }

                double lambart_value;
                double vect_u[] = {normal_at_intersection_point.x , normal_at_intersection_point.y , normal_at_intersection_point.z};
                double vect_v[] = {rayl.dir.x , rayl.dir.y , rayl.dir.z};
                lambart_value = dotProduct(vect_u , vect_v);
                //cout<<lambart_value<<endl;

                /*point rayr;
                rayr.x = 2*lambart_value*normal_at_intersection_point.x - rayl.dir.x;
                rayr.y = 2*lambart_value*normal_at_intersection_point.y - rayl.dir.y;
                rayr.z = 2*lambart_value*normal_at_intersection_point.z - rayl.dir.z;*/
                point rayr;
                rayr.x = rayl.dir.x - 2*lambart_value*normal_at_intersection_point.x;
                rayr.y = rayl.dir.y - 2*lambart_value*normal_at_intersection_point.y;
                rayr.z = rayl.dir.z - 2*lambart_value*normal_at_intersection_point.z;

                double phong_value;
                double vect_w[] = {rayr.x , rayr.y , rayr.z};
                double vect_x[] = {-r.dir.x , -r.dir.y , -r.dir.z};
                phong_value = dotProduct(vect_w , vect_x);

                if(lambart_value < 0)
                {
                    lambart_value = 0;
                }
                if(phong_value < 0)
                {
                    phong_value = 0;
                }

                catch_color[0] = catch_color[0] + lights[i].color[0]*coEfficients[1]*abs(lambart_value)*this->color[0];
                catch_color[1] = catch_color[1] + lights[i].color[1]*coEfficients[1]*abs(lambart_value)*this->color[1];
                catch_color[2] = catch_color[2] + lights[i].color[2]*coEfficients[1]*abs(lambart_value)*this->color[2];

                //cout<<this->color[0]<<" "<<lights[i].color[0]*coEfficients[1]*abs(lambart_value)*this->color[0]<<endl;

                /*catch_color[0] = catch_color[0] + lights[i].color[0]*coEfficients[2]* pow(abs(phong_value) , shine) *this->color[0];
                catch_color[1] = catch_color[1] + lights[i].color[1]*coEfficients[2]* pow(abs(phong_value) , shine) *this->color[1];
                catch_color[2] = catch_color[2] + lights[i].color[2]*coEfficients[2]* pow(abs(phong_value) , shine) *this->color[2];*/

                catch_color[0] = catch_color[0] + lights[i].color[0]*coEfficients[2]* pow(abs(phong_value) , shine);
                catch_color[1] = catch_color[1] + lights[i].color[1]*coEfficients[2]* pow(abs(phong_value) , shine);
                catch_color[2] = catch_color[2] + lights[i].color[2]*coEfficients[2]* pow(abs(phong_value) , shine);

            }

            //cout<<catch_color[0]<<" "<<catch_color[1]<<" "<<catch_color[2]<<" "<<endl;

            colored[0] = catch_color[0];
            colored[1] = catch_color[1];
            colored[2] = catch_color[2];

            //return t;
            if(level >= level_of_recursion)
            {
                return t;
            }

            point reflection_direction;
            double lambart_value_like_before;
            double vect_g[] = {normal_at_intersection_point.x , normal_at_intersection_point.y , normal_at_intersection_point.z};
            double vect_h[] = {r.dir.x , r.dir.y , r.dir.z};
            lambart_value_like_before = dotProduct(vect_g , vect_h);
            reflection_direction.x = r.dir.x - 2*lambart_value_like_before*normal_at_intersection_point.x;
            reflection_direction.y = r.dir.y - 2*lambart_value_like_before*normal_at_intersection_point.y;
            reflection_direction.z = r.dir.z - 2*lambart_value_like_before*normal_at_intersection_point.z;

            point advanced;
            advanced.x = intersection_point.x + reflection_direction.x;
            advanced.y = intersection_point.y + reflection_direction.y;
            advanced.z = intersection_point.z + reflection_direction.z;

            Ray ray_reflected(advanced.x , advanced.y , advanced.z , reflection_direction.x , reflection_direction.y , reflection_direction.z);

            int nearest = -1;
            double t_check , t_min = 100000;
            double* color_of_reflection = new double[3];
            double* dummy_yes_I_am = new double[3];

            for(int l=0 ; l<objects.size() ; l++)
            {
                t_check = objects[l]->intersect(ray_reflected , dummy_yes_I_am , 0);
                if(t_check > 0 && t_check < t_min)
                {
                    t_min = t_check;
                    nearest = l;
                    //cout<<"bingo"<<endl;
                }
            }

            if(nearest != -1)
            {
                t_min = objects[nearest]->intersect(ray_reflected , color_of_reflection , level+1);

                colored[0] = colored[0] + color_of_reflection[0]*coEfficients[3];
                colored[1] = colored[1] + color_of_reflection[1]*coEfficients[3];
                colored[2] = colored[2] + color_of_reflection[2]*coEfficients[3];
            }

            delete color_of_reflection;
            delete dummy_yes_I_am;

            return t;
        }
    }
};


//extern vector<Object*> objects;
//extern vector<Light> lights;


