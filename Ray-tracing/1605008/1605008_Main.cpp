/*
 * GLUT Shapes Demo
 *
 * Written by Nigel Stewart November 2003
 *
 * This program is test harness for the sphere, cone
 * and torus shapes in GLUT.
 *
 * Spinning wireframe and smooth shaded shapes are
 * displayed until the ESC or q key is pressed.  The
 * number of geometry stacks and slices can be adjusted
 * using the + and - keys.
 */
#include <windows.h>
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <bits/stdc++.h>
#include <string.h>
#include <list>
#include <fstream>

#include "classes.hpp"
#include "bitmap_image.hpp"

#define pi (2*acos(0.0))

using namespace std;

double cameraHeight;
double cameraAngle;
int drawgrid;
int drawaxes;
double angle;

double drawpoint = 0;

/*struct point
{
	double x,y,z;
};*/

struct point_vector
{
    double x,y,z;
};

struct point pos;
struct point up;
struct point c_right;
struct point look;

double wall = 100;

int level_of_recursion;
int pixel_along_both_dimensions;

int total_objects;

int light_sources;

vector<Object*> objects;
vector<Light> lights;

double floor_width = 1000;
double tile_width = 20;

//vector<Ray> rays;

/*double dotProduct(double vect_A[], double vect_B[])
{
    double product = 0;
    for (int i = 0; i < 3; i++)
        product = product + vect_A[i] * vect_B[i];
    return product;
}*/

/*point_vector crossProduct(double vect_A[] , double vect_B[])
{
    point_vector cross_P;

    cross_P.x = vect_A[1] * vect_B[2] - vect_A[2] * vect_B[1];
    cross_P.y = vect_A[2] * vect_B[0] - vect_A[0] * vect_B[2];
    cross_P.z = vect_A[0] * vect_B[1] - vect_A[1] * vect_B[0];

    return cross_P;
}*/

point rotation(point vect , point axes , double rotation_angle)
{
    point after_rotation;

    double cos_value = cos(rotation_angle);
    double one_minus_cos = 1 - cos_value;
    double sin_value = sin(rotation_angle);

    double vect_v[] = {vect.x , vect.y , vect.z};
    double vect_k[] = {axes.x , axes.y , axes.z};

    point k_cross_v = crossProduct(vect_k , vect_v);

    double k_dot_v = dotProduct(vect_k , vect_v);

    after_rotation.x = vect.x*cos_value + k_cross_v.x*sin_value + axes.x*k_dot_v*one_minus_cos;
    after_rotation.y = vect.y*cos_value + k_cross_v.y*sin_value + axes.y*k_dot_v*one_minus_cos;
    after_rotation.z = vect.z*cos_value + k_cross_v.z*sin_value + axes.z*k_dot_v*one_minus_cos;

    return after_rotation;

}

void rotate_about_u(double camAngle)
{
    point look_rotated = rotation(look , up , camAngle);
    look.x = look_rotated.x;
    look.y = look_rotated.y;
    look.z = look_rotated.z;

    point right_rotated = rotation(c_right , up , camAngle);
    c_right.x = right_rotated.x;
    c_right.y = right_rotated.y;
    c_right.z = right_rotated.z;
}

void rotate_about_l(double camAngle)
{
    point up_rotated = rotation(up , look , camAngle);
    up.x = up_rotated.x;
    up.y = up_rotated.y;
    up.z = up_rotated.z;

    point right_rotated = rotation(c_right , look , camAngle);
    c_right.x = right_rotated.x;
    c_right.y = right_rotated.y;
    c_right.z = right_rotated.z;
}

void rotate_about_r(double camAngle)
{
    point look_rotated = rotation(look , c_right , camAngle);
    look.x = look_rotated.x;
    look.y = look_rotated.y;
    look.z = look_rotated.z;

    point up_rotated = rotation(up , c_right , camAngle);
    up.x = up_rotated.x;
    up.y = up_rotated.y;
    up.z = up_rotated.z;
}

void move_forward()
{
    pos.x = pos.x + look.x*2;
    pos.y = pos.y + look.y*2;
    pos.z = pos.z + look.z*2;
}

void move_backwards()
{
    pos.x = pos.x - look.x*2;
    pos.y = pos.y - look.y*2;
    pos.z = pos.z - look.z*2;
}

void move_right()
{
    pos.x = pos.x + c_right.x*2;
    pos.y = pos.y + c_right.y*2;
    pos.z = pos.z + c_right.z*2;
}

void move_left()
{
    pos.x = pos.x - c_right.x*2;
    pos.y = pos.y - c_right.y*2;
    pos.z = pos.z - c_right.z*2;
}

void move_up()
{
    pos.x = pos.x + up.x*2;
    pos.y = pos.y + up.y*2;
    pos.z = pos.z + up.z*2;
}

void move_down()
{
    pos.x = pos.x - up.x*2;
    pos.y = pos.y - up.y*2;
    pos.z = pos.z - up.z*2;
}

void drawPoint()
{
    glPointSize(5);
    glColor3f(1 , 0 , 0);

    glBegin(GL_POINTS);

    glVertex3f(0 , 0 , 1);

    glEnd();
}

void drawAxes()
{
	if(drawaxes==1)
	{
		glColor3f(1.0, 1.0, 1.0);
		glBegin(GL_LINES);{
			glVertex3f( 500,0,0);
			glVertex3f(-500,0,0);

			glVertex3f(0,-500,0);
			glVertex3f(0, 500,0);

			glVertex3f(0,0, 500);
			glVertex3f(0,0,-500);
		}glEnd();
	}
}


void drawGrid()
{
	int i;
	if(drawgrid==1)
	{
		glColor3f(0.6, 0.6, 0.6);	//grey
		glBegin(GL_LINES);{
			for(i=-8;i<=8;i++){

				if(i==0)
					continue;	//SKIP the MAIN axes

				//lines parallel to Y-axis
				glVertex3f(i*10, -90, 0);
				glVertex3f(i*10,  90, 0);

				//lines parallel to X-axis
				glVertex3f(-90, i*10, 0);
				glVertex3f( 90, i*10, 0);
			}
		}glEnd();
	}
}

void drawSquare(double a)
{
    glColor3f(0.8,0.8,0.8);
	glBegin(GL_QUADS);{
		glVertex3f( a, 300,a);
		glVertex3f( a,300,-a);
		glVertex3f(-a,300,-a);
		glVertex3f(-a,300,a);
	}glEnd();
}

void drawSphere(double radius,int slices,int stacks)
{
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
        glColor3f((double)i/(double)stacks,(double)i/(double)stacks,(double)i/(double)stacks);
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

void load_data()
{
    ifstream input;
    input.open("F:\\MATERIALS\\4-1\\graphics\\draw\\scene.txt");
    //input.open("scene.txt");
    //freopen("scene.txt" , "r" , stdin);

    if(! input)
    {
        cout<<"what ?";
    }

    input>>level_of_recursion;
    input>>pixel_along_both_dimensions;
    input>>total_objects;

    //cout<<pixel_along_both_dimensions;

    string command;

    int i;
    for(i=0 ; i<total_objects ; i++)
    {
        input>>command;

        if(command.compare("sphere") == 0){
            double center_x , center_y , center_z;
            double radius;
            double color_r , color_g , color_b;
            double amb , dif , spec , rec_ref;
            int shine;

            input>>center_x;
            input>>center_y;
            input>>center_z;

            input>>radius;

            input>>color_r;
            input>>color_g;
            input>>color_b;

            input>>amb;
            input>>dif;
            input>>spec;
            input>>rec_ref;

            input>>shine;

            Object* sphere;
            sphere = new Sphere(center_x , center_y , center_z , radius);
            sphere->setColor(color_r , color_g , color_b);
            sphere->setCoEfficients(amb , dif , spec , rec_ref);
            sphere->setShine(shine);

            objects.push_back(sphere);

            //sphere->print();
        }

        else if(command.compare("triangle") == 0){
            double x1 , y1 , z1 , x2 , y2 , z2 , x3 , y3 , z3;
            double color_r , color_g , color_b;
            double amb , dif , spec , rec_ref;
            int shine;

            input>>x1;
            input>>y1;
            input>>z1;
            input>>x2;
            input>>y2;
            input>>z2;
            input>>x3;
            input>>y3;
            input>>z3;

            input>>color_r;
            input>>color_g;
            input>>color_b;

            input>>amb;
            input>>dif;
            input>>spec;
            input>>rec_ref;

            input>>shine;

            Object* triangle;
            triangle = new Triangle(x1 , y1 , z1 , x2 , y2 , z2 , x3 , y3 , z3);
            triangle->setColor(color_r , color_g , color_b);
            triangle->setCoEfficients(amb , dif , spec , rec_ref);
            triangle->setShine(shine);

            objects.push_back(triangle);

            //triangle->print();
        }

        else if(command.compare("general") == 0){
            double A,B,C,D,E,F,G,H,I,J;
            double cube_x , cube_y , cube_z , l ,  w ,  he;
            double color_r , color_g , color_b;
            double amb , dif , spec , rec_ref;
            int shine;

            input>>A;
            input>>B;
            input>>C;
            input>>D;
            input>>E;
            input>>F;
            input>>G;
            input>>H;
            input>>I;
            input>>J;

            input>>cube_x;
            input>>cube_y;
            input>>cube_z;

            input>>l;
            input>>w;
            input>>he;

            input>>color_r;
            input>>color_g;
            input>>color_b;

            input>>amb;
            input>>dif;
            input>>spec;
            input>>rec_ref;

            input>>shine;

            Object* general;
            general = new General(A , B , C , D , E , F , G , H , I , J , cube_x , cube_y , cube_z , l , w , he);
            general->setColor(color_r , color_g , color_b);
            general->setCoEfficients(amb , dif , spec , rec_ref);
            general->setShine(shine);

            objects.push_back(general);

            //general->print();
        }

    }

    input>>light_sources;

    for(i=0 ; i<light_sources ; i++)
    {
        double x , y , z , r , g , b;

        input>>x;
        input>>y;
        input>>z;

        input>>r;
        input>>g;
        input>>b;

        Light light(x , y , z);
        light.setColor(r , g , b);
        lights.push_back(light);

        //light.print();
    }

    input.close();

    Object* floor;
    floor = new Floor(floor_width , tile_width);

    objects.push_back(floor);

    //floor->print();
}

void capture()
{
    double image_width = pixel_along_both_dimensions;
    double image_height = pixel_along_both_dimensions;
    //double image_width = 500.0;
    //double image_height = 500.0;

    bitmap_image image(image_width,image_height);

    /*for(int i=0;i<image_height;i++)
    {
        for(int j=0;j<image_width;j++)
        {
            image.set_pixel(j,i,0,0,0);
        }
    }*/
    //rays.clear();
    image.clear();

    double plane_distance = (500.0/2.0) / (tan((80 * (pi / 180))/2.0));

    point top_left;
    top_left.x = pos.x + look.x*plane_distance - (c_right.x*500.0)/2.0 + (up.x*500.0)/2.0;
    top_left.y = pos.y + look.y*plane_distance - (c_right.y*500.0)/2.0 + (up.y*500.0)/2.0;
    top_left.z = pos.z + look.z*plane_distance - (c_right.z*500.0)/2.0 + (up.z*500.0)/2.0;

    double du = 500.0/image_width;
    double dv = 500.0/image_height;

    top_left.x = top_left.x + c_right.x*(0.5*du) - up.x*(0.5*dv);
    top_left.y = top_left.y + c_right.y*(0.5*du) - up.y*(0.5*dv);
    top_left.z = top_left.z + c_right.z*(0.5*du) - up.z*(0.5*dv);

    //cout<<top_left.x<<" "<<top_left.y<<" "<<top_left.z<<endl;
    //cout<<du<<" "<<dv<<endl;

    /*cout<<pos.x<<" "<<pos.y<<" "<<pos.z<<endl;
    cout<<look.x<<" "<<look.y<<" "<<look.z<<endl;
    cout<<c_right.x<<" "<<c_right.y<<" "<<c_right.z<<endl;
    cout<<up.x<<" "<<up.y<<" "<<up.z<<endl;*/

    for(int i=0;i<image_width;i++)
    {
        for(int j=0;j<image_height;j++)
        {
            int nearest = -1;
            double t , t_min = 100000;

            point cur_pixel;
            cur_pixel.x = top_left.x + c_right.x*(i*du) - up.x*(j*dv);
            cur_pixel.y = top_left.y + c_right.y*(i*du) - up.y*(j*dv);
            cur_pixel.z = top_left.z + c_right.z*(i*du) - up.z*(j*dv);

            //cout<<pos.x<<" "<<pos.y<<" "<<pos.z<<endl;
            Ray ray(pos.x , pos.y , pos.z , cur_pixel.x-pos.x , cur_pixel.y-pos.y , cur_pixel.z-pos.z);
            //rays.push_back(ray);

            double* color = new double[3];
            double* dummy = new double[3];

            for(int k=0 ; k<objects.size() ; k++)
            {
                t = objects[k]->intersect(ray , dummy , 0);
                if(t > 0 && t < t_min)
                {
                    t_min = t;
                    nearest = k;
                    //cout<<"bingo"<<endl;
                }
            }

            if(nearest != -1)
            {
                t_min = objects[nearest]->intersect(ray , color , 1);
                //cout<<color[0]<<" "<<color[1]<<" "<<color[2]<<endl;
                //cout<<nearest<<endl;
                if(color[0] > 1)
                {
                    color[0] = 1;
                }
                if(color[1] > 1)
                {
                    color[1] = 1;
                }
                if(color[2] > 1)
                {
                    color[2] = 1;
                }

                image.set_pixel(i,j,color[0]*255,color[1]*255,color[2]*255);
            }

            delete dummy;
            delete color;
        }
    }

    cout<<"bingo"<<endl;
    image.save_image("F:\\MATERIALS\\4-1\\graphics\\draw\\out.bmp");
}


void keyboardListener(unsigned char key, int x,int y){
	switch(key){

		case '1':
			rotate_about_u(cameraAngle);
			break;

        case '2':
            rotate_about_u(-cameraAngle);
            break;

        case '3':
            rotate_about_r(cameraAngle);
            break;

        case '4':
            rotate_about_r(-cameraAngle);
            break;

        case '5':
            rotate_about_l(cameraAngle);
            break;

        case '6':
            rotate_about_l(-cameraAngle);
            break;

        case '0':
            capture();
            break;

		default:
			break;
	}
}


void specialKeyListener(int key, int x,int y){
	switch(key){
		case GLUT_KEY_DOWN:		//down arrow key
            move_backwards();
			break;
		case GLUT_KEY_UP:		// up arrow key
			move_forward();
			break;

		case GLUT_KEY_RIGHT:
			move_right();
			break;
		case GLUT_KEY_LEFT:
			move_left();
			break;

		case GLUT_KEY_PAGE_UP:
		    move_up();
			break;
		case GLUT_KEY_PAGE_DOWN:
		    move_down();
			break;

		case GLUT_KEY_INSERT:
			break;

		case GLUT_KEY_HOME:
			break;
		case GLUT_KEY_END:
			break;

		default:
			break;
	}
}


void mouseListener(int button, int state, int x, int y){	//x, y is the x-y of the screen (2D)
	switch(button){
		case GLUT_LEFT_BUTTON:
            if(state == GLUT_DOWN){
            }
            break;

		case GLUT_RIGHT_BUTTON:
			if(state == GLUT_DOWN){		// 2 times?? in ONE click? -- solution is checking DOWN or UP
				drawaxes=1-drawaxes;
			}
			break;

		case GLUT_MIDDLE_BUTTON:
			//........
			break;

		default:
			break;
	}
}



void display(){

	//clear the display
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClearColor(0,0,0,0);	//color black
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	/********************
	/ set-up camera here
	********************/
	//load the correct matrix -- MODEL-VIEW matrix
	glMatrixMode(GL_MODELVIEW);

	//initialize the matrix
	glLoadIdentity();

    gluLookAt(pos.x , pos.y , pos.z ,       pos.x+look.x , pos.y+look.y , pos.z+look.z ,      up.x , up.y , up.z);

	//again select MODEL-VIEW
	glMatrixMode(GL_MODELVIEW);


	/****************************
	/ Add your objects from here
	****************************/
	//add objects

	drawAxes();
	//drawGrid();
	//drawSquare(wall);
	objects[total_objects]->draw();
	//objects[3]->draw();
	//objects[4]->draw();
	//objects[5]->draw();
	//objects[0]->draw();
	for(int i=0 ; i<total_objects ; i++)
    {
        glPushMatrix();
        {
            glTranslatef(objects[i]->reference_point.x , objects[i]->reference_point.y , objects[i]->reference_point.z);
            objects[i]->draw();
        }
        glPopMatrix();
    }

	//drawPoint();
	for(int i=0 ; i<light_sources ; i++)
    {
        lights[i].draw();
    }

    /*for(int i=0 ; i<rays.size() ; i++)
    {
        glBegin(GL_LINES);{
			glVertex3f( rays[i].start.x , rays[i].start.y , rays[i].start.z);
			glVertex3f(rays[i].start.x + rays[i].dir.x*200 ,rays[i].start.y + rays[i].dir.y*200 , rays[i].start.z + rays[i].dir.z*200);
		}glEnd();
    }*/

	//ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
	glutSwapBuffers();
}


void animate(){
	angle+=0.05;
	//codes for any changes in Models, Camera
	glutPostRedisplay();
}

void init(){
	//codes for initialization
	drawgrid=0;
	drawaxes=1;
	cameraHeight=150.0;
	cameraAngle=0.02;
	angle=0;

	pos.x = 100;
	pos.y = 100;
	pos.z = 0;

	up.x = 0;
	up.y = 0;
	up.z = 1;

	c_right.x = -(1/sqrt(2));
	c_right.y = 1/sqrt(2);
	c_right.z = 0;

	look.x = -(1/sqrt(2));
	look.y = -(1/sqrt(2));
	look.z = 0;

	//clear the screen
	glClearColor(0,0,0,0);

	/************************
	/ set-up projection here
	************************/
	//load the PROJECTION matrix
	glMatrixMode(GL_PROJECTION);

	//initialize the matrix
	glLoadIdentity();

	//give PERSPECTIVE parameters
	gluPerspective(80,	1,	1,	1000.0);
	//field of view in the Y (vertically)
	//aspect ratio that determines the field of view in the X direction (horizontally)
	//near distance
	//far distance
}

int main(int argc, char **argv){
	glutInit(&argc,argv);
	glutInitWindowSize(500, 500);
	glutInitWindowPosition(0, 0);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	//Depth, Double buffer, RGB color

	glutCreateWindow("Ray");

	init();
	load_data();
	/*for(int i=0 ; i<total_objects ; i++)
    {
        objects[i]->print();
    }*/
    /*Ray ray(100 , 100 , 0 , 0-100 , 0-100 , 0);
    double* dummy;
    double inter = objects[0]->intersect(ray , dummy , 1);
    cout<<inter<<endl;*/

	glEnable(GL_DEPTH_TEST);	//enable Depth Testing

	glutDisplayFunc(display);	//display callback function
	glutIdleFunc(animate);		//what you want to do in the idle time (when no drawing is occuring)

	glutKeyboardFunc(keyboardListener);
	glutSpecialFunc(specialKeyListener);
	glutMouseFunc(mouseListener);

	glutMainLoop();		//The main loop of OpenGL

	return 0;
}
