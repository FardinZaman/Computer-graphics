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
#include <time.h>

#define pi (2*acos(0.0))

double cameraHeight;
double cameraAngle;
int drawgrid;
int drawaxes;
double angle;

double drawpoint = 0;

double camera_x_mover;
double camera_y_mover;
double camera_z_mover;

struct point
{
	double x,y,z;
};

struct point_vector
{
    double x,y,z;
};



double red_circle_radius = 35;

point positions[5];
point_vector velocity[5];
point positions_changed[5];
double center_to_center_catch[5];
bool inside_red_circle[5];
point_vector previous[5];
point_vector crossed_previous[5][5];

double bubbles_radius = 5;

int bubbles_permitted = 0;
int bubble_tracker = 0;

double bubble_speed = 0.005;
double speed_changer = 0.001;
bool paused = false;
double remember_the_speed = 0;
bool already_crossed[5][5];
bool already_crossed_the_circle[5];

double dotProduct(double vect_A[], double vect_B[])
{
    double product = 0;
    for (int i = 0; i < 3; i++)
        product = product + vect_A[i] * vect_B[i];
    return product;
}


point_vector reflection(point_vector vect , point_vector normalized_axis)
{
    double velocity_now[] = {vect.x , vect.y , vect.z};
    double axis[] = {normalized_axis.x , normalized_axis.y , normalized_axis.z};

    double velocity_dot_axis = dotProduct(velocity_now , axis);

    point_vector reflected;

    reflected.x = vect.x - 2*velocity_dot_axis*normalized_axis.x;
    reflected.y = vect.y - 2*velocity_dot_axis*normalized_axis.y;
    reflected.z = 0;

    return reflected;
}


void initialize_bubbles(double radius)
{
    srand(time(0));
    int i;
    for(i=0 ; i<5 ; i++)
    {
        positions[i].x = radius;
        positions[i].y = radius;
        positions[i].z = 0;

        inside_red_circle[i] = false;

        //double x = i + 1;
        //double y = 5 - i;
        double x = rand();
        double y = rand();
        double value = sqrt(x*x + y*y);

        velocity[i].x = x/value;
        velocity[i].y = y/value;
        velocity[i].z = 0;

        already_crossed_the_circle[i] = false;

        int j;
        for(j=0 ; j<5 ; j++)
        {
            already_crossed[i][j] = false;
        }

        positions_changed[i].x = radius;
        positions_changed[i].y = radius;
        positions_changed[i].z = 0;
    }

}


void drawGreenSquare(double a)
{
    glColor3f(0 , 1 , 0);

    glBegin(GL_LINES);
    glVertex2f(0 , 0);
    glVertex2f(a , 0);
    glEnd();

    glBegin(GL_LINES);
    glVertex2f(a , a);
    glVertex2f(a , 0);
    glEnd();

    glBegin(GL_LINES);
    glVertex2f(0 , a);
    glVertex2f(a , a);
    glEnd();

    glBegin(GL_LINES);
    glVertex2f(0 , 0);
    glVertex2f(0 , a);
    glEnd();
}


void drawCircle(double radius,int segments)
{
    int i;
    struct point points[100];
    glColor3f(0.8,0.1,0.4);
    //generate points
    for(i=0;i<=segments;i++)
    {
        points[i].x=radius*cos(((double)i/(double)segments)*2*pi);
        points[i].y=radius*sin(((double)i/(double)segments)*2*pi);
    }
    //draw segments using generated points
    for(i=0;i<segments;i++)
    {
        glBegin(GL_LINES);
        {
			glVertex3f(points[i].x,points[i].y,0);
			glVertex3f(points[i+1].x,points[i+1].y,0);
        }
        glEnd();
    }
}

void drawRedCircle(double radius , int segments)
{
    int i;
    struct point points[100];
    glColor3f(1,0,0);
    //generate points
    for(i=0;i<=segments;i++)
    {
        points[i].x=radius*cos(((double)i/(double)segments)*2*pi);
        points[i].y=radius*sin(((double)i/(double)segments)*2*pi);
    }
    //draw segments using generated points
    for(i=0;i<segments;i++)
    {
        glBegin(GL_LINES);
        {
			glVertex3f(points[i].x,points[i].y,0);
			glVertex3f(points[i+1].x,points[i+1].y,0);
        }
        glEnd();
        //glFlush();
    }
}

void draw_bubbles(int how_many)
{
    int i;
    for(i=0 ; i<how_many ; i++)
    {
        double x_part = positions[i].x;
        double y_part = positions[i].y;

        glPushMatrix();
        glTranslatef(x_part , y_part , 0);
        drawCircle(bubbles_radius , 30);
        glPopMatrix();
    }
}


void keyboardListener(unsigned char key, int x,int y){
	switch(key){
        case 'p':
            if(paused == false)
            {
                remember_the_speed = bubble_speed;
                bubble_speed = 0;
                paused = true;
            }
            else
            {
                bubble_speed = remember_the_speed;
                paused = false;
            }
		default:
			break;
	}
}


void specialKeyListener(int key, int x,int y){
	switch(key){
		case GLUT_KEY_DOWN:		//down arrow key
            if(bubble_speed > 0.002)
            {
                bubble_speed -= speed_changer;
            }
            else if(paused == true)
            {
                if(remember_the_speed > 0.002)
                {
                    remember_the_speed -= speed_changer;
                }
            }
			break;
		case GLUT_KEY_UP:		// up arrow key
			if(bubble_speed < 0.2)
            {
                if(paused == false){
                    bubble_speed += speed_changer;
                }
                else
                {
                    remember_the_speed += speed_changer;
                }
            }
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

	//now give three info
	//1. where is the camera (viewer)?
	//2. where is the camera looking?
	//3. Which direction is the camera's UP direction?

	//gluLookAt(100,100,100,	0,0,0,	0,0,1);
	//gluLookAt(200*cos(cameraAngle), 200*sin(cameraAngle), cameraHeight,		0,0,0,		0,0,1);
	//gluLookAt(0+camera_x_mover, 0+camera_y_mover , 100+camera_z_mover,	0 + camera_x_mover , 0+camera_y_mover , 0+camera_z_mover ,	0,1,0);
    //gluLookAt(pos.x , pos.y , pos.z ,       pos.x+look.x , pos.y+look.y , pos.z+look.z ,      up.x , up.y , up.z);
    gluLookAt(50 , 50 , 80 ,      50 , 50 , 0 ,       0 , 1 , 0);

	//again select MODEL-VIEW
	glMatrixMode(GL_MODELVIEW);


	/****************************
	/ Add your objects from here
	****************************/
	//add objects

    drawGreenSquare(100);

    glPushMatrix();
    glTranslatef(50 , 50 , 0);
    drawRedCircle(red_circle_radius , 50);
    glPopMatrix();

    bubbles_permitted = bubble_tracker / 1000;
    if(bubble_tracker < 5000)
    {
        bubble_tracker++;
    }

    draw_bubbles(bubbles_permitted);

    //drawRedCircle(5 , 50);


	//ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
	glutSwapBuffers();
}


void animate(){
	//angle+=0.05;
	//codes for any changes in Models, Camera
	/*struct point_vector catch_the_result;
	struct point_vector catch_the_bubble;
	catch_the_result.x = 0;
	catch_the_result.y = 0;
	catch_the_result.z = 0;
	catch_the_bubble.x = 0;
	catch_the_bubble.y = 0;
	catch_the_bubble.z = 0;
	*/

	int i;
	for(i=0 ; i<bubbles_permitted ; i++)
    {
        //printf("%f" , positions[i].x);
        positions_changed[i].x = positions[i].x;
        positions_changed[i].y = positions[i].y;
        positions_changed[i].z = positions[i].z;
    }

	for(i=0 ; i<bubbles_permitted ; i++)
    {
        struct point_vector catch_the_result;
        point_vector catch_the_bubble[bubbles_permitted];
        catch_the_result.x = 0;
        catch_the_result.y = 0;
        catch_the_result.z = 0;
        //catch_the_bubble.x = 0;
        //catch_the_bubble.y = 0;
        //catch_the_bubble.z = 0;
        int bb;
        for(bb=0 ; bb<bubbles_permitted ; bb++)
        {
            catch_the_bubble[bb].x = 0;
            catch_the_bubble[bb].y = 0;
            catch_the_bubble[bb].z = 0;
        }

        bool done = false;

        double radius_difference_with_circle = red_circle_radius - bubbles_radius;
        double x_gap = 50 - positions[i].x;
        double y_gap = 50 - positions[i].y;
        double center_to_center = sqrt(x_gap*x_gap + y_gap*y_gap);

        center_to_center_catch[i] = center_to_center;

        if(center_to_center >= radius_difference_with_circle && inside_red_circle[i] == true)
        {
            //struct point_vector axis;
            //axis.x = x_gap;
            //axis.y = y_gap;

            if(already_crossed_the_circle[i] == false){

                struct point_vector normalized_axis;
                normalized_axis.x = x_gap / center_to_center;
                normalized_axis.y = y_gap / center_to_center;
                normalized_axis.z = 0;

                //struct point_vector catch_the_result;
                catch_the_result = reflection(velocity[i] , normalized_axis);
                previous[i] = catch_the_result;

                /*velocity[i].x = catch_the_result.x;
                velocity[i].y = catch_the_result.y;
                velocity[i].z = catch_the_result.z;

                positions[i].x += velocity[i].x*bubble_speed;
                positions[i].y += velocity[i].y*bubble_speed;*/

                done = true;


                already_crossed_the_circle[i] = true;
            }

            else if(already_crossed_the_circle[i] == true)
            {
                catch_the_result = previous[i];
                done = true;
            }
        }

        if(center_to_center < radius_difference_with_circle)
        {
            if(inside_red_circle[i] == false)
            {
                int k;
                for(k=0 ; k<bubbles_permitted ; k++)
                {
                    if(k != i && inside_red_circle[k] == true)
                    {
                        double two_bubbles_radius = 2 * bubbles_radius;
                        double x_difference = positions[i].x - positions[k].x;
                        double y_difference = positions[i].y - positions[k].y;
                        double two_bubbles_center_distance = sqrt(x_difference*x_difference + y_difference*y_difference);

                        if(two_bubbles_center_distance < two_bubbles_radius)
                        {
                            already_crossed[i][k] = true;
                            already_crossed[k][i] = true;
                        }
                    }
                }

                inside_red_circle[i] = true;
            }

            else if(center_to_center < radius_difference_with_circle && already_crossed_the_circle[i] == true)
            {
                already_crossed_the_circle[i] = false;

                //previous[i].x = 0;
                //previous[i].y = 0;
            }
        }


        if(inside_red_circle[i] == true){

            int j;
            int tracker = 0;
            for(j=0 ; j<bubbles_permitted ; j++)
            {
                //tracker++;
                if(j != i && inside_red_circle[j] == true)
                {
                    //tracker++;
                    double two_bubbles_radius = 2 * bubbles_radius;
                    double x_difference = positions[i].x - positions[j].x;
                    double y_difference = positions[i].y - positions[j].y;
                    double two_bubbles_center_distance = sqrt(x_difference*x_difference + y_difference*y_difference);

                    if(two_bubbles_center_distance <= two_bubbles_radius)
                    {
                        //tracker++;
                        if(already_crossed[i][j] == false){
                            struct point_vector normalized_vect;
                            normalized_vect.x = x_difference / two_bubbles_center_distance;
                            normalized_vect.y = y_difference / two_bubbles_center_distance;
                            normalized_vect.z = 0;

                            //struct point_vector catch_the_bubble;
                            catch_the_bubble[j] = reflection(velocity[i] , normalized_vect);
                            //catch_the_bubble[j].x = velocity[j].x;
                            //catch_the_bubble[j].y = velocity[j].y;


                            crossed_previous[i][j] = catch_the_bubble[j];

                            //velocity[i].x = catch_the_bubble.x;
                            //velocity[i].y = catch_the_bubble.y;
                            //velocity[i].z = catch_the_bubble.z;

                            //positions[i].x += velocity[i].x*bubble_speed;
                            //positions[i].y += velocity[i].y*bubble_speed;

                            done = true;


                            already_crossed[i][j] = true;
                            //already_crossed[j][i] = true;
                        }

                        /*else
                        {
                            catch_the_bubble[j] = crossed_previous[i][j];
                            done = true;
                        }*/
                    }

                    if(two_bubbles_center_distance > two_bubbles_radius && already_crossed[i][j] == true)
                    {
                        already_crossed[i][j] = false;
                        //already_crossed[j][i] = false;
                    }
                }
            }

            //printf("%d" , tracker);

        }

        if(done == true)
        {
            velocity[i].x = catch_the_result.x;
            velocity[i].y = catch_the_result.y;
            int b;
            for(b=0 ; b<bubbles_permitted ; b++)
            {
                velocity[i].x += catch_the_bubble[b].x;
                velocity[i].y += catch_the_bubble[b].y;

            }
            double velocity_value = sqrt(velocity[i].x*velocity[i].x + velocity[i].y*velocity[i].y);
            velocity[i].x = velocity[i].x / velocity_value;
            velocity[i].y = velocity[i].y / velocity_value;
            velocity[i].z = 0;

            positions_changed[i].x += velocity[i].x*bubble_speed;
            positions_changed[i].y += velocity[i].y*bubble_speed;
        }


        if(inside_red_circle[i] == false || done == false){
            //printf("30");

            positions_changed[i].x += velocity[i].x*bubble_speed;
            positions_changed[i].y += velocity[i].y*bubble_speed;

            double boundary_x_right = positions_changed[i].x + bubbles_radius;
            double boundary_x_left = positions_changed[i].x - bubbles_radius;
            double boundary_y_right = positions_changed[i].y + bubbles_radius;
            double boundary_y_left = positions_changed[i].y - bubbles_radius;

            if(boundary_x_right >= 100 || boundary_x_left <= 0)
            {
                velocity[i].x = -velocity[i].x;
            }

            if(boundary_y_right >= 100 || boundary_y_left <= 0)
            {
                velocity[i].y = -velocity[i].y;
            }
        }
    }

    //printf("go");

    for(i=0 ; i<bubbles_permitted ; i++)
    {
        double radius_difference_with_circle = red_circle_radius - bubbles_radius;
        double x_gap = 50 - positions_changed[i].x;
        double y_gap = 50 - positions_changed[i].y;
        double center_to_center = sqrt(x_gap*x_gap + y_gap*y_gap);

        if(center_to_center > radius_difference_with_circle && already_crossed_the_circle[i] == true && center_to_center > center_to_center_catch[i])
        {
            //already_crossed_the_circle[i] = false;
        }
        else
        {
            positions[i].x = positions_changed[i].x;
            positions[i].y = positions_changed[i].y;
            positions[i].z = 0;
        }
        //printf("go");
        //positions[i].x = positions_changed[i].x;
        //positions[i].y = positions_changed[i].y;
        //positions[i].z = 0;
    }

    //printf("%f" , bubble_speed);

	glutPostRedisplay();
}

void init(){
	//codes for initialization
	drawgrid=0;
	drawaxes=1;
	cameraHeight=150.0;
	cameraAngle=0.02;
	angle=0;

	initialize_bubbles(bubbles_radius);

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

	glutCreateWindow("Bubble");

	init();

	glEnable(GL_DEPTH_TEST);	//enable Depth Testing

	glutDisplayFunc(display);	//display callback function
	glutIdleFunc(animate);		//what you want to do in the idle time (when no drawing is occuring)

	glutKeyboardFunc(keyboardListener);
	glutSpecialFunc(specialKeyListener);
	glutMouseFunc(mouseListener);

	glutMainLoop();		//The main loop of OpenGL

	return 0;
}
