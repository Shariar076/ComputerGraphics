#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include <iostream>
#include <windows.h>
#include <glut.h>

#define pi (2*acos(0.0))


class point
{
public:
	double x, y, z;
};

class vector
{
public:
	double x, y, z;
};


class position
{
public:
	double x, y, z;
};

class up
{
public:
	double x, y, z;
};

class right
{
public:
	double x, y, z;
};

class look
{
public:
	double x, y, z;
};

double cameraAngle;

int drawaxes;
double rot_angle;
double base_angle_X;
double base_angle_Y;
double barrel_angle;
double base_r;
double base_l;
double barrel_r;
double barrel_l;
double target_a;
double target_d;
int num_bullets;
bool drawblood = false;
bool key_q = false;
bool key_w = false;
bool key_e = false;
bool key_r = false;
bool key_a = false;
bool key_s = false;
bool key_d = false;
bool key_f = false;

point bullet[100];
position pos;
up u;
right r;
look l;
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

void drawSquare(double a)
{
    //glColor3f(1.0,0.0,0.0);
	glBegin(GL_QUADS);{
		glVertex3f(a, a, 0);
		glVertex3f(a, -a, 0);
		glVertex3f(-a, -a, 0);
		glVertex3f(-a, a, 0);
	}glEnd();
}

void drawBlood(){
	
	if (drawblood){
//		std::cout << base_angle_X << " " << base_angle_Y << std::endl;
		double base_x = -(base_l + 2 * base_r)*sin(base_angle_Y*(pi/180));
		double base_y = (base_l + 2 * base_r)*cos(base_angle_Y*(pi / 180))*sin(base_angle_X*(pi / 180));
		double base_z = -(base_l + 2 * base_r)*cos(base_angle_Y*(pi / 180))*cos(base_angle_X*(pi / 180));

		double vec_x = -(barrel_l + 2 * barrel_r)*sin((base_angle_Y+barrel_angle)*(pi / 180));
		double vec_y = (barrel_l + 2 * barrel_r)*cos((base_angle_Y + barrel_angle)*(pi / 180))*sin(base_angle_X*(pi / 180));
		double vec_z = - (barrel_l + 2 * barrel_r)*cos((base_angle_Y + barrel_angle)*(pi / 180))*cos(base_angle_X*(pi / 180));
		
		double t = ((target_d+2) - base_z) / vec_z;
		double x = base_x + t*vec_x;
		double y = base_y + t*vec_y;
		double z = base_z + t*vec_z;
		if (x<target_a && x>-target_a && y<target_a && y>-target_a){
			bullet[num_bullets].x = x;
			bullet[num_bullets].y = y;
			bullet[num_bullets].z = z;
			num_bullets++;
		}
		drawblood = false;
	}
	for (int i = 0; i < num_bullets; i++){
		glPushMatrix(); {
			glTranslatef(bullet[i].x, bullet[i].y, bullet[i].z);
			glColor3f(1, 0, 0);
			drawSquare(5);
		}glPopMatrix();
	}
}

void drawCircle(double radius,int segments)
{
    int i;
    point points[100];
    glColor3f(0.7,0.7,0.7);
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

void drawCylinder(double radius, double height,int slices){
	point points[2][100];
	

	for (int i = 0; i < 2; i++)
	{
		int h = height*i;
		for (int j = 0; j <= slices; j++)
		{
			points[i][j].x = radius*cos(((double)j / (double)slices) * 2 * pi);
			points[i][j].y = radius*sin(((double)j / (double)slices) * 2 * pi);
			points[i][j].z = h;
		}
	}

	
	//glColor3f((double)i/(double)stacks,(double)i/(double)stacks,(double)i/(double)stacks);
	double shade = 1;
	for (int j = 0; j<slices; j++)
	{
		shade = 1 - shade;
		glColor3f(shade, shade, shade);
		glBegin(GL_QUADS); {
			//upper hemisphere
			glVertex3f(points[0][j].x, points[0][j].y, points[0][j].z);
			glVertex3f(points[0][j + 1].x, points[0][j + 1].y, points[0][j + 1].z);
			glVertex3f(points[1][j + 1].x, points[1][j + 1].y, points[1][j + 1].z);
			glVertex3f(points[1][j].x, points[1][j].y, points[1][j].z);
		}glEnd();
	}

}

void drawOpenHemisphere(double radius,int slices,int stacks)
{
	point points[100][100];
	int i,j;
	double h,r;
	//generate points
	for(i=0;i<=stacks;i++)
	{
		h = radius*sin(((double)i / (double)stacks)*(pi / 2));
		r = radius*cos(((double)i / (double)stacks)*(pi / 2));
		for(j=0;j<=slices;j++)
		{
			points[i][j].x = r*cos(((double)j / (double)slices) * 2 * pi);
			points[i][j].y = r*sin(((double)j / (double)slices) * 2 * pi);
			points[i][j].z = h;
		}
	}
	//draw quads using generated points
	double shade = 1;
	for(i=0;i<stacks;i++)
	{
        //glColor3f((double)i/(double)stacks,(double)i/(double)stacks,(double)i/(double)stacks);
		for(j=0;j<slices;j++)
		{
			shade = 1 - shade;
			glColor3f(shade, shade, shade);
			glBegin(GL_QUADS);{
			    //upper hemisphere
				glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
				glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);
			}glEnd();
		}
	}
}

void drawCloseHemisphere(double radius, int slices, int stacks)
{
	point points[100][100];
	int i, j;
	double h, r;
	//generate points
	for (i = 0; i <= stacks; i++)
	{
		h = radius*sin(((double)i / (double)stacks)*(pi / 2));
		r = radius*cos(((double)i / (double)stacks)*(pi / 2));
		for (j = 0; j <= slices; j++)
		{
			points[i][j].x = r*cos(((double)j / (double)slices) * 2 * pi);
			points[i][j].y = r*sin(((double)j / (double)slices) * 2 * pi);
			points[i][j].z = h;
		}
	}
	//draw quads using generated points
	double shade = 1;
	for (i = 0; i<stacks; i++)
	{
		//glColor3f((double)i/(double)stacks,(double)i/(double)stacks,(double)i/(double)stacks);
		for (j = 0; j<slices; j++)
		{
			shade = 1 - shade;
			glColor3f(shade, shade, shade);
			glBegin(GL_QUADS); {

				glVertex3f(points[i][j].x, points[i][j].y, -points[i][j].z);
				glVertex3f(points[i][j + 1].x, points[i][j + 1].y, -points[i][j + 1].z);
				glVertex3f(points[i + 1][j + 1].x, points[i + 1][j + 1].y, -points[i + 1][j + 1].z);
				glVertex3f(points[i + 1][j].x, points[i + 1][j].y, -points[i + 1][j].z);
			}glEnd();
		}
	}
}

void drawTipHemisphere(double radius, int slices, int stacks)
{
	point points[100][100];
	int i, j;
	double h, r;
	//generate points
	for (i = 0; i <= stacks; i++)
	{
		h = radius*sin(((double)i / (double)stacks)*(pi / 2));
		if (i == 0){
			r = radius*cos(((double)i / (double)stacks)*(pi / 2));
		}
		else{
			r += r-radius*cos(((double)i / (double)stacks)*(pi / 2));
		}
		for (j = 0; j <= slices; j++)
		{
			points[i][j].x = r*cos(((double)j / (double)slices) * 2 * pi);
			points[i][j].y = r*sin(((double)j / (double)slices) * 2 * pi);
			points[i][j].z = h;
		}
	}
	//draw quads using generated points
	double shade = 1;
	for (i = 0; i<stacks; i++)
	{
		//glColor3f((double)i/(double)stacks,(double)i/(double)stacks,(double)i/(double)stacks);
		for (j = 0; j<slices; j++)
		{
			shade = 1 - shade;
			glColor3f(shade, shade, shade);
			glBegin(GL_QUADS); {
				
				glVertex3f(points[i][j].x,points[i][j].y,-points[i][j].z);
				glVertex3f(points[i][j+1].x,points[i][j+1].y,-points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,-points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x,points[i+1][j].y,-points[i+1][j].z);
			}glEnd();
		}
	}
}

double value_(double x, double y, double z)
{
	return sqrt(x*x + y*y + z*z);
}

void keyboardListener(unsigned char key, int x,int y){
	switch(key){
		look l_temp, r_temp, u_temp;
	case '1':
		
		l_temp.x = l.x*cos(cameraAngle*(pi / 180)) - r.x*sin(cameraAngle*(pi / 180));
		l_temp.y = l.y*cos(cameraAngle*(pi / 180)) - r.y*sin(cameraAngle*(pi / 180));
		l_temp.z = l.z*cos(cameraAngle*(pi / 180)) - r.z*sin(cameraAngle*(pi / 180));
		
		
		r_temp.x = r.x*cos(cameraAngle*(pi / 180)) + l.x*sin(cameraAngle*(pi / 180));
		r_temp.y = r.y*cos(cameraAngle*(pi / 180)) + l.y*sin(cameraAngle*(pi / 180));
		r_temp.z = r.z*cos(cameraAngle*(pi / 180)) + l.z*sin(cameraAngle*(pi / 180));
		
		l.x = l_temp.x / value_(l_temp.x, l_temp.y, l_temp.z);
		l.y = l_temp.y / value_(l_temp.x, l_temp.y, l_temp.z);
		l.z = l_temp.z / value_(l_temp.x, l_temp.y, l_temp.z);

		r.x = r_temp.x / value_(r_temp.x, r_temp.y, r_temp.z);
		r.y = r_temp.y / value_(r_temp.x, r_temp.y, r_temp.z);
		r.z = r_temp.z / value_(r_temp.x, r_temp.y, r_temp.z);

		break;
	case '2':
		l_temp.x = l.x*cos(cameraAngle*(pi / 180)) + r.x*sin(cameraAngle*(pi / 180));
		l_temp.y = l.y*cos(cameraAngle*(pi / 180)) + r.y*sin(cameraAngle*(pi / 180));
		l_temp.z = l.z*cos(cameraAngle*(pi / 180)) + r.z*sin(cameraAngle*(pi / 180));
		
		r_temp.x = r.x*cos(cameraAngle*(pi / 180)) - l.x*sin(cameraAngle*(pi / 180));
		r_temp.y = r.y*cos(cameraAngle*(pi / 180)) - l.y*sin(cameraAngle*(pi / 180));
		r_temp.z = r.z*cos(cameraAngle*(pi / 180)) - l.z*sin(cameraAngle*(pi / 180));
		
		l.x = l_temp.x / value_(l_temp.x, l_temp.y, l_temp.z);
		l.y = l_temp.y / value_(l_temp.x, l_temp.y, l_temp.z);
		l.z = l_temp.z / value_(l_temp.x, l_temp.y, l_temp.z);

		r.x = r_temp.x / value_(r_temp.x, r_temp.y, r_temp.z);
		r.y = r_temp.y / value_(r_temp.x, r_temp.y, r_temp.z);
		r.z = r_temp.z / value_(r_temp.x, r_temp.y, r_temp.z);
		
		break;
	case '3':
		l_temp.x = l.x*cos(cameraAngle*(pi / 180)) + u.x*sin(cameraAngle*(pi / 180));
		l_temp.y = l.y*cos(cameraAngle*(pi / 180)) + u.y*sin(cameraAngle*(pi / 180));
		l_temp.z = l.z*cos(cameraAngle*(pi / 180)) + u.z*sin(cameraAngle*(pi / 180));

		u_temp.x = u.x*cos(cameraAngle*(pi / 180)) - l.x*sin(cameraAngle*(pi / 180));
		u_temp.y = u.y*cos(cameraAngle*(pi / 180)) - l.y*sin(cameraAngle*(pi / 180));
		u_temp.z = u.z*cos(cameraAngle*(pi / 180)) - l.z*sin(cameraAngle*(pi / 180));

		l.x = l_temp.x / value_(l_temp.x, l_temp.y, l_temp.z);
		l.y = l_temp.y / value_(l_temp.x, l_temp.y, l_temp.z);
		l.z = l_temp.z / value_(l_temp.x, l_temp.y, l_temp.z);

		u.x = u_temp.x / value_(u_temp.x, u_temp.y, u_temp.z);
		u.y = u_temp.y / value_(u_temp.x, u_temp.y, u_temp.z);
		u.z = u_temp.z / value_(u_temp.x, u_temp.y, u_temp.z);
		
		break;
	case '4':
		l_temp.x = l.x*cos(cameraAngle*(pi / 180)) - u.x*sin(cameraAngle*(pi / 180));
		l_temp.y = l.y*cos(cameraAngle*(pi / 180)) - u.y*sin(cameraAngle*(pi / 180));
		l_temp.z = l.z*cos(cameraAngle*(pi / 180)) - u.z*sin(cameraAngle*(pi / 180));

		u_temp.x = u.x*cos(cameraAngle*(pi / 180)) + l.x*sin(cameraAngle*(pi / 180));
		u_temp.y = u.y*cos(cameraAngle*(pi / 180)) + l.y*sin(cameraAngle*(pi / 180));
		u_temp.z = u.z*cos(cameraAngle*(pi / 180)) + l.z*sin(cameraAngle*(pi / 180));
		
		l.x = l_temp.x / value_(l_temp.x, l_temp.y, l_temp.z);
		l.y = l_temp.y / value_(l_temp.x, l_temp.y, l_temp.z);
		l.z = l_temp.z / value_(l_temp.x, l_temp.y, l_temp.z);

		u.x = u_temp.x / value_(u_temp.x, u_temp.y, u_temp.z);
		u.y = u_temp.y / value_(u_temp.x, u_temp.y, u_temp.z);
		u.z = u_temp.z / value_(u_temp.x, u_temp.y, u_temp.z);

		break;
	case '5':
		u_temp.x = u.x*cos(cameraAngle*(pi / 180)) + r.x*sin(cameraAngle*(pi / 180));
		u_temp.y = u.y*cos(cameraAngle*(pi / 180)) + r.y*sin(cameraAngle*(pi / 180));
		u_temp.z = u.z*cos(cameraAngle*(pi / 180)) + r.z*sin(cameraAngle*(pi / 180));
		
		r_temp.x = r.x*cos(cameraAngle*(pi / 180)) - u.x*sin(cameraAngle*(pi / 180));
		r_temp.y = r.y*cos(cameraAngle*(pi / 180)) - u.y*sin(cameraAngle*(pi / 180));
		r_temp.z = r.z*cos(cameraAngle*(pi / 180)) - u.z*sin(cameraAngle*(pi / 180));

		r.x = r_temp.x / value_(r_temp.x, r_temp.y, r_temp.z);
		r.y = r_temp.y / value_(r_temp.x, r_temp.y, r_temp.z);
		r.z = r_temp.z / value_(r_temp.x, r_temp.y, r_temp.z);
		u.x = u_temp.x / value_(u_temp.x, u_temp.y, u_temp.z);
		u.y = u_temp.y / value_(u_temp.x, u_temp.y, u_temp.z);
		u.z = u_temp.z / value_(u_temp.x, u_temp.y, u_temp.z);

		break;
	case '6':
		u_temp.x = u.x*cos(cameraAngle*(pi / 180)) - r.x*sin(cameraAngle*(pi / 180));
		u_temp.y = u.y*cos(cameraAngle*(pi / 180)) - r.y*sin(cameraAngle*(pi / 180));
		u_temp.z = u.z*cos(cameraAngle*(pi / 180)) - r.z*sin(cameraAngle*(pi / 180));

		r_temp.x = r.x*cos(cameraAngle*(pi / 180)) + u.x*sin(cameraAngle*(pi / 180));
		r_temp.y = r.y*cos(cameraAngle*(pi / 180)) + u.y*sin(cameraAngle*(pi / 180));
		r_temp.z = r.z*cos(cameraAngle*(pi / 180)) + u.z*sin(cameraAngle*(pi / 180));

		u.x = u_temp.x / value_(u_temp.x, u_temp.y, u_temp.z);
		u.y = u_temp.y / value_(u_temp.x, u_temp.y, u_temp.z);
		u.z = u_temp.z / value_(u_temp.x, u_temp.y, u_temp.z);

		r.x = r_temp.x / value_(r_temp.x, r_temp.y, r_temp.z);
		r.y = r_temp.y / value_(r_temp.x, r_temp.y, r_temp.z);
		r.z = r_temp.z / value_(r_temp.x, r_temp.y, r_temp.z);

		break;
	case 'q':
		key_q = true;
		break;
	case 'w':
		key_w = true;
		break;
	case 'e':
		key_e = true;
		break;
	case 'r':
		key_r = true;
		break;
	case 'a':
		key_a = true;
		break;
	case 's':
		key_s = true;
		break;
	case 'd':
		key_d = true;
		break;
	case 'f':
		key_f = true;
		break;
	default:
		break;
	}
}

void specialKeyListener(int key, int x,int y){
	switch(key){
		case GLUT_KEY_DOWN:		//down arrow key
			pos.x -= l.x;
			pos.y -= l.y;
			pos.z -= l.z;
			break;
		case GLUT_KEY_UP:		// up arrow key
			pos.x += l.x;
			pos.y += l.y;
			pos.z += l.z;
			break;

		case GLUT_KEY_RIGHT:
			pos.x += r.x;
			pos.y += r.y;
			pos.z += r.z;
			break;
		case GLUT_KEY_LEFT:
			pos.x -= r.x;
			pos.y -= r.y;
			pos.z -= r.z;
			break;

		case GLUT_KEY_PAGE_UP:
			pos.x += u.x;
			pos.y += u.y;
			pos.z += u.z;
			break;
		case GLUT_KEY_PAGE_DOWN:
			pos.x -= u.x;
			pos.y -= u.y;
			pos.z -= u.z;
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
			if (state == GLUT_DOWN){
				drawblood =true;
			}

			break;

		case GLUT_RIGHT_BUTTON:
			if (state == GLUT_DOWN){		// 2 times?? in ONE click? -- solution is checking DOWN or UP
				drawaxes = 1 - drawaxes;
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

	//now give three info
	//1. where is the camera (viewer)?
	//2. where is the camera looking?
	//3. which direction is the camera's up direction?

	//glulookat(100,100,100,	0,0,0,	0,0,1);
	//gluLookAt(200*cos(cameraAngle), 200*sin(cameraAngle), cameraHeight,		0,0,0,		0,0,1);
	//gluLookAt(0,0,200,	0,0,0,	0,1,0);
	gluLookAt(pos.x, pos.y, pos.z, pos.x + l.x, pos.y + l.y, pos.z + l.z, u.x, u.y, u.z);

	//again select MODEL-VIEW
	glMatrixMode(GL_MODELVIEW);


	/****************************
	/ Add your objects from here
	****************************/
	//add objects

	drawAxes();

    //glColor3f(1,0,0);
	

    //drawCircle(30,24);

    //drawCone(20,50,24);
	
	glPushMatrix(); {
		glRotatef(base_angle_X, 1, 0, 0);
		glRotatef(base_angle_Y, 0, 1, 0);
		glTranslatef(0, 0, -base_r);
		drawOpenHemisphere(base_r, 26, 35);
		glTranslatef(0, 0, -base_l);
		drawCylinder(base_r, base_l, 26);
		drawCloseHemisphere(base_r, 26, 35);
		glTranslatef(0, 0, -base_r);
		glRotatef(barrel_angle, 0, 1, 0);
		glTranslatef(0, 0, -barrel_r);
		glRotatef(rot_angle, 0, 0, 1);
		drawOpenHemisphere(barrel_r, 26, 35);
		glTranslatef(0, 0, -barrel_l);
		drawCylinder(barrel_r, barrel_l, 26);
		drawTipHemisphere(barrel_r, 26, 1);
	}glPopMatrix();

	glPushMatrix(); {
		glColor3f(0.5, 0.5, 0.5);
		glTranslatef(0, 0, target_d);
		drawSquare(target_a);
	}glPopMatrix();
	glPushMatrix(); {
		drawBlood();
	}glPopMatrix();


	//ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
	glutSwapBuffers();
}

void animate(){
	if (key_q){
		if (base_angle_X<90)base_angle_X += 1;
		key_q = false;
	}
	if (key_w){
		if (base_angle_X>-25)base_angle_X -= 1;
		key_w = false;
	}
	if (key_e){
		if (base_angle_Y>-25)base_angle_Y -= 1;
		key_e = false;
	}
	if (key_r){
		if (base_angle_Y<25)base_angle_Y += 1;
		key_r = false;
	}
	if (key_a){
		if (barrel_angle < 25)barrel_angle += 1;
		key_a = false;
	}
	if (key_s){
		if (barrel_angle>-25)barrel_angle -= 1;
		key_s = false;
	}
	if (key_d){
		if (rot_angle<90)rot_angle += 1;
		key_d = false;
	}
	if (key_f){
		if (rot_angle>-90)rot_angle -= 1;
		key_f = false;
	}
	//codes for any changes in Models, Camera
	glutPostRedisplay();
}

void init(){
	//codes for initialization

	drawaxes=1;
	cameraAngle=1;
	rot_angle = 0;
	base_angle_X = 0;
	base_angle_Y = 0;
	barrel_angle = 0;
	num_bullets = 0;
	target_a = 150;
	target_d = -400;
	pos.x = 25;
	pos.y = 45;
	pos.z = 80;

	u.x = 1;
	u.y = 0;
	u.z = 0;

	r.x = 0;
	r.y = -1 / sqrt(2.0);
	r.z = 1 / sqrt(2.0);

	l.x = 0;
	l.y = -1 / sqrt(2.0);
	l.z = -1 / sqrt(2.0);

	base_l = 30;
	base_r = 18;
	barrel_l = 150;
	barrel_r = 12;
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
	gluPerspective(120,	1,	1,	1000.0);
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

	glutCreateWindow("My OpenGL Program");

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
