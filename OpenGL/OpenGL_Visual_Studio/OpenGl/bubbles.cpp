#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include <windows.h>
#include <glut.h>

#define pi (2*acos(0.0))



struct point
{
	double x, y;
};

class vector
{
public:
	double x, y;
};



int drawaxes;
double speed;
point p1, p2;
vector v1, v2;
double r_circle;

double r_bubble;

void drawAxes()
{
	if (drawaxes == 1)
	{
		glColor3f(1.0, 1.0, 1.0);
		glBegin(GL_LINES); {
			glVertex3f(200, 0, 0);
			glVertex3f(-200, 0, 0);

			glVertex3f(0, -200, 0);
			glVertex3f(0, 200, 0);

			glVertex3f(0, 0, 200);
			glVertex3f(0, 0, -200);
		}glEnd();
	}
}



void drawbubble(point p){
	int i;
	int segments = 25;
	struct point points[100];
	//generate points
	for (i = 0; i <= segments; i++)
	{
		points[i].x = r_bubble*cos(((double)i / (double)segments) * 2 * pi) + p.x;
		points[i].y = r_bubble*sin(((double)i / (double)segments) * 2 * pi) + p.y;
	}
	//draw segments using generated points
	for (i = 0; i<segments; i++)
	{
		glBegin(GL_LINES);
		{
			glVertex3f(points[i].x, points[i].y, 0);
			glVertex3f(points[i + 1].x, points[i + 1].y, 0);
		}
		glEnd();
	}
}

void drawCircle(double radius, int segments)
{
	int i;
	struct point points[100];
	glColor3f(0.7, 0.7, 0.7);
	//generate points
	for (i = 0; i <= segments; i++)
	{
		points[i].x = radius*cos(((double)i / (double)segments) * 2 * pi);
		points[i].y = radius*sin(((double)i / (double)segments) * 2 * pi);
	}
	//draw segments using generated points
	for (i = 0; i<segments; i++)
	{
		glBegin(GL_LINES);
		{
			glVertex3f(points[i].x, points[i].y, 0);
			glVertex3f(points[i + 1].x, points[i + 1].y, 0);
		}
		glEnd();
	}
}

void drawArrow(){
	double p1_x = p1.x + v1.x*(r_bubble - 5);
	double p1_y = p1.y + v1.y*(r_bubble - 5);
	
	double a1_x = p1_x + v1.x * 5;
	double a1_y = p1_y + v1.y * 5;
	double b1_x = p1_x - v1.y * 2;
	double b1_y = p1_y + v1.x * 2;
	double c1_x = p1_x + v1.y * 2;
	double c1_y = p1_y - v1.x * 2;

	
	double p2_x = p2.x + v2.x*(r_bubble - 5);
	double p2_y = p2.y + v2.y*(r_bubble - 5);
	
	double a2_x = p2_x + v2.x * 5;
	double a2_y = p2_y + v2.y * 5;
	double b2_x = p2_x - v2.y * 2;
	double b2_y = p2_y + v2.x * 2;
	double c2_x = p2_x + v2.y * 2;
	double c2_y = p2_y - v2.x * 2;
	
	glBegin(GL_LINES);
	{
		glColor3f(1, 1, 0);
		glVertex3f(p1.x, p1.y, 0);
		glVertex3f(p1_x, p1_y, 0);
		glColor3f(0, 1, 0);
		glVertex3f(p2.x, p2.y, 0);
		glVertex3f(p2_x, p2_y, 0);
	}
	glEnd();
	glBegin(GL_TRIANGLES);
	{
		glColor3f(1, 0, 0);
		glVertex3f(a1_x, a1_y, 0);
		glVertex3f(b1_x, b1_y, 0);
		glVertex3f(c1_x, c1_y, 0);
		glColor3f(1, 0, 0);
		glVertex3f(a2_x, a2_y, 0);
		glVertex3f(b2_x, b2_y, 0);
		glVertex3f(c2_x, c2_y, 0);
	}
	glEnd();
}

double value_(double x, double y, double z)
{
	return sqrt(x*x + y*y + z*z);
}

void keyboardListener(unsigned char key, int x, int y){
	switch (key){
	
	default:
		break;
	}
}

void specialKeyListener(int key, int x, int y){
	switch (key){
		vector v_t;
	case GLUT_KEY_DOWN:		//down arrow key
		break;
	case GLUT_KEY_UP:		// up arrow key

		break;

	case GLUT_KEY_RIGHT:
		v_t.x = v1.x*cos(3 * (pi / 180)) + v1.y*sin(3 * (pi / 180));
		v_t.y = v1.y*cos(3 * (pi / 180)) - v1.x*sin(3 * (pi / 180));
		v1.x = v_t.x / sqrt(v_t.x*v_t.x + v_t.y*v_t.y);
		v1.y = v_t.y / sqrt(v_t.x*v_t.x + v_t.y*v_t.y);
		break;
	case GLUT_KEY_LEFT:
		v_t.x = v1.x*cos(3 * (pi / 180)) - v1.y*sin(3 * (pi / 180));
		v_t.y = v1.y*cos(3 * (pi / 180)) + v1.x*sin(3 * (pi / 180));
		v1.x = v_t.x / sqrt(v_t.x*v_t.x + v_t.y*v_t.y);
		v1.y = v_t.y / sqrt(v_t.x*v_t.x + v_t.y*v_t.y);
		break;

	case GLUT_KEY_PAGE_UP:
		break;
	case GLUT_KEY_PAGE_DOWN:

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
	switch (button){
	case GLUT_LEFT_BUTTON:
		if (state == GLUT_DOWN){		// 2 times?? in ONE click? -- solution is checking DOWN or UP
			drawaxes = 1 - drawaxes;
		}
		break;

	case GLUT_RIGHT_BUTTON:
		//........
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
	glClearColor(0, 0, 0, 0);	//color black
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
	gluLookAt(0,0,200,	0,0,0,	0,1,0);
	//gluLookAt(pos.x, pos.y, pos.z, pos.x + l.x, pos.y + l.y, pos.z + l.z, u.x, u.y, u.z);

	//again select MODEL-VIEW
	glMatrixMode(GL_MODELVIEW);


	/****************************
	/ Add your objects from here
	****************************/
	//add objects

	drawAxes();
	
	drawCircle(r_circle,50);
	glColor3f(1, 1, 0);
	drawbubble(p1);
	glColor3f(0, 1, 0);
	drawbubble(p2);
	drawArrow();
	
	//ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
	glutSwapBuffers();
}

double distance(point p){
	return sqrt(p.x*p.x + p.y*p.y);
}

void animate(){
	double dist = sqrt((p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y));
	if (dist <= 2 * r_bubble){
		//p1
		double p_x = (p2.x - p1.x) / sqrt((p2.x - p1.x)*(p2.x - p1.x) + (p2.y - p1.y)*(p2.y - p1.y));
		double p_y = (p2.y - p1.y) / sqrt((p2.x - p1.x)*(p2.x - p1.x) + (p2.y - p1.y)*(p2.y - p1.y));
		double a_n = p_x*v1.x + p_y*v1.y;
		v1.x = v1.x - 2 * a_n*p_x;
		v1.y = v1.y - 2 * a_n*p_y;
		p_x = -p_x;
		p_y = -p_y;
		a_n = p_x*v2.x + p_y*v2.y;
		v2.x = v2.x - 2 * a_n*p_x;
		v2.y = v2.y - 2 * a_n*p_y;
	}
	if (distance(p1)+r_bubble>=r_circle){
		
		double p_x = -p1.x / sqrt(p1.x*p1.x + p1.y*p1.y);
		double p_y = -p1.y / sqrt(p1.x*p1.x + p1.y*p1.y);
		double a_n = p_x*v1.x + p_y*v1.y;
		v1.x = v1.x - 2 * a_n*p_x;
		v1.y = v1.y - 2 * a_n*p_y;
	}
	if (distance(p2) + r_bubble >= r_circle){
		double p_x = -p2.x / sqrt(p2.x*p2.x + p2.y*p2.y);
		double p_y = -p2.y / sqrt(p2.x*p2.x + p2.y*p2.y);
		double a_n = p_x*v2.x + p_y*v2.y;
		v2.x = v2.x - 2 * a_n*p_x;
		v2.y = v2.y - 2 * a_n*p_y;
	}
	p1.x += v1.x*speed;
	p1.y += v1.y*speed;
	p2.x += v2.x*speed;
	p2.y += v2.y*speed;

	glutPostRedisplay();
}

void init(){
	//codes for initialization
	drawaxes = 1;
	r_circle = 150;
	r_bubble = 15;
	speed = 0.03;

	p1.x = 25;
	p1.y = 25;
	p2.x = -25;
	p2.y = 25;


	v1.x = 1 / sqrt(2.0);
	v1.y = -1 / sqrt(2.0);
	
	v2.x = -1 / sqrt(2.0);
	v2.y = -1 / sqrt(2.0);

	//clear the screen
	glClearColor(0, 0, 0, 0);

	/************************
	/ set-up projection here
	************************/
	//load the PROJECTION matrix
	glMatrixMode(GL_PROJECTION);

	//initialize the matrix
	glLoadIdentity();

	//give PERSPECTIVE parameters
	gluPerspective(80, 1, 1, 1000.0);
	//field of view in the Y (vertically)
	//aspect ratio that determines the field of view in the X direction (horizontally)
	//near distance
	//far distance
}

int main(int argc, char **argv){
	glutInit(&argc, argv);
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
