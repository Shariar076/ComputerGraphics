#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<iostream>

#include <windows.h>
#include <glut.h>

#define pi (2*acos(0.0))

class point2d
{
public:
	double x, y;
};


point2d cp[200];
point2d crvPoint[10000];


int cpidx;
int state_;
bool drawCrv;
bool drawgm;
bool follow;
int targetidx;
int crvidx;
double t;

void drawSquare()
{
	glBegin(GL_QUADS);
	{
		glVertex3d(3, 3, 0);
		glVertex3d(3, -3, 0);
		glVertex3d(-3, -3, 0);
		glVertex3d(-3, 3, 0);
	}
	glEnd();
}
void drawCircle(double radius, int segments)
{
	int i;
	point2d points[100];
	glColor3f(0, 1, 1);
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
			glBegin(GL_TRIANGLES);
			{
				glVertex3f(0, 0,0);
				glVertex3f(points[i].x, points[i].y, 0);
				glVertex3f(points[i + 1].x, points[i + 1].y, 0);

			}
		}
		glEnd();
	}
}
void drawGeometry(){
	for (int i = 0; i < cpidx; i++)
	{
		
		//glColor3f(1, 1, 0);
		glPushMatrix();
		{
			if (i % 2)glColor3f(1, 1, 0);
			else glColor3f(0, 1, 0);
			glTranslatef(cp[i].x, cp[i].y, 0);
			drawSquare();
		}glPopMatrix();
		if (i % 2){
			glColor3f(1, 1, 1);
			double x = cp[i].x - cp[i - 1].x;
			double y = cp[i].y - cp[i - 1].y;
			double len = sqrt(x*x + y*y);
			double v_x = x / len;
			double v_y = y / len;
			double l_x = cp[i - 1].x + v_x*0.9*len;
			double l_y = cp[i - 1].y + v_y*0.9*len;
			double b_x = l_x - v_y * 5;
			double b_y = l_y + v_x * 5;

			double c_x = l_x + v_y * 5;
			double c_y = l_y - v_x * 5;
			glBegin(GL_LINES); {
				glVertex3f(cp[i - 1].x, cp[i - 1].y, 0);
				glVertex3f(l_x, l_y, 0);
			}glEnd();
			glColor3f(1, 0, 0);
			glBegin(GL_TRIANGLES); {
				glVertex3f(cp[i].x, cp[i].y, 0);
				glVertex3f(b_x, b_y, 0);
				glVertex3f(c_x, c_y, 0);
			}glEnd();
		}
	}
}
void drawCurve(){
	double p1_x = cp[0].x;
	double p1_y = cp[0].y;
	double r1_x = (cp[1].x - cp[0].x);
	double r1_y = (cp[1].y - cp[0].y);

	crvidx = 0;
	for (int i = 2; i < cpidx + 1; i += 2)
	{
		double p4_x = cp[i % cpidx].x;
		double p4_y = cp[i % cpidx].y;
		double r4_x = (cp[(i + 1) % cpidx].x - cp[i % cpidx].x);
		double r4_y = (cp[(i + 1) % cpidx].y - cp[i % cpidx].y);

		double a_x = 2 * p1_x - 2 * p4_x + r1_x + r4_x;
		double b_x = -3 * p1_x + 3 * p4_x - 2 * r1_x - r4_x;
		double c_x = r1_x;
		double d_x = p1_x;

		double a_y = 2 * p1_y - 2 * p4_y + r1_y + r4_y;
		double b_y = -3 * p1_y + 3 * p4_y - 2 * r1_y - r4_y;
		double c_y = r1_y;
		double d_y = p1_y;

		p1_x = p4_x;
		p1_y = p4_y;
		r1_x = r4_x;
		r1_y = r4_y;

		point2d point[201];
		double t = 0;

		for (int j = 0; j <= 200; j++){
			t = (float)j / 200;
			point[j].x = a_x*t*t*t + b_x*t*t + c_x*t + d_x;
			point[j].y = a_y*t*t*t + b_y*t*t + c_y*t + d_y;

			crvPoint[crvidx].x = point[j].x;
			crvPoint[crvidx].y = point[j].y;
			crvidx++;

		}

		glColor3f(1, 1, 1);
		glLineWidth(2.0);
		for (int j = 0; j < 200; j++)
		{
			glBegin(GL_LINES); {

				glVertex3f(point[j].x, point[j].y, 0);
				glVertex3f(point[j + 1].x, point[j + 1].y, 0);
			}
			glEnd();
		}
		glLineWidth(1.0);
	}
	//glColor3f(1, 1, 1);
	////glLineWidth(2.0);
	//for (int i = 0; i < crvidx - 1; i++)
	//{
	//	glBegin(GL_LINES); {

	//		glVertex3f(crvPoint[i].x, crvPoint[i].y, 0);
	//		glVertex3f(crvPoint[i + 1].x, crvPoint[i + 1].y, 0);
	//	}
	//	glEnd();
	//}
	
}
void followCurve(){	
	glPushMatrix(); {
		glTranslatef(crvPoint[(int)t].x, crvPoint[(int)t].y, 0.1);
		drawCircle(7, 10);
	}glPopMatrix();

}

void keyboardListener(unsigned char key, int x, int y){
	switch (key){

	case '1':
		break;
	case 'u':
		state_ = 2;
		break;
	case 'g':
		drawgm = !drawgm;
		break;
	case 'a':
		follow = !follow;
		t = 0; // start following from the start
		break;
	default:
		break;
	}
}

void specialKeyListener(int key, int x, int y){
	switch (key){
	case GLUT_KEY_DOWN:		//down arrow key
		break;
	case GLUT_KEY_UP:		// up arrow key
		break;

	case GLUT_KEY_RIGHT:
		break;
	case GLUT_KEY_LEFT:
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

double distance(double x1, double y1, double x2, double y2){
	return sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));
}

void mouseListener(int button, int state, int x, int y){	//x, y is the x-y of the screen (2D)
	switch (button){
	case GLUT_LEFT_BUTTON:
		if (state == GLUT_DOWN){		// 2 times?? in ONE click? -- solution is checking DOWN or UP
			if (state_ == 0){
				std::cout << x << " " << y << std::endl;
				cp[cpidx].x = (double)x;
				cp[cpidx].y = (double)(600 - y);
				cpidx++;
			}
			else if (state_==2){
				double min = 99999;
				for (int i = 0; i < cpidx; i++)
				{
					double d = distance(x, 600-y, cp[i].x, cp[i].y);
					if (d < min){
						min = d;
						targetidx = i;
					}
				}
				state_ = 3;
			}
			else if (state_ == 3){
				cp[targetidx].x = (double)x;
				cp[targetidx].y = (double)(600 - y);
				state_ = 1;
			}
		}
		break;

	case GLUT_RIGHT_BUTTON:

		if (state == GLUT_DOWN && cpidx > 2 && cpidx % 2 == 0 && state_== 0){
			drawCrv = true;
			state_ = 1;
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
	//3. Which direction is the camera's UP direction?

	//gluLookAt(100,100,100,	0,0,0,	0,0,1);
	//gluLookAt(150*cos(cameraAngle), 150*sin(cameraAngle), cameraHeight,		0,0,0,		0,0,1);
	gluLookAt(0, 0, 0, 0, 0, -1, 0, 1, 0);


	//again select MODEL-VIEW
	glMatrixMode(GL_MODELVIEW);


	/****************************
	/ Add your objects from here
	****************************/
	//add objects

	if(drawgm)drawGeometry();
	if (drawCrv)drawCurve();
	if(follow && drawCrv)followCurve();
	if (state_ == 3){
		glPushMatrix(); {
			glTranslatef(cp[targetidx].x, cp[targetidx].y, 0.1);
			drawCircle(7, 15);
		}glPopMatrix();
	}

	//ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
	glutSwapBuffers();
}

void animate(){

	
	if (t >= crvidx)t = 0;
	else t += 0.5;
	//codes for any changes in Models, Camera
	glutPostRedisplay();
}

void init(){
	//codes for initialization
	t = 0;
	cpidx = 0;
	state_ = 0;
	drawCrv = false;
	drawgm = true;
	follow = false;
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
	gluOrtho2D(0, 800, 0, 600);
	//gluPerspective(80,	1,	1,	1000.0);
	//field of view in the Y (vertically)
	//aspect ratio that determines the field of view in the X direction (horizontally)
	//near distance
	//far distance
}

int main(int argc, char **argv){
	glutInit(&argc, argv);
	glutInitWindowSize(800, 600);
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
