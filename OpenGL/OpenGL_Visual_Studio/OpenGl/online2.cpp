#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<iostream>
#include<fstream>
#include <windows.h>
#include <glut.h>

#define pi (2*acos(0.0))

using namespace std;
class point3d
{
public:
	double x, y, z;
};


point3d cp[200];
point3d crvPoint[100000];
point3d crvVec[100000];


int cpidx;
int state_;
/*
0 - input
1 - draw
2 - update mode
3 - select
4 - animate
*/
bool drawgm;
bool follow;
int targetidx;
int crvidx;
double t;
double camH;
ifstream file;

void drawAxes()
{

	glColor3f(1.0, 1.0, 1.0);
	glBegin(GL_LINES); {
		glVertex3f(500, 0, 0);
		glVertex3f(-500, 0, 0);

		glVertex3f(0, -500, 0);
		glVertex3f(0, 500, 0);

		glVertex3f(0, 0, 500);
		glVertex3f(0, 0, -500);
	}glEnd();
	
}


void drawCircle(double radius, int segments)
{
	int i;
	point3d points[100];
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
				glVertex3f(0, 0, 0);
				glVertex3f(points[i].x, points[i].y, 0);
				glVertex3f(points[i + 1].x, points[i + 1].y, 0);

			}
		}
		glEnd();
	}
}

double value_(double x, double y, double z)
{
	return sqrt(x*x + y*y + z*z);
}
void drawCurve(){
	double p1_x = cp[0].x;
	double p1_y = cp[0].y;
	double p1_z = cp[0].z;
	double r1_x = (cp[1].x - cp[0].x);
	double r1_y = (cp[1].y - cp[0].y);
	double r1_z = (cp[1].z - cp[0].z);

	crvidx = 0;
	for (int i = 2; i < cpidx + 1; i += 2)
	{
		double p4_x = cp[i % cpidx].x;
		double p4_y = cp[i % cpidx].y;
		double p4_z = cp[i % cpidx].z;

		double r4_x = (cp[(i + 1) % cpidx].x - cp[i % cpidx].x);
		double r4_y = (cp[(i + 1) % cpidx].y - cp[i % cpidx].y);
		double r4_z = (cp[(i + 1) % cpidx].z - cp[i % cpidx].z);

		double a_x = 2 * p1_x - 2 * p4_x + r1_x + r4_x;
		double b_x = -3 * p1_x + 3 * p4_x - 2 * r1_x - r4_x;
		double c_x = r1_x;
		double d_x = p1_x;

		double a_y = 2 * p1_y - 2 * p4_y + r1_y + r4_y;
		double b_y = -3 * p1_y + 3 * p4_y - 2 * r1_y - r4_y;
		double c_y = r1_y;
		double d_y = p1_y;

		double a_z = 2 * p1_z - 2 * p4_z + r1_z + r4_z;
		double b_z = -3 * p1_z + 3 * p4_z - 2 * r1_z - r4_z;
		double c_z = r1_z;
		double d_z = p1_z;

		p1_x = p4_x;
		p1_y = p4_y;
		p1_z = p4_z;

		r1_x = r4_x;
		r1_y = r4_y;
		r1_z = r4_z;

		point3d point[201];
		double t = 1.0 / 200.0;
		double f_x = d_x;
		double f_y = d_y;
		double f_z = d_z;

		double df_x = a_x*t*t*t + b_x*t*t + c_x*t;
		double df_y = a_y*t*t*t + b_y*t*t + c_y*t;
		double df_z = a_z*t*t*t + b_z*t*t + c_z*t;

		double df2_x = 6 * a_x*t*t*t + 2 * b_x*t*t;
		double df2_y = 6 * a_y*t*t*t + 2 * b_y*t*t;
		double df2_z = 6 * a_z*t*t*t + 2 * b_z*t*t;
		
		double df3_x = 6 * a_x*t*t*t;
		double df3_y = 6 * a_y*t*t*t;
		double df3_z = 6 * a_z*t*t*t;

		for (int j = 0; j <= 200; j++){
			/*t = (float)j / 200;*/
			point[j].x = f_x;
			point[j].y = f_y;
			point[j].z = f_z;

			crvPoint[crvidx].x = point[j].x;
			crvPoint[crvidx].y = point[j].y;
			crvPoint[crvidx].z = point[j].z;

			crvVec[crvidx].x = df_x / value_(df_x, df_y, df_z);
			crvVec[crvidx].y = df_y / value_(df_x, df_y, df_z);
			crvVec[crvidx].z = df_z / value_(df_x, df_y, df_z);

			crvidx++;

			f_x += df_x;
			f_y += df_y;
			f_z += df_z;

			df_x += df2_x;
			df_y += df2_y;
			df_z += df2_z;

			df2_x += df3_x;
			df2_y += df3_y;
			df2_z += df3_z;
		}

		glColor3f(1, 0, 0);
		glLineWidth(2.0);
		for (int j = 0; j < 200; j++)
		{
			glBegin(GL_LINES); {

				glVertex3f(point[j].x, point[j].y, point[j].z);
				glVertex3f(point[j + 1].x, point[j + 1].y, point[j + 1].z);
			}
			glEnd();
		}
		glLineWidth(1.0);
	}


}
void followCurve(){
	glPushMatrix(); {

		glTranslatef(crvPoint[(int)t].x, crvPoint[(int)t].y, crvPoint[(int)t].z);
		glRotatef(asin(crvVec[(int)t].x), 1, 0, 0);
		glRotatef(asin(crvVec[(int)t].y), 0, 1, 0);
		glRotatef(asin(crvVec[(int)t].z), 0, 0, 1);

		//glRotatef(crvVec[(int)t].x, crvVec[(int)t].y, crvVec[(int)t].z);
		drawCircle(10, 15);
	}glPopMatrix();

}

void keyboardListener(unsigned char key, int x, int y){
	switch (key){

	case '1':
		break;

	case 'a':
		if (state_ == 1)state_ = 0;
		else if (state_ == 0)state_ = 1;
		t = 0; // start following from the start
		break;
	default:
		break;
	}
}

void specialKeyListener(int key, int x, int y){
	switch (key){
	case GLUT_KEY_DOWN:		//down arrow key
		camH-=2;
		break;
	case GLUT_KEY_UP:		// up arrow key
		camH+=2;
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

double distance(double x1, double y1, double z1,double x2, double y2, double z2){
	return sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2) + (z1 - z2)*(z1 - z2));
}

void mouseListener(int button, int state, int x, int y){	//x, y is the x-y of the screen (2D)
	switch (button){
	case GLUT_LEFT_BUTTON:
		if (state == GLUT_DOWN){		// 2 times?? in ONE click? -- solution is checking DOWN or UP
		}
		break;

	case GLUT_RIGHT_BUTTON:

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
	gluLookAt(500, 500, camH, 0, 0, 0, 0, 0, 1);


	//again select MODEL-VIEW
	glMatrixMode(GL_MODELVIEW);


	/****************************
	/ Add your objects from here
	****************************/
	//add objects
	drawAxes();

	drawCurve();
	if (state_ == 1)followCurve();

	//ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
	glutSwapBuffers();
}


void animate(){

	if (t >= crvidx)t = 0;

	else t += 0.05;
	//else{
	//	int ci = (int)t;
	//	point2d p1 = crvPoint[ci % crvidx];
	//	point2d p2 = crvPoint[(ci + 1) % crvidx];
	//	double div=distance(p1.x, p1.y, p2.x, p2.y);
	//	t += 0.3/div;
	//}
	//codes for any changes in Models, Camera
	glutPostRedisplay();
}

void init(){
	//codes for initialization
	t = 0;
	cpidx = 0;
	state_ = 0;
	drawgm = true;
	follow = false;
	camH = 200;

	file.open("model.txt");
	int n;
	file >> n;
	cout << n << endl;
	cpidx = n;
	for (int i = 0; i < n; i++){
		file >> cp[i].x >> cp[i].y >> cp[i].z;
	}

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
	//gluOrtho2D(0, 800, 0, 600);
	gluPerspective(80,	1,	1,	2000.0);
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
