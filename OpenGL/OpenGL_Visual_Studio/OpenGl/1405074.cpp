#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include <iostream>
#include <map>
#include <glut.h>
#include <fstream>
#include <vector>
#include <string>
#include "bitmap_image.hpp"
#define pi (2*acos(0.0))
#define epsilon (1.0e-6)
using namespace std;

struct point
{
	double x,y,z;
};

class color {
public:
    double r, g, b;
    color(double r, double g, double b) {
        this->r = r;
        this->g = g;
        this->b = b;
    }
    color() {
    }
};

class objects {
public:
    string type;
    point center;
    double radius, width, height;
    double ambient, diffuse, specular, reflection;
    color objColor;
    double shininess;
    objects(){
    }

};

class source{
public:
    point position;
    double falloff;
    point lookAt;
    double cutOff;
};

/*
A simple class to hold the color components, r, g, b of a certain shade.
*/
int spotlightCount;
vector <source> spotLightList;
color **textureBuffer;
int textHeight, textWidth;
bitmap_image b_img ("texture.bmp");
int objectCount, sourceCount;
vector <objects> objectList;
vector <source> sourceList;
color Background(0, 0, 0);
double FoV;
double aspect_ratio;
double near_;
double far_;
double cameraHeight;
double cameraAngle;
int drawgrid;
int drawaxes;
double angle;
int recursionLevel;
int pixel;
double checkerBoardWidth;
point pos, u, r, l;
ifstream description;
double ambientC, diffuseC, reflectionC;
bool bonus;

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


void drawAxes()
{
	if(drawaxes==1)
	{

	glColor3f(1, 0, 0);
		glBegin(GL_LINES);{
			glVertex3f( 2000,0,0);
			glVertex3f(-2000,0,0);

			glVertex3f(0,-2000,0);
			glVertex3f(0, 2000,0);

			glVertex3f(0,0, 2000);
			glVertex3f(0,0,-2000);
		}glEnd();
	}
}

void drawSquare(double a)
{
	glBegin(GL_QUADS);{
		glVertex3f( a/2, a/2,0);
		glVertex3f( a/2,-a/2,0);
		glVertex3f(-a/2,-a/2,0);
		glVertex3f(-a/2, a/2,0);
	}glEnd();
}

void drawCheckerBoardX(){
    for(int i=0;i<=3*far_/checkerBoardWidth; i++){
        glColor3f(1, 1, 1);
        glPushMatrix();
        {
            glTranslatef(2*i*checkerBoardWidth, 0, 0);
            drawSquare(checkerBoardWidth);
        }
        glPopMatrix();
        glPushMatrix();
        {
            glTranslatef(-2*i*checkerBoardWidth, 0, 0);
            drawSquare(checkerBoardWidth);
        }
        glPopMatrix();
        glColor3f(0, 0, 0);
        glPushMatrix();
        {
            glTranslatef((2*i + 1)*checkerBoardWidth, 0, 0);
            drawSquare(checkerBoardWidth);
        }
        glPopMatrix();
        glPushMatrix();
        {
            glTranslatef(-(2*i + 1)*checkerBoardWidth, 0, 0);
            drawSquare(checkerBoardWidth);
        }
        glPopMatrix();
    }
}

void drawCheckerBoard(){

    glColor3f(1, 1, 1);
    glPushMatrix();
    glTranslatef(.5*checkerBoardWidth, .5*checkerBoardWidth, 0);
    for(int i=0;i<=3*far_/checkerBoardWidth; i++){
        glPushMatrix();
        if(i%2){
            glTranslatef(checkerBoardWidth, 0, 0);
        }
        glPushMatrix();
        {
            glTranslatef(0, i*checkerBoardWidth, 0);
            drawCheckerBoardX();
        }
        glPopMatrix();
        glPushMatrix();
        {
            glTranslatef(0, -i*checkerBoardWidth, 0);
            drawCheckerBoardX();
        }
        glPopMatrix();
        glPopMatrix();
    }
    glPopMatrix();
}

void drawPyramid(objects pyramid){
    drawSquare(pyramid.width);

    glBegin(GL_TRIANGLES);
    {
        glVertex3f(pyramid.width/2,pyramid.width/2,0);
        glVertex3f(pyramid.width/2,-pyramid.width/2,0);
        glVertex3f(0,0,pyramid.height);
    }
    glEnd();
    glBegin(GL_TRIANGLES);
    {
        glVertex3f(pyramid.width/2,-pyramid.width/2,0);
        glVertex3f(-pyramid.width/2,-pyramid.width/2,0);
        glVertex3f(0,0,pyramid.height);
    }
    glEnd();glBegin(GL_TRIANGLES);
    {
        glVertex3f(-pyramid.width/2,-pyramid.width/2,0);
        glVertex3f(-pyramid.width/2,pyramid.width/2,0);
        glVertex3f(0,0,pyramid.height);
    }
    glEnd();glBegin(GL_TRIANGLES);
    {
        glVertex3f(-pyramid.width/2,pyramid.width/2,0);
        glVertex3f(pyramid.width/2,pyramid.width/2,0);
        glVertex3f(0,0,pyramid.height);
    }
    glEnd();
}

void drawObjects(){
    for(int i=0;i<objectList.size();i++){
        glColor3f(objectList[i].objColor.r, objectList[i].objColor.g, objectList[i].objColor.b);
        //cout << objectList[i].objColor.r << " " << objectList[i].objColor.g << " " << objectList[i].objColor.b << endl;
        if(objectList[i].type == "sphere"){
            glPushMatrix();
            glTranslatef(objectList[i].center.x, objectList[i].center.y, objectList[i].center.z);
            drawSphere(objectList[i].radius, 50, 50);
            glPopMatrix();
        }
        else{
            glPushMatrix();
            glTranslatef(objectList[i].center.x, objectList[i].center.y, objectList[i].center.z);
            drawPyramid(objectList[i]);
            glPopMatrix();
        }
    }
}

void drawSources(){
    for(int i =0;i<sourceCount;i++){
        glPushMatrix();
        glColor3f(1, 1, 1);
        glTranslatef(sourceList[i].position.x, sourceList[i].position.y, sourceList[i].position.z);
        drawSphere(5, 50, 50);
        glPopMatrix();
    }
    for(int i =0;i<spotlightCount;i++){
        glPushMatrix();
        glColor3f(1, 1, 1);
        glTranslatef(spotLightList[i].position.x, spotLightList[i].position.y, spotLightList[i].position.z);
        drawSphere(5, 50, 50);
        glPopMatrix();
    }
}

void move_right(){
    pos.x += 3*r.x;
    pos.y += 3*r.y;
    pos.z += 3*r.z;
}

void move_left(){
    pos.x -= 3*r.x;
    pos.y -= 3*r.y;
    pos.z -= 3*r.z;
}

void move_up(){
    pos.x += 3*u.x;
    pos.y += 3*u.y;
    pos.z += 3*u.z;
}

void move_down(){
    pos.x -= 3*u.x;
    pos.y -= 3*u.y;
    pos.z -= 3*u.z;
}

void move_forward(){
    pos.x += 3*l.x;
    pos.y += 3*l.y;
    pos.z += 3*l.z;
}

void move_backward_(){
    pos.x -= 3*l.x;
    pos.y -= 3*l.y;
    pos.z -= 3*l.z;
}

void look_right(){
    point temp;
    temp.x = l.x*cos((pi*1)/180) + r.x*sin((pi*1)/180);
    temp.y = l.y*cos((pi*1)/180) + r.y*sin((pi*1)/180);
    temp.z = l.z*cos((pi*1)/180) + r.z*sin((pi*1)/180);
    r.x = r.x*cos((pi*1)/180) - l.x*sin((pi*1)/180);
    r.y = r.y*cos((pi*1)/180) - l.y*sin((pi*1)/180);
    r.z = r.z*cos((pi*1)/180) - l.z*sin((pi*1)/180);
    l = temp;
}

void look_left(){
    point temp;
    temp.x = l.x*cos((pi*1)/180) - r.x*sin((pi*1)/180);
    temp.y = l.y*cos((pi*1)/180) - r.y*sin((pi*1)/180);
    temp.z = l.z*cos((pi*1)/180) - r.z*sin((pi*1)/180);
    r.x = r.x*cos((pi*1)/180) + l.x*sin((pi*1)/180);
    r.y = r.y*cos((pi*1)/180) + l.y*sin((pi*1)/180);
    r.z = r.z*cos((pi*1)/180) + l.z*sin((pi*1)/180);
    l = temp;
}

void look_up(){
    point temp;
    temp.x = l.x*cos((pi*1)/180) + u.x*sin((pi*1)/180);
    temp.y = l.y*cos((pi*1)/180) + u.y*sin((pi*1)/180);
    temp.z = l.z*cos((pi*1)/180) + u.z*sin((pi*1)/180);
    u.x = u.x*cos((pi*1)/180) - l.x*sin((pi*1)/180);
    u.y = u.y*cos((pi*1)/180) - l.y*sin((pi*1)/180);
    u.z = u.z*cos((pi*1)/180) - l.z*sin((pi*1)/180);
    l = temp;
}

void look_down(){
    point temp;
    temp.x = l.x*cos((pi*1)/180) - u.x*sin((pi*1)/180);
    temp.y = l.y*cos((pi*1)/180) - u.y*sin((pi*1)/180);
    temp.z = l.z*cos((pi*1)/180) - u.z*sin((pi*1)/180);
    u.x = u.x*cos((pi*1)/180) + l.x*sin((pi*1)/180);
    u.y = u.y*cos((pi*1)/180) + l.y*sin((pi*1)/180);
    u.z = u.z*cos((pi*1)/180) + l.z*sin((pi*1)/180);
    l = temp;
}

void tilt_clockwise(){
    point temp;
    temp.x = u.x*cos((pi*1)/180) + r.x*sin((pi*1)/180);
    temp.y = u.y*cos((pi*1)/180) + r.y*sin((pi*1)/180);
    temp.z = u.z*cos((pi*1)/180) + r.z*sin((pi*1)/180);
    r.x = r.x*cos((pi*1)/180) - u.x*sin((pi*1)/180);
    r.y = r.y*cos((pi*1)/180) - u.y*sin((pi*1)/180);
    r.z = r.z*cos((pi*1)/180) - u.z*sin((pi*1)/180);
    u = temp;
}

void tilt_counterClockwisw(){
    point temp;
    temp.x = u.x*cos((pi*1)/180) - r.x*sin((pi*1)/180);
    temp.y = u.y*cos((pi*1)/180) - r.y*sin((pi*1)/180);
    temp.z = u.z*cos((pi*1)/180) - r.z*sin((pi*1)/180);
    r.x = r.x*cos((pi*1)/180) + u.x*sin((pi*1)/180);
    r.y = r.y*cos((pi*1)/180) + u.y*sin((pi*1)/180);
    r.z = r.z*cos((pi*1)/180) + u.z*sin((pi*1)/180);
    u = temp;
}

point cross(point a, point b)
    {
        point temp;
        temp.x = a.y*b.z - a.z*b.y;
        temp.y = b.x*a.z - b.z*a.x;
        temp.z = a.x*b.y - a.y*b.x;
        return temp;
    }

double dot(point a, point b){
    return a.x*b.x + a.y*b.y + a.z*b.z;
}

point negateV(point p){
    point temp;
    temp.x = - p.x;
    temp.y = - p.y;
    temp.z = - p.z;
    return temp;
}

point scaleV(point p, double s){
    point temp;
    temp.x = s*p.x;
    temp.y = s*p.y;
    temp.z = s*p.z;
    return temp;
}

point normalize(point p){
    double val = sqrt(p.x*p.x + p.y*p.y + p.z*p.z);
    point temp;
    temp.x = p.x/val;
    temp.y = p.y/val;
    temp.z = p.z/val;
    return temp;
}

point minusV(point a, point b){
    point temp;
    temp.x = a.x - b.x;
    temp.y = a.y - b.y;
    temp.z = a.z - b.z;
    return temp;
}

double length(point a, point b){
    return sqrt((a.x-b.x)*(a.x-b.x) + (a.y-b.y)*(a.y-b.y) + (a.z-b.z)*(a.z-b.z));
}

point getPointAtT(point p, point v, double t){
    point temp;
    temp.x = p.x + t*v.x;
    temp.y = p.y + t*v.y;
    temp.z = p.z + t*v.z;
    return temp;
}


double getFarT(point ray, point startPoint){
    point temp = getPointAtT(startPoint, l, far_-near_);
    double farD = -dot(l, temp);
    return -(farD + dot(l, startPoint))/dot(l, ray);
}



bool checkInside(point t1, point t2, point t3, point p){
    double l1, l2, l3;
    l1 = length(t1, t2);
    l2 = length(t1, t3);
    l3 = length(t3, t2);
    double s = (l1 + l2 + l3)/2;
    double area = sqrt(s*(s-l1)*(s-l2)*(s-l3));
    l1 = length(p, t2);
    l2 = length(p, t3);
    l3 = length(t3, t2);
    s = (l1 + l2 + l3)/2;
    double area1 = sqrt(s*(s-l1)*(s-l2)*(s-l3));
    l1 = length(t1, t2);
    l2 = length(t1, p);
    l3 = length(p, t2);
    s = (l1 + l2 + l3)/2;
    double area2 = sqrt(s*(s-l1)*(s-l2)*(s-l3));
    l1 = length(t1, p);
    l2 = length(t1, t3);
    l3 = length(t3, p);
    s = (l1 + l2 + l3)/2;
    double area3 = sqrt(s*(s-l1)*(s-l2)*(s-l3));
    if (abs(area - (area1 + area2 + area3)) < 0.00001){
        return true;
    }
    else return false;
}

double pyramidT(objects pyramid, point ray, point start){
    double dir[4][2] = {{.5, .5}, {.5, -.5}, {-.5, -.5}, {-.5, .5}};
    double t = 100000;
    for(int i=0;i<4;i++){
        int j = (i+1)%4;
        point t1, t2, t3;
        t1.x = pyramid.center.x + pyramid.width*dir[i][0];
        t1.y = pyramid.center.y + pyramid.width*dir[i][1];
        t1.z = pyramid.center.z;
        t2.x = pyramid.center.x + pyramid.width*dir[j][0];
        t2.y = pyramid.center.y + pyramid.width*dir[j][1];
        t2.z = pyramid.center.z;
        t3.x = pyramid.center.x;
        t3.y = pyramid.center.y;
        t3.z = pyramid.center.z + pyramid.height;
        point a1 = minusV(t1, t2);
        point a2 = minusV(t1, t3);
        point normal = normalize(cross(a1, a2));
        double d = - (normal.x*pyramid.center.x + normal.y*pyramid.center.y + normal.z*(pyramid.center.z+pyramid.height));
        double tps = -(d + dot(normal, start))/dot(normal,ray);
        point intersectPoint = getPointAtT(start, ray, tps);
        if (!checkInside(t1, t2, t3, intersectPoint)) continue;
        t = min(tps, t);
    }
    double baseT = -(-pyramid.center.z + start.z)/ray.z;
    point basePoint = getPointAtT(start, ray, baseT);
    if(basePoint.x>pyramid.center.x+.5*pyramid.width || basePoint.x<pyramid.center.x-.5*pyramid.width)
        return t;
    if(basePoint.y>pyramid.center.y+.5*pyramid.width || basePoint.y<pyramid.center.y-.5*pyramid.width)
        return t;

    return min(t, baseT);
}

double getInterSection(point ray, point startPoint){
    double t = 100000;
    for(int i=0;i<objectList.size();i++){
        if (objectList[i].type == "sphere"){
            double b = 2*dot(startPoint, ray);
            b -= 2*dot(objectList[i].center, ray);
            double c = dot(startPoint, startPoint);
            c = c - objectList[i].radius*objectList[i].radius;
            c -= 2*dot(objectList[i].center, startPoint);
            c += dot(objectList[i].center, objectList[i].center);
            double d = b*b - 4*c;
            if(d>=0){
                d = sqrt(d);
                double ts = min(.5*(-b-d),.5*(-b+d));
                if (ts>0 && ts<t){
                    t = ts;
                    //cout << "sphere" << endl;

                }
            }
        }
        else{
            double ts = pyramidT(objectList[i], ray, startPoint);
            if (ts>0 && ts<t){
                    t = ts;
                    //cout << "pyramid" << endl;

                }
        }
    }
    double tcb = -startPoint.z/ray.z;
    if(tcb>0) t = min(t, tcb);
    return t;
}

color checkColor(color c){
    color temp;
    if(c.r>1) temp.r = 1;
    else temp.r = c.r;
    if(c.g>1) temp.g = 1;
    else temp.g = c.g;
    if(c.b>1) temp.b = 1;
    else temp.b = c.b;
    return temp;
}


color trace(point ray, point startPoint, int level){
    if(level>recursionLevel){
        color temp(0, 0, 0);
        return temp;
    }
    double t = 100000;
    int type = -1;
    objects top;
    point Normal;
    double farT = getFarT(ray, startPoint);
    if(level>0) farT = 1000;
    color retCol = Background;
    for(int i=0;i<objectList.size();i++){
        if (objectList[i].type == "sphere"){
            double b = 2*dot(startPoint, ray);
            b -= 2*dot(objectList[i].center, ray);
            double c = dot(startPoint, startPoint);
            c = c - objectList[i].radius*objectList[i].radius;
            c -= 2*dot(objectList[i].center, startPoint);
            c += dot(objectList[i].center, objectList[i].center);
            double d = b*b - 4*c;
            if(d>=0){
                d = sqrt(d);
                double ts = min(.5*(-b-d),.5*(-b+d));
                if (ts<t && ts>0 && ts<farT){
                    retCol = objectList[i].objColor;
                    t = ts;
                    type = 1;
                    top = objectList[i];
                    point inter = getPointAtT(startPoint, ray, ts);
                    Normal = normalize(minusV(inter, objectList[i].center));
                }
            }
        }
        else{
            objects pyramid = objectList[i];
            double dir[4][2] = {{.5, .5}, {.5, -.5}, {-.5, -.5}, {-.5, .5}};
            for(int k=0;k<4;k++){
                int j = (k+1)%4;
                point t1, t2, t3;
                t1.x = pyramid.center.x + pyramid.width*dir[k][0];
                t1.y = pyramid.center.y + pyramid.width*dir[k][1];
                t1.z = pyramid.center.z;
                t2.x = pyramid.center.x + pyramid.width*dir[j][0];
                t2.y = pyramid.center.y + pyramid.width*dir[j][1];
                t2.z = pyramid.center.z;
                t3.x = pyramid.center.x;
                t3.y = pyramid.center.y;
                t3.z = pyramid.center.z + pyramid.height;
                point a1 = minusV(t1, t2);
                point a2 = minusV(t1, t3);
                point normal = normalize(cross(a1, a2));
                double d = - (normal.x*pyramid.center.x + normal.y*pyramid.center.y + normal.z*(pyramid.center.z+pyramid.height));
                double tps = -(d + dot(normal, startPoint))/dot(normal,ray);
                point intersectPoint = getPointAtT(startPoint, ray, tps);
                if (!checkInside(t1, t2, t3, intersectPoint)) continue;
                if (tps<t && tps>0 && tps<farT){
                    retCol = objectList[i].objColor;
                    t = tps;
                    type = 1;
                    top = objectList[i];
                    Normal = negateV(normal);
                }
            }
            double baseT = -(-pyramid.center.z + startPoint.z)/ray.z;
            point basePoint = getPointAtT(startPoint, ray, baseT);
            if(basePoint.x>pyramid.center.x+.5*pyramid.width || basePoint.x<pyramid.center.x-.5*pyramid.width)
                continue;
            if(basePoint.y>pyramid.center.y+.5*pyramid.width || basePoint.y<pyramid.center.y-.5*pyramid.width)
                continue;
            if (baseT<t && baseT>0 && baseT<farT){
                    retCol = objectList[i].objColor;
                    t = baseT;
                    type = 1;
                    top = objectList[i];
                    point temp;
                    temp.x = 0; temp.y = 0; temp.z = 1;
                    Normal = temp;
                }

        }
    }


    double tcb = -startPoint.z/ray.z;
    //cout << getNFT(ray, near_) << " " << getNFT(ray, far_) << endl;
    if (tcb < t  && tcb>0 && (tcb<farT)){
        t = tcb;
        double x_loc = startPoint.x + t*ray.x;
        double y_loc = startPoint.y + t*ray.y;
        int x_c = (int)(x_loc/checkerBoardWidth);
        int y_c = (int)(y_loc/checkerBoardWidth);
        int pixel_x, pixel_y;
        if(bonus){
            double textX = x_loc - x_c*checkerBoardWidth;
            double textY = y_loc - y_c*checkerBoardWidth;
            if(textX<0) textX = checkerBoardWidth + textX;
            if(textY<0) textY = checkerBoardWidth + textY;
            double pix_x = checkerBoardWidth/textWidth;
            double pix_y = checkerBoardWidth/textHeight;
            pixel_x = (int) textX/pix_x;
            pixel_y = (int) textY/pix_y;
        }
        if (x_loc < 0) x_c = 1 - x_c;
        if (y_loc < 0) y_c = 1 - y_c;
        if (x_c%2 == y_c%2){
            type = 2;
            if(bonus) retCol = textureBuffer[pixel_x][pixel_y];
            else{
                color temp(1, 1, 1);
                retCol = temp;
            }
            point temp1;
            temp1.x = 0; temp1.y = 0; temp1.z = 1;
            Normal = temp1;
        }
        else{
            type = 2;
            if(bonus) retCol = textureBuffer[pixel_x][pixel_y];
            else{
                color temp(0, 0, 0);
                retCol = temp;
            }
            point temp1;
            temp1.x = 0; temp1.y = 0; temp1.z = 1;
            Normal = temp1;
        }
    }
    color finalCol(0, 0, 0);
    if(type != -1){
        double amb, dif, shi, spec, refl;
        point intersect = getPointAtT(startPoint, ray, t);
        if(type == 1){
            amb = top.ambient;
            dif = top.diffuse;
            shi = top.shininess;
            spec = top.specular;
            refl = top.reflection;
        }
        else{
            amb = ambientC;
            dif = diffuseC;
            shi = 0;
            spec = 0;
            refl = reflectionC;
        }
        finalCol.r += retCol.r*amb;
        finalCol.g += retCol.g*amb;
        finalCol.b += retCol.b*amb;
        double lambert = 0, phong = 0;
        for(int i=0;i<sourceCount;i++){
            point sourcePos = sourceList[i].position;
            point sourceDir = normalize(minusV(intersect, sourcePos));
            double interT = getInterSection(sourceDir, sourcePos);
            double checkT = (intersect.x - sourcePos.x)/sourceDir.x;
            if(interT>0 && interT < checkT - .01) continue;
            point toSource = negateV(sourceDir);
            Normal = normalize(Normal);
            double distance = length(sourcePos, intersect);
            double scalingFactor = exp(-distance*distance*sourceList[i].falloff);
            lambert += dot(toSource, Normal)*scalingFactor;
            point reflectedLight = scaleV(Normal, 2*dot(negateV(toSource), Normal));
            reflectedLight = normalize(minusV(negateV(toSource), reflectedLight));
            phong += pow(dot(negateV(reflectedLight), ray), shi) * scalingFactor;
        }
        for(int i=0;i<spotlightCount;i++){
            point sourcePos = spotLightList[i].position;
            point sourceDir = normalize(minusV(intersect, sourcePos));
            double interT = getInterSection(sourceDir, sourcePos);
            double checkT = (intersect.x - sourcePos.x)/sourceDir.x;
            if(interT>0 && interT < checkT - .01) continue;
            point spotDir = minusV(spotLightList[i].lookAt, spotLightList[i].position);
            spotDir = normalize(spotDir);
            double angle = acos(dot(spotDir, normalize(sourceDir)));
            if (angle*180/pi > spotLightList[i].cutOff) continue;
            point toSource = negateV(sourceDir);
            Normal = normalize(Normal);
            double distance = length(sourcePos, intersect);
            double scalingFactor = exp(-distance*distance*spotLightList[i].falloff);
            lambert += dot(toSource, Normal)*scalingFactor;
            point reflectedLight = scaleV(Normal, 2*dot(negateV(toSource), Normal));
            reflectedLight = normalize(minusV(negateV(toSource), reflectedLight));
            phong += pow(dot(negateV(reflectedLight), ray), shi) * scalingFactor;
        }
        point reflectedRay = scaleV(Normal, 2*dot(ray, Normal));
        reflectedRay = normalize(minusV(ray, reflectedRay));
        finalCol.r += retCol.r*dif*lambert + retCol.r*spec*phong;
        finalCol.g += retCol.g*dif*lambert + retCol.g*spec*phong;
        finalCol.b += retCol.b*dif*lambert + retCol.b*spec*phong;
		color reflection(0, 0, 0);
        //color reflection = trace(reflectedRay, getPointAtT(intersect, reflectedRay, .0001), level + 1);
        finalCol.r += refl*reflection.r;
        finalCol.g += refl*reflection.g;
        finalCol.b += refl*reflection.b;

    }
    return checkColor(finalCol);
}

void render(){
    l = normalize(l);
    r = normalize(r);
    u = normalize(u);
    cout << l.x << " " << l.y << " " << l.z << endl;
    cout << r.x << " " << r.y << " " << r.z << endl;
    cout << u.x << " " << u.y << " " << u.z << endl;
    point midPoint = getPointAtT(pos, l, near_);
    double halfD = near_*tan((pi*(FoV/2))/180);
    double incr = (2*halfD)/pixel;
    point cornerPoint;
    cornerPoint.x = midPoint.x + u.x*(halfD-.5*incr) - r.x*(halfD-.5*incr);
    cornerPoint.y = midPoint.y + u.y*(halfD-.5*incr) - r.y*(halfD-.5*incr);
    cornerPoint.z = midPoint.z + u.z*(halfD-.5*incr) - r.z*(halfD-.5*incr);
    color** pixels = new color*[pixel];
    for (int i = 0; i < pixel; i++) {
        pixels[i] = new color [pixel];
        for (int j = 0; j < pixel; j++) {
            pixels[i][j] = Background;
        }
    }
    cout << "Rendering image." << endl;
    for(int i = 0; i<pixel; i++){
        for(int j = 0; j<pixel; j++){
            point startPoint;
            startPoint.x = cornerPoint.x + i*r.x*incr - j*u.x*incr;
            startPoint.y = cornerPoint.y + i*r.y*incr - j*u.y*incr;
            startPoint.z = cornerPoint.z + i*r.z*incr - j*u.z*incr;
            point ray;
            ray.x = startPoint.x - pos.x;
            ray.y = startPoint.y - pos.y;
            ray.z = startPoint.z - pos.z;
            ray = normalize(ray);
            //cout << "(" << i << ", " << "): " << ray.x << " " << ray.y << " " << ray.z << endl;
            pixels[i][j] = trace(ray, startPoint, 0);
        }
        if ((i+1)%(int)(pixel/10) == 0){
            cout << "Rendering " << 100*(i+1)/(pixel) << "% complete." << endl;
        }
    }
    cout << "Rendering complete." << endl;
    bitmap_image image(pixel, pixel);
    for (int x = 0; x < pixel; x++) {
        for (int y = 0; y < pixel; y++) {
            image.set_pixel(x, y, pixels[x][y].r*255, pixels[x][y].g*255, pixels[x][y].b*255);
        }
    }
    image.save_image("out.bmp");
}


void keyboardListener(unsigned char key, int x,int y){
	switch(key){
        case '0':
            render();
            break;
		case '1':
			look_left();
			break;
        case '2':
			look_right();
			break;
        case '3':
			look_up();
			break;
        case '4':
			look_down();
			break;
        case '5':
			tilt_clockwise();
			break;
        case '6':
			tilt_counterClockwisw();
			break;
        case ' ':
            bonus = !bonus;
		default:
			break;
	}
}


void specialKeyListener(int key, int x,int y){
	switch(key){
		case GLUT_KEY_DOWN:		//down arrow key
			move_backward_();
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


		case 32:
		    bonus = !bonus;

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
			if(state == GLUT_DOWN){		// 2 times?? in ONE click? -- solution is checking DOWN or UP

			}
			break;


		default:
			break;
	}
}



void display(){

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClearColor(0,0,0,0);	//color black
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	/********************
	/ set-up camera here
	********************/
	glMatrixMode(GL_MODELVIEW);

	glLoadIdentity();
    gluLookAt(pos.x,pos.y,pos.z,	pos.x + l.x,pos.y + l.y,pos.z + l.z,	u.x,u.y,u.z);


	glMatrixMode(GL_MODELVIEW);


	/****************************
	/ Add your objects from here
	****************************/
    //drawAxes();
    if(!bonus) drawCheckerBoard();
    drawObjects();
    drawSources();
glutSwapBuffers();
}


void animate(){
	glutPostRedisplay();
}

void init(){
    description.open("description.txt");
    bonus = false;
	drawgrid=0;
	drawaxes=1;
	cameraHeight=150.0;
	cameraAngle=1.0;
	angle=0;
	pos.x = 150;
	pos.y = 150;
	pos.z = 60;
	u.x = 0;
	u.y = 0;
	u.z = 1;
	r.x = -1/sqrt(2);
	r.y = 1/sqrt(2);
	r.z = 0;
	l.x = -1/sqrt(2);
	l.y = -1/sqrt(2);
	l.z = 0;

	glClearColor(0,0,0,0);

	/************************
	/ set-up projection here
	************************/
	glMatrixMode(GL_PROJECTION);

	glLoadIdentity();

    description >> near_ >> far_ >> FoV >> aspect_ratio;
    description >> recursionLevel >> pixel >> checkerBoardWidth;
    description >> ambientC >> diffuseC >> reflectionC;
	gluPerspective(FoV,	aspect_ratio, near_, far_);
	description >> objectCount;
	for(int i=0; i<objectCount; i++){
        objects newObject;
        description >> newObject.type;
        description >> newObject.center.x >> newObject.center.y >> newObject.center.z;
        if(newObject.type == "sphere"){
            description >> newObject.radius;
        }
        else{
            description >> newObject.width >> newObject.height;
        }
        description >> newObject.objColor.r >> newObject.objColor.g >> newObject.objColor.b;
        description >> newObject.ambient >> newObject.diffuse >> newObject.specular >> newObject.reflection;
        description >> newObject.shininess;
        objectList.push_back(newObject);
	}
	description >> sourceCount;
	for(int i=0; i<sourceCount; i++){
        source newSource;
        description >> newSource.position.x >> newSource.position.y >> newSource.position.z;
        description >> newSource.falloff;
        sourceList.push_back(newSource);
	}
	description >> spotlightCount;
	for(int i=0; i<spotlightCount; i++){
        source newSource;
        description >> newSource.position.x >> newSource.position.y >> newSource.position.z;
        description >> newSource.falloff;
        description >> newSource.lookAt.x >> newSource.lookAt.y >> newSource.lookAt.z;
        description >> newSource.cutOff;
        spotLightList.push_back(newSource);
	}
	textHeight = b_img.height();
    textWidth = b_img.width();
    textureBuffer = new color* [textWidth];
    for (int i = 0; i < textWidth; i++) {
        textureBuffer[i] = new color [textHeight];
        for (int j = 0; j < textHeight; j++) {
            unsigned char r, g, b;
            b_img.get_pixel(i, j, r, g, b);
            color c(r/255.0, g/255.0, b/255.0);
            textureBuffer[i][j] = c;
        }
    }
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

