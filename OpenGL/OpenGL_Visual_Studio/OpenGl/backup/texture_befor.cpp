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
bitmap_image b_img("texture.bmp");
#define pi (2*acos(0.0))
#define epsilon (1.0e-6)
using namespace std;
class point
{
public:
	double x, y, z;
	point(){

	}
	point(double x, double y, double z) {
		this->x = x;
		this->y = y;
		this->z = z;
	}
	point operator+(const point& p) {
		point p1(x + p.x, y + p.y, z + p.z);
		return p1;
	}

	point operator-(const point& p) {
		point p1(x - p.x, y - p.y, z - p.z);
		return p1;
	}

	point operator*(double m) {
		point p(x*m, y*m, z * m);
		return p;
	}

	static double dot(point a, point b) {
		return a.x * b.x + a.y * b.y + a.z * b.z;
	}

	static double distance(point a, point b) {
		return sqrt((a.x - b.x)*(a.x - b.x) + (a.y - b.y)*(a.y - b.y) + (a.z - b.z)*(a.z - b.z));
	}
};

class Color {
public:
	double r, g, b;
	Color(double r, double g, double b) {
		this->r = r;
		this->g = g;
		this->b = b;
	}
	Color() {
	}
	void print(){
		cout << r << " " << g << " " << b << endl;
	}
	void normalize(){
		double s = r + g + b;
		r /= r / s;
		g /= g / s;
		b /= b / s;
	}
};

class Vector {
public:
	double x, y, z;

	Vector() {
	}
	// constructs a vector with given components

	Vector(double x, double y, double z) {
		this->x = x;
		this->y = y;
		this->z = z;
	}

	// keeps the direction same. recalculates the vector to be unit.

	void normalize() {
		double r = sqrt(x * x + y * y + z * z);
		x = x / r;
		y = y / r;
		z = z / r;
	}

	double value(){
		return sqrt(x*x + y*y + z*z);
	}

	Vector rotate(Vector X, Vector a, double angle) {
		double rad = angle * (pi / 180.0);
		Vector res = X * cos(rad) + a * (X.dot(X, a)*(1 - cos(rad))) + X.cross(a, X) * sin(rad);
		return res;
	}

	// add two vectors

	Vector operator+(const Vector& v) {
		Vector v1(x + v.x, y + v.y, z + v.z);
		return v1;
	}

	// subtract one vector from another

	Vector operator-(const Vector& v) {
		Vector v1(x - v.x, y - v.y, z - v.z);
		return v1;
	}

	// scale a vector with a given coefficient

	Vector operator*(double m) {
		Vector v(x*m, y*m, z * m);
		return v;
	}


	static double dot(Vector a, point b) {
		return a.x * b.x + a.y * b.y + a.z * b.z;
	}

	// get the dot product of two vectors

	static double dot(Vector a, Vector b) {
		return a.x * b.x + a.y * b.y + a.z * b.z;
	}

	// get the cross product of two vectors

	static Vector cross(Vector a, Vector b) {
		Vector v(a.y * b.z - a.z * b.y, b.x * a.z - b.z * a.x, a.x * b.y - a.y * b.x);
		return v;
	}

	// print a vector. only for testing purposes.

	void print() {
		std::cout << "Vector" << std::endl;
		std::cout << x << " " << y << " " << z << std::endl;
	}
};

class position
{
public:
	double x, y, z;
};

class Up
{
public:
	double x, y, z;
};

class Right
{
public:
	double x, y, z;
};

class Look
{
public:
	double x, y, z;
};

class checkerboard{
public:
	double a;
	double ambient;
	double diffuse;
	double reflection;
};

class sphere{
public:
	point center;
	double radius;
	Color color;
	double ambient;
	double diffuse;
	double specular;
	double reflection;
	double shininess;
};

class pyramid{
public:
	point lowest;
	double height;
	double width;
	Color color;
	double ambient;
	double diffuse;
	double specular;
	double reflection;
	double shininess;
};

class spotlight{
public:
	point position;
	double falloff; 
	point look;
	double cutoff;
};

class normlight{
public:
	point position;
	double falloff;
};
double cameraAngle;

int drawaxes;

position pos;
Up u;
Right r;
Look l;
double near_;
double far_;
double fovY;
double aspRatio;
int recLevel;
int nPixels;
int nObjets;
int nNorms;
int nSpots;
int gridDim;
Color backgroud;
checkerboard infchecker;
vector<sphere> spheres;
vector<pyramid> pyramids;
vector<spotlight> spotlights;
vector<normlight> normlights;
map<pair<int, int>, Color>gridColMap;
bool texture_rend;

point intersection;
Vector Norm;
Color surf_color;
double surf_amb;
double surf_diff;
double surf_spec;
double surf_ref;
double surf_shine;

Color **textureBuffer;
int text_height, text_width;

void setTextureBuffer(){
	int height, width;
	text_height= height = b_img.height();
	text_width = width = b_img.width();
	textureBuffer = new Color*[width];
	for (int i = 0; i < width; i++) {
		textureBuffer[i] = new Color[height];
		for (int j = 0; j < height; j++) {
			unsigned char r, g, b;
			b_img.get_pixel(i, j, r, g, b);
			Color c(r / 255.0, g / 255.0, b / 255.0);
			textureBuffer[i][j] = c;
		}
	}
}

Color textureColor(double img_height, double img_width, int pos_y, int pos_x){
	int text_y = (text_height / img_height)*pos_y;
	int text_x = (text_width / img_width)*pos_x;
	
	if (text_y > text_height || text_x > text_width)cout << "something's wrong";
	Color ret_col = textureBuffer[text_y][text_x];

	return ret_col;
}
double triangleArea(point p1, point p2, point p3){
	Vector a(p2.x - p1.x, p2.y - p1.y, p2.z - p1.z);
	Vector b(p3.x - p1.x, p3.y - p1.y, p3.z - p1.z);
	double area = a.cross(a, b).value() / 2;
	return area;
}


bool isInside(point p1, point p2, point p3, point p){
	double area1 = triangleArea(p, p1, p2);
	double area2 = triangleArea(p, p2, p3);
	double area3 = triangleArea(p, p3, p1);
	double area = triangleArea(p1, p2, p3);
	double total_area = area1 + area2 + area3;
	return abs(total_area-area) < epsilon ? true : false;
}

bool check_intersect(point P, Vector toSource, double distance){
	double t;
	double D = 0;
	//point intersection; Vector Norm; Color surf_color;
	Vector n(0, 0, 1);
	t = -(D + n.dot(n, P)) / (n.dot(n, toSource));
	if (t > 0){
		intersection.x = P.x + t*toSource.x;
		intersection.y = P.y + t*toSource.y;
		intersection.z = P.z + t*toSource.z;
		Norm = n;
		
		if (P.distance(P, intersection) < distance){
			int grid_i, grid_j;
			if (intersection.x > 0)grid_i = intersection.x / infchecker.a;
			else grid_i = gridDim / 2 - intersection.x / infchecker.a + 1;

			if (intersection.y > 0)grid_j = intersection.y / infchecker.a;
			else grid_j = gridDim / 2 - intersection.y / infchecker.a + 1;

			surf_color = gridColMap[make_pair(grid_i, grid_j)];
			surf_amb=infchecker.ambient;
			surf_diff = infchecker.diffuse;
			surf_spec = 0;
			surf_ref = infchecker.reflection;
			surf_shine = 0;
			return true;
		}
	}

	
	for (int i = 0; i < spheres.size(); i++){		
		
		point origin_t;
		origin_t.x = P.x - spheres[i].center.x;
		origin_t.y = P.y - spheres[i].center.y;
		origin_t.z = P.z - spheres[i].center.z;

		double Rdd = toSource.dot(toSource, toSource);
		double Roo = P.dot(origin_t, origin_t);
		double Rdo = toSource.dot(toSource, origin_t);

		double a = 1;
		double b = 2 * Rdo;
		double c = Roo - spheres[i].radius*spheres[i].radius;
		if ((b*b - 4 * a*c) > 0){
			double t1 = (-b + sqrt(b*b - 4 * a*c)) / (2 * a);
			double t2 = (-b - sqrt(b*b - 4 * a*c)) / (2 * a);
			t = min(t1, t2);
			if (t > 0){
				intersection.x = P.x + t*toSource.x;
				intersection.y = P.y + t*toSource.y;
				intersection.z = P.z + t*toSource.z;
				Norm = Vector(intersection.x, intersection.y, intersection.z);
				Norm.normalize();
				surf_color = spheres[i].color;
				surf_amb = spheres[i].ambient;
				surf_diff = spheres[i].diffuse;
				surf_spec = spheres[i].specular;
				surf_ref = spheres[i].reflection;
				surf_shine = spheres[i].shininess;
				if (P.distance(P, intersection) < distance)return true;
			}
		}
	}

	for (int i = 0; i < pyramids.size(); i++){
		//double D;
		point top(pyramids[i].width / 2, pyramids[i].width / 2, pyramids[i].height);
		point p1(0, 0, 0);
		point p2(0, pyramids[i].width, 0);
		point p3(pyramids[i].width, pyramids[i].width, 0);
		point p4(pyramids[i].width, 0, 0);
		Vector v1(p1.x - top.x, p1.y - top.y, p1.z - top.z);
		Vector v2(p2.x - top.x, p2.y - top.y, p2.z - top.z);
		Vector v3(p3.x - top.x, p3.y - top.y, p3.z - top.z);
		Vector v4(p4.x - top.x, p4.y - top.y, p4.z - top.z);
		point origin_t;
		origin_t.x = P.x - pyramids[i].lowest.x;
		origin_t.y = P.y - pyramids[i].lowest.y;
		origin_t.z = P.z - pyramids[i].lowest.z;

		surf_color = pyramids[i].color;
		surf_color = pyramids[i].color;
		surf_amb = pyramids[i].ambient;
		surf_diff = pyramids[i].diffuse;
		surf_spec = pyramids[i].specular;
		surf_ref = pyramids[i].reflection;
		surf_shine = pyramids[i].shininess;
		//side-1
		Vector n1 = v1.cross(v1, v2);
		n1.normalize();
		D = -n1.dot(n1, top);
		t = -(D + n1.dot(n1, P)) / n1.dot(n1, toSource);
		intersection.x = origin_t.x + t*toSource.x;
		intersection.y = origin_t.y + t*toSource.y;
		intersection.z = origin_t.z + t*toSource.z;
		Norm = n1;
		if (isInside(top, p1, p2, intersection))
		if (P.distance(P, intersection) < distance)return true;
		//side-2
		Vector n2 = v2.cross(v2, v3);
		n2.normalize();
		D = -n2.dot(n2, top);
		t = -(D + n2.dot(n2, P)) / n2.dot(n2, toSource);
		if (t>0){
			intersection.x = origin_t.x + t*toSource.x;
			intersection.y = origin_t.y + t*toSource.y;
			intersection.z = origin_t.z + t*toSource.z;
			Norm = n2;
			if (isInside(top, p2, p3, intersection))
			if (P.distance(P, intersection) < distance)return true;
		}
		//side-3
		Vector n3 = v3.cross(v3, v4);
		n3.normalize();
		D = -n3.dot(n3, top);
		t = -(D + n3.dot(n3, P)) / n3.dot(n3, toSource);
		if (t>0){
			intersection.x = origin_t.x + t*toSource.x;
			intersection.y = origin_t.y + t*toSource.y;
			intersection.z = origin_t.z + t*toSource.z;
			Norm = n3;
			if (isInside(top, p3, p4, intersection))
			if (P.distance(P, intersection) < distance)return true;
		}
		//side-4
		Vector n4 = v4.cross(v4, v1);
		n4.normalize();
		D = -n4.dot(n4, top);
		t = -(D + n4.dot(n4, P)) / n4.dot(n4, toSource);
		if (t>0){
			intersection.x = origin_t.x + t*toSource.x;
			intersection.y = origin_t.y + t*toSource.y;
			intersection.z = origin_t.z + t*toSource.z;
			Norm = n4;
			if (isInside(top, p4, p1, intersection))
			if (P.distance(P, intersection) < distance)return true;
		}
		//base
		Vector n5 = Vector(0, 0, 1);
		D = 0;
		t = -(D + n5.dot(n5, P)) / n5.dot(n5, toSource);
		if (t>0){
			intersection.x = origin_t.x + t*toSource.x;
			intersection.y = origin_t.y + t*toSource.y;
			intersection.z = origin_t.z + t*toSource.z;
			Norm = n5;
			if ((intersection.x > 0 && intersection.x < pyramids[i].width&&intersection.y>0 && intersection.y < pyramids[i].width))
			if (P.distance(P, intersection) < distance)return true;
		}
	}

	return false;
}

bool check_cutoff(point P, point S, point look, double cutoff){
	Vector v1(P.x - S.x, P.y - S.y, P.z - S.z);
	v1.normalize();
	Vector v2(look.x - S.x, look.y - S.y, look.z - S.z);
	v2.normalize();
	double angle = acos(v1.dot(v1, v2))*(180 / pi);
	//cout << angle << endl;
	return angle > cutoff ? true : false;
}


Color setColor(Color current, double amb, double diff, double spec, double ref,  double shininess, point P, Vector fromEye, Vector N, int depthOfRecursion){
	double lambert = 0, phong = 0;
	Vector R = fromEye - N * 2 * N.dot(fromEye, N);
	R.normalize();
	Color reflected;
	//reflectance
	if (depthOfRecursion == 0)return Color(0, 0, 0);
	else {
		
		if (check_intersect(P, R, far_)){
			reflected = setColor(surf_color, surf_amb, surf_diff, surf_spec, surf_ref, surf_shine, intersection, R, Norm, depthOfRecursion - 1);
		}
		else return Color(0, 0, 0);
		//reflected = setColor(current, amb, diff, spec, ref, shininess, P, R, N, depthOfRecursion - 1);
	}

	//diffuse
	for (int i = 0; i < normlights.size(); i++){
		point S = normlights[i].position;
		Vector toSource(P.x - S.x, P.y - S.y, P.z - S.z);
		double distance = toSource.value();
		toSource.normalize();
		if (check_intersect(P, toSource, distance))continue;
		double scaling_factor = exp(-distance*distance*normlights[i].falloff);
		lambert += max(N.dot(toSource, N)*scaling_factor, 0.0);
		phong += max(pow(R.dot(R, N), shininess)*scaling_factor, 0.0);
	}
	//specular
	for (int i = 0; i < spotlights.size(); i++){
		point S = spotlights[i].position;
		Vector toSource(P.x - S.x, P.y - S.y, P.z - S.z);
		double distance = toSource.value();
		toSource.normalize();
		if (check_intersect(P, toSource, distance) || check_cutoff(P, S, spotlights[i].look, spotlights[i].cutoff))continue;
		double scaling_factor = exp(-distance*distance*normlights[i].falloff);
		//cout << N.dot(toSource, N) << " " << scaling_factor << endl;
		lambert += N.dot(toSource, N)*scaling_factor, 0.0;
		phong += pow(R.dot(R, N), shininess)*scaling_factor, 0.0;
	}
	lambert = max(epsilon, lambert);
	phong = max(epsilon, phong);
	Color color;
	color.r = diff*lambert*current.r + diff*phong*current.r + amb*current.r + ref*reflected.r;
	color.g = diff*lambert*current.g + diff*phong*current.g + amb*current.g + ref*reflected.g;
	color.b = diff*lambert*current.b + diff*phong*current.b + amb*current.b + ref*reflected.b;
	//cout << phong << " " << lambert << endl;
	return color;
}

void render_image(){
	Color** pixels = new Color*[nPixels];
	double** ts = new double*[nPixels];
	for (int i = 0; i < nPixels; i++) {
		pixels[i] = new Color[nPixels];
		for (int j = 0; j < nPixels; j++) {
			pixels[i][j] = backgroud;
		}
		ts[i] = new double[nPixels];
		for (int j = 0; j < nPixels; j++) {
			ts[i][j] = +9999999; // a very large value intended as +INFINITY
		}
	}

	vector<point> pointBuffer;
	double fovX = fovY*aspRatio;
	double screen_height = near_*tan(fovY*pi / 360);
	double screen_width = near_*tan(fovX*pi / 360);
	double dy = (2 * screen_height) / nPixels;
	double dx = (2 * screen_width) / nPixels;
	point toppoint;
	toppoint.x = pos.x + l.x*near_ + screen_height*u.x - screen_width*r.x;
	toppoint.y = pos.y + l.y*near_ + screen_height*u.y - screen_width*r.y;
	toppoint.z = pos.z + l.z*near_ + screen_height*u.z - screen_width*r.z;
	for (int i = 0; i < nPixels; i++){
		for (int j = 0; j < nPixels; j++){
			point new_point;
			new_point.x = toppoint.x - i*dy*u.x + j*dx*r.x;
			new_point.y = toppoint.y - i*dy*u.y + j*dx*r.y;
			new_point.z = toppoint.z - i*dy*u.z + j*dx*r.z;
			pointBuffer.push_back(new_point);
		}
	}
	for (int i = 0; i < pointBuffer.size(); i++){
		int pixel_x = i%nPixels;
		int pixel_y = (int)i / nPixels;

		point origin;
		origin.x = pointBuffer[i].x;
		origin.y = pointBuffer[i].y;
		origin.z = pointBuffer[i].z;
		Vector Rd;
		Rd.x = pointBuffer[i].x - pos.x;
		Rd.y = pointBuffer[i].y - pos.y;
		Rd.z = pointBuffer[i].z - pos.z;
		Rd.normalize();
		double D = 0;
		Vector n(0,0,1);
		double t = -(D + n.dot(n, origin)) / (n.dot(n, Rd));
		if (t > 0 && t <= far_){
			point intersect;
			intersect.x = origin.x + t*Rd.x;
			intersect.y = origin.y + t*Rd.y;
			intersect.z = origin.z + t*Rd.z;
			int grid_i, grid_j;
			int img_y, img_x;
			if (intersect.x > 0){
				grid_i = intersect.x / infchecker.a;
				img_x = int(intersect.x) % int(infchecker.a);
			}
			else{
				grid_i = gridDim / 2 - intersect.x / infchecker.a + 1;
				img_x = int(gridDim / 2 - intersect.x) % int(infchecker.a);
			}

			if (intersect.y > 0){
				grid_j = intersect.y / infchecker.a;
				img_y = int(intersect.y) % int(infchecker.a);
			}
			else{
				grid_j = gridDim / 2 - intersect.y / infchecker.a + 1;
				img_y = int(gridDim / 2 - intersect.y) % int(infchecker.a);
			}


			if (t < ts[pixel_x][pixel_y]){
				Color curr_col;
				if (texture_rend)curr_col = textureColor(infchecker.a, infchecker.a, img_y, img_x);
				else curr_col = gridColMap[make_pair(grid_i, grid_j)];;
				double amb = infchecker.ambient;
				double diff = infchecker.diffuse;
				double spec = 0;
				double ref = 0;
				double shine = 0;
				point P = intersect;
				Vector fromEye(P.x - origin.x, P.y - origin.y, P.z - origin.z);
				fromEye.normalize();
				Color pixel_color = setColor(curr_col, amb, diff, spec, ref, shine, P, fromEye, n, recLevel);
				pixels[pixel_x][pixel_y] = pixel_color;
				ts[pixel_x][pixel_y] = t;
			}

		}
		for (int j = 0; j<pyramids.size(); j++){
			double D, t;
			point top(pyramids[j].width / 2, pyramids[j].width / 2, pyramids[j].height);
			point p1(0, 0, 0);
			point p2(0, pyramids[j].width, 0);
			point p3(pyramids[j].width, pyramids[j].width, 0);
			point p4(pyramids[j].width, 0, 0);

			Vector v1(p1.x - top.x, p1.y - top.y, p1.z - top.z);
			Vector v2(p2.x - top.x, p2.y - top.y, p2.z - top.z);
			Vector v3(p3.x - top.x, p3.y - top.y, p3.z - top.z);
			Vector v4(p4.x - top.x, p4.y - top.y, p4.z - top.z);

			point origin_t;
			origin_t.x = origin.x - pyramids[j].lowest.x;
			origin_t.y = origin.y - pyramids[j].lowest.y;
			origin_t.z = origin.z - pyramids[j].lowest.z;
			//side-1
			Vector n1 = v1.cross(v1, v2);
			n1.normalize();
			D = -n1.dot(n1, top);
			t = -(D + n1.dot(n1, origin_t)) / n1.dot(n1, Rd);
			if (t > 0 && t <= far_){
				point intersect;
				intersect.x = origin_t.x + t*Rd.x;
				intersect.y = origin_t.y + t*Rd.y;
				intersect.z = origin_t.z + t*Rd.z;

				if (t < ts[pixel_x][pixel_y] && isInside(top, p1, p2, intersect)){
					Color curr_col = pyramids[j].color;
					double amb = pyramids[j].ambient;
					double diff = pyramids[j].diffuse;
					double spec = pyramids[j].specular;
					double ref = pyramids[j].reflection;
					double shine = pyramids[j].shininess;
					point P = intersect;
					Vector fromEye(P.x - origin_t.x, P.y - origin_t.y, P.z - origin_t.z);
					fromEye.normalize();
					Color pixel_color = setColor(curr_col, amb, diff, spec, ref, shine, P, fromEye, n, recLevel);
					pixels[pixel_x][pixel_y] = pixel_color;
					//pixels[pixel_x][pixel_y] = pyramids[j].color;
					ts[pixel_x][pixel_y] = t;
				}
			}
			//side-2
			Vector n2 = v2.cross(v2, v3);
			n2.normalize();
			D = -n2.dot(n2, top);
			t = -(D + n2.dot(n2, origin_t)) / (n2.dot(n2, Rd));
			if (t > 0 && t <= far_){
				point intersect;
				intersect.x = origin_t.x + t*Rd.x;
				intersect.y = origin_t.y + t*Rd.y;
				intersect.z = origin_t.z + t*Rd.z;

				if (t < ts[pixel_x][pixel_y] && isInside(top, p2, p3, intersect)){
					Color curr_col = pyramids[j].color;
					double amb = pyramids[j].ambient;
					double diff = pyramids[j].diffuse;
					double spec = pyramids[j].specular;
					double ref = pyramids[j].reflection;
					double shine = pyramids[j].shininess;
					point P = intersect;
					Vector fromEye(P.x - origin_t.x, P.y - origin_t.y, P.z - origin_t.z);
					fromEye.normalize();
					Color pixel_color = setColor(curr_col, amb, diff, spec, ref, shine, P, fromEye, n, recLevel);
					pixels[pixel_x][pixel_y] = pixel_color;
					//pixels[pixel_x][pixel_y] = pyramids[j].color;
					ts[pixel_x][pixel_y] = t;
				}
			}
			//side-3
			Vector n3 = v3.cross(v3, v4);
			n3.normalize();
			D = -n3.dot(n3, top);
			t = -(D + n3.dot(n3, origin_t)) / (n3.dot(n3, Rd));
			if (t > 0 && t <= far_){
				point intersect;
				intersect.x = origin_t.x + t*Rd.x;
				intersect.y = origin_t.y + t*Rd.y;
				intersect.z = origin_t.z + t*Rd.z;
				if (t < ts[pixel_x][pixel_y] && isInside(top, p3, p4, intersect)){
					Color curr_col = pyramids[j].color;
					double amb = pyramids[j].ambient;
					double diff = pyramids[j].diffuse;
					double spec = pyramids[j].specular;
					double ref = pyramids[j].reflection;
					double shine = pyramids[j].shininess;
					point P = intersect;
					Vector fromEye(P.x - origin_t.x, P.y - origin_t.y, P.z - origin_t.z);
					fromEye.normalize();
					Color pixel_color = setColor(curr_col, amb, diff, spec, ref, shine, P, fromEye, n, recLevel);
					pixels[pixel_x][pixel_y] = pixel_color;
					//pixels[pixel_x][pixel_y] = pyramids[j].color;
					ts[pixel_x][pixel_y] = t;
				}
			}

			//side-4
			Vector n4 = v4.cross(v4, v1);
			n4.normalize();
			D = -n4.dot(n4, top);
			t = -(D + n4.dot(n4, origin_t)) / (n4.dot(n4, Rd));
			if (t > 0 && t <= far_){
				point intersect;
				intersect.x = origin_t.x + t*Rd.x;
				intersect.y = origin_t.y + t*Rd.y;
				intersect.z = origin_t.z + t*Rd.z;

				if (t < ts[pixel_x][pixel_y] && isInside(top, p4, p1, intersect)){
					Color curr_col = pyramids[j].color;
					double amb = pyramids[j].ambient;
					double diff = pyramids[j].diffuse;
					double spec = pyramids[j].specular;
					double ref = pyramids[j].reflection;
					double shine = pyramids[j].shininess;
					point P = intersect;
					Vector fromEye(P.x - origin_t.x, P.y - origin_t.y, P.z - origin_t.z);
					fromEye.normalize();
					Color pixel_color = setColor(curr_col, amb, diff, spec, ref, shine, P, fromEye, n, recLevel);
					pixels[pixel_x][pixel_y] = pixel_color;
					//pixels[pixel_x][pixel_y] = pyramids[j].color;
					ts[pixel_x][pixel_y] = t;
				}
			}
			//base
			Vector n5 = Vector(0, 0, 1);
			D = 0;
			t = -(D + n5.dot(n5, origin_t)) / (n5.dot(n5, Rd));
			if (t > 0 && t <= far_){
				point intersect;
				intersect.x = origin_t.x + t*Rd.x;
				intersect.y = origin_t.y + t*Rd.y;
				intersect.z = origin_t.z + t*Rd.z;

				if (t < ts[pixel_x][pixel_y] && (intersect.x>0 && intersect.x<pyramids[j].width&&intersect.y>0 && intersect.y<pyramids[j].width)){
					Color curr_col = pyramids[j].color;
					double amb = pyramids[j].ambient;
					double diff = pyramids[j].diffuse;
					double spec = pyramids[j].specular;
					double ref = pyramids[j].reflection;
					double shine = pyramids[j].shininess;
					point P = intersect;
					Vector fromEye(P.x - origin_t.x, P.y - origin_t.y, P.z - origin_t.z);
					fromEye.normalize();
					Color pixel_color = setColor(curr_col, amb, diff, spec, ref, shine, P, fromEye, n, recLevel);
					pixels[pixel_x][pixel_y] = pixel_color;
					//pixels[pixel_x][pixel_y] = pyramids[j].color;
					ts[pixel_x][pixel_y] = t;
				}
			}
		}
		for (int j = 0; j < spheres.size(); j++){
			point origin_t;
			origin_t.x = origin.x - spheres[j].center.x;
			origin_t.y = origin.y - spheres[j].center.y;
			origin_t.z = origin.z - spheres[j].center.z;

			double Rdd = Rd.dot(Rd, Rd);
			double Roo = origin_t.dot(origin_t, origin_t);
			double Rdo = Rd.dot(Rd, origin_t);

			double a = 1;
			double b = 2 * Rdo;
			double c = Roo - spheres[j].radius*spheres[j].radius;
			if ((b*b - 4 * a*c) > 0){
				double t1 = (-b + sqrt(b*b - 4 * a*c)) / (2 * a);
				double t2 = (-b - sqrt(b*b - 4 * a*c)) / (2 * a);
				double t = min(t1, t2);
				point intersect;
				intersect.x = origin_t.x + t*Rd.x;
				intersect.y = origin_t.y + t*Rd.y;
				intersect.z = origin_t.z + t*Rd.z;
				if (t < ts[pixel_x][pixel_y] && t>0 && t <= far_){
					Color curr_col = spheres[j].color;
					double amb = spheres[j].ambient;
					double diff = spheres[j].diffuse;
					double spec = spheres[j].specular;
					double ref = spheres[j].reflection;
					double shine = spheres[j].shininess;
					point P = intersect;
					Vector n = Vector(intersect.x, intersect.y, intersect.z);
					n.normalize();
					Vector fromEye(P.x - origin_t.x, P.y - origin_t.y, P.z - origin_t.z);
					fromEye.normalize();
					Color pixel_color = setColor(curr_col, amb, diff, spec, ref, shine, P, fromEye, n, recLevel);
					pixels[pixel_x][pixel_y] = pixel_color;
					//pixels[pixel_x][pixel_y] = spheres[j].color;
					ts[pixel_x][pixel_y] = t;
				}
			}
		}

	}
	cout << "done" << endl;
	
	bitmap_image image(nPixels, nPixels);
	for (int x = 0; x < nPixels; x++) {
		for (int y = 0; y < nPixels; y++) {
			image.set_pixel(x, y, pixels[x][y].r * 255, pixels[x][y].g * 255, pixels[x][y].b * 255);
		}
	}
	image.save_image("out.bmp");

}

void drawSquare(double a)
{
	glBegin(GL_QUADS); {

		glVertex3f(0, 0, 0);
		glVertex3f(0, a, 0);		
		glVertex3f(a, a, 0);
		glVertex3f(a, 0, 0);
	}glEnd();
}

void drawPyramid(double width, double height){
	glBegin(GL_QUADS); {
		glVertex3f( 0, 0, 0);
		glVertex3f(0, width, 0);
		glVertex3f(width, width, 0);
		glVertex3f(width, 0, 0);
	}glEnd();
	glBegin(GL_TRIANGLES); {
		glVertex3f(width / 2, width / 2, height);
		glVertex3f(0, 0, 0);
		glVertex3f(width, 0, 0);

		glVertex3f(width / 2, width / 2, height);
		glVertex3f(0, width, 0);
		glVertex3f(width, width, 0);
		
		glVertex3f(width / 2, width / 2, height);
		glVertex3f(width, width, 0);
		glVertex3f(width, 0, 0);

		glVertex3f(width / 2, width / 2, height);
		glVertex3f(width, 0, 0);
		glVertex3f(0, 0, 0);

	}glEnd();


}


void drawSphere(double radius, int slices, int stacks)
{
	struct point points[100][100];
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
	for (i = 0; i<stacks; i++)
	{
		//glColor3f((double)i / (double)stacks, (double)i / (double)stacks, (double)i / (double)stacks);
		for (j = 0; j<slices; j++)
		{
			glBegin(GL_QUADS); {
				//upper hemisphere
				glVertex3f(points[i][j].x, points[i][j].y, points[i][j].z);
				glVertex3f(points[i][j + 1].x, points[i][j + 1].y, points[i][j + 1].z);
				glVertex3f(points[i + 1][j + 1].x, points[i + 1][j + 1].y, points[i + 1][j + 1].z);
				glVertex3f(points[i + 1][j].x, points[i + 1][j].y, points[i + 1][j].z);
				//lower hemisphere
				glVertex3f(points[i][j].x, points[i][j].y, -points[i][j].z);
				glVertex3f(points[i][j + 1].x, points[i][j + 1].y, -points[i][j + 1].z);
				glVertex3f(points[i + 1][j + 1].x, points[i + 1][j + 1].y, -points[i + 1][j + 1].z);
				glVertex3f(points[i + 1][j].x, points[i + 1][j].y, -points[i + 1][j].z);
			}glEnd();
		}
	}
}
void drawCheckerboard(double a){
	int k = ceil(log2(far_/a));
	gridDim = (int)pow(2,k);
	double startX = (gridDim / 2)*a;
	int col;
	for (int i = 0; i < gridDim; i++){
		double startY = (gridDim / 2)*a;
		if (i % 2)col = 0;
		else col =1;

		for (int j = 0; j < gridDim; j++){
			glPushMatrix(); 
			{
				Color color;
				if (col){
					color.r = 1;
					color.g = 1;
					color.b = 1;
				}
				else{
					color.r = 0;
					color.g = 0;
					color.b = 0;
				}

				glColor3f(color.r * 255, color.g * 255, color.b * 255);
				glTranslatef(startX, startY, 0);
				gridColMap[make_pair(i,j)] = color;
				drawSquare(a);
			}
			glPopMatrix();
			
			startY -= a;
			col = 1 - col;
		}
		startX -= a;
	}
	
}
double value_(double x, double y, double z)
{
	return sqrt(x*x + y*y + z*z);
}

void keyboardListener(unsigned char key, int x, int y){
	switch (key){
		Look l_temp, r_temp, u_temp;
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
	case '0':
		render_image();
		break;
	case ' ':
		setTextureBuffer();
		cout << "texture loaded" << endl;
		texture_rend = 1 - texture_rend;
		break;
	default:
		break;
	}

}

void specialKeyListener(int key, int x, int y){
	switch (key){
	case GLUT_KEY_DOWN:		//down arrow key
		pos.x -= l.x * 2;
		pos.y -= l.y * 2;
		pos.z -= l.z * 2;
		break;
	case GLUT_KEY_UP:		// up arrow key
		pos.x += l.x * 2;
		pos.y += l.y * 2;
		pos.z += l.z * 2;
		break;

	case GLUT_KEY_RIGHT:
		pos.x += r.x * 2;
		pos.y += r.y * 2;
		pos.z += r.z * 2;
		break;
	case GLUT_KEY_LEFT:
		pos.x -= r.x * 2;
		pos.y -= r.y * 2;
		pos.z -= r.z * 2;
		break;

	case GLUT_KEY_PAGE_UP:
		pos.x += u.x * 2;
		pos.y += u.y * 2;
		pos.z += u.z * 2;
		break;
	case GLUT_KEY_PAGE_DOWN:
		pos.x -= u.x * 2;
		pos.y -= u.y * 2;
		pos.z -= u.z * 2;
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
		if (state == GLUT_DOWN){
		}

		break;

	case GLUT_RIGHT_BUTTON:
		if (state == GLUT_DOWN){		// 2 times?? in ONE click? -- solution is checking DOWN or UP
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
	drawCheckerboard(infchecker.a);

	
	for (int i = 0; i < spheres.size(); i++){
		glPushMatrix();
		{
			glColor3f(spheres[i].color.r, spheres[i].color.g, spheres[i].color.b);
			glTranslatef(spheres[i].center.x, spheres[i].center.y, spheres[i].center.z);
			drawSphere(spheres[i].radius,40,40);
		}
		glPopMatrix();
	}
	for (int i = 0; i < pyramids.size(); i++){
		glPushMatrix();
		{
			glColor3f(pyramids[i].color.r, pyramids[i].color.g, pyramids[i].color.b);
			glTranslatef(pyramids[i].lowest.x, pyramids[i].lowest.y, pyramids[i].lowest.z);
			drawPyramid(pyramids[i].width, pyramids[i].height);
		}
		glPopMatrix();
	}

	for (int i = 0; i < normlights.size(); i++){
		glPushMatrix();
		{
			glColor3f(255, 255, 255);
			glTranslatef(normlights[i].position.x, normlights[i].position.y, normlights[i].position.z);
			drawSphere(5, 10, 10);
		}
		glPopMatrix();
	}

	for (int i = 0; i < spotlights.size(); i++){
		glPushMatrix();
		{
			glColor3f(255, 255, 255);
			glTranslatef(spotlights[i].position.x, spotlights[i].position.y, spotlights[i].position.z);
			drawSphere(5, 10, 10);
		}
		glPopMatrix();
	}

	//ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
	glutSwapBuffers();
}

void animate(){
	//codes for any changes in Models, Camera
	glutPostRedisplay();
}

void init(){
	//codes for initialization

	drawaxes = 1;
	cameraAngle = 1;
	pos.x = 200;
	pos.y = 200;
	pos.z = 80;

	u.x = 0;
	u.y = 0;
	u.z = 1;

	r.x = -1 / sqrt(2.0);
	r.y = 1 / sqrt(2.0);
	r.z = 0;

	l.x = -1/sqrt(2.0);
	l.y = -1 / sqrt(2.0);
	l.z = 0;
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
	//gluPerspective(120, 1, 1, 1000.0);
	gluPerspective(fovY, aspRatio, near_, far_);
	//field of view in the Y (vertically)
	//aspect ratio that determines the field of view in the X direction (horizontally)
	//near distance
	//far distance
}

void inputFile(){
	ifstream input;
	input.open("description.txt");
	input >> near_ >> far_;
	input >> fovY;
	input >> aspRatio;
	input >>  recLevel;
	input >> nPixels;
	input >> infchecker.a;
	input >> infchecker.ambient >> infchecker.diffuse >> infchecker.reflection;
	input >> nObjets;
	string command;
	for (int i = 0; i < nObjets; i++){
		input >> command;
		if (command == "sphere"){
			sphere tempSph;
			input >> tempSph.center.x >> tempSph.center.y >> tempSph.center.z;
			input >> tempSph.radius;
			input >> tempSph.color.r >> tempSph.color.g >> tempSph.color.b;
			input >> tempSph.ambient >> tempSph.diffuse >> tempSph.specular>>tempSph.reflection;
			input >> tempSph.shininess;
			spheres.push_back(tempSph);
		}
		if (command == "pyramid"){
			pyramid tempPyr;
			input >> tempPyr.lowest.x >> tempPyr.lowest.y >> tempPyr.lowest.z;
			input >> tempPyr.width >> tempPyr.height;
			input >> tempPyr.color.r >> tempPyr.color.g >> tempPyr.color.b;
			input >> tempPyr.ambient >> tempPyr.diffuse >> tempPyr.specular >> tempPyr.reflection;
			input >> tempPyr.shininess;
			pyramids.push_back(tempPyr);
		}
	}
	input >> nNorms;
	for (int i = 0; i < nNorms; i++){
		normlight tempNorm;
		input >> tempNorm.position.x >> tempNorm.position.y >> tempNorm.position.z >> tempNorm.falloff;
		normlights.push_back(tempNorm);
	}
	input >> nSpots;
	for (int i = 0; i < nSpots; i++){
		spotlight tempSpot;
		input >> tempSpot.position.x >> tempSpot.position.y >> tempSpot.position.z >> tempSpot.falloff;
		input >> tempSpot.look.x >> tempSpot.look.y >> tempSpot.look.z;
		input >> tempSpot.cutoff;
		spotlights.push_back(tempSpot);
	}
	input.close();
}

int main(int argc, char **argv){
	glutInit(&argc, argv);
	glutInitWindowSize(500, 500);
	glutInitWindowPosition(0, 0);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	//Depth, Double buffer, RGB color

	glutCreateWindow("My OpenGL Program");
	inputFile();
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