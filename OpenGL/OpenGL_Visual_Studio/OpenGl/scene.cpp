#include <iostream>
#include <fstream>
#include <cassert>
#include <iomanip>
#include <cmath>
#include <stack>
#include <queue>
#include <vector>
#include "bitmap_image.hpp"
using namespace std;

#define pi (2*acos(0.0))
#define epsilon (1.0e-6)

class homogeneous_point {
public:
	double x, y, z, w;

	// set the three coordinates, set w to 1

	homogeneous_point(double x, double y, double z) {
		this->x = x;
		this->y = y;
		this->z = z;
		this->w = 1;
	}

	/*
	default constructor. does nothing. allows declarations like below:
	matrix m;
	therefore, usage is dangerous
	*/
	homogeneous_point() {
		this->w = 1;
	}

	// constructs a homogeneous point with given coordinates. forces w to be 1.0
	// if w is zero, raises error

	homogeneous_point(double x, double y, double z, double w) {
		assert(w != 0);
		this->x = x / w;
		this->y = y / w;
		this->z = z / w;
		this->w = 1;
	}

	// adds two points. returns a point forcing w to be 1.0

	homogeneous_point operator+(const homogeneous_point& point) {
		double x = this->x + point.x;
		double y = this->y + point.y;
		double z = this->z + point.z;
		double w = this->w + point.w;
		homogeneous_point p(x, y, z, w);
		return p;
	}

	// subtracts one point from another. returns a point forcing w to be 1.0

	homogeneous_point operator-(const homogeneous_point& point) {
		double x = this->x - point.x;
		double y = this->y - point.y;
		double z = this->z - point.z;
		double w = this->w - point.w;
		homogeneous_point p(x, y, z, w);
	}

	// Print the coordinates of a point. exists for testing purpose.

	void print() {
		cout << "Point: " << endl;
		cout << x << " " << y << " " << z << " " << w << endl;
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
		cout << "Vector" << endl;
		cout << x << " " << y << " " << z << endl;
	}
};

/*
The matrices are forced to be 4x4. This is because in this assignment, we will deal with points in triangles.
Maximum # of points that we will deal with at once is 3. And all the standard matrices are 4x4 (i.e. scale, translation, rotation etc.)
*/
class matrix {
public:
	double values[4][4];
	int num_rows, num_cols;

	// only set the number of rows and cols

	matrix(int rows, int cols) {
		assert(rows <= 4 && cols <= 4);
		num_rows = rows;
		num_cols = cols;
	}

	// prepare an nxn square matrix

	matrix(int n) {
		assert(n <= 4);
		num_rows = num_cols = n;
	}

	// prepare and return an identity matrix of size nxn

	static matrix make_identity(int n) {
		assert(n <= 4);
		matrix m(n);
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				if (i == j)
					m.values[i][j] = 1;
				else
					m.values[i][j] = 0;
			}
		}
		return m;
	}

	// print the matrix. exists for testing purposes

	void print() {
		cout << "Matrix:" << endl;
		for (int i = 0; i < num_rows; i++) {
			for (int j = 0; j < num_cols; j++) {
				cout << values[i][j] << "\t";
			}
			cout << endl;
		}
	}

	// add the two matrices. Raise error if dimension mismatches

	matrix operator+(const matrix& m) {
		assert(this->num_rows == m.num_rows);
		assert(this->num_cols == m.num_cols);

		matrix m1(num_rows, num_cols);
		for (int i = 0; i < num_rows; i++) {
			for (int j = 0; j < num_cols; j++) {
				m1.values[i][j] = values[i][j] + m.values[i][j];
			}
		}
		return m1;
	}

	// subtract a matrix from another. raise error if dimension mismatches

	matrix operator-(const matrix& m) {
		assert(this->num_rows == m.num_rows);
		assert(this->num_cols == m.num_cols);

		matrix m1(num_rows, num_cols);
		for (int i = 0; i < num_rows; i++) {
			for (int j = 0; j < num_cols; j++) {
				m1.values[i][j] = values[i][j] - m.values[i][j];
			}
		}
		return m1;
	}

	// multiply two matrices. allows statements like m1 = m2 * m3; raises error is dimension mismatches

	matrix operator*(const matrix& m) {
		assert(this->num_cols == m.num_rows);
		matrix m1(this->num_rows, m.num_cols);

		for (int i = 0; i < m1.num_rows; i++) {
			for (int j = 0; j < m1.num_cols; j++) {
				double val = 0;
				for (int k = 0; k < this->num_cols; k++) {
					val += this->values[i][k] * m.values[k][j];
				}
				m1.values[i][j] = val;
			}
		}
		return m1;
	}

	// multiply a matrix with a constant

	matrix operator*(double m) {
		matrix m1(this->num_rows, this->num_cols);
		for (int i = 0; i < num_rows; i++) {
			for (int j = 0; j < num_cols; j++) {
				m1.values[i][j] = m * this->values[i][j];
			}
		}
		return m1;
	}

	// multiply a 4x4 matrix with a homogeneous point and return the resulting point.
	// usage: homogeneous_point p = m * p1;
	// here, m is a 4x4 matrix, intended to be the transformation matrix
	// p1 is the point on which the transformation is being made
	// p is the resulting homogeneous point

	homogeneous_point operator*(const homogeneous_point& p) {
		assert(this->num_rows == this->num_cols && this->num_rows == 4);

		matrix m(4, 1);
		m.values[0][0] = p.x;
		m.values[1][0] = p.y;
		m.values[2][0] = p.z;
		m.values[3][0] = p.w;

		matrix m1 = (*this) * m;
		homogeneous_point p1(m1.values[0][0], m1.values[1][0], m1.values[2][0], m1.values[3][0]);
		return p1;
	}

	// return the transpose of a matrix

	matrix transpose() {
		matrix m(num_cols, num_rows);
		for (int i = 0; i < num_rows; i++) {
			for (int j = 0; j < num_cols; j++) {
				m.values[j][i] = values[i][j];
			}
		}
		return m;
	}

};

/*
A simple class to hold the color components, r, g, b of a certain shade.
*/
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


double eye_x, eye_y, eye_z;
double look_x, look_y, look_z;
double up_x, up_y, up_z;
double fov_x, fov_y, aspectRatio, near, far;
color backgroud;
int screen_x, screen_y;
vector<color>colors;
int n_inp = 0;
int n_triangle;

void scan_convert() {
	ifstream stage3;
	stage3.open("stage3.txt");

	color** pixels = new color*[screen_x];
	double** zs = new double*[screen_x];
	for (int i = 0; i < screen_x; i++) {
		pixels[i] = new color[screen_y];
		for (int j = 0; j < screen_y; j++) {
			pixels[i][j] = backgroud;
		}
		zs[i] = new double[screen_y];
		for (int j = 0; j < screen_y; j++) {
			zs[i][j] = +20; // a very large value intended as +INFINITY
		}
	}

	// perform scan conversion, populate the 2D array pixels
	// the array zs is the z-buffer.
	for (int i = 0; i < n_inp; i += 3) {
		homogeneous_point points[3];
		color cur_color = colors[i / 3];
		stage3 >> points[0].x >> points[0].y >> points[0].z;
		stage3 >> points[1].x >> points[1].y >> points[1].z;
		stage3 >> points[2].x >> points[2].y >> points[2].z;
		homogeneous_point sorted[3];

		if (points[0].y > points[1].y && points[0].y > points[2].y) {
			sorted[0] = points[0];
			if (points[1].y > points[2].y) {
				sorted[1] = points[1];
				sorted[2] = points[2];
			}
			else {
				sorted[1] = points[2];
				sorted[2] = points[1];
			}
		}
		else if (points[1].y > points[0].y && points[1].y > points[2].y) {
			sorted[0] = points[1];
			if (points[0].y > points[2].y) {
				sorted[1] = points[0];
				sorted[2] = points[2];
			}
			else {
				sorted[1] = points[2];
				sorted[2] = points[0];
			}
		}
		else if (points[2].y > points[0].y && points[2].y > points[1].y) {
			sorted[0] = points[2];
			if (points[0].y > points[1].y) {
				sorted[1] = points[0];
				sorted[2] = points[1];
			}
			else {
				sorted[1] = points[1];
				sorted[2] = points[0];
			}
		}

		double delta_y = 2.0 / screen_y;
		double delta_x = 2.0 / screen_x;
		double y = sorted[0].y > 1.0 ? 1.0 : sorted[0].y;
		while (y >= sorted[1].y && y >= -1) {
			double xa = sorted[0].x - (sorted[0].x - sorted[1].x)*(sorted[0].y - y) / (sorted[0].y - sorted[1].y);
			double xb = sorted[0].x - (sorted[0].x - sorted[2].x)*(sorted[0].y - y) / (sorted[0].y - sorted[2].y);
			double za = (sorted[0].z - (sorted[0].z - sorted[1].z)*(sorted[0].y - y) / (sorted[0].y - sorted[1].y));
			double zb = (sorted[0].z - (sorted[0].z - sorted[2].z)*(sorted[0].y - y) / (sorted[0].y - sorted[2].y));

			if (xb < xa) {
				double temp;
				temp = xa;
				xa = xb;
				xb = temp;
				temp = za;
				za = zb;
				zb = temp;

			}

			if (xb >= -1 && xa <= 1) {
				if (xa < -1)xa = -1;
				if (xb > 1)xb = 1;
				for (double x = xa; x <= xb; x += delta_x) {
					double zp = zb - (zb - za)*(xb - x) / (xb - xa);
					int p_x = (int)((x + 1) * screen_x / 2);
					int p_y = (int)(screen_y - (y + 1) * screen_y / 2);
					if (zs[p_x][p_y] - epsilon >= zp) {
						zs[p_x][p_y] = zp;
						pixels[p_x][p_y] = cur_color;
					}
				}
				y -= delta_y;
			}
		}

		while (y >= sorted[2].y && y >= -1) {
			double xa = sorted[1].x - (sorted[1].x - sorted[2].x)*(sorted[1].y - y) / (sorted[1].y - sorted[2].y);
			double xb = sorted[0].x - (sorted[0].x - sorted[2].x)*(sorted[0].y - y) / (sorted[0].y - sorted[2].y);
			double za = (sorted[1].z - (sorted[1].z - sorted[2].z)*(sorted[1].y - y) / (sorted[1].y - sorted[2].y));
			double zb = (sorted[0].z - (sorted[0].z - sorted[2].z)*(sorted[0].y - y) / (sorted[0].y - sorted[2].y));
			if (xb < xa) {
				double temp;
				temp = xa;
				xa = xb;
				xb = temp;
				temp = za;
				za = zb;
				zb = temp;
			}

			if (xb >= -1 && xa <= 1) {
				if (xa < -1)xa = -1;
				if (xb > 1) xb = 1;
				for (double x = xa; x <= xb; x += delta_x) {
					double zp = zb - (zb - za)*(xb - x) / (xb - xa);
					int p_x = (int)((x + 1) * screen_x / 2);
					int p_y = (int)(screen_y - (y + 1) * screen_y / 2);
					if (zs[p_x][p_y] - epsilon >= zp) {
						zs[p_x][p_y] = zp;
						pixels[p_x][p_y] = cur_color;
					}
				}
				y -= delta_y;
			}
		}



	}

	// the following code generates a bmp image. do not change this.
	bitmap_image image(screen_x, screen_y);
	for (int x = 0; x < screen_x; x++) {
		for (int y = 0; y < screen_y; y++) {
			image.set_pixel(x, y, pixels[x][y].r, pixels[x][y].g, pixels[x][y].b);
		}
	}
	image.save_image("out.bmp");

	// free the dynamically allocated memory

}

void stage3() {
	if (near == far) return;
	ifstream stage2;
	ofstream stage3;
	stage2.open("stage2.txt");
	stage3.open("stage3.txt");
	stage3 << std::fixed;
	stage3 << std::setprecision(7);

	// process input from stage2 and write to stage3

	double fov_x = fov_y * aspectRatio;
	double t = near * tan(fov_y * pi / 360);
	double r = near * tan(fov_x * pi / 360);
	matrix mat(4);
	matrix P = mat.make_identity(4);
	P.values[0][0] = near / r;
	P.values[1][1] = near / t;
	P.values[2][2] = -(far + near) / (far - near);
	P.values[3][2] = -1;
	P.values[2][3] = -(2 * far * near) / (far - near);
	P.values[3][3] = 0;

	vector<color>colors_new;
	for (int t = 0; t < n_inp; t += 3) {
		color cur_color = colors[t / 3];
		homogeneous_point points[3];
		stage2 >> points[0].x >> points[0].y >> points[0].z;
		stage2 >> points[1].x >> points[1].y >> points[1].z;
		stage2 >> points[2].x >> points[2].y >> points[2].z;

		vector<homogeneous_point>P1;
		vector<homogeneous_point>P2;
		for (int i = 0; i < 3; i++) {
			int pos = (2 + i) % 3;
			homogeneous_point s = points[pos];
			homogeneous_point p = points[(pos + 1) % 3];
			if (s.z < -near && p.z < -near) {
				P1.push_back(p);
			}
			else if (s.z >= -near && p.z < -near) {
				//				cout << "here1" << endl;
				double intersect_x = s.x + (p.x - s.x)*((-near - s.z) / (p.z - s.z));
				double intersect_y = s.y + (p.y - s.y)*((-near - s.z) / (p.z - s.z));
				P1.push_back(homogeneous_point(intersect_x, intersect_y, -near));
				P1.push_back(p);
			}
			else if (s.z < -near && p.z >= -near) {
				//				cout << "here2" << endl;
				double intersect_x = s.x + (p.x - s.x)*((-near - s.z) / (p.z - s.z));
				double intersect_y = s.y + (p.y - s.y)*((-near - s.z) / (p.z - s.z));
				P1.push_back(homogeneous_point(intersect_x, intersect_y, -near));
			}
		}

		for (int i = 0; i < P1.size(); i++) {
			int pos = (P1.size() - 1 + i) % P1.size();
			homogeneous_point s = P1[pos];
			homogeneous_point p = P1[(pos + 1) % P1.size()];
			if (s.z > -far && p.z > -far) {
				P2.push_back(p);
			}
			else if (s.z <= -far && p.z > -far) {
				//				cout << "here3" << endl;
				double intersect_x = s.x + (p.x - s.x)*((-far - s.z) / (p.z - s.z));
				double intersect_y = s.y + (p.y - s.y)*((-far - s.z) / (p.z - s.z));
				P2.push_back(homogeneous_point(intersect_x, intersect_y, -far));
				P2.push_back(p);
			}
			else if (s.z > -far && p.z <= -far) {
				//				cout << "here4" << endl;
				double intersect_x = s.x + (p.x - s.x)*((-far - s.z) / (p.z - s.z));
				double intersect_y = s.y + (p.y - s.y)*((-far - s.z) / (p.z - s.z));
				P2.push_back(homogeneous_point(intersect_x, intersect_y, -far));
			}
		}

		if (P2.size() != 0) {
			for (int i = 0; i < P2.size() - 2; i++) {
				homogeneous_point p_new;
				p_new = P * P2[0];
				stage3 << p_new.x << " " << p_new.y << " " << p_new.z << endl;
				p_new = P * P2[i + 1];
				stage3 << p_new.x << " " << p_new.y << " " << p_new.z << endl;
				p_new = P * P2[i + 2];
				stage3 << p_new.x << " " << p_new.y << " " << p_new.z << endl;
				stage3 << endl;
				colors_new.push_back(cur_color);
			}
		}
	}

	colors = colors_new;
	n_inp = colors_new.size() * 3;
	stage3.close();
	stage2.close();

}

void stage2() {
	ifstream stage1;
	ofstream stage2;
	stage1.open("stage1.txt");
	stage2.open("stage2.txt");
	stage2 << std::fixed;
	stage2 << std::setprecision(7);
	// collect input from stage1 and process, write output to stage2
	matrix mat(4);
	Vector l(look_x - eye_x, look_y - eye_y, look_z - eye_z);
	l.normalize();
	Vector r = l.cross(l, Vector(up_x, up_y, up_z));
	r.normalize();
	Vector u = l.cross(r, l);
	matrix T = mat.make_identity(4);
	T.values[0][3] = -eye_x;
	T.values[1][3] = -eye_y;
	T.values[2][3] = -eye_z;

	matrix R = mat.make_identity(4);
	R.values[0][0] = r.x;
	R.values[0][1] = r.y;
	R.values[0][2] = r.z;

	R.values[1][0] = u.x;
	R.values[1][1] = u.y;
	R.values[1][2] = u.z;

	R.values[2][0] = -l.x;
	R.values[2][1] = -l.y;
	R.values[2][2] = -l.z;

	matrix V = R*T;
	for (int i = 0; i < n_inp; i++) {
		homogeneous_point p(0, 0, 0, 1);
		stage1 >> p.x >> p.y >> p.z;
		homogeneous_point p_new = V*p;

		stage2 << p_new.x << " " << p_new.y << " " << p_new.z << endl;
		if ((i + 1) % 3 == 0)stage2 << endl;
	}

	stage1.close();
	stage2.close();

}

void stage1() {
	ifstream scene;
	ofstream stage1;
	scene.open("scene.txt");
	stage1.open("stage1.txt");
	stage1 << std::fixed;
	stage1 << std::setprecision(7);

	string command;

	scene >> eye_x >> eye_y >> eye_z;
	scene >> look_x >> look_y >> look_z;
	scene >> up_x >> up_y >> up_z;
	scene >> fov_y >> aspectRatio >> near >> far;
	scene >> screen_x >> screen_y;
	scene >> backgroud.r >> backgroud.g >> backgroud.b;

	// take other commands as input from scene in a loop
	// process accordingly
	// write to stage1
	stack<matrix> S;
	stack<int>scope;
	matrix mat(4);
	int curr_scope = 0;
	S.push(mat.make_identity(4));
	while (true) {
		scene >> command;
		if (command == "triangle") {
			//cout << "triangle" << endl;
			homogeneous_point point(0, 0, 0, 1);
			color col;
			for (int i = 0; i < 3; i++) {
				scene >> point.x >> point.y >> point.z;
				homogeneous_point p = S.top() * point;
				stage1 << p.x << " " << p.y << " " << p.z << endl;
				n_inp++;
			}
			stage1 << endl;
			scene >> col.r >> col.g >> col.b;
			colors.push_back(col);
		}
		else if (command == "translate") {
			//cout << "translate" << endl;
			matrix T = mat.make_identity(4);
			scene >> T.values[0][3] >> T.values[1][3] >> T.values[2][3];
			S.push(S.top() * T);
			curr_scope++;
		}
		else if (command == "scale") {
			//cout << "scale" << endl;
			matrix T = mat.make_identity(4);
			scene >> T.values[0][0] >> T.values[1][1] >> T.values[2][2];
			S.push(S.top() * T);
			curr_scope++;
		}
		else if (command == "rotate") {
			//cout << "rotate" << endl;
			matrix T = mat.make_identity(4);
			Vector a;
			double angle;
			scene >> angle >> a.x >> a.y >> a.z;
			a.normalize();
			Vector c1 = a.rotate(Vector(1, 0, 0), a, angle);
			Vector c2 = a.rotate(Vector(0, 1, 0), a, angle);
			Vector c3 = a.rotate(Vector(0, 0, 1), a, angle);
			T.values[0][0] = c1.x;
			T.values[0][1] = c2.x;
			T.values[0][2] = c3.x;

			T.values[1][0] = c1.y;
			T.values[1][1] = c2.y;
			T.values[1][2] = c3.y;

			T.values[2][0] = c1.z;
			T.values[2][1] = c2.z;
			T.values[2][2] = c3.z;
			S.push(S.top() * T);
			curr_scope++;
		}
		else if (command == "push") {
			//cout << "push" << endl;
			scope.push(curr_scope);
			curr_scope = 0;
		}
		else if (command == "pop") {
			//cout << "pop" << endl;
			for (int i = 0; i < curr_scope; i++)S.pop();
			curr_scope = scope.top();
			scope.pop();
		}
		else if (command == "end") {
			//cout << "end" << endl;
			break;
		}
	}

	scene.close();
	stage1.close();

}

int main() {
	cout << std::fixed;
	cout << std::setprecision(4);

	stage1();
	stage2();
	stage3();
	scan_convert();
	//	while (true);
	return 0;
}