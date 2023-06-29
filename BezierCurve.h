#ifndef BEZIERCURVE_H
#define BEZIERCURVE_H

#include <igl/opengl/glfw/Viewer.h>
#include <Eigen/Dense> // This line includes the Eigen library.

struct BezierCurve
{
	BezierCurve(Eigen::MatrixXd P);

	// evaluate the coordinate of curvePoint with parameter t  
	Eigen::MatrixXd eval(double t);

	// display by libigl
	void show(igl::opengl::glfw::Viewer& viewer, double resolution = 0.01);
	
	Eigen::MatrixXd controlP;
};
#endif // !BEZIERCURVE_H





