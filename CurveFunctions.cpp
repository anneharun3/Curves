#include "CurveFunctions.h"
using std::vector;
#include <Eigen/LU>
using Eigen::MatrixXd;
using Eigen::Matrix4d;
using Eigen::Vector4d;
#include <cassert>
#include <cmath>
#include <iostream>
#include "jsassert.h"
#undef assert
#define assert(cond) jsAssert(cond)
// Call these to raise a dialog box or log to the javascript console for debugging.
// NOTE: You can pass either a const char* or an std::string.
// void jsAlert( const std::string& msg );
// void jsLog( const std::string& msg );
// void jsWarn( const std::string& msg );
// void jsError( const std::string& msg );
namespace Curve{
// Bezier helper functions.
namespace{
// You may add helper functions here.
}
// Evaluate a cubic Bezier curve at location 't'.
Point EvaluateCubicBezierCurve( const Point& p0, const Point& p1, const Point& p2, const Point& p3, const real_t t, EvaluateCubicBezierCurveApproach approach ){
    if( BernsteinApproach == approach ) return EvaluateCubicBezierCurveBernstein( p0, p1, p2, p3, t );
    else if( MatrixApproach == approach ) return EvaluateCubicBezierCurveMatrix( p0, p1, p2, p3, t );
    else if( CasteljauApproach == approach ) return EvaluateCubicBezierCurveCasteljau( p0, p1, p2, p3, t );
    else {
        assert( !"Unknown EvaluateCubicBezierCurveApproach" );
        return Point(-31337,-31337);
    }
}
Point EvaluateCubicBezierCurveBernstein( const Point& p0, const Point& p1, const Point& p2, const Point& p3, const real_t t ){
    // ADD YOUR CODE HERE
	const Point& k=(pow((1-t),3)*p0+(3*t*pow((1-t),2)*p1)+(3*pow(t,2)*(1-t)*p2)+(pow(t,3)*p3));
    return k;
}
Point EvaluateCubicBezierCurveMatrix( const Point& p0, const Point& p1, const Point& p2, const Point& p3, const real_t t ){
    // ADD YOUR CODE HERE	
	Point l,m,n,o,p;
	l = (-p0+3*p1-3*p2+p3);
	m = (3*p0-6*p1+3*p2);
	n = (-3*p0+3*p1);
	o = p0;
	p = l*pow(t,3) +m*pow(t,2) +n*pow(t,1) +o ;
    return p;
    }
Point EvaluateCubicBezierCurveCasteljau( const Point& p0, const Point& p1, const Point& p2, const Point& p3, const real_t t ){
    // ADD YOUR CODE HERE
	const Point& a=((1-t)*((1-t)*((1-t)*p0+t*p1)+t*((1-t)*p1+t*p2))+t*((1-t)*((1-t)*p1+t*p2)+t*((1-t)*p2+t*p3)));
    return a;
}
// Evaluate a cubic Bezier spline with control points 'controlPoints' arranged
//     on_curve ( off_curve off_curve on_curve )+
// at positive integer 'samplesPerCurve' locations along each curve.
// Returns the sampled points.
std::vector< Point > EvaluateCubicBezierSpline( const std::vector< Point >& controlPoints, const int samplesPerCurve, EvaluateCubicBezierCurveApproach approach ){
    assert( controlPoints.size() >= 4 );
    assert( samplesPerCurve > 0 );
    // ADD YOUR CODE HERE
    std::vector< Point > result;
	int l=controlPoints.size();
    real_t a=0.0;
	for(int i=0;i<l-3;i=i+3){
	a=0.0;
	for(int j=0;j<samplesPerCurve;j++){
	while(a<=1.0){
	result.push_back(EvaluateCubicBezierCurve(controlPoints[i],controlPoints[i+1],controlPoints[i+2],controlPoints[i+3],a,approach));
	a=a+(1.00/samplesPerCurve);
	}
	}
	}
    return result;
}

/// ======================================================================================
// Evaluate a cubic Hermite curve at location 't'.
Point EvaluateCubicHermiteCurve( const Point& p0, const Point& dp0, const Point& p1, const Point& dp1, const real_t t ){
    // ADD YOUR CODE HERE
    const Point& a=p0*(2*pow(t,3)-3*pow(t,2)+1)+p1*(-2*pow(t,3)+3*pow(t,2))+dp0*(pow(t,3)-2*pow(t,2)+t)+dp1*(pow(t,3)-pow(t,2));
    return a;
}
// Evaluate a cubic Hermite spline with control points 'controlPoints' arranged:
// p0 derivative_at_p0 ( p1 derivative_at_p1 )+
// at positive integer 'samplesPerCurve' locations along each curve.
// Upon return, 'curvePointsOut' is cleared and replaced with the sampled points.
std::vector< Point > EvaluateCubicHermiteSpline(const std::vector< Point >& controlPoints, const int samplesPerCurve ){
    assert( controlPoints.size() >= 4 );
    assert( samplesPerCurve > 0 );
    std::vector< Point > result;
	int k=controlPoints.size();
    real_t a=0.0;
    for(int i=0;i<=k-2;i=i+2){
    a = 0.0;
    for(int j=0;j<samplesPerCurve;j++){
    result.push_back(EvaluateCubicHermiteCurve(controlPoints[i],controlPoints[i+1],controlPoints[i+2],controlPoints[i+3],a));
    a=a+(1.00/samplesPerCurve);
    }
    }
    return result;
}
// Given a cubic Hermite spline with control points 'controlPoints' arranged:
//     p0 derivative_at_p0 ( p1 derivative_at_p1 )+
// replaces the derivative entries with values that result in a C2 continuous
// Hermite spline.
// NOTE: 'controlPoints' is an input and output parameter. The derivative entries are replaced.
// Hint: To implement C2 continuous Hermite splines, you need to solve a system of equations.
//       This is best done with a matrix.
//       Then you solve a linear system Ac = p, where c are the unknown derivatives.
//       I have included the Eigen linear algebra package that can solve linear systems.
//       Below is an example showing how to use the linear system solver.
//  MatrixXd A(3,3);
//  MatrixXd c(3,1);
//  MatrixXd p(3,1);
//  A(0,0) = 1.0; A(0,1) = 0.0; A(0,2) = 0.0;
//  A(1,0) = 0.0; A(1,1) = 1.0; A(1,2) = 0.0;
//  A(2,0) = 0.0; A(2,1) = 0.0; A(2,2) = 1.0;
//  p(0,0) = 1.0; p(1,0) = 2.0; p(2,0) = 3.0;
//  c = A.fullPivLu().solve(p);
//  The result will be stored in c as follows: c(0,0) = 1.0; c(1,0) = 2.0; c(3,0) = 3.0, which satisfies Ac = p.
void CalculateHermiteSplineDerivativesForC2Continuity( std::vector< Point >& controlPoints ){
    if( controlPoints.size() < 4 ) return;   
    assert( controlPoints.size() >= 4 );
    assert( controlPoints.size() % 2 == 0 );
     int i,j,n =(controlPoints.size())/2.0;
     MatrixXd A(n,n);
     MatrixXd dp(n,2);
     MatrixXd p(n,2);
     for(i=0;i<n;i++)
      for(j=0;j<n;j++)
		 A(i,j)=0.0;
    A(0,0)=2;A(0,1)=1;
    A((n-1),(n-1))=1;A((n-1),(n-2))=2;
    for(i=1,j=1;i<=n-2;i++){
        A(i,j-1)=1;A(i,j)=4;A(i,j+1)=1;
          j++;  
      }
    p(0,0)=3*(controlPoints.at(2).x()-controlPoints.at(0).x());
    p(0,1)=3*(controlPoints.at(2).y()-controlPoints.at(0).y());
    for(i=1,j=0;i<=n-2;i++){
        p(i,0)=3*(controlPoints.at(j+4).x()-controlPoints.at(j).x());
        p(i,1)=3*(controlPoints.at(j+4).y()-controlPoints.at(j).y());
        j=j+2;
    } 
    p(n-1,0) = 3 * (controlPoints.at(2*n-2).x() - controlPoints.at(2*n-4).x());
    p(n-1,1) = 3 * (controlPoints.at(2*n-2).y() - controlPoints.at(2*n-4).y());
    dp = A.fullPivLu().solve(p);
    for(i=0,j=1;j<n;i++){
        controlPoints.at(j).x()= dp(i,0); 
        controlPoints.at(j).y()= dp(i,1) ;
        j = j+2;
    }  
	}
/// ======================================================================================

// Evaluate a Catmull-Rom Spline with control points 'controlPoints' arranged:
//     p0 p1 p2 ( p3 )+
// at positive integer 'samplesPerCurve' locations along each curve.
// Upon return, 'curvePointsOut' is cleared and replaced with the sampled points.
std::vector< Point > EvaluateCatmullRomSpline( const std::vector< Point >& controlPoints, const int samplesPerCurve, const real_t alpha ){
    assert( controlPoints.size() >= 4);
    assert( samplesPerCurve>0);
    std::vector< Point > result;
    real_t a=0.0;int i;
    for(i=0;i<samplesPerCurve;i++){
        result.push_back(EvaluateCatmullRomCurve(2*controlPoints[0]-controlPoints[1],controlPoints[0],controlPoints[1],controlPoints[2],a,alpha));
        a=a+(1/samplesPerCurve);
    }
	int n=controlPoints.size();
    for(i=0;i<=n-2;i=i+1){
        a = 0.0;
        for(int j=0;j<samplesPerCurve;j++){
            a=a+(1.00/samplesPerCurve);
            result.push_back(EvaluateCatmullRomCurve(controlPoints[i-1],controlPoints[i],controlPoints[i+1],controlPoints[i+2],a,alpha));            
        }      
    }
    a =0.0;
    // Plot the last curve
    for(i=0;i<samplesPerCurve;i++){
        result.push_back(EvaluateCatmullRomCurve(controlPoints[n-3],controlPoints[n-2],controlPoints[n-1],2*controlPoints[n-1] - controlPoints[n-2],a,alpha));
       a=a+(1/samplesPerCurve);
    }     
    return result;
}
// Evaluate a cubic Catmull-Rom Spline curve at location 't'.
Point EvaluateCatmullRomCurve( const Point& p0, const Point& p1, const Point& p2, const Point& p3, const real_t t, const real_t alpha ){
    // ADD YOUR CODE HERE
    MatrixXd a(4,4);MatrixXd ta(1,4);
	MatrixXd pr(1,4);MatrixXd px(4,1);MatrixXd py(4,1);
	a(0,0)=-0.5;a(0,1)=1.5;a(0,2)=-1.5;a(0,3)=0.5;a(1,0)=1.0;a(1,1)=-2.5;a(1,2)=2.0;a(1,3)=-0.5;a(2,0)=-0.5;a(2,1)=0.0;a(2,2)=0.5;a(2,3)=0.0;a(3,0)=0.0;a(3,1)=1.0;a(3,2)=0.0;a(3,3)=0.0;
	ta(0,0)=pow(t,3);ta(0,1)=pow(t,2);ta(0,2)=pow(t,1);ta(0,3)=1.0;
	pr=ta*a;   
	const Point& p=pr(0,0)*p0+pr(0,1)*p1+pr(0,2)*p2+pr(0,3)*p3;
    return p;
}
/// ======================================================================================
// B-Spline helper functions
namespace{
// You may add helper functions here.
}
// Evaluate a cubic B-Spline with control points 'controlPoints' arranged:
//     p0 p1 p2 ( p3 )+
// at positive integer 'samplesPerCurve' locations along each curve.
// Upon return, 'curvePointsOut' is cleared and replaced with the sampled points.
std::vector< Point > EvaluateCubicBSpline( const std::vector< Point >& controlPoints, const int samplesPerCurve ){
    assert( controlPoints.size() >= 4 );
    assert( samplesPerCurve > 0 );
    // ADD YOUR CODE HERE
    std::vector< Point > result;
    return result;
}
// Evaluate a cubic B-Spline curve at location 't'.
Point EvaluateCubicBSplineCurve( const Point& p0, const Point& p1, const Point& p2, const Point& p3, const real_t t ){
    // ADD YOUR CODE HERE
    return Point(0,0);
}

// Compute cubic BSpline control points that interpolate the given points.
std::vector< Point > ComputeBSplineFromInterpolatingPoints( const std::vector< Point >& interpPoints ){
    assert( interpPoints.size() >= 2 );
    // ADD YOUR CODE HERE
    std::vector< Point > controlPoints;
    return controlPoints;
}
// Given a sequence of cubic BSpline control points, returns the interpolating points
// that could have been used to create them via ComputeBSplineFromInterpolatingPoints().
std::vector< Point > ComputeInterpolatingPointsFromBSpline( const std::vector< Point >& controlPoints ){
    assert( controlPoints.size() >= 4 );
    const std::vector< Point >& C = controlPoints;
    std::vector< Point > result;
    // The interpolated points are at the start and end of each cubic B-Spline
    // (they are continuous).
    // So let's just sample the t=1 point on every cubic BSpline,
    // as well as the t=0 point of the first one.
    result.push_back( EvaluateCubicBSplineCurve( C[0], C[1], C[2], C[3], 0. ) );
    for( int i = 0; i+3 < C.size(); ++i ){
        result.push_back( EvaluateCubicBSplineCurve( C[i], C[i+1], C[i+2], C[i+3], 1. ) );
    }
    return result;
}
}
