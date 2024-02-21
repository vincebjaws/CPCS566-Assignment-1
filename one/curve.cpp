#include "curve.h"
#include "extra.h"
#ifdef WIN32
#include <windows.h>
#endif
#include <GL/gl.h>
using namespace std;

namespace
{
    // Approximately equal to.  We don't want to use == because of
    // precision issues with floating point.
    inline bool approx( const Vector3f& lhs, const Vector3f& rhs )
    {
        const float eps = 1e-8f;
        return ( lhs - rhs ).absSquared() < eps;
    }

    
}
    

Curve evalBezier( const vector< Vector3f >& P, unsigned steps )
{
    // Check
    if( P.size() < 4 || P.size() % 3 != 1 )
    {
        cerr << "evalBezier must be called with 3n+1 control points." << endl;
        exit( 0 );
    }

    // TODO:
    // You should implement this function so that it returns a Curve
    // (e.g., a vector< CurvePoint >).  The variable "steps" tells you
    // the number of points to generate on each piece of the spline.
    // At least, that's how the sample solution is implemented and how
    // the SWP files are written.  But you are free to interpret this
    // variable however you want, so long as you can control the
    // "resolution" of the discretized spline curve with it.

    // Make sure that this function computes all the appropriate
    // Vector3fs for each CurvePoint: V,T,N,B.
    // [NBT] should be unit and orthogonal.

    // Also note that you may assume that all Bezier curves that you
    // receive have G1 continuity.  Otherwise, the TNB will not be
    // be defined at points where this does not hold.

    cerr << "\t>>> evalBezier has been called with the following input:" << endl;

    cerr << "\t>>> Control points (type vector< Vector3f >): "<< endl;
    for( unsigned i = 0; i < P.size(); ++i )
    {
        cerr << "\t>>> " << P[i] << endl;
    }

    cerr << "\t>>> Steps (type steps): " << steps << endl;
    //cerr << "\t>>> Returning empty curve." << endl;

    // Right now this will just return this empty curve.
    
    int subCurves = 0;
    int controlPnts = 0;
    
    for (unsigned i=0; i< P.size(); i++){	//identify subcurves here
    	
    	controlPnts++;
    	
    	if (controlPnts == 4){		//every 4 control points is 1 subcurve
    		controlPnts = 1;	//set back to 1, end of 1 subcurve is the start of next subcurve
    		subCurves++;
    	}
    }
    
    // Preallocate a curve with # of subCurves, and # of steps in each subCurve
    Curve R(subCurves * steps + 1);
    
    //spline matrix, for 4 control points, monomial basis
    Matrix4f bernstein (1.0f, -3.0f, 3.0f, -1.0f,
    			0.0f, 3.0f, -6.0f, 3.0f,
    			0.0f, 0.0f, 3.0f, -3.0f,
    			0.0f, 0.0f, 0.0f, 1.0f);
    			
    //spline matrix derivative
    Matrix4f bernsteinderiv (0.0f, -3.0f, 6.0f, -3.0f,
    			     0.0f, 3.0f, -12.0f, 9.0f,
    			     0.0f, 0.0f, 6.0f, -9.0f,
    			     0.0f, 0.0f, 0.0f, 3.0f);
    			
    Vector3f binormal;
    Vector3f prev_binormal;
    
    int index = 0;
    for (unsigned i = 0; i < P.size() - 3; i += 3){
    	
    	binormal = (i == 0) ? Vector3f(0.0f, 0.0f, 1.0f) : prev_binormal;	//if binormal = (i == 0), then Vector3f(0.0f, 0.0f, 1.0f), else prev_binormal
    	
    	//control points matrix
    	Matrix4f controlPoints(P[i+0][0], P[i+1][0], P[i+2][0], P[i+3][0], 
			       P[i+0][1], P[i+1][1], P[i+2][1], P[i+3][1], 
			       P[i+0][2], P[i+1][2], P[i+2][2], P[i+3][2], 
			       0.f, 0.f, 0.f, 0.f);
			       
    	// Fill it in counterclockwise
   	for( unsigned j = 0; j <= steps; ++j ){

		// step
		float t = float( j ) / steps;
		
		//monomial basis matrix
    	Vector4f monomialbasis(1, t, t*t, t*t*t);

		// Initialize position
		// calculate V
		R[index].V = Vector3f( (controlPoints * bernstein * monomialbasis)[0], 
				       (controlPoints * bernstein * monomialbasis)[1],
				       (controlPoints * bernstein * monomialbasis)[2] );	//V = position
		
		// Tangent vector is first derivative
		R[index].T = Vector3f( (controlPoints * bernsteinderiv * monomialbasis)[0], 
				       (controlPoints * bernsteinderiv * monomialbasis)[1],
				       (controlPoints * bernsteinderiv * monomialbasis)[2] ).normalized();	//T = tangent
		
		// Normal vector is second derivative
		R[index].N = Vector3f::cross(binormal, R[index].T).normalized();	//N = vector 90 degrees w/ respect to T

		// Finally, binormal is facing up.
		R[index].B = Vector3f::cross(R[index].T, R[index].N).normalized();
		
		//keep track of binormal
		prev_binormal = R[index].B;
		
		index++;
  	}
    
    }
    			
    
    
    //calculate position: control points matrix * basis * monomial basis = position
    
    //calculate tangent: control points matrix * basis * monomial basis = tangent 
    
    //calculate normal: previous binormal and tangent cross product, then normalize
    
    //calculate binormal: tangent and normal cross product, then save as previous binormal for next recursive call
   
    
    return R;
}

Curve evalBspline( const vector< Vector3f >& P, unsigned steps )
{
    // Check
    if( P.size() < 4 )
    {
        cerr << "evalBspline must be called with 4 or more control points." << endl;
        exit( 0 );
    }

    // TODO:
    // It is suggested that you implement this function by changing
    // basis from B-spline to Bezier.  That way, you can just call
    // your evalBezier function
    
    cerr << "\t>>> evalBSpline has been called with the following input:" << endl;

    cerr << "\t>>> Control points (type vector< Vector3f >): "<< endl;
    
    for( unsigned i = 0; i < P.size(); ++i )
    {
        cerr << "\t>>> " << P[i] << endl;
    }

    cerr << "\t>>> Steps (type steps): " << steps << endl;

    Curve BSplineCurve;
    

    Matrix4f Matbernstein (1.0f, -3.0f, 3.0f, -1.0f,
    			0.0f, 3.0f, -6.0f, 3.0f,
    			0.0f, 0.0f, 3.0f, -3.0f,
    			0.0f, 0.0f, 0.0f, 1.0f);    

    // Bspline transform
    
    Matrix4f MatBSpline(1.0f/6, -3.0f/6, 3.0f/6, -1.0f/6,
    			4.0f/6, 0.0f/6, -6.0f/6, 3.0f/6,
    			1.0f/6, 3.0f/6, 3.0f/6, -3.0f/6,
    			0.0f/6, 0.0f/6, 0.0f/6, 1.0f/6);
    			
    Matrix4f Matbernstein_inverse = Matbernstein.inverse(false, 0.0);
    
    // Plot control points
    
    for (unsigned i = 0; i < P.size()-3; i++){
    	vector<Vector3f> newControlPoints;
    	
        Matrix4f controlPoints(P[i+0][0], P[i+1][0], P[i+2][0], P[i+3][0], 
		       P[i+0][1], P[i+1][1], P[i+2][1], P[i+3][1], 
		       P[i+0][2], P[i+1][2], P[i+2][2], P[i+3][2], 
		       0.f, 0.f, 0.f, 0.f);
		       
	//y(t) = Geomerty Basis * Spline basis * Monoial basis
		       
	Matrix4f Geometery_Matrix = controlPoints * MatBSpline * Matbernstein_inverse;
	
	// Transform vector4f from matrix to Vector3f
	
	for (unsigned iter = 0; iter < 4; iter++){
		Vector4f col = Geometery_Matrix.getCol(iter);
		Vector3f transform = Vector3f(col[0], col[1], col[2]);
		newControlPoints.push_back(transform);
	}
    	
    	// New Sub Curve from control points and steps
    	
    	Curve subR = evalBezier(newControlPoints,steps);

	// assign new curve point

    	for(unsigned j = 1; j <subR.size();j++)
    	{
    		CurvePoint P_t;
    		P_t.V = subR[j].V;
    		P_t.T = subR[j].T;
    		P_t.N = subR[j].N;
    		P_t.B = subR[j].B;
    		BSplineCurve.push_back(P_t);
    	}
    }
    return BSplineCurve;
}

Curve evalCircle( float radius, unsigned steps )
{
    // This is a sample function on how to properly initialize a Curve
    // (which is a vector< CurvePoint >).
    
    // Preallocate a curve with steps+1 CurvePoints
    Curve R( steps+1 );

    // Fill it in counterclockwise
    for( unsigned i = 0; i <= steps; ++i )
    {
        // step from 0 to 2pi
        float t = 2.0f * M_PI * float( i ) / steps;

        // Initialize position
        // We're pivoting counterclockwise around the y-axis
        R[i].V = radius * Vector3f( cos(t), sin(t), 0 );
        
        // Tangent vector is first derivative
        R[i].T = Vector3f( -sin(t), cos(t), 0 );
        
        // Normal vector is second derivative
        R[i].N = Vector3f( -cos(t), -sin(t), 0 );

        // Finally, binormal is facing up.
        R[i].B = Vector3f( 0, 0, 1 );
    }

    return R;
}

void drawCurve( const Curve& curve, float framesize )
{
    // Save current state of OpenGL
    glPushAttrib( GL_ALL_ATTRIB_BITS );

    // Setup for line drawing
    glDisable( GL_LIGHTING ); 
    glColor4f( 1, 1, 1, 1 );
    glLineWidth( 1 );
    
    // Draw curve
    glBegin( GL_LINE_STRIP );
    for( unsigned i = 0; i < curve.size(); ++i )
    {
        glVertex( curve[ i ].V );
    }
    glEnd();

    glLineWidth( 1 );

    // Draw coordinate frames if framesize nonzero
    if( framesize != 0.0f )
    {
        Matrix4f M;

        for( unsigned i = 0; i < curve.size(); ++i )
        {
            M.setCol( 0, Vector4f( curve[i].N, 0 ) );
            M.setCol( 1, Vector4f( curve[i].B, 0 ) );
            M.setCol( 2, Vector4f( curve[i].T, 0 ) );
            M.setCol( 3, Vector4f( curve[i].V, 1 ) );

            glPushMatrix();
            glMultMatrixf( M );
            glScaled( framesize, framesize, framesize );
            glBegin( GL_LINES );
            glColor3f( 1, 0, 0 ); glVertex3d( 0, 0, 0 ); glVertex3d( 1, 0, 0 );
            glColor3f( 0, 1, 0 ); glVertex3d( 0, 0, 0 ); glVertex3d( 0, 1, 0 );
            glColor3f( 0, 0, 1 ); glVertex3d( 0, 0, 0 ); glVertex3d( 0, 0, 1 );
            glEnd();
            glPopMatrix();
        }
    }
    
    // Pop state
    glPopAttrib();
}

