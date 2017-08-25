#include "GL/glut.h"
#include <math.h>
#include <vector>


class CylLine
{
public:

  class CylPoint3D
  {
  public:
    float mX, mY, mZ;

	float getNorm()
	{
		return sqrt( mX * mX + mY * mY + mZ * mZ );
	}

	void normalize()
	{
		float norm = getNorm();

		if ( norm != 0 )
		{	
			mX /= norm;
			//mX = -mX;
			mY /= norm;
			//mY = -mY;
			mZ /= norm;
			//mZ = -mZ;
		}
	}

	float multVec( const float x, const float y, const float z )
	{
		return x * mX + y * mY + z * mZ;
	}
  };

  float                    mRadius;
  int                      mNumOfPointsOnCyl;
  std::vector<CylPoint3D>  mPointsOnCircle;
  std::vector<CylPoint3D>  mNormalsOnCircle;
  const float              mPi;
  const float			   mInvalidAngle;


  CylLine(float radius, float numOfPointsOnCyl)
    : mRadius(radius), mNumOfPointsOnCyl(numOfPointsOnCyl),
      mPi(3.14159265358979323846), mInvalidAngle(400.0f)
  {
    calcCircle();
  }

  void operator()(float firstX, float firstY, float firstZ,
             float secondX, float secondY, float secondZ)
  {
	  CylPoint3D my_start, my_xAx, my_zAx, my_yAx, my_norm;

	  my_start.mX = firstX;
	  my_start.mY = firstY;
	  my_start.mZ = firstZ;

	  // the normal direction
	  my_norm.mX = secondX - firstX;
	  my_norm.mY = secondY - firstY;
	  my_norm.mZ = secondZ - firstZ;

	  // the y axis
	  my_yAx.mX = secondX - firstX;
	  my_yAx.mY = secondY - firstY;
	  my_yAx.mZ = secondZ - firstZ;
	  my_yAx.normalize();
	  
	  if ( secondX - firstX == 0 && secondZ - firstZ == 0 )
	  {
		  my_xAx.mX = 1;
		  my_xAx.mY = 0;
		  my_xAx.mZ = 0;
		  
		  my_zAx.mX = 0;
		  my_zAx.mY = 0;
		  my_zAx.mZ = 1;
	  }
	  else
	  {
		 if ( secondX - firstX == 0 )
		 {
			  my_xAx.mX = -1;
			  my_xAx.mY = 0;
			  my_xAx.mZ = 0;
			  
			  my_zAx.mX = 0;
			  my_zAx.mY = 1;
			  my_zAx.mZ = 0;
		 }
		 else
		 {
			 if ( secondZ - firstZ == 0 )
			 {
				  my_xAx.mX = 0;
				  my_xAx.mY = 0;
				  my_xAx.mZ = 1;
				  
				  my_zAx.mX = 0;
				  my_zAx.mY = 1;
				  my_zAx.mZ = 0;
			 }
			 else
			 {
				 if ( my_yAx.mY != 0 )
				 {
					 calcOrthoVec(	my_xAx.mX, my_xAx.mY, my_xAx.mZ,	// result
						 my_yAx.mX, my_yAx.mY, my_yAx.mZ, 
						 my_yAx.mX, 0,		  my_yAx.mZ  );
					 
					 my_xAx.normalize();
					 
					 calcOrthoVec(	my_zAx.mX, my_zAx.mY, my_zAx.mZ, 
						 my_yAx.mX, my_yAx.mY, my_yAx.mZ, 
						 my_xAx.mX, my_xAx.mY, my_xAx.mZ  );
					 
					 my_zAx.normalize();
				 }
				 else
				 {
					my_xAx.mX = -my_yAx.mZ;
					my_xAx.mY = 0;
					my_xAx.mZ = my_yAx.mX;
					my_xAx.normalize();
					
					calcOrthoVec(	my_zAx.mX, my_zAx.mY, my_zAx.mZ, 
									my_yAx.mX, my_yAx.mY, my_yAx.mZ, 
									my_xAx.mX, my_xAx.mY, my_xAx.mZ  );
					
					my_zAx.normalize();
				 }
			 }
		 }
	  }

	  renderCylQuads(my_xAx, my_yAx, my_zAx, my_norm, my_start);
	  
	  // RENDER ENDPOINTS SPHERES
	  glPushMatrix();
	  {
		  glTranslatef(firstX, firstY, firstZ);
		  glutSolidSphere(mRadius, mNumOfPointsOnCyl, mNumOfPointsOnCyl);
		  
		  glTranslatef(secondX - firstX, secondY - firstY, secondZ - firstZ);
		  glutSolidSphere(mRadius, mNumOfPointsOnCyl, mNumOfPointsOnCyl);
	  }
	  glPopMatrix();
  }  

  void renderCylQuads(const CylPoint3D &my_xAx, const CylPoint3D &my_yAx, const CylPoint3D &my_zAx, const CylPoint3D &my_norm, const CylPoint3D &my_start)
  {
	CylPoint3D tmpPoint;

    int nextI;
    glBegin(GL_QUADS);
    {
      for (int i = 0; i < mNumOfPointsOnCyl; ++i)
      {
        nextI = (i + 1) %  mNumOfPointsOnCyl;
  
		tmpPoint.mX = mNormalsOnCircle[nextI].multVec( my_xAx.mX, my_yAx.mX, my_zAx.mX);
		tmpPoint.mY = mNormalsOnCircle[nextI].multVec( my_xAx.mY, my_yAx.mY, my_zAx.mY);
		tmpPoint.mZ = mNormalsOnCircle[nextI].multVec( my_xAx.mZ, my_yAx.mZ, my_zAx.mZ);
		glNormal3f(-tmpPoint.mX, -tmpPoint.mY, -tmpPoint.mZ);

		tmpPoint.mX = mPointsOnCircle[nextI].multVec( my_xAx.mX, my_yAx.mX, my_zAx.mX) + my_start.mX;
		tmpPoint.mY = mPointsOnCircle[nextI].multVec( my_xAx.mY, my_yAx.mY, my_zAx.mY) + my_start.mY;
		tmpPoint.mZ = mPointsOnCircle[nextI].multVec( my_xAx.mZ, my_yAx.mZ, my_zAx.mZ) + my_start.mZ;
		glVertex3f(tmpPoint.mX, tmpPoint.mY, tmpPoint.mZ);
        
		tmpPoint.mX = mNormalsOnCircle[nextI].multVec( my_xAx.mX, my_yAx.mX, my_zAx.mX);
		tmpPoint.mY = mNormalsOnCircle[nextI].multVec( my_xAx.mY, my_yAx.mY, my_zAx.mY);
		tmpPoint.mZ = mNormalsOnCircle[nextI].multVec( my_xAx.mZ, my_yAx.mZ, my_zAx.mZ);
		glNormal3f(-tmpPoint.mX, -tmpPoint.mY, -tmpPoint.mZ);

		tmpPoint.mX = mPointsOnCircle[nextI].multVec( my_xAx.mX, my_yAx.mX, my_zAx.mX) + my_norm.mX + my_start.mX;
		tmpPoint.mY = mPointsOnCircle[nextI].multVec( my_xAx.mY, my_yAx.mY, my_zAx.mY) + my_norm.mY + my_start.mY;
		tmpPoint.mZ = mPointsOnCircle[nextI].multVec( my_xAx.mZ, my_yAx.mZ, my_zAx.mZ) + my_norm.mZ + my_start.mZ;
		glVertex3f(tmpPoint.mX, tmpPoint.mY, tmpPoint.mZ);

		tmpPoint.mX = mNormalsOnCircle[i].multVec( my_xAx.mX, my_yAx.mX, my_zAx.mX);
		tmpPoint.mY = mNormalsOnCircle[i].multVec( my_xAx.mY, my_yAx.mY, my_zAx.mY);
		tmpPoint.mZ = mNormalsOnCircle[i].multVec( my_xAx.mZ, my_yAx.mZ, my_zAx.mZ);
		glNormal3f(-tmpPoint.mX, -tmpPoint.mY, -tmpPoint.mZ);

		tmpPoint.mX = mPointsOnCircle[i].multVec( my_xAx.mX, my_yAx.mX, my_zAx.mX) + my_norm.mX + my_start.mX;
		tmpPoint.mY = mPointsOnCircle[i].multVec( my_xAx.mY, my_yAx.mY, my_zAx.mY) + my_norm.mY + my_start.mY;
		tmpPoint.mZ = mPointsOnCircle[i].multVec( my_xAx.mZ, my_yAx.mZ, my_zAx.mZ) + my_norm.mZ + my_start.mZ;
		glVertex3f(tmpPoint.mX, tmpPoint.mY, tmpPoint.mZ);

		tmpPoint.mX = mNormalsOnCircle[i].multVec( my_xAx.mX, my_yAx.mX, my_zAx.mX);
		tmpPoint.mY = mNormalsOnCircle[i].multVec( my_xAx.mY, my_yAx.mY, my_zAx.mY);
		tmpPoint.mZ = mNormalsOnCircle[i].multVec( my_xAx.mZ, my_yAx.mZ, my_zAx.mZ);
		glNormal3f(-tmpPoint.mX, -tmpPoint.mY, -tmpPoint.mZ);

		tmpPoint.mX = mPointsOnCircle[i].multVec( my_xAx.mX, my_yAx.mX, my_zAx.mX) + my_start.mX;
		tmpPoint.mY = mPointsOnCircle[i].multVec( my_xAx.mY, my_yAx.mY, my_zAx.mY) + my_start.mY;
		tmpPoint.mZ = mPointsOnCircle[i].multVec( my_xAx.mZ, my_yAx.mZ, my_zAx.mZ) + my_start.mZ;
		glVertex3f(tmpPoint.mX, tmpPoint.mY, tmpPoint.mZ);
      }
    }
    glEnd();
  }

  float calcAngle( float vectX, float vectY )
  {
    const float  norm  = sqrt(vectX * vectX + vectY * vectY);
    float        angle;

	if ( norm == 0 )
		return mInvalidAngle;

    angle = 90.0 - acos( vectX / norm ) * 180.0 / mPi;
    if (vectY > 0)
      angle = -angle;

    return angle;
  }

  void calcOrthoVec(	float &orthoX, float &orthoY, float &orthoZ, 
						float VecOneX, float VecOneY, float VecOneZ, 
						float VecTwoX, float VecTwoY, float VecTwoZ )
  {
	orthoX = VecOneY * VecTwoZ - VecOneZ * VecTwoY;
	orthoY = VecOneZ * VecTwoX - VecOneX * VecTwoZ;
	orthoZ = VecOneX * VecTwoY - VecTwoX * VecOneY;
  }

  void calcCircle()
  {
    float deltaAngle = 2*mPi / mNumOfPointsOnCyl;
    float angle;
    CylPoint3D      currPoint;
    currPoint.mY = 0;

    for(angle = 0; angle < 2*mPi; angle += deltaAngle)
    {
      currPoint.mX = sin(angle);  
      currPoint.mZ = cos(angle);  
      mNormalsOnCircle.push_back(currPoint);

      currPoint.mX *= mRadius;  
      currPoint.mZ *= mRadius;  
      mPointsOnCircle.push_back(currPoint);
    }
  }

};