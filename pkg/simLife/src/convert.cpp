/**
 * @file convert.cpp
 * @date 04/26/2016
 *
 * @brief methods for conversion of R objects
 *
 * @author M.Baaske
 */

#include <R.h>
#include <Rdefines.h>

#include "GeometricPrimitives.h"

#define COPY_C2R_MATRIX(M,R,DIM)                  \
do {                                              \
    int _i, _j;                                   \
    for (_i = 0; _i < DIM; _i++)                  \
      for (_j = 0; _j < DIM; _j++)                \
      REAL((R))[_i+DIM*_j] = (M)[_i][_j];         \
} while(0)


#define COPY_R2C_MATRIX(M,R,DIM)                  \
do {                                              \
    int _i, _j;                                   \
    for (_i = 0; _i < DIM; _i++)                  \
      for (_j = 0; _j < DIM; _j++)                \
      (M)[_i][_j] = REAL((R))[_i+DIM*_j];         \
} while(0)

#define SET_REAL_VECTOR(Rv, vec)  			\
do {										\
  int _i;		                            \
  for(_i = 0; _i < LENGTH( (Rv) ); _i++) {	\
    REAL((Rv))[_i] = (vec)[_i];				\
  }											\
} while(0)

SEXP getVar(SEXP name, SEXP rho)
{
    SEXP ans;
    if(!isString(name) || length(name) != 1)
        error("name is not a single string");
    if(!isEnvironment(rho))
        error("rho should be an environment");
    ans = findVar(install(CHAR(STRING_ELT(name, 0))), rho);
    return ans;
}


///* get the elements of a list */
SEXP getListElement (SEXP list, const char *str)
{
     SEXP elmt = R_NilValue;
     SEXP names = getAttrib(list, R_NamesSymbol);

     for (int i = 0; i < length(list); i++)
         if(std::strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
             elmt = VECTOR_ELT(list, i);
             break;
         }
     return elmt;

}

SEXP convert_R_Ellipse2(STGM::CEllipse2 &ellipse) {
  SEXP R_tmp = R_NilValue;
  SEXP R_center, R_minor, R_major, R_ab, R_A;
  const char *nms[] = {"id", "type", "center", "A", "ab", "minor", "major", "phi", "S", ""};

  PROTECT(R_tmp = mkNamed(VECSXP,nms));
  PROTECT(R_center = allocVector(REALSXP, 2));
  PROTECT(R_ab = allocVector(REALSXP, 2));
  PROTECT(R_A = allocMatrix(REALSXP,2,2));
  PROTECT(R_minor = allocVector(REALSXP, 2));
  PROTECT(R_major = allocVector(REALSXP, 2));

  STGM::CVector2d &center = ellipse.center();
  SET_REAL_VECTOR(R_center,center);

  STGM::CVector2d &minor = ellipse.minorAxis();
  SET_REAL_VECTOR(R_minor,minor);

  STGM::CVector2d &major = ellipse.majorAxis();
  SET_REAL_VECTOR(R_major,major);

  REAL(R_ab)[0] = ellipse.a();    // major semi-axis (for both prolate/oblate)
  REAL(R_ab)[1] = ellipse.b();	  // minor semi-axis (for both prolate/oblate)

  const STGM::CMatrix2d &A = ellipse.MatrixA();
  COPY_C2R_MATRIX(A,R_A,2);

  SET_VECTOR_ELT(R_tmp,0,ScalarInteger(ellipse.Id()));
  SET_VECTOR_ELT(R_tmp,1,ScalarInteger(STGM::ELLIPSE_2D));
  SET_VECTOR_ELT(R_tmp,2,R_center);
  SET_VECTOR_ELT(R_tmp,3,R_A);
  SET_VECTOR_ELT(R_tmp,4,R_ab);
  SET_VECTOR_ELT(R_tmp,5,R_minor);
  SET_VECTOR_ELT(R_tmp,6,R_major);

  /* return angle between [0,2pi] */
  SET_VECTOR_ELT(R_tmp,7,ScalarReal(ellipse.phi()));
  SET_VECTOR_ELT(R_tmp,8,ScalarReal(ellipse.b()/ellipse.a()));

  UNPROTECT(6);
  return(R_tmp);

}


STGM::CSphere convert_C_Sphere(SEXP R_sphere)
{
	STGM::CVector3d ctr(REAL(getListElement( R_sphere, "center")));

	return STGM::CSphere(ctr,REAL(getListElement(R_sphere, "r"))[0],
			INTEGER(getListElement( R_sphere, "id"))[0],
			translateChar(asChar(getAttrib(R_sphere, install("label")))),
			INTEGER(getAttrib(R_sphere, install("interior")))[0]);
}

STGM::Spheres convert_C_Spheres(SEXP R_spheres)
{
  STGM::Spheres spheres;
  size_t num = (size_t) LENGTH(R_spheres);
  spheres.reserve(num);

  for(size_t i=0; i < num; i++) {
      spheres.push_back(convert_C_Sphere( VECTOR_ELT(R_spheres,i) ) );
  }

  return spheres;
}

STGM::CCylinder convert_C_Cylinder(SEXP R_cyl)
{
  int interior = 1;
  const char *label = "N";

  if(!isNull(getAttrib(R_cyl, install("label"))))
    label = translateChar(asChar(getAttrib(R_cyl, install("label"))));
  else{
	 label = "N";
	 warning(("Undefined attribute `label`. Set to 'N'."));
  }
  if(!isNull(getAttrib(R_cyl, install("interior"))))
    interior = INTEGER(getAttrib(R_cyl, install("interior")))[0];
  else {
	  interior = 0;
	  warning(("Undefined attribute `interior`. Set to zero."));
  }

  if(!std::strcmp(label,"F")) {		/* Ferrit - actually a sphere as a cylinder because of application of FBA */
      SEXP R_ctr;
      PROTECT( R_ctr = getListElement( R_cyl, "center"));
      STGM::CVector3d ctr(REAL(R_ctr)), u(0,0,1);

      UNPROTECT(1);
      return STGM::CCylinder(ctr,u,0,REAL(getListElement(R_cyl, "r"))[0],0,0,
               INTEGER(getListElement(R_cyl, "id"))[0], label, interior);

  } else {
      SEXP R_ctr, R_u, R_angles;
      PROTECT( R_ctr    =  getListElement( R_cyl, "center"));
      PROTECT( R_u      =  getListElement( R_cyl, "u"));
      PROTECT( R_angles =  getListElement( R_cyl, "angles"));
      STGM::CVector3d ctr(REAL(R_ctr)), u(REAL(R_u));

      UNPROTECT(3);
      return STGM::CCylinder(ctr,u,
    		  	REAL(getListElement(R_cyl, "h"))[0],
                REAL(getListElement(R_cyl, "r"))[0],
				REAL(R_angles)[0],REAL(R_angles)[1],
                INTEGER(getListElement(R_cyl, "id"))[0],
				label, interior);
  }

}

STGM::Cylinders convert_C_Cylinders(SEXP R_cyls)
{
  SEXP R_cyl;
  STGM::Cylinders cylinders;
  cylinders.reserve(LENGTH(R_cyls));

  for(int i=0; i < LENGTH(R_cyls); i++) {
      PROTECT(R_cyl = VECTOR_ELT(R_cyls,i));
      STGM::CVector3d ctr(REAL(VECTOR_ELT( R_cyl, 1)));
      STGM::CVector3d u(REAL(VECTOR_ELT( R_cyl, 5)));

      cylinders.push_back(
    	STGM::CCylinder(ctr,u,
    		  REAL(VECTOR_ELT(R_cyl, 4))[0],
    		  REAL(VECTOR_ELT(R_cyl, 6))[0],
    		  REAL(VECTOR_ELT(R_cyl, 7))[0],
			  REAL(VECTOR_ELT(R_cyl, 7))[1],
			  INTEGER(VECTOR_ELT(R_cyl, 0))[0],
			  translateChar(asChar(getAttrib(R_cyl, install("label")))),
			  INTEGER(getAttrib(R_cyl, install("interior")))[0]));

      UNPROTECT(1);

  }

  return cylinders;
}

STGM::CSpheroid convert_C_Spheroid(SEXP R_spheroid)
{
  STGM::CVector3d ctr(REAL(VECTOR_ELT( R_spheroid, 1)));
  STGM::CVector3d   u(REAL(VECTOR_ELT( R_spheroid, 2)));

  SEXP R_acb, R_angles;
  PROTECT( R_acb    = VECTOR_ELT(R_spheroid, 3));
  PROTECT( R_angles = VECTOR_ELT(R_spheroid, 4));
  UNPROTECT(2);

  return STGM::CSpheroid(ctr,REAL(R_acb)[0],REAL(R_acb)[1],REAL(R_acb)[2],u,
			  REAL(R_angles)[0],REAL(R_angles)[1],
			  INTEGER(getListElement( R_spheroid, "id"))[0],
			  translateChar(asChar(getAttrib(R_spheroid, install("label")))),
			  INTEGER(getAttrib(R_spheroid, install("interior")))[0]);
}


STGM::Spheroids convert_C_Spheroids(SEXP R_spheroids)
{
  STGM::Spheroids spheroids;
  size_t num = (size_t) LENGTH(R_spheroids);
  spheroids.reserve(num);

  for(size_t i=0; i < num; i++)  {
      spheroids.push_back( convert_C_Spheroid( VECTOR_ELT(R_spheroids,i) ) );
  }

  return spheroids;
}


STGM::CEllipse2 convert_C_Ellipse2(SEXP R_E)
{
   STGM::CVector2d ctr(REAL(VECTOR_ELT(R_E,2)));
   STGM::CMatrix2d A(REAL(VECTOR_ELT(R_E,3)));

   return STGM::CEllipse2(A,ctr,INTEGER(VECTOR_ELT(R_E,0))[0]);

}

STGM::Ellipses2 convert_C_Ellipses2(SEXP R_E) {
	STGM::Ellipses2 ellipses2;
	size_t num = (size_t) LENGTH(R_E);

	ellipses2.reserve(num);
	for(size_t i=0; i < num; i++)  {
		  ellipses2.push_back( convert_C_Ellipse2( VECTOR_ELT(R_E,i) ) );
	}

	return ellipses2;
}



SEXP convert_R_SphereSystem(STGM::Spheres& spheres) {
  int ncomps=3;

  SEXP R_resultlist = R_NilValue;
  PROTECT(R_resultlist = allocVector(VECSXP, spheres.size()) );

  SEXP R_tmp, R_names, R_center;
  for(size_t k=0;k<spheres.size();k++)
  {
    STGM::CSphere &sphere = spheres[k];

    PROTECT(R_tmp = allocVector(VECSXP,ncomps));
    PROTECT(R_center = allocVector(REALSXP, 3));

    REAL(R_center)[0]=sphere.center()[0];
    REAL(R_center)[1]=sphere.center()[1];
    REAL(R_center)[2]=sphere.center()[2];

    SET_VECTOR_ELT(R_tmp,0,ScalarInteger(sphere.Id()));
    SET_VECTOR_ELT(R_tmp,1,R_center);
    SET_VECTOR_ELT(R_tmp,2,ScalarReal(sphere.r()));

    PROTECT(R_names = allocVector(STRSXP, ncomps));
    SET_STRING_ELT(R_names, 0, mkChar("id"));
    SET_STRING_ELT(R_names, 1, mkChar("center"));
    SET_STRING_ELT(R_names, 2, mkChar("r"));

    setAttrib(R_tmp, R_NamesSymbol, R_names);
    SET_VECTOR_ELT(R_resultlist,k,R_tmp);
    UNPROTECT(3);
  }

  UNPROTECT(1);
  return R_resultlist;
}
