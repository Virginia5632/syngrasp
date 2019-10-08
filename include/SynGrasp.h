// Syngrasp Header file
#include <iostream>
#include <math.h>  
#include <Eigen/Dense>
#include <Eigen/SVD>
using namespace std;

using Eigen::Matrix;
using Eigen::MatrixXd;
using Eigen::Matrix4d;
using Eigen::VectorXd;
using Eigen::Vector3d;

#define _pi 3.14159265

// Structures
struct Finger { //finger = struct('n',[],'DHpars',[],'base',[],'q',[],'qin',[]);
  int n;
  MatrixXd DHpars;
  MatrixXd base;
  VectorXd q;
  VectorXd qin;  
  MatrixXd joints;
  } ;  

struct Hand { //hand = struct('F',[],'n',[],'m',[],'q',[],'qin',[],'qinf',[],'ctype',[], ...
            // 'ftips',[],'S',[],'cp',[],'Kq',[],'Kz',[],'H',[],'J',[],'JS',[],'Wr',[]);
  string type;
  MatrixXd T;
  vector<Finger> F;
  int n; // Number of fingers
  int m; // number of degrees of freedom of the hand
  VectorXd q;
  VectorXd qin; //Finger index=[1 1 1 1  2 2 2 2  3 3 3 3 4 4 4 4]
  VectorXd qinf;//Joint index=[ 1 2 3 4  1 2 3 4  1 2 3 4  1 2 3 4]
  vector<int> qtype;  //ROT:0 ; TRANSL:1
  int ctype;
  MatrixXd ftips; // contact points location matrix
  MatrixXd cp;     
  MatrixXd J;
  MatrixXd G;
  } ;

struct Cube {   
  string type;
  Vector3d center, dim; 
  Matrix4d Htr, base;  
  MatrixXd G;
  typedef Matrix<double, 3, 4> MatrixFace;
  vector<MatrixFace> faces;
  typedef Matrix<double, 3, 6> MatrixNormals;
  MatrixNormals normals, means;
  MatrixXd normals_cp;
  double weight;
} ;  

// UTILS
MatrixXd SGrotx (double x);
MatrixXd SGroty (double x);
MatrixXd SGrotz (double x);
MatrixXd SGskew(Vector3d t);
Matrix4d SGtransl(double a, double b, double c);
MatrixXd SGDHMatrix(VectorXd v);
Hand SGjoints(Hand hand);
MatrixXd SGgTildeMatrix(MatrixXd Cp, Vector3d cm);
Cube SGcube(Matrix4d Htr,double lx, double ly, double lz, double w);
int SGfaceDetector (MatrixXd CP, Cube Obj);
Eigen::Matrix3d SGeul2rotmZYX (Eigen::Vector3d eulZYX);
Eigen::Vector3d SGrotmZYX2eul (Eigen::Matrix3d rotm);

// RESOURCES
double distancePointPlane (Vector3d point, VectorXd plane);

// GRASP MODELLING
Finger SGmakeFinger(MatrixXd DHpars,MatrixXd base,VectorXd q);
Hand SGmakeHand(vector<Finger> F);
MatrixXd SGfingertips(Hand hand);
Hand SGmoveHand_simp(Hand hand,VectorXd q);
Hand SGmoveHand_Gen(Hand hand,VectorXd q);

// GRASP ANALYSIS;
VectorXd SG_Simple_Stab_ManualDx(VectorXd w, MatrixXd Kp, Cube obj, VectorXd fmax, VectorXd &flims, MatrixXd &cp_estab, int maxiter, double mu);

// GRASP DEFINITION
void SG_CubeContact_Simp (Hand &Hand,  Cube &Cube);

// Big functions
int addition (int a, int b);
Hand AllegroRight(MatrixXd T);


class Pinv{
  public: 
    MatrixXd pinv(MatrixXd mat){          
      int xsize=mat.rows();
      int ysize=mat.cols();
      MatrixXd pmat(ysize,xsize);
	  
	  svd = Eigen::JacobiSVD<Eigen::MatrixXd>(ysize,ysize,Eigen::ComputeThinU | Eigen::ComputeThinV);
      #if EIGEN_MINOR_VERSION <= 0
	  FPL = Eigen::FullPivLU<Eigen::MatrixXd>(ysize,ysize);
      #endif  	  
	  	  
	  svd.compute(mat);
      #if EIGEN_MINOR_VERSION <= 0
	  FPL.compute(mat);
	  pmat = getDampedPinv(mat, svd, FPL);
      #else
	  pmat = getDampedPinv(mat, svd);
      #endif
      return pmat;  
    }
    Eigen::JacobiSVD<Eigen::MatrixXd> svd;
    #if EIGEN_MINOR_VERSION <= 0
    Eigen::FullPivLU<Eigen::MatrixXd> FPL;
    #endif  
    
    #if EIGEN_MINOR_VERSION <= 0
    Eigen::MatrixXd getDampedPinv(  const Eigen::MatrixXd& J,
				    const Eigen::JacobiSVD<Eigen::MatrixXd>& svd,
				    const Eigen::FullPivLU<Eigen::MatrixXd>& fpl) const
    {
	      int rank = fpl.rank();
	      Eigen::MatrixXd singularValuesInv(J.cols(), J.cols());
	      singularValuesInv.setZero(J.rows(), J.rows());

	      double lambda = svd.singularValues().minCoeff();


	      if(svd.singularValues().minCoeff() >= 1e-3)
	      {
		  for(unsigned int i = 0; i < rank; ++i)
		      singularValuesInv(i,i) = 1./svd.singularValues()[i];
	      } else {
		  //double lambda = std::pow(lambda_max,2) * (1. -  std::pow(svd.singularValues()[rank-1]/sigma_min,2));
		  for(unsigned int i = 0; i < rank; ++i)
		      singularValuesInv(i,i) =
			  svd.singularValues()[i]/(std::pow(svd.singularValues()[i],2)+lambda*lambda);
	    }


	    return svd.matrixV()* singularValuesInv *svd.matrixU().transpose();
    }
    #else
    Eigen::MatrixXd getDampedPinv(  const Eigen::MatrixXd& J,
				    const Eigen::JacobiSVD<Eigen::MatrixXd>& svd) const
    {
	int rank = svd.rank();
	Eigen::MatrixXd singularValuesInv(J.cols(), J.cols());
	//singularValuesInv.setZero(J.rows(), J.rows());
	singularValuesInv.setZero(6,6);  //I FORCED IT BECASUE IF NOT IT TAKES TOO BIG SIZE, FILL WITH 0S

	double lambda = svd.singularValues().minCoeff();

	if(svd.singularValues().minCoeff() >= 1e-3)
	{
	    for(unsigned int i = 0; i < rank; ++i)
		singularValuesInv(i,i) = 1./svd.singularValues()[i];
	} else {
	    //double lambda = std::pow(lambda_max,2) * (1. -  std::pow(svd.singularValues()[rank-1]/sigma_min,2));
	    for(unsigned int i = 0; i < rank; ++i)
		singularValuesInv(i,i) =
			svd.singularValues()[i]/(std::pow(svd.singularValues()[i],2)+lambda*lambda);
	}	
	return svd.matrixV()* singularValuesInv *svd.matrixU().transpose();	
    }
    #endif
   
};
