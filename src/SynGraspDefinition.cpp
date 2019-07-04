#include <SynGrasp.h>

void SG_CubeContact_Simp (Hand &hand,  Cube &cube){
  Vector3d center=cube.center;
  MatrixXd CPs=hand.cp.block(0,0,3,hand.cp.cols());
  MatrixXd normals(CPs.rows(), CPs.cols());  
  int face_contact;
  for(int i=0; i<CPs.cols(); i++){
    face_contact=SGfaceDetector (CPs.block<3,1>(0,i), cube);   
    normals.block<3,1>(0,i)=cube.normals.block<3,1>(0,face_contact);
  } 
  cube.normals_cp=normals;
  cube.base=Matrix4d::Identity();
  cube.base.block<3,1>(0,3)=cube.center;
}
