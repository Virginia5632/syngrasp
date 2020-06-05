#include <SynGrasp.h>

int addition (int a, int b)
{
  int r;
  r=a+b;
  return r;
}

MatrixXd SGrotx(double x){
  MatrixXd R = MatrixXd::Identity(3,3);
  R<< 1,  0, 0,
      0,  cos(x),  -sin(x),
      0,  sin(x),  cos(x);    
   return R;
}

MatrixXd SGroty(double x){
  MatrixXd R = MatrixXd::Identity(3,3); 
  R<< cos(x),  0,  sin(x),
      0,       1,       0,
      -sin(x), 0,       cos(x);
  return R;
}

MatrixXd SGrotz(double x){
  MatrixXd R = MatrixXd::Identity(3,3);
  R<<cos(x),  -sin(x),  0,   
     sin(x),  cos(x),   0,
     0,       0,        1;
  return R;
}

// MatrixXd SGrot(Axis axis, double x)
// {
//     MatrixXd R = MatrixXd::Identity(3,3);
//     switch (axis)
//     {
//     case X:
//         R << 1, 0,         0,
//              0, cos(x), -sin(x),
//              0, sin(x),  cos(x);
//         break;
//     case Y:
//         R  << cos(x), 0, sin(x),
//               0,        1, 0,
//              -sin(x), 0, cos(x);
//         break;
//     case Z:
//         R << cos(x), -sin(x), 0,
//              sin(x),  cos(x), 0,
//              0,         0,        1;
//         break;  
//     }
//     return R;
// }


MatrixXd SGskew(Vector3d t){
    MatrixXd S(3,3);
    S << 0, -t(2), t(1),
         t(2), 0, -t(0),
        -t(1), t(0), 0;
   return S;
}

Matrix4d SGtransl(double a, double b, double c){
//returns the homogeneous transformation matrix translation by P
  Vector3d P(a,b,c);
  Matrix4d H = Matrix4d::Identity(4,4);
  H.block<3,1>(0,3) = P;
  return H;
}

MatrixXd SGDHMatrix(VectorXd v){
  double alpha=v(0);
  double A=v(1);
  double theta=v(2);
  double D=v(3);
// alpha,theta -> radians
// Denavit-Hartenberg homogeneous transformation matrix
  MatrixXd H(4,4);
  H<<cos(theta), -sin(theta)*cos(alpha), sin(theta)*sin(alpha), A*cos(theta),
    sin(theta), cos(theta)*cos(alpha), -cos(theta)*sin(alpha), A*sin(theta),
    0,          sin(alpha),             cos(alpha),            D,
    0,          0,                      0,                       1;
  return H;
}

Hand SGjoints(Hand hand){
  Hand newHand;
  newHand=hand;
  MatrixXd referenceJoint(4,4);
  MatrixXd localTransf(4,4);    
  for(unsigned i=0; i<hand.n; i++){
      newHand.F[i].joints = MatrixXd::Zero(3,hand.F[i].n+1);  
      newHand.F[i].joints.block<3,1>(0,0) = hand.F[i].base.block<3,1>(0,3); //(1:3,4);
      referenceJoint = hand.F[i].base;         
      for (unsigned j=0; j<hand.F[i].n;j++){
       localTransf=SGDHMatrix(hand.F[i].DHpars.row(j));            
       referenceJoint = referenceJoint*localTransf;   
       newHand.F[i].joints.block<3,1>(0,j+1) = referenceJoint.block<3,1>(0,3); //(1:3,4)];        
      }  
  }      
  return newHand;
  
}

MatrixXd SGgTildeMatrix(MatrixXd Cp, Vector3d cm){
// Builds the grasp matrix G for a configuration with n contact point c_i and normals n_i. Contact point coordinates in base frame are passed as rows of argument Cp (nx3). cm is the center of mass vetor.  
  int ndim=Cp.rows();
  int ncp=Cp.cols();
  MatrixXd G = MatrixXd::Zero(ndim,ncp*ndim); 
  MatrixXd Gtilde;
  int k=1;
  for(unsigned i=0; i<ncp; i++){
    G.block(0,k-1,ndim,ndim)=MatrixXd::Identity(ndim,ndim);    //matrix.block(i,j,p,q);matrix.block<p,q>(i,j);
    k = k + ndim;
  }
  if(ndim == 3) {
    MatrixXd Gl= MatrixXd::Zero(ndim,ncp*ndim);   
    k=1;
    for(unsigned i=0; i<ncp; i++){      
      Gl.block(0,k-1,ndim,ndim)<< 0, -(Cp(2,i)-cm(2)), (Cp(1,i)-cm(1)),
				  (Cp(2,i)-cm(2)), 0, -(Cp(0,i)-cm(0)),
				  -(Cp(1,i)-cm(1)), (Cp(0,i)-cm(0)), 0;				  
      k = k + ndim;
    }
    Gtilde = MatrixXd::Zero(2*ndim,2*ncp*ndim);
    // Gtilde = [G zeros(ndim,ndim*n);Gl G];
    Gtilde.block(0,0,ndim,ncp*ndim)=G; 
    Gtilde.block(0,ncp*ndim,ndim,ncp*ndim)=MatrixXd::Zero(ndim,ncp*ndim);
    Gtilde.block(ndim,0,ndim,ncp*ndim)=Gl;
    Gtilde.block(ndim,ncp*ndim,ndim,ncp*ndim)=G;
  } 
  else {
    MatrixXd Gl= MatrixXd::Zero(1,ncp*ndim);
    k=1;
    for(unsigned i=0; i<ncp; i++){
      Gl.block(0,k-1,1,ndim)<< -(Cp(1,i)-cm(1)), (Cp(0,i)-cm(0)); 
      k = k + ndim;
    }
    Gtilde = MatrixXd::Zero((ndim+1),(ncp*ndim+ncp));
    //Gtilde = [G zeros(ndim,n);Gl ones(1,n)];
    Gtilde.block(0,0,ndim,ncp*ndim)=G; 
    Gtilde.block(0,ncp*ndim,ndim,ncp)=MatrixXd::Zero(ndim,ncp);
    Gtilde.block(ndim,0,1,ncp*ndim)=Gl;
    Gtilde.block(ndim,ncp*ndim,1,ncp)<<1,1,1;
  }
  return Gtilde;

}
/*function Gtilde = SGgTildeMatrix(Cp,cm)
% Builds the grasp matrix G for a configuration
% with n contact point c_i and normals n_i.
% Contact point coordinates in base frame are passed
% as rows of argument Cp (nx3). 
% cm is the center of mass vetor.

[n ndim] = size(Cp);
[n ndim] = size(Cp);
G = zeros(ndim, n*ndim);
I = eye(ndim);
k = 1; % first column index
for i = 1:n
    % Same as G = []; G = [G I]; - dimensions statically determined
    G(:,k:(k+ndim-1)) = I;
    k = k + ndim;
end

if ndim == 3
    Gl = zeros(ndim, n*ndim);
    k = 1;
    for i=1:n        
        Cx = [ 0 -(Cp(i,3)-cm(3)) (Cp(i,2)-cm(2));
               (Cp(i,3)-cm(3)) 0 -(Cp(i,1)-cm(1));
              -(Cp(i,2)-cm(2)) (Cp(i,1)-cm(1)) 0];
        Gl(:,k:(k+ndim-1)) = Cx;    
    k = k + ndim;
    end
    
elseif ndim == 2
    Gl = zeros(1, n*ndim);
    k = 1;
    for i=1:n        
        Cx = [-(Cp(i,2)-cm(2)) (Cp(i,1)-cm(1))]; 
        Gl(:,k:(k+ndim-1)) = Cx;        
    k = k + ndim;
    end    
end

if ndim==3
    Gtilde = [G zeros(ndim,ndim*n);Gl G];
else
    Gtilde = [G zeros(ndim,n);Gl ones(1,n)];
end
*/


Cube SGcube(Matrix4d Htr,double lx, double ly, double lz, double w){
  Cube newCube;
  newCube.type = "cube";
  newCube.center = Htr.block<3,1>(0,3);
  newCube.weight=w;
  
  Matrix4d vertex;
  Vector3d v0,v1,v2,v3,v4,v5,v6,v7;
  //vertices:
  vertex=Htr*SGtransl(lx/2,-ly/2,-lz/2) ;
    v0=vertex.block<3,1>(0,3); //vertex 1
  vertex=Htr*SGtransl(lx/2,ly/2,-lz/2) ;
    v1=vertex.block<3,1>(0,3); //vertex 2
  vertex=Htr*SGtransl(-lx/2,ly/2,-lz/2);
    v2=vertex.block<3,1>(0,3); //vertex 3
  vertex=Htr*SGtransl(-lx/2,-ly/2,-lz/2);
    v3=vertex.block<3,1>(0,3); //vertex 4
  vertex=Htr*SGtransl(-lx/2,-ly/2,lz/2);
    v4=vertex.block<3,1>(0,3); //vertex 5
  vertex=Htr*SGtransl(lx/2,-ly/2,lz/2);
    v5=vertex.block<3,1>(0,3); //vertex 6
  vertex=Htr*SGtransl(lx/2,ly/2,lz/2);
    v6=vertex.block<3,1>(0,3); //vertex 7
  vertex=Htr*SGtransl(-lx/2,ly/2,lz/2);
    v7=vertex.block<3,1>(0,3); //vertex 8

  //faces
  newCube.Htr = Htr;
  newCube.dim << lx,ly,lz;
  for (unsigned i=0; i<6; i++){	
    newCube.faces.push_back(MatrixXd::Zero(3,4));
  }
  newCube.faces[0]<<v0,v1,v6,v5;
  newCube.faces[1]<<v1,v2,v7,v6;
  newCube.faces[2]<<v2,v3,v4,v7;
  newCube.faces[3]<<v3,v0,v5,v4;
  newCube.faces[4]<<v0,v3,v2,v1;
  newCube.faces[5]<<v6,v7,v4,v5;
 
  //normals and mean points
  for (unsigned i=0; i<6; i++){	
      Vector3d w0,w1,w2,m;
      w0 = newCube.faces[i].block<3,1>(0,0);
      w1 = newCube.faces[i].block<3,1>(0,1);
      w2 = newCube.faces[i].block<3,1>(0,2);
      m=(w0-w1).cross((w2-w1));
      newCube.normals.block<3,1>(0,i)=m/m.norm();      
      newCube.means.block<3,1>(0,i)=newCube.faces[i].rowwise().mean();
  }    
  return newCube;
}

int SGfaceDetector (MatrixXd CP, Cube Obj){
  int face;
  Vector3d v1,v2,v3;  
  VectorXd plane(9), d(6);
  double dmin;
  for (unsigned i=0; i<6; i++){
    v1=Obj.faces[i].block<3,1>(0,0);
    v2=Obj.faces[i].block<3,1>(0,1);
    v3=Obj.faces[i].block<3,1>(0,2);
    plane.segment<3>(0)=v1;
    plane.segment<3>(3)=v2-v1;
    plane.segment<3>(6)=v3-v1;
    d(i)=fabs(distancePointPlane (CP, plane));        
  }  
  std::ptrdiff_t facerow, facecol;
  dmin = d.minCoeff(&facerow,&facecol); 
  return facerow;
}

Eigen::Matrix3d SGeul2rotmZYX (Eigen::Vector3d eulZYX){
  Eigen::Matrix3d rotm;
  
  rotm(0,0)=cos(eulZYX(1))*cos(eulZYX(0));
  rotm(0,1)=sin(eulZYX(2))*sin(eulZYX(1))*cos(eulZYX(0))-cos(eulZYX(2))*sin(eulZYX(0));
  rotm(0,2)=cos(eulZYX(2))*sin(eulZYX(1))*cos(eulZYX(0))+sin(eulZYX(2))*sin(eulZYX(0));
  
  rotm(1,0)=cos(eulZYX(1))*sin(eulZYX(0));
  rotm(1,1)=sin(eulZYX(2))*sin(eulZYX(1))*sin(eulZYX(0))+cos(eulZYX(2))*cos(eulZYX(0));
  rotm(1,2)=cos(eulZYX(2))*sin(eulZYX(1))*sin(eulZYX(0))-sin(eulZYX(2))*cos(eulZYX(0));
  
  rotm(2,0)=-sin(eulZYX(1));
  rotm(2,1)=sin(eulZYX(2))*cos(eulZYX(1));
  rotm(2,2)=cos(eulZYX(2))*cos(eulZYX(1));
  
  return rotm;
}

Eigen::Vector3d SGrotmZYX2eul (Eigen::Matrix3d rotm){
  Eigen::Vector3d eulZYX;
  double sy;
  sy=sqrt(rotm(0,0) * rotm(0,0) +  rotm(1,0) * rotm(1,0));
  eulZYX(0) = atan2(rotm(1,0), rotm(0,0));
  eulZYX(1) = atan2(-rotm(2,0), sy);
  eulZYX(2) = atan2(rotm(2,1) , rotm(2,2)); 
  
  
  return eulZYX;
}



void KpComp(Eigen::MatrixXd& Kpout, Eigen::MatrixXd& KobjDes, const Eigen::MatrixXd& KobjDesPer, const Hand arm, const Cube Obj){   
  
  Eigen::VectorXd Kmax(6);
  Eigen::VectorXd Kmin(6);
  Kmax << 500.0,500.0,500.0,50.0,50.0,50.0;  
  Kmin<<30.0,30.0,30.0,5.0,5.0,5.0;
  
  Eigen::MatrixXd Kp;
  Eigen::MatrixXd KpComp(Eigen::MatrixXd);
  
  Eigen::Matrix3d KT=Eigen::Matrix3d::Identity();
  Eigen::Matrix3d KR=Eigen::Matrix3d::Identity();
  Eigen::Matrix3d KTR=Eigen::Matrix3d::Identity();
  Eigen::Matrix3d Kaux=Eigen::Matrix3d::Identity();
  
  Eigen::MatrixXd wKobjDesPer=Eigen::MatrixXd::Identity(6,6);
  Eigen::MatrixXd wKobjDes=Eigen::MatrixXd::Identity(6,6);
  Eigen::Matrix3d wKobjTmin=Eigen::MatrixXd::Identity(3,3);
  Eigen::Matrix3d wKobjTmax=Eigen::MatrixXd::Identity(3,3);
  Eigen::Matrix3d wKobjRmin=Eigen::MatrixXd::Identity(3,3);
  Eigen::Matrix3d wKobjRmax=Eigen::MatrixXd::Identity(3,3);
  Eigen::Matrix3d wKobjRmin_act=Eigen::Matrix3d::Identity();
  Eigen::Matrix3d wKobjRmax_act=Eigen::Matrix3d::Identity();

  Eigen::MatrixXd KmaxM=Eigen::MatrixXd::Identity(6,6);
  Eigen::MatrixXd KminM=Eigen::MatrixXd::Identity(6,6);
  Eigen::MatrixXd Rot_K=Eigen::MatrixXd::Identity(6,6);
  Eigen::Matrix3d SkewC1, SkewC2;
  Eigen::MatrixXd wKmax1=Eigen::MatrixXd::Identity(6,6);
  Eigen::MatrixXd wKmax2=Eigen::MatrixXd::Identity(6,6);
  Eigen::MatrixXd wKmin1=Eigen::MatrixXd::Identity(6,6);
  Eigen::MatrixXd wKmin2=Eigen::MatrixXd::Identity(6,6);
  Eigen::MatrixXd wK1=Eigen::MatrixXd::Identity(6,6);
  Eigen::MatrixXd wK2=Eigen::MatrixXd::Identity(6,6);
  Eigen::MatrixXd K1=Eigen::MatrixXd::Identity(6,6);
  Eigen::MatrixXd K2=Eigen::MatrixXd::Identity(6,6);      
  
 
  // Max and Min stiffness of each robot in the world frame
  Kpout=Eigen::MatrixXd::Identity(arm.n*6,arm.n*6);
  KobjDes=Eigen::MatrixXd::Identity(6,6);
  

  // Max and Min stiffness of each robot in the world frame
  for (int i=0; i<3; i++){
    KmaxM(i,i)=Kmax(i);
  }

  for (int i=0; i<3; i++){
    KminM(i,i)=Kmin(i);
  }
  
  
  if (arm.n==1){     
  // From Obj Percentage to actual Values 
  for (int i=0; i<3; i++){
    for (int j=0; j<3; j++){
      KobjDes(i,j)=KminM(i,j)+0.01*KobjDesPer(i,j)*(KmaxM(i,j)-KminM(i,j));
    }
  }   
  
  // Build output
  Kpout.block<3,3>(0,0)=KobjDes.block<3,3>(0,0);
  Kpout.block<3,3>(3,3)=KobjDes.block<3,3>(3,3);   
  }
    
  if (arm.n==2){    
    Rot_K.block<3,3>(0,0)=arm.F[0].base.block<3,3>(0,0);
    Rot_K.block<3,3>(3,3)=arm.F[0].base.block<3,3>(0,0);;        
    wKmax1=Rot_K.transpose().inverse()*KmaxM*Rot_K.inverse();
    wKmin1=Rot_K.transpose().inverse()*KminM*Rot_K.inverse();
          
    Rot_K.block<3,3>(0,0)=arm.F[1].base.block<3,3>(0,0);;
    Rot_K.block<3,3>(3,3)=arm.F[1].base.block<3,3>(0,0);;        
    wKmax2=Rot_K.transpose().inverse()*KmaxM*Rot_K.inverse();
    wKmin2=Rot_K.transpose().inverse()*KminM*Rot_K.inverse();
    
    // Max and Min stiffness of the object in the world frame
    Rot_K.block<3,3>(0,0)=Obj.Htr.block<3,3>(0,0);
    Rot_K.block<3,3>(3,3)=Obj.Htr.block<3,3>(0,0);
    wKobjDesPer=Rot_K.transpose().inverse()*KobjDesPer*Rot_K.inverse();
    SkewC1=Obj.G.block<3,3>(3,0);
    SkewC2=Obj.G.block<3,3>(3,3);
    
    wKobjTmin=wKmin1.block<3,3>(0,0)+wKmin2.block<3,3>(0,0);
    wKobjTmax=wKmax1.block<3,3>(0,0)+wKmax2.block<3,3>(0,0);
    
    wKobjRmin=(wKmin1.block<3,3>(3,3)+wKmin2.block<3,3>(3,3))+SkewC1*wKmin1.block<3,3>(0,0)*SkewC1.transpose()+SkewC2*wKmin2.block<3,3>(0,0)*SkewC2.transpose();
    wKobjRmax=(wKmax1.block<3,3>(3,3)+wKmax2.block<3,3>(3,3))+SkewC1*wKmax1.block<3,3>(0,0)*SkewC1.transpose()+SkewC2*wKmax2.block<3,3>(0,0)*SkewC2.transpose();  
    
    // From Obj Percentage to actual Values
    /////////// Translation
    for (int i=0; i<3; i++){
      for (int j=0; j<3; j++){
        wKobjDes(i,j)=wKobjTmin(i,j)+0.01*wKobjDesPer(i,j)*(wKobjTmax(i,j)-wKobjTmin(i,j));
      }
    }  
    KT=0.5*wKobjDes.block<3,3>(0,0);
    
    /////////// Couple terms
    KTR=Eigen::Matrix3d::Zero();
    wKobjDes.block<3,3>(0,3)=2.0*KTR+KT*(SkewC1.transpose()+SkewC2.transpose());
    wKobjDes.block<3,3>(3,0)=wKobjDes.block<3,3>(0,3);
    
    ////////Rotation
    Kaux=SkewC1*KT*SkewC1.transpose()+SkewC2*KT*SkewC2.transpose(); //+KTR*(SkewC1.transpose()+SkewC2.transpose())+(SkewC1+SkewC2)*KTR;
    wKobjRmin_act=(wKmin1.block<3,3>(3,3)+wKmin2.block<3,3>(3,3))+Kaux;
    wKobjRmax_act=(wKmax1.block<3,3>(3,3)+wKmax2.block<3,3>(3,3))+Kaux;

    for (int i=0; i<3; i++){
      for (int j=0; j<3; j++){
        wKobjDes(i+3,j+3)=wKobjRmin_act(i,j)+0.01*wKobjDesPer(i+3,j+3)*(wKobjRmax_act(i,j)-wKobjRmin_act(i,j));
      }
    }
  //   cout<<"Kaux: \n"<<Kaux<<endl;
    KR=0.5*(wKobjDes.block<3,3>(3,3)-Kaux);
    wKobjDes.block<3,3>(3,3)=2.0*KR+Kaux;
    
    //From World to Local frame  
    wK1.block<3,3>(0,0)=KT;
    wK1.block<3,3>(3,3)=KR;
    Rot_K.block<3,3>(0,0)=arm.F[0].base.block<3,3>(0,0);;
    Rot_K.block<3,3>(3,3)=arm.F[0].base.block<3,3>(0,0);;  
    K1=Rot_K.transpose()*wK1*Rot_K;  
    
    wK2.block<3,3>(0,0)=KT;
    wK2.block<3,3>(3,3)=KR;
    Rot_K.block<3,3>(0,0)=arm.F[1].base.block<3,3>(0,0);;
    Rot_K.block<3,3>(3,3)=arm.F[1].base.block<3,3>(0,0);;  
    K2=Rot_K.transpose()*wK2*Rot_K;  
  
    
    Rot_K.block<3,3>(0,0)=Obj.Htr.block<3,3>(0,0);
    Rot_K.block<3,3>(3,3)=Obj.Htr.block<3,3>(0,0);
    KobjDes=Rot_K.transpose()*wKobjDes*Rot_K; 
    
    // Build output
    Kpout.block<3,3>(0,0)=K1.block<3,3>(0,0);
    Kpout.block<3,3>(3,3)=K2.block<3,3>(0,0);
    Kpout.block<3,3>(0,6)=KTR;
    Kpout.block<3,3>(6,0)=KTR;
    Kpout.block<3,3>(3,9)=KTR;
    Kpout.block<3,3>(9,3)=KTR;
    Kpout.block<3,3>(6,6)=K1.block<3,3>(3,3);
    Kpout.block<3,3>(9,9)=K2.block<3,3>(3,3);

    std::cout<<"wKobjTmin: \n"<<wKobjTmin<<std::endl;
    std::cout<<"wKobjTmax: \n"<<wKobjTmax<<std::endl;
    std::cout<<"wKobjRmin: \n"<<wKobjRmin<<std::endl;
    std::cout<<"wKobjRmin_act: \n"<<wKobjRmin_act<<std::endl;
    std::cout<<"wKobjRmax: \n"<<wKobjRmax<<std::endl;
    std::cout<<"wKobjRmax_act: \n"<<wKobjRmax_act<<std::endl;    
        
    std::cout<<"wK1: \n"<<wK1<<std::endl;
    std::cout<<"K1: \n"<<K1<<std::endl;
    std::cout<<"wK2: \n"<<wK1<<std::endl;
    std::cout<<"K2: \n"<<K1<<std::endl;
    
    std::cout<<"wKobjDesPer: \n"<<wKobjDesPer<<std::endl;
    std::cout<<"wKobjDes: \n"<<wKobjDes<<std::endl;
  }
  
  std::cout<<"KobjDesPer: \n"<<KobjDesPer<<std::endl; 
  std::cout<<"KobjDes: \n"<< KobjDesPer<<std::endl;
} 
