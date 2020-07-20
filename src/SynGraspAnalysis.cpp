#include <SynGrasp.h>

VectorXd SG_Simple_Stab_ManualDx(VectorXd w, MatrixXd Kp, Cube obj, VectorXd fmax, VectorXd &flims, MatrixXd &cp_estab, int maxiter, double mu){
  int iter=0;
  int nc=obj.normals_cp.cols();
  
  double deltaN=0.02;
  double alpha=1.0/sqrt(1.0+pow(mu,2));
  double kint=-0.01;
  
  VectorXd fc=VectorXd::Zero(2); 
  VectorXd cestab=VectorXd::Zero(nc);
  VectorXd deltaX=VectorXd::Zero(3*nc);
  VectorXd lambda=VectorXd::Zero(3*nc);  
  VectorXd sigma=VectorXd::Zero(nc); 
  MatrixXd Gcontact;
  //Pinv pGcontact;  
  MatrixXd pG=MatrixXd::Zero(3,6); ;
  flims=VectorXd::Zero(nc);  
  
  while (flims.sum()==0 && cestab.sum()<nc && iter<maxiter){
     iter=iter+1;    
     for(int i=0; i<nc; i++){   
       deltaX.segment(i*3,3)=deltaN*obj.normals_cp.block<3,1>(0,i);       
       Gcontact=obj.G.block<6,3>(0,i*3);       
       //lambda.segment<3>(i*3)=-pGcontact.pinv(Gcontact)*w+Kp.block<3,3>(i*3,i*3)*deltaX.segment<3>(i*3);         
	 pG=Gcontact.completeOrthogonalDecomposition().pseudoInverse();
	 lambda.segment<3>(i*3)=-pG*w+Kp.block<3,3>(i*3,i*3)*deltaX.segment<3>(i*3);	
       sigma(i)=alpha*lambda.segment<3>(i*3).norm()-lambda.segment<3>(i*3).transpose()*obj.normals_cp.block<3,1>(0,i);      
       cp_estab(1,i)=sigma(i);          
       if (sigma(i)<kint){         
        cp_estab(0,i)=1;
        cestab(i)=1;
       }
       fc(i)=lambda.segment<3>(i*3).norm();
    }               
   for (int j=0; j<flims.size(); j++){
       if (fc(j)>fmax(j)){
         flims(j)=1;
       }
       else{
         flims(j)=0;
       }     
    }
   if (cestab.sum()<nc){
     deltaN=deltaN+0.01; //Increment desired position in normal direction by 1cm
    }    
  }   
  return deltaX;
}


void KpComp(Eigen::MatrixXd& Kpout, Eigen::MatrixXd& KobjDes, const Eigen::MatrixXd& KobjDesPer, const Hand arm, const Cube Obj){   
  
  Eigen::VectorXd Kmax(6);
  Eigen::VectorXd Kmin(6);
  Kmax << 1000.0,1000.0,1000.0,75.0,75.0,75.0;
  Kmin <<50.0,50.0,50.0,1.0,1.0,1.0;
  
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
//   Kpout=Eigen::MatrixXd::Identity(6,6);
  int n=arm.n;
//   Kpout=Eigen::MatrixXd::Identity(arm.n*6,arm.n*6);
  Kpout=Eigen::MatrixXd::Identity(n*6,n*6);
  KobjDes=Eigen::MatrixXd::Identity(6,6);
  

  // Max and Min stiffness of each robot in the world frame
  for (int i=0; i<6; i++){
    KmaxM(i,i)=Kmax(i);
  }
  for (int i=0; i<6; i++){
    KminM(i,i)=Kmin(i);
  }
  
  
  if (arm.n==1){     
  // From Obj Percentage to actual Values 
        
        KobjDes=KminM+0.01*KobjDesPer*(KmaxM-KminM);
   std::cout << "HERE I AM : " << KobjDes << std::endl << KminM << std::endl << KmaxM << std::endl << KminM << std::endl << KobjDesPer << std::endl;
  // Build output
        Kpout=KobjDes; 
  }
    
  else if (arm.n==2){ 
    Rot_K.block<3,3>(0,0)=arm.F[0].base.block<3,3>(0,0);
    Rot_K.block<3,3>(3,3)=arm.F[0].base.block<3,3>(0,0);        
    wKmax1=Rot_K.transpose().inverse()*KmaxM*Rot_K.inverse();
    wKmin1=Rot_K.transpose().inverse()*KminM*Rot_K.inverse();
          
    Rot_K.block<3,3>(0,0)=arm.F[1].base.block<3,3>(0,0);
    Rot_K.block<3,3>(3,3)=arm.F[1].base.block<3,3>(0,0);        
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
  std::cout<<"KobjDes: \n"<< KobjDes<<std::endl;
} 



void ObjForce(Eigen::VectorXd& FMobj, const Eigen::MatrixXd& ee_FM_robots, const Hand arm, const Cube Obj){
  //Frobots (6xn) End-effector frame --> rows (Fx, Fy, Fx, Mx, My, Mz); columns (robot1, robot2, ...robotn)
   FMobj=Eigen::VectorXd::Zero(6);
   Eigen::VectorXd obj_F_ee=Eigen::VectorXd::Zero(3);
   Eigen::VectorXd obj_M_ee=Eigen::VectorXd::Zero(3);
   
   Eigen::Matrix3d w_Rot_obj=Eigen::Matrix3d::Identity();
   Eigen::Matrix3d w_Rot_ee=Eigen::Matrix3d::Identity();
   vector<Eigen::Matrix3d> w_Rots_ee;
   vector<Eigen::Matrix3d> obj_Rots_ee;
   
   w_Rot_obj=Obj.Htr.block<3,3>(0,0);
   
   for (int i=0; i<arm.n; i++){
     w_Rots_ee.push_back(Eigen::Matrix3d::Identity());
     obj_Rots_ee.push_back(Eigen::Matrix3d::Identity());
   }     
     
   Eigen::Matrix4d referencejoint=Eigen::Matrix4d::Identity();
   Eigen::Matrix4d localTransf=Eigen::Matrix4d::Identity();
   for(int i=0; i<arm.n; i++){
     referencejoint=arm.F[i].base;
     for (int j=0; j<arm.F[i].n; j++){
       localTransf=SGDHMatrix(arm.F[i].DHpars.row(j));
       referencejoint=referencejoint*localTransf;       
      }
      w_Rots_ee[i]=referencejoint.block<3,3>(0,0);
      obj_Rots_ee[i]=w_Rot_obj.inverse() * w_Rots_ee[i];
      
      obj_F_ee=obj_Rots_ee[i]*ee_FM_robots.block<3,1>(0,i);
      obj_M_ee=obj_Rots_ee[i]*ee_FM_robots.block<3,1>(3,i);
      
      FMobj.block<3,1>(0,0)=FMobj.block<3,1>(0,0)+obj_F_ee;
      FMobj.block<3,1>(3,0)=FMobj.block<3,1>(3,0)+obj_M_ee;
   }                 
}
