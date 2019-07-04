#include <SynGrasp.h>

VectorXd SG_Simple_Stab_ManualDx(VectorXd w, MatrixXd Kp, Cube obj, VectorXd fmax, VectorXd &flims, MatrixXd &cp_estab, int maxiter, double mu){
  int iter=0;
  int nc=obj.normals_cp.cols();
  
  double deltaN=0.01;
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
