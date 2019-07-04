#include <SynGrasp.h>

Finger SGmakeFinger(MatrixXd DHpars, MatrixXd base,VectorXd q){
  Finger newFinger;
  newFinger.DHpars = DHpars;
  newFinger.n = DHpars.rows();
  newFinger.base = base;
  
  newFinger.q = q;
  int sizeq;
  sizeq=q.size();
  VectorXd tmp = VectorXd::Zero(sizeq,1);  
  newFinger.qin=tmp; 
  for(unsigned i=0; i<sizeq; i++){   
    newFinger.qin(i)=i; //index indicating the position of each joint in the finger
    //newFinger.qin(i)=newFinger.qin(i-1)+1;
  } 
//   cout << "Finger \n";
//   cout << newFinger.joints<< '\n';
  return newFinger;
}

Hand SGmakeHand(vector<Finger> F){
  Hand newHand;  
  newHand.F=F;
  newHand.n=F.size(); // Number of fingers
    
  int m=0;
  //cout << "whatever: \n";
  for(unsigned i=0; i<newHand.n; i++){
    m=m+F[i].n;
    newHand.q.resize(newHand.q.size()+F[i].n);
    newHand.qin.resize(newHand.qin.size()+F[i].n);
    newHand.qinf.resize(newHand.qinf.size()+F[i].n);     
    for (unsigned j=0; j<F[i].n; j++){
      newHand.q(i*4+j)=F[i].q[j];      
      newHand.qin(i*4+j)=i;      
      newHand.qinf(i*4+j)=j;
//       cout << newHand.qinf(i*4+j) << endl;
    }
  }
  newHand.m=m; // number of degrees of freedom of the hand
  newHand.ctype = 1;
  newHand.type = "User Defined Hand";
   
  newHand.ftips = SGfingertips(newHand);
  newHand = SGjoints(newHand);
  
  return newHand;
}
// %%% Pre-allocation
// hand = struct([],[],..
//             [],'S',[],'Kq',[],'Kz',[],'H',[],'J',[],'JS',[],'Wr',[]);
//
// hand.S = eye(size(hand.q,1),size(hand.q,1)); 
// 
// nj = length(hand.q);
// hand.Kq = eye(nj);
// hand.Kz = eye(nj);
// hand.Kw = eye(6);
// 
//            
// hand.limit = zeros(length(hand.q),2);
// hand.limit(:,2) = pi/2*ones(length(hand.q),1);
// 
// hand.active = ones(length(hand.q),1);
// 
// hand.Wr = zeros(6,1);
// 

MatrixXd SGfingertips(Hand hand){
  vector<int> cwhere(hand.n);   
   for(unsigned i=0; i<hand.n; i++){
    cwhere[i]=i; 
  }  
  MatrixXd cp(3, hand.n);  
  MatrixXd referenceJoint(4,4);
  MatrixXd localTransf(4,4);
  
  cp << MatrixXd::Zero(3, hand.n);
  //Calculate the position of the contact point with respect to the base reference system
   for(unsigned i=0; i<cwhere.size(); i++){    
     referenceJoint=hand.F[cwhere[i]].base; //base of the joint         
     for (unsigned j=1; j<hand.F[cwhere[i]].n+1;j++){	
       localTransf=SGDHMatrix(hand.F[cwhere[i]].DHpars.row(j-1));  
       referenceJoint = referenceJoint*localTransf;   
       
    }
       
     //the contact point is defined on the tip of each finger
     cp.block<3,1>(0,i) = referenceJoint.block<3,1>(0,3);
//       cout << "matrix is: \n";
//       cout << cp << endl;
//     the following numbers will be necessary when the position of the
//     contact point will be set arbitrarily on a generic link
//     fourth row: finger number where the contact point is 
  }    
  return cp;
}

Hand SGmoveHand_simp(Hand hand,VectorXd q){
// JUST moves the hand and places the new Fingertips  
  Hand newHand;
  newHand = hand;   
  int n=0;
  Finger F_old;
  MatrixXd DHpars(4,4);
  VectorXd q_finger(4);
  for(unsigned i=0; i<hand.n; i++){   
    F_old = hand.F[i];
    DHpars=F_old.DHpars;
    n=F_old.n;
     for(unsigned j=0; j<4; j++){   
      DHpars(j,2)=q(i*4+j)+F_old.DHpars(j,2);  
      q_finger(j)=q(i*4+j);
     }
     newHand.F[i] = SGmakeFinger(DHpars,F_old.base,q_finger);
  } 
  newHand.q = q;
  newHand.ftips = SGfingertips(newHand);
  newHand = SGjoints(newHand);
  
  return newHand;
}


Hand SGmoveHand_Gen(Hand hand,VectorXd q){
//JUST moves the hand and places the new Fingertips with fingers of different DOF
  // Also valid for translation joints
    Hand newHand;
    newHand = hand;   
    int n=0;
    Finger F_old;
    int	k = 1;
    MatrixXd DHpars(4,4);
    VectorXd q_finger(10);    
   for(unsigned i=0; i<hand.n; i++){   
      F_old = hand.F[i];
      DHpars = F_old.DHpars;
      n = F_old.n;
      for (unsigned j=0; j<n; j++){	
	if (hand.qtype[k-1+j]==0){
	  DHpars(j,2) = q(k-1+j)+DHpars(j,2);
	}
	else{
	  DHpars(j,3) = q(k-1+j)+DHpars(j,3);
	}	
      }
      for (unsigned ind=0; ind<n; ind++){
	q_finger(ind)=q(k-1+ind);
      }      
      newHand.F[i] = SGmakeFinger(DHpars,F_old.base,q_finger.segment(0,n));      	     
      k = k + n;      
    }
    newHand.q = q;
    newHand.ftips = SGfingertips(newHand);    
    newHand = SGjoints(newHand);
    
    return newHand;
}
