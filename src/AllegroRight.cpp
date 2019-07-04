#include <SynGrasp.h>

Hand AllegroRight(MatrixXd T) {
  ////////////////////////////////////////////////////////////////
  // P A R A D I G M A T I C   H A N D    P A  R A M E T E R S
  ////////////////////////////////////////////////////////////////
  double FinTipDim=28/1000; 
  
  vector<MatrixXd> DHpars; //Create an empty vector of matrices
  vector<MatrixXd> base; //Create an empty vector of matrices
  vector<Finger> F; //Create an empty vector of matrices
  
  MatrixXd tmp(4,4); tmp.setZero(); //Create a matrix and I fill with zeros
  for(int i = 0; i < 4; ++i){ //Put a matrix inside the vector (every push_back a copy is created)
    DHpars.push_back(tmp); 
    base.push_back(tmp);
  }
//   double a = DHpars[0](1,1); //I request the elment (1,1) of first matrix//   
//   Matrix4f M = DHpars[0]*DHpars[1];
//   cout << "Matrix multiplication";
//   cout << M<< '\n';
    
     
// INDEX FINGER
  double x11 = -45.098/1000.0;
  double y11 = 14.293/1000.0;

  double a12 = 54.0/1000.0;
  double a13 = 38.4/1000.0;
  double a14 = 43.7/1000.0;
  
  double rotbaseindex = 5.0*_pi/180.0;

  base[0]<< 1.0, 0.0, 0.0, x11,
	    0.0, 1.0, 0.0, y11,
	    0.0, 0.0, 1.0, 0.0,
	    0.0, 0.0, 0.0, 1.0;  
  base[0].block<3,3>(0,0)=SGrotz(rotbaseindex)*SGrotz(-_pi/2)*SGroty(-_pi/2);
//   cout << "Base is: \n";
//   cout << base[0] << endl; 
  DHpars[0]<< _pi/2, 0.0, 0.0, 0.0,
	      0.0, a12, _pi/2, 0.0,            
	      0.0, a13, 0.0, 0.0,
	      0.0, a14, 0.0, 0.0; 
	      
  // MIDDLE FINGER
  double x21 = 0.0;
  double y21 = 16.6/1000.0;

  double a22 = a12;
  double a23 = a13;
  double a24 = a14;

  base[1]<< 1.0, 0.0, 0.0, x21,
 	    0.0, 1.0, 0.0, y21,
 	    0.0, 0.0, 1.0, 0.0,
 	    0.0, 0.0, 0.0, 1.0;
  base[1].block<3,3>(0,0)=SGrotz(-_pi/2)*SGroty(-_pi/2); 
  DHpars[1]<< _pi/2, 0.0, 0.0, 0.0,
	      0.0, a22, _pi/2, 0.0,            
	      0.0, a23, 0.0, 0.0,
	      0.0, a24, 0.0, 0.0;  

  // LAST FINGER
  double x31 = 45.098/1000.0;
  double y31 = 14.293/1000.0;

  double a32 = a12;
  double a33 = a13;
  double a34 = a14;
  
  double rotbaselast = -5.0*_pi/180.0; 
  base[2]<< 1.0, 0.0, 0.0, x31,
	    0.0, 1.0, 0.0, y31,
	    0.0, 0.0, 1.0, 0.0,
	    0.0, 0.0, 0.0, 1.0;
  base[2].block<3,3>(0,0)= SGrotz(rotbaselast)*SGrotz(-_pi/2)*SGroty(-_pi/2); 
  
  DHpars[2]<< _pi/2, 0.0, 0.0, 0.0,
	      0.0, a32, _pi/2, 0.0,            
	      0.0, a33, 0.0, 0.0,
	      0.0, a34, 0.0, 0.0;  
  // THUMB
  double x41 = -16.958/1000.0;
  double y41 = -73.288/1000.0;  //(-73.288+14.293)/1000;
  double z41 = 18.2/1000.0;

  double x42 = -72.147/1000.0;
  double y42 = -78.116/1000.0;
  double z42 = 13.2/1000.0;

  double x43 = -72.147/1000.0;
  double y43 = -78.116/1000.0;
  double z43 = 13.2/1000.0;

  double x44 = -123.351/1000.0;
  double y44 = -82.596/1000.0;
  double z44 = 13.2/1000.0;
  
  double rotthumb = 5.0*_pi/180.0;
  base[3]<< 1.0, 0.0, 0.0, x41,
	    0.0, 1.0, 0.0, y41,
	    0.0, 0.0, 1.0, z41,
	    0.0, 0.0, 0.0, 1.0;
  base[3].block<3,3>(0,0)= SGrotz(rotthumb)*SGrotz(_pi/2)*SGroty(-_pi/2);
 
  double d41 = -sqrt(pow(z41-z42,2)+pow(y41-y42,2)); // d41 = -sqrt((x41-x42)^2+(y41-y42)^2);
  double a41 = z41-z42;
  double x42dh = x41 -d41*sin(rotthumb);
  double y42dh = y41 +d41*cos(rotthumb);
  double a42 = -sqrt(pow(x42-x41,2)+pow(z42-z41,2)); //a42 = -sqrt((x43-x42dh)^2+(y43-y42dh)^2);  
  double a43 = sqrt(pow(x44-x43,2)+pow(y44-y43,2));
  double a44 = 59.3/1000.0; // to check...

  DHpars[3]<< -_pi/2, -a41, 0.0,  0.0,   // -_pi/2 a41 0  d41 ;    
	      -_pi/2, 0.0, _pi/2, -a42, //  a42+40.8/1000
	      0.0, -a43, _pi/2, 0.0,
	      0.0, -a44, 0.0, 0.0;  

  
// // BUILD THE HAND
  int nfingers=DHpars.size(); 
  double joints;

  for(unsigned i=0; i<nfingers; i++){
     joints = DHpars[i].rows();
     VectorXd q = VectorXd::Zero(joints, 1);  
     F.push_back(Finger());
     F[i]=SGmakeFinger(DHpars[i], T*base[i],q);  
  }
  
  Hand newHand;
  newHand = SGmakeHand(F);
  newHand.type ="Paradigmatic";
  newHand.T=T;
  
//   cout << "Ftips \n";
//   cout << newHand.ftips<< '\n';
 
  return newHand;
 }  


