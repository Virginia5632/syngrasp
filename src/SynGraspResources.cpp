#include <SynGrasp.h>

double distancePointPlane (Vector3d point, VectorXd plane){
  //Point (x y x)
  // Plane (x y z, dx1 dy1 dz1, dx2, dy2, dz2)
  double d=0.0;
  Vector3d crossPlane, normPlane, dPointPlane, delement ;
  crossPlane=plane.segment<3>(3).cross(plane.segment<3>(6));
  
  dPointPlane=plane.segment<3>(0)-point;
  for (unsigned i=0; i<3; i++){
    normPlane(i)=crossPlane(i)/crossPlane.norm();    
    delement(i)=normPlane(i)*dPointPlane(i);    
  }   
  d=-delement.sum();
  
  return d;
  
}
