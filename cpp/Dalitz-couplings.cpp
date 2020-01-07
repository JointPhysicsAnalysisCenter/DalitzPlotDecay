//
//  Dalitz-couplings.cpp
//  
//
//  Created by Vincent Mathieu on 06/10/2019.
//

#include "Dalitz-couplings.hpp"

// *********************************************************************************
double recoupling_coef(int k, double sk, double IsoHel[4], double CMFhel[4],
                       double spin[4], double mass2[4], double inv[4]){
  // k = #chain
  // sk = isobar spin
  // IsoHel = helicities of {0,1,2,3} in the isobar rest frame
  // CMFhel = heliciticies of {0,1,2,3} in the CMF
  
  int i=1,j=2;
  if( k<1 || k>3){return 0.0;}
  if(k==1){i = 2 ; j = 3;} if(k==2){i = 3 ; j = 1;} if(k==3){i = 1 ; j = 2;}
  
  double cos_theta_star_k1 = cos_center_of_mass_angle(k, 1, mass2, inv);
  double cos_theta_ij = cos_isobar_frame_angle(i, j, mass2, inv);
  double cos_omega_ik = cos_wigner_angle(i, k, mass2, inv);
  double cos_omega_jk = cos_wigner_angle(j, k, mass2, inv);
  
  
  double d1 = wignerD(spin[0], CMFhel[0], IsoHel[k]-CMFhel[k], cos_theta_star_k1 );
  double d2 = wignerD(sk     , IsoHel[k], CMFhel[i]-CMFhel[j], cos_theta_ij );
  double d3 = wignerD(spin[i], IsoHel[i], CMFhel[i], cos_omega_ik );
  double d4 = wignerD(spin[j], IsoHel[j], CMFhel[j], cos_omega_jk );
  d3 = d3*pow(-1,IsoHel[i]-CMFhel[i]);  // d(-w)
  
  /*
  // to debug
  if(k==2){
    printf("in the dalitz: %6.4f ; %6.4f ; %6.4f ; %6.4f \n", d1, d2, d3, d4);
  } */
  
  
  return d1*d2*d3*d4;
}

// *********************************************************************************
double cos_center_of_mass_angle(int i, int j, double mass2[4], double inv[4]){
  // return the COSINE of the angle between i and j in the (ijk) = (0) rest frame
  
  if (i==j){return 1.0;}    // i,j are different
  if ( i<1 || j<1 || i>3 || j>3 ){return 1.0;}
  int k = 6-(i+j);        // i+j+k = 6;
  if(!check_invariants(mass2, inv)){ return 1.0;}
  
  // denominator
  double lam1 = lambda(mass2[0], mass2[i], inv[i]);
  double lam2 = lambda(mass2[0], mass2[j], inv[j]);
  double deno = sqrt(lam1*lam2);
  
  // numerator
  double num = (mass2[0]+mass2[i]-inv[i])*(mass2[0]+mass2[j]-inv[j]);
  num =  num + 2.*mass2[0]*(mass2[i] + mass2[j] - inv[k]);
  
  return num/deno;
}

// *********************************************************************************
double cos_isobar_frame_angle(int i, int j, double mass2[4], double inv[4]){
  // return the COSINE of the angle between i and -k in the (ij) = (k0) rest frame
  
  // basic checks
  if (i==j){return 1.0;}    // i,j are different
  if ( i<1 || j<1 || i>3 || j>3 ){return 1.0;}
  int k = 6-(i+j);        // i+j+k = 6;
  if(!check_invariants(mass2, inv)){ return 1.0;}
  
  // denominator
  double lam1 = lambda(mass2[i], mass2[j], inv[k]);
  double lam2 = lambda(mass2[0], mass2[k], inv[k]);
  double deno = sqrt(lam1*lam2);
  
  // numerator
  double num = inv[k]*(inv[i]-inv[j]) + (mass2[i]-mass2[j])*(mass2[0]-mass2[k]);
  
  return num/deno;
}

// *********************************************************************************
double cos_wigner_angle(int i, int j, double mass2[4], double inv[4]){
  // return the COSINE of the angle of the Wigner rotation for
  // particle i boosted from the (ij) frame to the (ijk) rest frame
 
  // basic checks
  if (i==j){return 1.0;}    // i,j are different
  if ( i<1 || j<1 || i>3 || j>3 ){return 1.0;}
  int k = 6-(i+j);        // i+j+k = 6;
  if(!check_invariants(mass2, inv)){ return 1.0;}
  
  // denominator
  double lam1 = lambda(mass2[i], mass2[j], inv[k]);
  double lam2 = lambda(mass2[0], mass2[i], inv[i]);
  double deno = sqrt(lam1*lam2);
  
  // numerator
  double num = (mass2[0]+mass2[i]-inv[i])*(mass2[i]-mass2[j]+inv[k]);
  num =  num - 2.*mass2[i]*(mass2[0] - mass2[k] + inv[k]);
  
  return num/deno;
}


// *********************************************************************************
bool check_invariants(double mass2[4], double inv[4]){
  // check kinematical relations

  double m0 = sqrt(mass2[0]), m1 = sqrt(mass2[1]), m2 = sqrt(mass2[2]), m3 = sqrt(mass2[3]);
  if (inv[1]<0 || inv[2]<0 || inv[3] <0){return false;}
  if ( sqrt(inv[1])< m2+m3 || sqrt(inv[2])< m3+m1 || sqrt(inv[3])< m1+m2 ) {return false;}
  if ( sqrt(inv[1])> m0-m1 || sqrt(inv[2])> m0-m2 || sqrt(inv[3])> m0-m3 ) {return false;}
  
  // need to be inside the Dalitz plot
  double E2 = (inv[3] - mass2[1] + mass2[2])/(2.*sqrt(inv[3]));
  double E3 = (mass2[0] - inv[3] - mass2[3])/(2.*sqrt(inv[3]));
  double s1min = pow(E2+E3,2) - pow( sqrt(E2*E2-mass2[2]) + sqrt(E3*E3-mass2[3]),2);
  double s1max = pow(E2+E3,2) - pow( sqrt(E2*E2-mass2[2]) - sqrt(E3*E3-mass2[3]),2);
  //printf("boundaries: %6.4f , %6.4f\n", s1min, s1max);
  if ( inv[1] >= s1max || inv[1] <= s1min ) {return false;}
  
  return true;
}



// *********************************************************************************
double lambda(double a, double b, double c){
  // triangle function
  return a*a + b*b + c*c - 2*(a*b + b*c + c*a);
}

// *********************************************************************************
double wignerD(double J, double M1, double M2, double z){
  /*
   *  compute the Wigner D^j_{M1,M2}( z = cos[theta] )
   *  using the summation formula by Wigner
   */
  
  double logfact[51] = {0.0, 0.0,
    6.93147180559945309e-1, 1.79175946922805500e00, 3.17805383034794562e00,
    4.78749174278204599e00, 6.57925121201010100e00, 8.52516136106541430e00,
    1.06046029027452502e01, 1.28018274800814696e01, 1.51044125730755153e01,
    1.75023078458738858e01, 1.99872144956618861e01, 2.25521638531234229e01,
    2.51912211827386815e01, 2.78992713838408916e01, 3.06718601060806728e01,
    3.35050734501368889e01, 3.63954452080330536e01, 3.93398841871994940e01,
    4.23356164607534850e01, 4.53801388984769080e01, 4.84711813518352239e01,
    5.16066755677643736e01, 5.47847293981123192e01, 5.80036052229805199e01,
    6.12617017610020020e01, 6.45575386270063311e01, 6.78897431371815350e01,
    7.12570389671680090e01, 74.6582363488301644, 78.0922235533153106,
    81.5579594561150372, 85.0544670175815174, 88.5808275421976788,
    92.1361756036870925, 95.7196945421432025, 99.3306124547874269,
    102.968198614513813, 106.631760260643459, 110.320639714757395,
    114.034211781461703, 117.771881399745072, 121.533081515438634,
    125.317271149356895, 129.123933639127215, 132.952575035616310,
    136.802722637326368, 140.673923648234259, 144.565743946344886,
    148.477766951773032
  };
  
  // USE THE PHYSICS CONVENTION!!!
  M1 = -M1; M2 = -M2;
  
  double wd = 0.0, fac = 0.0, half = 0.0, deno = 0.0, num = 0.0;
  double cosh = sqrt((1+z)/2.), sinh = sqrt((1-z)/2.);
  int lb = 0;     // lower bound for s
  int ub = J+M2;  // upper bound for s
  if(M2-M1>0)   lb = M2-M1;
  if(J-M1<J+M2) ub = J-M1;
  
  int i1= J+M1,i2 = J-M1,i3 = J+M2,i4 = J-M2, i5 = M1-M2;
  if( J+abs(M1) > 50 || J+abs(M2) > 50 ){ printf("not enough factorials stored!\n");}
  
  num = (logfact[i1] + logfact[i2] + logfact[i3] + logfact[i4])/2.;
  
  for(int s = lb ; s<=ub ; s++){
    half = pow(cosh, 2*J+M2-M1-2*s)*pow(sinh, M1-M2+2*s) ;
    deno = logfact[i3-s] + logfact[s] + logfact[i5+s] + logfact[i2-s];
    fac = exp(num-deno);
    wd = wd + half*pow(-1,s)*fac;
  }
  
  return wd;
}


// ***************************************************
int isfrac(double x){
  int val = 1;
  double intpart;
  if( modf(x, &intpart) == 0.0){ val = 0;}
  return val;
}

// ***************************************************
long int factorial(int n)
{
  if(n > 1)
    return n * factorial(n - 1);
  else
    return 1;
}

