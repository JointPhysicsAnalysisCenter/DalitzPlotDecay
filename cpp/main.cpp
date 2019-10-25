//
//  main.cpp
//  Dalitz-plot
//
//  Created by Vincent Mathieu on 06/10/2019.
//  Copyright © 2019 Vincent Mathieu. All rights reserved.
//

#include "main.hpp"
#include "Dalitz-couplings.hpp"


int main(int argc, const char * argv[]) {
  printf("Code starts:\n");
  clock_t start = clock(), diff;
  
  double recoup[4] = {0.0};   // recoupling coefficients
  double inv[4]    = {0};     // the Mandelstam invariants sigma_{0,1,2,3}
  double IsoHel[4] = {0};
  
  // comment/uncomment accordingly to change the example
  /*
   * mass2[0-3]    are the masses squared of particle 0-3
   * spin[0-3]     are the spin of particle 0-3
   * spinIso[0-3]  are the spin of the isobar in the chain 1-3
   *               (only the k component is used in the k-chain - 0 is never used)
   * CMFhel[0-3]   are the CMF helicities of particle 0-3
   * IsoHel[0-3]   are the helicities of particle i and j in the k-isobar rest frame
   *               the k-component is the helicity of the isobar in the k-chain
   */
  
  //--------------------------- CHOOSE THE EXAMPLE -----------------------------------
  
  // Example A: Lambda_c(0) --> p(1) K(2) Pi(3)
  // isobars are: K* (1-chain) ; Delta (2-chain) ; Lambda (3-chain)
  // /*
  printf("Example A: Lambda_c --> p K Pi\n");
  printf("isobars are: K* (1-chain) ; Delta (2-chain) ; Lambda (3-chain)\n\n");
  double mass2[4]   = { MLAMC*MLAMC, MP*MP, MK*MK, MPI*MPI };
  double spin[4]    = {1./2, 1./2, 0.  , 0.};
  double spinIso[4] = {0   , 1.  , 1./2, 1./2};
  double CMFhel[4]  = {1./2, 1./2, 0.  , 0.};
  // */
  
  /*
  // Example B: B_0(0) --> Psi(1) Pi(2) K(3)
  // isobars are: K* (1-chain) ; no isobar (2-chain) ; Z (3-chain)
  printf("Example B: B_0(0) --> Psi(1) Pi(2) K(3)\n");
  printf("isobars are: K* (1-chain) ; no isobar (2-chain) ; Z (3-chain)\n\n");
  double mass2[4]   = { MB*MB, MJPSI2S*MJPSI2S, MPI*MPI, MK*MK};
  double spin[4]    = {0., 1., 0., 0.};
  double spinIso[4] = {0., 1., 0., 1.};
  double CMFhel[4]  = {0., 1., 0., 0.};
  */
  
  /*
  // Example C: Lambda_b(0) --> p(1) K(2) J/psi(3)
  // isobars are: no isobar (1-chain) ; Pentaquark (2-chain) ; Lambda (3-chain)
  printf("Example C: Lambda_b(0) --> p(1) K(2) J/psi(3)\n");
  printf("isobars are: no isobar (1-chain) ; Pentaquark (2-chain) ; Lambda (3-chain)\n\n");
  double mass2[4]   = { MLAMB*MLAMB, MP*MP, MK*MK, MJPSI*MJPSI };
  double spin[4]    = {1./2, 1./2, 0.  , 1.  };
  double spinIso[4] = {0.  , 0.  , 3./2, 1./2};
  double CMFhel[4]  = {1./2,-1./2, 0.  , 1.  };
  */

  //------------------------ END OF EXAMPLE DATA   --------------------------------------
  
  // loop over the Dalitz plot
  int nb = 5; // number of bins in each directions.
  double s12_0 = pow(sqrt(mass2[1]) + sqrt(mass2[2]),2);
  double s23_0 = pow(sqrt(mass2[2]) + sqrt(mass2[3]),2);
  double s12_1 = pow(sqrt(mass2[0]) - sqrt(mass2[3]),2);
  double s23_1 = pow(sqrt(mass2[0]) - sqrt(mass2[1]),2);
  double ds12 = (s12_1-s12_0)/(nb+1.), ds23 = (s23_1-s23_0)/(nb+1.);
  double eps = 0.001; // need to start just after the boundary otherwise denominator = 0
  
  printf("boudaries in 12: %6.4f - %6.4f\n", s12_0, s12_1);
  printf("boudaries in 23: %6.4f - %6.4f\n\n", s23_0, s23_1);
  
  
  for (double s12 = s12_0 + eps; s12 <= s12_1; s12 += ds12) {
  for (double s23 = s23_0 + eps; s23 <= s23_1; s23 += ds23) {
    inv[3] = s12; inv[1] = s23;
    inv[2] = mass2[0] + mass2[1] + mass2[2] + mass2[3] - s12 - s23;
    if(!check_invariants(mass2, inv)){ continue;}
    
    /*
    //--------------------------------------------------------------------------
    // Example A: Lambda_c(0) --> p(1) K(2) Pi(3)
    // 1-chain K*(+1)--> K(0) + pi(0)
    // in the 1-chain, theta_{k1} = 0 and thus CMFhel[0] == IsoHel[k]-CMFhel[k]
    IsoHel[1] = 1.; IsoHel[2] = 0.; IsoHel[3] = 0.;
    recoup[1] = recoupling_coef(1, spinIso[1], IsoHel, CMFhel, spin, mass2, inv);
    
    // 2-chain Delta(-1/2)--> p(-1/2) + pi(0)
    IsoHel[2] = -1./2; IsoHel[1] = -1./2; IsoHel[3] = 0.;
    recoup[2] = recoupling_coef(2, spinIso[2], IsoHel, CMFhel, spin, mass2, inv);
    
    // 3-chain Lambda(+1/2)--> p(+1/2) + K(0)
    IsoHel[3] = -1./2; IsoHel[1] = +1./2; IsoHel[2] = 0.;
    recoup[3] = recoupling_coef(3, spinIso[3], IsoHel, CMFhel, spin, mass2, inv);
    */
    
    /*
    //--------------------------------------------------------------------------
    // Example B: B_0(0) --> Psi(1) Pi(2) K(3)
    // 1-chain K*(+1)--> K(0) + pi(0)
    // in the 1-chain, theta_{k1} = 0 and thus CMFhel[0] == IsoHel[k]-CMFhel[k]
    IsoHel[1] = 1.; IsoHel[2] = 0.; IsoHel[3] = 0.;
    recoup[1] = recoupling_coef(1, spinIso[1], IsoHel, CMFhel, spin, mass2, inv);
    
    // 2-chain: nothing
    recoup[2] = 0.0;
    
    // 3-chain Z(0)--> psi(+1) + pi(0)
    // spin 0 for particle 0 and 3 so IsoHel[3] == 0
    IsoHel[3] = +0; IsoHel[1] = 1; IsoHel[2] = 0.;
    recoup[3] = recoupling_coef(3, spinIso[3], IsoHel, CMFhel, spin, mass2, inv);
    */
    
    /*
    //--------------------------------------------------------------------------
    // Example C: Lambda_b(0) --> p(1) K(2) J/psi(3)
    // isobars are: no isobar (1-chain) ; Pentaquark (2-chain) ; Lambda (3-chain)
    // 1-chain: nothing
    recoup[1] = 0.0;
    
    // in the 2-chain Pentaquark(-1/2) --> p(1/2) + J/psi(-1)
    IsoHel[2] = -1./2; IsoHel[1] = 1./2; IsoHel[3] = -0.;
    recoup[2] = recoupling_coef(2, spinIso[2], IsoHel, CMFhel, spin, mass2, inv);
    
    // 3-chain Lambda(1/2)--> p(+1/2) + K(0)
    IsoHel[3] = 1./2; IsoHel[1] = 1./2; IsoHel[2] = 0.;
    recoup[3] = recoupling_coef(3, spinIso[3], IsoHel, CMFhel, spin, mass2, inv);
    */
    
    //printf("(s12,s23) = (%6.4f, %6.4f) ; recouplings[1-3] =  (%6.4f,%6.4f,%6.4f)\n", s12, s23, recoup[1], recoup[2], recoup[3]);
    
  }}
  
  
  // Check with Misha
  // double center_of_mass_angle(int i, int j, double mass2[4], double inv[4]);
  // double isobar_frame_angle(int i, int j, double mass2[4], double inv[4]);
  // double wigner_angle(int i, int j, double mass2[4], double inv[4]);
  double s12 = 4.0, s23 = 1.0;
  inv[3] = s12; inv[1] = s23;
  inv[2] = mass2[0] + mass2[1] + mass2[2] + mass2[3] - s12 - s23;
  if(!check_invariants(mass2, inv)){ printf("Not inside the Dalit plot\n");}
  printf("s12 = %6.4f ; s23 = %6.4f\n", s12, s23);
  
  
  printf("CMF : 12 = %6.4f ; 23 = %6.4f ; 31 = %6.4f ;\n",
         cos(center_of_mass_angle(1, 2, mass2, inv)), cos(center_of_mass_angle(2, 3, mass2, inv)),
         cos(center_of_mass_angle(3, 1, mass2, inv)) );
  
  printf("ISO : 12 = %6.4f ; 23 = %6.4f ; 31 = %6.4f ;\n",
        cos(isobar_frame_angle(1, 2, mass2, inv)), cos(isobar_frame_angle(2, 3, mass2, inv)),
        cos(isobar_frame_angle(3, 1, mass2, inv)) );
  printf("ISO : 21 = %6.4f ; 32 = %6.4f ; 13 = %6.4f ;\n",
         cos(isobar_frame_angle(2, 1, mass2, inv)), cos(isobar_frame_angle(3, 2, mass2, inv)),
         cos(isobar_frame_angle(1, 3, mass2, inv)) );
  
  printf("WIG : 12 = %6.4f ; 23 = %6.4f ; 31 = %6.4f ;\n",
         cos(wigner_angle(1, 2, mass2, inv)), cos(wigner_angle(2, 3, mass2, inv)),
         cos(wigner_angle(3, 1, mass2, inv)) );
  
  printf("WIG : 21 = %6.4f ; 32 = %6.4f ; 13 = %6.4f ;\n",
         cos(wigner_angle(2, 1, mass2, inv)), cos(wigner_angle(3, 2, mass2, inv)),
         cos(wigner_angle(1, 3, mass2, inv)) );
  
  
  
  diff = clock() - start;
  long int msec = diff * 1000 / CLOCKS_PER_SEC;
  long int min = (msec/1000-(msec/1000)%60)/60;

  printf("\nTime taken: %ld minutes %ld seconds %ld milliseconds\n",min, (msec/1000)%60, msec%1000);
  return 0;
}
