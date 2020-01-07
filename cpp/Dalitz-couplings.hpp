//
//  Dalitz-couplings.hpp
//  
//
//  Created by Vincent Mathieu on 06/10/2019.
//

#ifndef Dalitz_couplings_hpp
#define Dalitz_couplings_hpp

//#include "main.hpp"

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <complex.h>
#include <unistd.h>

double lambda(double a, double b, double c);
double wignerD(double J, double M1, double M2, double z);
int isfrac(double x);
long int factorial(int n);

double recoupling_coef(int k, double sk, double IsoHel[4], double CMFhel[4],
                       double spin[4], double mass2[4], double inv[4]);
double cos_center_of_mass_angle(int i, int j, double mass2[4], double inv[4]);
double cos_isobar_frame_angle(int i, int j, double mass2[4], double inv[4]);
double cos_wigner_angle(int i, int j, double mass2[4], double inv[4]);
bool check_invariants(double mass2[4], double inv[4]);

#endif /* Dalitz_couplings_hpp */
