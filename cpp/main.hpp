//
//  main.hpp
//  Dalitz-plot
//
//  Created by Vincent Mathieu on 06/10/2019.
//  Copyright © 2019 Vincent Mathieu. All rights reserved.
//

#ifndef main_h
#define main_h

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <complex.h>
#include <unistd.h>

//#define PATH "/Users/vincent/Library/Mobile Documents/com~apple~CloudDocs/Projects/Gluex/Moments/EtaPi_Moments/"
//#define PATH ""


// masses in GeV
// mesons
#define MPI       0.13957 //061
#define MK        0.49367 //7
#define META      0.547682
#define METAP     0.95778
#define MJPSI     3.096900
#define MJPSI2S   3.686097
#define MB        5.27963

// baryons
#define MP        0.938 //27203
#define MLAMC     2.28646
#define MLAMB     5.61960


#if !defined( M_PI )
  #define M_PI 3.141592654    // otherwise it doesn't compile on FuturGrid
#endif



#endif /* main_h */
