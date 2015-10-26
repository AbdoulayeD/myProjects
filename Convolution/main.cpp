//
//  main.cpp
//  Convolution
//
//  Created by DIOP Abdoulaye  on 20/10/2015.
//  Copyright © 2015 DIOP Abdoulaye . All rights reserved.
//

//typedef convolType int;
#include <iostream>
#include "convolution.h"
int main(int argc, const char * argv[]) {
    Image<int> img("lena.ascii.pgm");
    double fact=1.0;
    double bia=0.0;
    //Convolution<int> convol(img);
    //Convolution<float> convol(img,fact, bia);
    //Convolution<int> convol(img,edges<int>(9));
    Convolution<float> convol(img,boundaries<float>(3),fact,bia);
    convol.save("output.pgm");
    return 0;
}
