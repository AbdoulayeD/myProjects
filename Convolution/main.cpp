//
//  main.cpp
//  Convolution
//
//  Created by DIOP Abdoulaye  on 20/10/2015.
//  Copyright © 2015 DIOP Abdoulaye . All rights reserved.
//


#include <iostream>
#include "convolution.h"
int main(int argc, const char * argv[]) {
    Image<int> img("lena.ascii.pgm");
    Convolution<int> convol(img);
    return 0;
}
