//
//  Pixel.h
//  Convolution
//
//  Created by DIOP Abdoulaye  on 20/10/2015.
//  Copyright © 2015 DIOP Abdoulaye . All rights reserved.
//

#ifndef Pixel_h
#define Pixel_h

class Pixel
{
public:
    int r;
    int g;
    int b;

    Pixel (int ri, int rg ,int rb) { r=ri; g=rg; b=rb;}

    Pixel (const Pixel& a){r=a.r; g=a.g; b=a.b;}

    Pixel operator=( const Pixel& a) { r=a.r; g=a.g; b=a.b;}
    //T& operator(Pixel b)
    ~Pixel(){};

};

#endif /* Pixel_h */
