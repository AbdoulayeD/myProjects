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
protected:
    int r;
    int g;
    int b;
public:
    Pixel(int ri, int rg ,int rb): r(ri), g(rg), b(rb){}
    ~Pixel(){};
    
};

#endif /* Pixel_h */
