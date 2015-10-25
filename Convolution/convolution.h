//
//  Convolution.h
//  Convolution
//
//  Created by DIOP Abdoulaye  on 20/10/2015.
//  Copyright © 2015 DIOP Abdoulaye . All rights reserved.
//

#ifndef Convolution_h
#define Convolution_h
#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <iterator>
#include  "image.h"


template<typename T>
class Convolution
{
    
    
protected:
    Image<T> cible;
    int     nkernel;
public:
    Convolution(Image<T> a);//: cible(a){}
    void save(std::string);
    ~Convolution(){};
};

int index(int i, int j,int n){return i*n+j;}

template<typename T>
std::vector<T> decalage()
{
    std::vector<T> kernel(9);
    kernel[index(0,0,3)]=0;
    kernel[index(0,1,3)]=1;
    kernel[index(0,2,3)]=0;
    
    kernel[index(1,0,3)]=0;
    kernel[index(1,1,3)]=0;
    kernel[index(1,2,3)]=0;
    
    kernel[index(2,0,3)]=0;
    kernel[index(2,1,3)]=0;
    kernel[index(2,2,3)]=0;
    
    
    return kernel;
}


template<typename T>
std::vector<T> sharpen()
{
            std::vector<T> kernel(9);
            kernel[index(0,0,3)]=-1;
            kernel[index(0,1,3)]=-1;
            kernel[index(0,2,3)]=-1;
    
            kernel[index(1,0,3)]=-1;
            kernel[index(1,1,3)]=9;
            kernel[index(1,2,3)]=-1;
    
            kernel[index(2,0,3)]=-1;
            kernel[index(2,1,3)]=-1;
            kernel[index(2,2,3)]=-1;
    
    
    return kernel;
}
template<typename T>
std::vector<T> motionblur()
{
    std::vector<T> kernel(9);
    kernel[index(0,0,3)]=1;
    kernel[index(0,1,3)]=0;
    kernel[index(0,2,3)]=0;
    
    kernel[index(1,0,3)]=0;
    kernel[index(1,1,3)]=1;
    kernel[index(1,2,3)]=0;
    
    kernel[index(2,0,3)]=0;
    kernel[index(2,1,3)]=0;
    kernel[index(2,2,3)]=1;
    
    
    return kernel;
}
template<typename T>
std::vector<T> blur(int a)
{
    std::vector<T> kernel(a);
    
    return kernel;
}



template<typename T>
Convolution <T>::Convolution(Image<T> a): cible(a)
{
    int h=cible.geth();
    int w=cible.getw();
    int mv=cible.getmv();
    nkernel=3;
    std::vector<T> kernel(nkernel*nkernel);
    kernel=motionblur<T>();
    
    std::cout<<std::endl<<"---->Try test:"<<std::endl;
    std::cout<<"Kernel:"<<nkernel<<std::endl;
    for(int i=0;i<nkernel;i++){
        for(int j=0;j<nkernel;j++)
            std::cout<<kernel[index(i,j,nkernel)]<<" ";
        std::cout<<std::endl;
        
    }
    std::cout<<"Image:"<<std::endl;
    std::cout<<"w:"<<w<<std::endl;
    std::cout<<"h:"<<h<<std::endl;
    std::cout<<"mv:"<<mv<<std::endl;
    for(int i=0;i<h;i++){
        for(int j=0;j<w;j++)
            std::cout<<a(i,j)<<" ";
        std::cout<<std::endl;
        
    }
    
    /*for(int x = 0; x < w; x++)
        for(int y = 0; y < h; y++)
        {
            T newnc;
            
            for(int filterX = 0; filterX <nkernel; filterX++)
                for(int filterY = 0;filterY<nkernel; filterY++)
                {
                    int imageX = (x-nkernel /2 + filterX + w) % w;
                    int imageY = (y-nkernel /2 + filterY + h) % h;
                    newnc += a(imageX,imageY)* kernel[index(filterX,filterY,nkernel)];

                }
            
            cible(x,y) = (newnc > mv) ?  mv : newnc ;
            
          
        }*/
    //int x=0;
    //int y=0;
            for(int x = 0; x < h; x++)
                for(int y = 0; y < w; y++)
                {
                	if(x==0 || (x == h-1) || y==0 || (y == w-1))
                    {   }
                    else
                    {
                        cible(x,y)=a(x,y)*kernel[4] + a(x-1,y)*kernel[3] + a(x+1,y)*kernel[5] + a(x,y+1)*kernel[7] + a(x,y-1)*kernel[1]+
                                        a(x-1,y-1)*kernel[0] + a(x-1,y+1)*kernel[6] + a(x+1,y-1)*kernel[2] + a(x+1,y+1)*kernel[8];
                    
                        cible(x,y) = (cible(x,y) > mv) ?  mv : cible(x,y);
                    }
                }
    
    
    std::cout<<std::endl<<"---->Try test 2:"<<std::endl;
    std::cout<<"Image:"<<std::endl;
    std::cout<<"w:"<<w<<std::endl;
    std::cout<<"h:"<<h<<std::endl;
    for(int i=0;i<h;i++){
        for(int j=0;j<w;j++)
            std::cout<<cible(i,j)<<" ";
        std::cout<<std::endl;
    }
    std::string output="output.pgm";
    
    std::ofstream file(output,std::ios::out);
    if(file)
    {
        file <<cible.getmn()<<std::endl;
        file <<cible.getw()<<" "<<cible.geth()<<std::endl;
        file <<cible.getmv()<<std::endl;
        for(int i=0;i<h;i++){
            for(int j=0;j<w;j++)
                file <<cible(i,j)<<" ";
            file<<std::endl;
        
        }
        
        file.close();
    }
    //std::ifstream file(imgname, std::ios::in);
        

    
}



template<typename T>
void Convolution<T>::save(std::string a)
{
    
}




#endif /* Convolution_h */
