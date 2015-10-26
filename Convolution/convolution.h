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
#include <cmath>
#include  "image.h"


template<typename T>
class Convolution
{


protected:
    Image<int> cible;
    int    nkernel;
    double  factor;
    double    bias;
public:
    Convolution(Image<int> a);
    Convolution(Image<int> a, double fact, double bia);
    Convolution(Image<int> a, std::vector<T> kernelt);
    Convolution(Image<int> a, std::vector<T> kernelt, double fact, double bia);

    void save(std::string output);
    /*
    friend std::vector<T>  decalage(int n);
    friend std::vector<T>  sharpen(int n);
    friend std::vector<T>  motionblur(int n);
    friend std::vector<T>  edges(int n);*/

    ~Convolution(){};
};

int index(int i, int j,int n){return i*n+j;}

template<typename T>
std::vector<T> decalage(int n)
{
    std::vector<T> kernel(n*n);
    if(n%2==0)
        std::cerr<<"Bad kernel size"<<std::endl;
    else
    {
        for(int i=0;i<n*n;i++)
        {
            kernel[i]=0;
        }
        kernel[n/2]=2;
    }
    return kernel;
}

template<typename T>
std::vector<T> edges(int n)
{
    std::vector<T> kernel(n*n);
    if(n%2==0)
        std::cerr<<"Bad kernel size"<<std::endl;
    else
    {
        for(int i=0;i<n*n;i++)
        {
            kernel[i]=1;
        }
        kernel[n*n/2]=-7;
    }
    return kernel;
}

template<typename T>
std::vector<T> sharpen(int n)
{
    std::vector<T> kernel(n*n);
    if(n%2==0)
        std::cerr<<"Bad kernel size"<<std::endl;
    else
    {
        for(int i=0;i<n*n;i++)
        {
            kernel[i]=-1;
        }
        kernel[n*n/2]=9;
    }
    return kernel;
}

template<typename T>
std::vector<T> motionblur(int n)
{
    std::vector<T> kernel(n*n);
    if(n%2==0)
        std::cerr<<"Bad kernel size"<<std::endl;
    else
    {

        for(int i=0;i<n*n;i++)
        {
            kernel[i]=0;
        }

        for(int i=0;i<n*n;i+=n+1)
        {
            kernel[i]=1;
        }
    }
    return kernel;
}
template<typename T>
std::vector<T> blur(int n)
{
    std::vector<T> kernel(n*n);
    if(n%2==0)
        std::cerr<<"Bad kernel size"<<std::endl;
    else
    {
        if(n==3)
        {
            for(int i=0;i<n*n;i++)
            {
                kernel[i]=0.2;
            }
            kernel[0]=0.0;
            kernel[1]=0.0;
            kernel[6]=0.0;
            kernel[8]=0.0;
        }
        else if(n==5)
        {
            for(int i=0;i<25;i++)
            {
                kernel[i]=1;
            }
            kernel[0]=0;
            kernel[1]=0;
            kernel[3]=0;
            kernel[4]=0;
            kernel[5]=0;
            kernel[9]=0;
            kernel[15]=0;
            kernel[19]=0;
            kernel[20]=0;
            kernel[21]=0;
            kernel[23]=0;
            kernel[24]=0;
        }
        else
        {
            for(int i=0;i<n*n;i++)
            {
                kernel[i]=1;
            }
        }
    }
    return kernel;
}
template<typename T>
std::vector<T>boundaries(int n)
{
    std::vector<T> kernel(n*n);
    if(n%2==0)
        std::cerr<<"Bad kernel size"<<std::endl;
    else if(n>3)
        std::cerr<<"Bad kernel size, chose size(3) for this boundaries filter"<<std::endl;
    else
    {
      for(int i=0;i<n*n;i++)
      {
        kernel[i]=1;
      }
        kernel[4]=-4;
        kernel[0]=0;
        kernel[1]=0;
        kernel[6]=0;
        kernel[8]=0;
    }
        return kernel;
}


template<typename T>
void Convolution<T>::save(std::string output)
{
    if(cible.getw()==0){
        std::cerr<<"None operations where made in the images, so no need to be saved"<<std::endl;
    }
    else{
        std::ofstream file(output,std::ios::out);
        if(file)
        {
            file <<cible.getmn()<<std::endl;
            file <<cible.getw()<<" "<<cible.geth()<<std::endl;
            file <<cible.getmv()<<std::endl;
            for(int i=0;i<cible.geth();i++){
                for(int j=0;j<cible.getw();j++)
                    file <<cible(i,j)<<" ";
                file<<std::endl;
            }
            file.close();
        }
        else
            std::cerr <<"Imposible d'écrire dan le fichier image cyble.\n"<<std::endl;
    }
}


template<typename T>
Convolution <T>::Convolution(Image<int> a, double fact, double bia): cible(a), factor(fact), bias(bia)
{
    int n;
    int h=cible.geth();
    int w=cible.getw();
    int mv=cible.getmv();
    nkernel=3;
    n=nkernel*nkernel;
    std::vector<T> kernel(n);
    kernel=sharpen<T>(nkernel);
    std::cout<<"Sharpen Filter uses by Default"<<std::endl;
    std::cout<<std::endl<<"---->Try test:"<<std::endl;
    std::cout<<"Kernel:"<<nkernel<<std::endl;
    for(int i=0;i<nkernel;i++){
        for(int j=0;j<nkernel;j++)
            std::cout<<kernel[index(i,j,nkernel)]<<" ";
        std::cout<<std::endl;

    }
    std::cout<<"Bias:"<<bias<<std::endl;
    std::cout<<"factor:"<<factor<<std::endl;

    std::cout<<"Image:"<<std::endl;
    std::cout<<"w:"<<w<<std::endl;
    std::cout<<"h:"<<h<<std::endl;
    std::cout<<"mv:"<<mv<<std::endl;
    for(int i=0;i<h;i++){
        for(int j=0;j<w;j++)
            std::cout<<a(i,j)<<" ";
        std::cout<<std::endl;

    }
    #pragma omp parallel for
    for(int x = 0; x < h ; x++)
        for(int y = 0; y < w; y++)
        {
            double newval= 0.0;
            for(int kx=-1;kx<2;kx++)
            for(int ky=-1;ky<2;ky++)
            {
                int imagex=(x+kx)%h;
                int imagey=(y+ky)%w;

                newval+=a(imagex,imagey)*kernel[n/2 + kx + nkernel*ky];
            }
            newval=newval*factor +bias;

            if (newval<0.0) newval=0.0;
            if (newval > mv) newval=mv;
            cible(x,y)=static_cast<T>(newval);
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

}
template<typename T>
Convolution <T>::Convolution(Image<int> a): cible(a),nkernel(3),factor(1.0),bias(0.0)
{
    int n;
    int h=cible.geth();
    int w=cible.getw();
    int mv=cible.getmv();
    n=nkernel*nkernel;
    std::vector<T> kernel(n);
    kernel=sharpen<T>(nkernel);
    std::cout<<"Sharpen Filter uses by Default, with factor = 1.0 and Bias = 0.0"<<std::endl;

    std::cout<<std::endl<<"---->Try test:"<<std::endl;
    std::cout<<"Kernel:"<<nkernel<<std::endl;
    for(int i=0;i<nkernel;i++){
        for(int j=0;j<nkernel;j++)
            std::cout<<kernel[index(i,j,nkernel)]<<" ";
        std::cout<<std::endl;

    }
    std::cout<<"Bias:"<<bias<<std::endl;
    std::cout<<"factor:"<<factor<<std::endl;

    std::cout<<"Image:"<<std::endl;
    std::cout<<"w:"<<w<<std::endl;
    std::cout<<"h:"<<h<<std::endl;
    std::cout<<"mv:"<<mv<<std::endl;
    for(int i=0;i<h;i++){
        for(int j=0;j<w;j++)
            std::cout<<a(i,j)<<" ";
        std::cout<<std::endl;

    }
    #pragma omp parallel for
    for(int x = 0; x < h ; x++)
        for(int y = 0; y < w; y++)
        {
            double newval=0;
            for(int kx=-1;kx<2;kx++)
            for(int ky=-1;ky<2;ky++)
            {
                int imagex=(x+kx)%h;
                int imagey=(y+ky)%w;

                newval+=a(imagex,imagey)*kernel[n/2 + kx + nkernel*ky];
            }
            newval=newval*factor +bias;


            /*newval = (newval > 0) ? newval: std::fabs(cible(x,y));
            cible(x,y) = (newval > mv) ?  mv : static_cast<T>(newval);*/
            if (newval<0.0) newval=0.0;
            if (newval > mv) newval=mv;
            cible(x,y)=static_cast<T>(newval);
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

}

template<typename T>
Convolution <T>::Convolution(Image<int> a,std::vector<T> kernelt):cible(a), factor(1.0), bias(0.0)
{
    int n;
    int h=cible.geth();
    int w=cible.getw();
    int mv=cible.getmv();
    n=kernelt.size();
    nkernel=sqrt(n);
    std::cout<<"Default factor = 1.0 and Bias = 0.0"<<std::endl;
    std::cout<<std::endl<<"---->Try test:"<<std::endl;
    std::cout<<"Kernel:"<<nkernel<<std::endl;
    for(int i=0;i<nkernel;i++){
        for(int j=0;j<nkernel;j++)
            std::cout<<kernelt[index(i,j,nkernel)]<<" ";
        std::cout<<std::endl;

    }
    std::cout<<"Bias:"<<bias<<std::endl;
    std::cout<<"factor:"<<factor<<std::endl;

    std::cout<<"Image:"<<std::endl;
    std::cout<<"w:"<<w<<std::endl;
    std::cout<<"h:"<<h<<std::endl;
    std::cout<<"mv:"<<mv<<std::endl;
    for(int i=0;i<h;i++){
        for(int j=0;j<w;j++)
            std::cout<<a(i,j)<<" ";
        std::cout<<std::endl;

    }
    #pragma omp parallel for
    for(int x = 0; x < h ; x++)
        for(int y = 0; y < w; y++)
        {
            double newval=0;
            for(int kx=-1;kx<2;kx++)
            for(int ky=-1;ky<2;ky++)
            {
                int imagex=(x+kx)%h;
                int imagey=(y+ky)%w;

                newval+=a(imagex,imagey)*kernelt[n/2 + kx + nkernel*ky];
            }
            newval=newval*factor +bias;

            /*newval = (newval > 0) ? newval: std::fabs(cible(x,y));
            cible(x,y) = (newval > mv) ?  mv : static_cast<T>(newval) ;*/
              if (newval<0.0) newval=0.0;
            if (newval > mv) newval=mv;
            cible(x,y)=static_cast<T>(newval);
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

}

template<typename T>
Convolution <T>::Convolution(Image<int> a,std::vector<T> kernelt, double fact , double bia):cible(a), factor(fact), bias(bia)
{
    int n;
    int h=cible.geth();
    int w=cible.getw();
    int mv=cible.getmv();
    n=kernelt.size();
    nkernel=sqrt(n);
    std::cout<<std::endl<<"---->Try test:"<<std::endl;
    std::cout<<"Kernel:"<<nkernel<<std::endl;
    for(int i=0;i<nkernel;i++){
        for(int j=0;j<nkernel;j++)
            std::cout<<kernelt[index(i,j,nkernel)]<<" ";
        std::cout<<std::endl;

    }
    std::cout<<"Bias:"<<bias<<std::endl;
    std::cout<<"factor:"<<factor<<std::endl;

    std::cout<<"Image:"<<std::endl;
    std::cout<<"w:"<<w<<std::endl;
    std::cout<<"h:"<<h<<std::endl;
    std::cout<<"mv:"<<mv<<std::endl;
    for(int i=0;i<h;i++){
        for(int j=0;j<w;j++)
            std::cout<<a(i,j)<<" ";
        std::cout<<std::endl;

    }
    #pragma omp parallel for
    for(int x = 0; x < h ; x++)
        for(int y = 0; y < w; y++)
        {
            double newval=0;
            for(int kx=-1;kx<2;kx++)
            for(int ky=-1;ky<2;ky++)
            {
                int imagex=(x+kx)%h;
                int imagey=(y+ky)%w;

                newval+=a(imagex,imagey)*kernelt[n/2 + kx + nkernel*ky];
            }
            newval=newval*factor +bias;

            /*newval = (newval > 0) ? newval: std::fabs(cible(x,y));
            cible(x,y) = (newval > mv) ?  mv : static_cast<T>(newval);*/
             if (newval<0.0) newval=0.0;
            if (newval > mv) newval=mv;
            cible(x,y)=static_cast<T>(newval);
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

}



#endif /* Convolution_h */
