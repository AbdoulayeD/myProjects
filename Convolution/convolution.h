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
#include <mpi.h>
#include  "image.h"

struct infop{
    int  rank;
    int  nproc;
    int ideb;
    int ifin;
    int nloc;
};
template<typename T>
struct kstruct
{
    std::vector<T> k;
    float factor;
    float bias;

};

template<typename T>
class Convolution
{


protected:
    Image<int> cible;
    int    nkernel;
    double  factor;
    double    bias;
public:
    Convolution(Image<int> a, struct kstruct<T> kernelt);
    Convolution(Image<int> a, struct kstruct<T> kernelt,struct infop pinfo);

    void save(std::string output,struct infop pinfo);

    ~Convolution(){};
};

int index(int i, int j,int n){return i*n+j;}

template<typename T>
struct kstruct<T> sharpen(int level)
{
    struct kstruct<T> kedges;
    kedges.factor=1.0;
    kedges.bias=0.0;
    if (level==2)
    {
        for(auto i=0;i<9;i++)
            kedges.k.push_back(1);
        kedges.k[4]=-7;
    }
    if(level==1)
    {
        for(auto i=0;i<9;i++)
            kedges.k.push_back(-1);
        kedges.k[4]=9;
    }
    std::cout<<"Filter : sharpen"<<std::endl;
    return kedges;

}



template<typename T>
struct kstruct<T> motionblur()
{
    struct kstruct<T> kmblur;
    kmblur.factor=1.0/9.0;
    kmblur.bias=0.0;

        for(auto i=0;i<9*9;i++)
            kmblur.k.push_back(0);

        for(auto i=0;i<9*9;i+=9+1)
            kmblur.k[i]=1;
    std::cout<<"Filter : Motion Blur"<<std::endl;
    return kmblur;
}

template<typename T>
struct kstruct<T> blur(int level)
{
    struct kstruct<T> kblur;
    kblur.factor=1.0;
    kblur.bias=0.0;
        if(level==1)
        {
          for(auto i=0;i<9;i++)
                kblur.k.push_back(0.2);

            kblur.k[0]=0.0;
            kblur.k[2]=0.0;
            kblur.k[6]=0.0;
            kblur.k[8]=0.0;

        }
        if (level==2)
        {
            kblur.factor=1.0/13.0;
            for(auto i=0;i<25;i++)
                kblur.k.push_back(1);

            kblur.k[0]=0;
            kblur.k[1]=0;
            kblur.k[3]=0;
            kblur.k[4]=0;
            kblur.k[5]=0;
            kblur.k[9]=0;
            kblur.k[15]=0;
            kblur.k[19]=0;
            kblur.k[20]=0;
            kblur.k[21]=0;
            kblur.k[23]=0;
            kblur.k[24]=0;
        }
        if (level >2)
        {
            kblur.factor=1.0/(level*level);
            for(auto i=0;i<level*level;i++)
                kblur.k.push_back(1);

        }
     std::cout<<"Filter : Blur"<<std::endl;
    return kblur;
}
template<typename T>
struct kstruct<T> boundaries()
{
    struct kstruct<T> kboundaries;
    kboundaries.factor=1.0;
    kboundaries.bias=0.0;

      for(auto i=0;i<9;i++)
          kboundaries.k.push_back(1);

        kboundaries.k[4]=-4;
        kboundaries.k[0]=0;
        kboundaries.k[1]=0;
        kboundaries.k[6]=0;
        kboundaries.k[8]=0;
         std::cout<<"Filter : Boundaries"<<std::endl;
        return kboundaries;
}


template<typename T>
void Convolution<T>::save(std::string output,struct infop pinfo)
{

        auto h =cible.geth();
        auto w = (cible.getmn()=="P3") ? cible.getw()*3 : cible.getw();

            std::ofstream file(output,std::ios::out);
            if(file)
            {
                file <<cible.getmn()<<std::endl;
                file <<cible.getw()<<" "<<cible.geth()<<std::endl;
                file <<cible.getmv()<<std::endl;

                for(auto i=0;i<h;i++)
                {
                    for(auto j=0;j<w;j++)
                    {
                        file<<static_cast<int>(cible(i,j))<<" ";
                    }
                    file<<std::endl;
                }
                file.close();
            }
            else
                std::cerr <<"Imposible d'écrire dan le fichier image cyble.\n"<<std::endl;


}



template<typename T>
Convolution <T>::Convolution(Image<int> a,struct kstruct<T> kernelt ):cible(a)
{
    auto mv=cible.getmv();
    auto n=kernelt.k.size();

    nkernel=sqrt(n);

    std::cout<<"Kernel:"<<nkernel<<std::endl;
    for(auto i=0;i<nkernel;i++){
        for(auto j=0;j<nkernel;j++)
            std::cout<<kernelt.k[index(i,j,nkernel)]<<" ";
        std::cout<<std::endl;

    }
    factor=kernelt.factor;
    bias=kernelt.bias;
    std::cout<<"mv:"<<mv<<std::endl;
    std::cout<<"Bias:"<<bias<<std::endl;
    std::cout<<"factor:"<<factor<<std::endl;
    std::cout<<std::endl;


    auto h =cible.geth();
    auto w = (cible.getmn()=="P3") ? cible.getw()*3 : cible.getw();
    auto stride = (cible.getmn()=="P3") ? 3 : 1;
    auto mod = (cible.getmn()=="P3") ? 3 : 0;
//#pragma omp parallel for
    for(auto x = 0; x < h; x++)
        for(auto y = 0; y < w ; y++)
        {
            auto newval=0.0;


            for(auto kx=-1*(nkernel/2);kx<=nkernel/2;kx++)
                for(auto ky=-1*(nkernel/2);ky<=nkernel/2;ky++)
                {
                     unsigned int imagex=(x+kx)%h;
                     unsigned int imagey=(y+ky*stride)%(w+mod);
                    newval+=a(imagex,imagey)*kernelt.k[n/2 + kx + nkernel*ky];
                }
            newval=factor*newval + bias;

            if (newval< 0.0) newval= 0.0;
            if (newval > mv) newval=mv;
            cible(x,y)=newval;

        }
}

template<typename T>
Convolution <T>::Convolution(Image<int> a, struct kstruct<T> kernelt, struct infop pinfo ):cible(a)
{
    MPI_Status sta;
    auto mv=cible.getmv();
    auto tag=20;
    auto n=kernelt.k.size();
    nkernel=sqrt(n);
    factor=kernelt.factor;
    bias=kernelt.bias;

    if(pinfo.rank==0){
        std::cout<<"Kernel:"<<nkernel<<std::endl;
        for(auto i=0;i<nkernel;i++){
            for(auto j=0;j<nkernel;j++)
                std::cout<<kernelt.k[index(i,j,nkernel)]<<" ";
            std::cout<<std::endl;

        }
        factor=kernelt.factor;
        bias=kernelt.bias;
        std::cout<<"mv:"<<mv<<std::endl;
        std::cout<<"Bias:"<<bias<<std::endl;
        std::cout<<"factor:"<<factor<<std::endl;
        std::cout<<std::endl;
    }


    auto h =cible.geth();
    auto w = (cible.getmn()=="P3") ? cible.getw()*3 : cible.getw();
    auto stride = (cible.getmn()=="P3") ? 3 : 1;
    auto mod = (cible.getmn()=="P3") ? 3 : 0 ;
    //std::cout<<"Convolution start:"<<pinfo.rank<<std::endl;

    #pragma omp parallel for collapse(2)
        for(auto x = pinfo.ideb; x < pinfo.ifin; x++)
            for(auto y = 0; y < w ; y++)
            {
                auto newval=0.0;
                for(auto kx=-1*(nkernel/2);kx<=nkernel/2;kx++)
                    for(auto ky=-1*(nkernel/2);ky<=nkernel/2;ky++)
                    {
                        unsigned int imagex=(x+kx)%h;
                        unsigned int imagey=(y+ky*stride)%(w+mod);
                        newval+=a(imagex,imagey)*kernelt.k[n/2 + kx + nkernel*ky];
                    }
                newval=factor*newval + bias;

                if (newval< 0.0) newval= 0.0;
                if (newval > mv) newval=mv;
                cible(x,y)=newval;

            }
        //std::cout<<"Convolution stop:"<<pinfo.rank<<std::endl;
        if (pinfo.rank != 0)
            MPI_Send(&cible(pinfo.ideb,0),pinfo.nloc*w,MPI_INT,0,tag,MPI_COMM_WORLD);


        if(pinfo.rank == 0)
            for(auto i=1;i<pinfo.nproc;i++)
                 MPI_Recv(&cible(i*pinfo.nloc,0), pinfo.nloc*w ,MPI_INT,i,tag,MPI_COMM_WORLD,&sta);

    //std::cout<<"send done:"<<pinfo.rank<<std::endl;
    MPI_Barrier(MPI_COMM_WORLD);
}

#endif /* Convolution_h */
