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
    Convolution(Image<int> a);
    Convolution(Image<int> a, double fact, double bia);
    Convolution(Image<int> a, std::vector<T> kernelt);
    Convolution(Image<int> a, std::vector<T> kernelt, double fact, double bia);
    Convolution(Image<int> a, std::vector<T> kernelt, struct infop pinfo);
    Convolution(Image<int> a, struct kstruct<T> kernelt);
    Convolution(Image<int> a, std::vector<T> kernelt, double fact, double bia, struct infop pinfo);
    Convolution(Image<int> a, struct kstruct<T> kernelt,struct infop pinfo);
    void save(std::string output,struct infop pinfo);
    /*
    friend std::vector<T>  decalage(int n);
    friend std::vector<T>  sharpen(int n);
    friend std::vector<T>  motionblur(int n);
    std::vector<T>  edges(int n);*/

    ~Convolution(){};
};

int index(int i, int j,int n){return i*n+j;}

/*
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
        kernel[n*n/2]=1;
    }
    return kernel;
}

template<typename T>
std::vector<T> Convolution<T>::edges(int n)
{
    factor=1.0;
    bias=0.0;
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
}*/

/*
template<typename T>
std::vector<T> edges(int n)
{
    //factor=1.0;
    //bias=0.0;
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
*/

/*
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
 */

template<typename T>
struct kstruct<T> sharpen(int level)
{
    struct kstruct<T> kedges;
    kedges.factor=1.0;
    kedges.bias=0.0;
    if (level==2)
    {
        for(int i=0;i<9;i++)
            kedges.k.push_back(1);
        kedges.k[4]=-7;
    }
    if(level==1)
    {
        for(int i=0;i<9;i++)
            kedges.k.push_back(-1);
        kedges.k[4]=9;
    }

    return kedges;
    
}



template<typename T>
struct kstruct<T> motionblur()
{
    struct kstruct<T> kmblur;
    kmblur.factor=1.0/9.0;
    kmblur.bias=0.0;

        for(int i=0;i<9*9;i++)
            kmblur.k.push_back(0);
    
        for(int i=0;i<9*9;i+=9+1)
            kmblur.k[i]=1;
    
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
          for(int i=0;i<9;i++)
                kblur.k.push_back(0.2);
            
            kblur.k[0]=0.0;
            kblur.k[2]=0.0;
            kblur.k[6]=0.0;
            kblur.k[8]=0.0;

        }
        if (level==2)
        {
            kblur.factor=1.0/13.0;
            for(int i=0;i<25;i++)
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
            for(int i=0;i<level*level;i++)
                kblur.k.push_back(1);
            
        }
    
    return kblur;
}
template<typename T>
struct kstruct<T> boundaries()
{
    struct kstruct<T> kboundaries;
    kboundaries.factor=1.0;
    kboundaries.bias=0.0;

      for(int i=0;i<9;i++)
          kboundaries.k.push_back(1);

        kboundaries.k[4]=-4;
        kboundaries.k[0]=0;
        kboundaries.k[1]=0;
        kboundaries.k[6]=0;
        kboundaries.k[8]=0;
    
        return kboundaries;
}


template<typename T>
void Convolution<T>::save(std::string output,struct infop pinfo)
{
    MPI_Barrier(MPI_COMM_WORLD);

    if(pinfo.rank==0)
    {
        int h =cible.geth();
        int w = (cible.getmn()=="P3") ? cible.getw()*3 : cible.getw();
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

                for(int i=0;i<h;i++)
                {
                    for(int j=0;j<w;j++)
                    {
                        file <<cible(i,j);
                        if( i < w-1)
                            file<<" ";
                        //if(j>0 && j%2 == 0)
                            //file<<" ";
                    }
                    file<<std::endl;
                }
                file.close();
            }
            else
                std::cerr <<"Imposible d'écrire dan le fichier image cyble.\n"<<std::endl;
        }
    }
}

/*
template<typename T>
Convolution <T>::Convolution(Image<int> a, double fact, double bia): cible(a), factor(fact), bias(bia)
{
    int n;

    int mv=cible.getmv();
    nkernel=3;
    n=nkernel*nkernel;
    std::vector<T> kernel(n);
    kernel=sharpen<T>(nkernel);
    std::cout<<"Sharpen Filter uses by Default"<<std::endl;
    std::cout<<"Kernel:"<<nkernel<<std::endl;
    for(int i=0;i<nkernel;i++){
        for(int j=0;j<nkernel;j++)
            std::cout<<kernel[index(i,j,nkernel)]<<" ";
        std::cout<<std::endl;

    }
    std::cout<<"Bias:"<<bias<<std::endl;
    std::cout<<"factor:"<<factor<<std::endl;
    int h =cible.geth();
    int w = (cible.getmn()=="P3") ? cible.getw()*3 : cible.getw();
    #pragma omp parallel for
    for(int x = 0; x < h ; x++)
        for(int y = 0; y < w; y++)
        {
            double newval= 0.0;
            for(int kx=-1;kx<2;kx++)
            for(int ky=-1;ky<2;ky++)
            {
                int imagex=(x+kx+1)%h;
                int imagey=(y+ky+1)%w;

                newval+=a(imagex,imagey)*kernel[n/2 + (kx+1) + nkernel*(ky+1)];
            }
            newval=newval*factor +bias;
            if (newval< 0.0) newval= -1*newval;
            //if (newval<0.0) newval=0.0;
            if (newval > mv) newval=mv;
            cible(x,y)=static_cast<int>(newval);
        }

}

template<typename T>

 Convolution <T>::Convolution(Image<int> a): cible(a)
{
 
    int mv=cible.getmv();
    

    
    std::cout<<"Kernel:"<<nkernel<<std::endl;

        
    }
    std::cout<<"mv:"<<mv<<std::endl;
    std::cout<<"Bias:"<<bias<<std::endl;
    std::cout<<"factor:"<<factor<<std::endl;
    std::cout<<"n:"<<n<<std::endl;
    std::cout<<"nkernel:"<<nkernel<<std::endl;
    
    
    int h =cible.geth();
    int w = (cible.getmn()=="P3") ? cible.getw()*3 : cible.getw();
    int stride = (cible.getmn()=="P3") ? 3 : 1;
    int mod = (cible.getmn()=="P3") ? 3 : 0;
#pragma omp parallel for
    for(int x = 0; x < h; x++)
        for(int y = 0; y < w ; y++)
        {
            double newval=0.0;
            
            
            for(int kx=-1*(nkernel/2);kx<=nkernel/2 ;kx++)
                for(int ky=-1*(nkernel/2);ky<=nkernel/2;ky++)
                {
                    
                    
                    int imagex=(x+kx)%h;
                    int imagey=(y+ky*stride)%(w+mod);
                    
                    
                    
                    newval+=static_cast<double>(a(imagex,imagey))*kernelt[n/2 + kx + nkernel*ky];
                }
            newval=factor*newval + bias;
            
            if (newval< 0.0) newval= 0.0;
            if (newval > mv) newval=mv;
            cible(x,y)=static_cast<int>(newval);
            
        }
}*/



/*template<typename T>
Convolution <T>::Convolution(Image<int> a,std::vector<T> kernelt, double fact , double bia):cible(a), factor(fact), bias(bia)
{
    int n;
    int mv=cible.getmv();
    n=kernelt.size();

    nkernel=sqrt(n);

    std::cout<<"Kernel:"<<nkernel<<std::endl;
    for(int i=0;i<nkernel;i++){
        for(int j=0;j<nkernel;j++)
            std::cout<<kernelt[index(i,j,nkernel)]<<" ";
        std::cout<<std::endl;

    }
    std::cout<<"mv:"<<mv<<std::endl;
    std::cout<<"Bias:"<<bias<<std::endl;
    std::cout<<"factor:"<<factor<<std::endl;
    std::cout<<"n:"<<n<<std::endl;
    std::cout<<"nkernel:"<<nkernel<<std::endl;

    
    int h =cible.geth();
    int w = (cible.getmn()=="P3") ? cible.getw()*3 : cible.getw();
    int stride = (cible.getmn()=="P3") ? 3 : 1;
    int mod = (cible.getmn()=="P3") ? 3 : 0;
    #pragma omp parallel for
    for(int x = 0; x < h; x++)
        for(int y = 0; y < w ; y++)
        {
            double newval=0.0;
            
            
            for(int kx=-1*(nkernel/2);kx<=nkernel/2;kx++)
                for(int ky=-1*(nkernel/2);ky<=nkernel/2;ky++)
                {
                    
                    
                    int imagex=(x+kx)%h;
                    int imagey=(y+ky*stride)%(w+mod);
                    
                    
                    
                    newval+=static_cast<double>(a(imagex,imagey))*kernelt[n/2 + kx + nkernel*ky];
                }
            newval=factor*newval + bias;
            
            if (newval< 0.0) newval= 0.0;
            if (newval > mv) newval=mv;
            cible(x,y)=static_cast<int>(newval);
            
        }



        std::cout<<"Convolution end"<<std::endl;

}
*/

/*template<typename T>
Convolution <T>::Convolution(Image<int> a,std::vector<T> kernelt, double fact , double bia,struct infop pinfo ):cible(a), factor(fact), bias(bia)
{
    int n;
    int mv=cible.getmv();
    int tag=20;
    MPI_Status sta;
    MPI_Request req;
    MPI_Request treq[pinfo.nproc-1];
    MPI_Status tsta[pinfo.nproc-1];
    int nbreq=pinfo.nproc-1;

    n=kernelt.size();

    nkernel=sqrt(n);
    std::cout<<"Kernel:"<<nkernel<<std::endl;
    for(int i=0;i<nkernel;i++){
        for(int j=0;j<nkernel;j++)
            std::cout<<kernelt[index(i,j,nkernel)]<<" ";
        std::cout<<std::endl;

    }
    std::cout<<"Bias:"<<bias<<std::endl;
    std::cout<<"factor:"<<factor<<std::endl;


    int h =cible.geth();
    int w = (cible.getmn()=="P3") ? cible.getw()*3 : cible.getw();
    int stride = (cible.getmn()=="P3") ? 3 : 1;
    int mod = (cible.getmn()=="P3") ? 3 : 0 ;
    #pragma omp parallel for
    for(int x = pinfo.ideb; x < pinfo.ifin; x++)
        for(int y = 0; y < w ; y++)
        {
            double newval=0.0;


            for(int kx=-1*(nkernel/2);kx<=nkernel/2;kx++)
                for(int ky=-1*(nkernel/2);ky<=nkernel/2;ky++)
                {


                int imagex=(x+kx)%h;
                int imagey=(y+ky*stride)%(w+mod);



                newval+=static_cast<double>(a(imagex,imagey))*kernelt[n/2 + kx + nkernel*ky];
                }
                newval=factor*newval + bias;
            
                if (newval< 0.0) newval= 0.0;
                if (newval > mv) newval=mv;
                cible(x,y)=static_cast<int>(newval);
            
        }
    if (pinfo.rank > 0){
        MPI_Isend(&cible(pinfo.ideb,0),pinfo.nloc*w,MPI_INT,0,tag,MPI_COMM_WORLD,&req);
        MPI_Wait(&req,&sta);
    }
    if(pinfo.rank == 0){
        for(int i=1;i<pinfo.nproc;i++)
            MPI_Irecv(&cible(i*pinfo.nloc,0),pinfo.nloc*w,MPI_INT,i,tag,MPI_COMM_WORLD,&treq[i-1]);
        MPI_Waitall(nbreq,treq,tsta);
    }
    std::cout<<"Convolution done"<<std::endl;
}*/

template<typename T>
Convolution <T>::Convolution(Image<int> a,struct kstruct<T> kernelt ):cible(a)
{
    int n;
    int mv=cible.getmv();
    n=kernelt.k.size();
    
    nkernel=sqrt(n);
    
    std::cout<<"Kernel:"<<nkernel<<std::endl;
    for(int i=0;i<nkernel;i++){
        for(int j=0;j<nkernel;j++)
            std::cout<<kernelt.k[index(i,j,nkernel)]<<" ";
        std::cout<<std::endl;
        
    }
    factor=kernelt.factor;
    bias=kernelt.bias;
    std::cout<<"mv:"<<mv<<std::endl;
    std::cout<<"Bias:"<<bias<<std::endl;
    std::cout<<"factor:"<<factor<<std::endl;
    std::cout<<"n:"<<n<<std::endl;
    std::cout<<"nkernel:"<<nkernel<<std::endl;
    
    
    int h =cible.geth();
    int w = (cible.getmn()=="P3") ? cible.getw()*3 : cible.getw();
    int stride = (cible.getmn()=="P3") ? 3 : 1;
    int mod = (cible.getmn()=="P3") ? 3 : 0;
#pragma omp parallel for
    for(int x = 0; x < h; x++)
        for(int y = 0; y < w ; y++)
        {
            double newval=0.0;
            
            
            for(int kx=-1*(nkernel/2);kx<=nkernel/2;kx++)
                for(int ky=-1*(nkernel/2);ky<=nkernel/2;ky++)
                {
                    
                    
                    int imagex=(x+kx)%h;
                    int imagey=(y+ky*stride)%(w+mod);
                    
                    
                    
                    newval+=static_cast<double>(a(imagex,imagey))*kernelt.k[n/2 + kx + nkernel*ky];
                }
            newval=factor*newval + bias;
            
            if (newval< 0.0) newval= 0.0;
            if (newval > mv) newval=mv;
            cible(x,y)=static_cast<int>(newval);
            
        }
}

template<typename T>
Convolution <T>::Convolution(Image<int> a, struct kstruct<T> kernelt, struct infop pinfo ):cible(a)
{
    int n;
    int mv=cible.getmv();
    int tag=20;
    MPI_Status sta;
    MPI_Request req;
    MPI_Request treq[pinfo.nproc-1];
    MPI_Status tsta[pinfo.nproc-1];
    int nbreq=pinfo.nproc-1;
    
    n=kernelt.k.size();
    nkernel=sqrt(n);
    
    std::cout<<"Kernel:"<<nkernel<<std::endl;
    for(int i=0;i<nkernel;i++){
        for(int j=0;j<nkernel;j++)
            std::cout<<kernelt.k[index(i,j,nkernel)]<<" ";
        std::cout<<std::endl;
        
    }
    factor=kernelt.factor;
    bias=kernelt.bias;    std::cout<<"Bias:"<<bias<<std::endl;
    std::cout<<"factor:"<<factor<<std::endl;
    
    
    int h =cible.geth();
    int w = (cible.getmn()=="P3") ? cible.getw()*3 : cible.getw();
    int stride = (cible.getmn()=="P3") ? 3 : 1;
    int mod = (cible.getmn()=="P3") ? 3 : 0 ;
#pragma omp parallel for
    for(int x = pinfo.ideb; x < pinfo.ifin; x++)
        for(int y = 0; y < w ; y++)
        {
            double newval=0.0;
            
            
            for(int kx=-1*(nkernel/2);kx<=nkernel/2/*-1*/;kx++)
                for(int ky=-1*(nkernel/2);ky<=nkernel/2/*-1*/;ky++)
                {
                    
                    
                    int imagex=(x+kx)%h;
                    int imagey=(y+ky*stride)%(w+mod);
                    
                    
                    
                    newval+=static_cast<double>(a(imagex,imagey))*kernelt.k[n/2 + kx + nkernel*ky];
                }
            newval=factor*newval + bias;
            
            if (newval< 0.0) newval= 0.0;
            if (newval > mv) newval=mv;
            cible(x,y)=static_cast<int>(newval);
            
        }
    if (pinfo.rank > 0){
        MPI_Isend(&cible(pinfo.ideb,0),pinfo.nloc*w,MPI_INT,0,tag,MPI_COMM_WORLD,&req);
        MPI_Wait(&req,&sta);
    }
    if(pinfo.rank == 0){
        for(int i=1;i<pinfo.nproc;i++)
            MPI_Irecv(&cible(i*pinfo.nloc,0),pinfo.nloc*w,MPI_INT,i,tag,MPI_COMM_WORLD,&treq[i-1]);
        MPI_Waitall(nbreq,treq,tsta);
    }
    std::cout<<"Convolution done"<<std::endl;
}



#endif /* Convolution_h */
