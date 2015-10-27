//
//  Image.h
//  Convolution
//
//  Created by DIOP Abdoulaye  on 20/10/2015.
//  Copyright © 2015 DIOP Abdoulaye . All rights reserved.
//

#ifndef Image_h
#define Image_h
#include <cstdio>
#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <iterator>

namespace std{

    template<typename T>
    std::ostream& operator<< (std::ostream&  os ,const std::vector<T>&  vec)
    {
		using std::begin;
		using std::end;
        std::copy(begin(vec),end(vec), std::ostream_iterator<T>(os," "));
        return os;
    }
    template<typename T>
    std::istream& operator>> (std::istream&  is , std::vector<T>& vec)
    {
        std::copy(std::istream_iterator<T>(is),std::istream_iterator<T>(),std::back_inserter(vec));
        return is;
    }

}
template<typename T>
class Image
{
    // friend class Convolution;
protected:
    int width;
    int height;
    int maxval;
    std::string mNumber;
    std::vector<T> data;
    //std::vector<int> data;
public:
    Image(std::string imgname);
    int index(int i,int j){return i*width+j;}
    friend std::istream& operator>> (std::istream&  is , Image<T>& a);
    friend std::ostream& operator<< (std::ostream&  os , const Image<T>& a);
    int geth()const{return height;}
    int getw()const{return width;}
    int getmv()const{return maxval;}
    std::string getmn()const{return mNumber;}
    T& operator()(int i,int j){ return data[index(i,j)]; }
};

template<typename T>

Image<T>::Image(std::string imgname)
{
    std::ifstream file(imgname, std::ios::in);
    if(file)
    {
        file >>mNumber>>width>>height>>maxval;
        file >>data;
        std::cout<<"Magic Number:"<<mNumber<<std::endl;
        std::cout<<"Width:"<<width<<std::endl;
        std::cout<<"Height:"<<height<<std::endl;
        std::cout<<"Maxval:"<<maxval<<std::endl;
        //std::cout<< "data:\n"<<data<<std::endl;
        file.close();
    }
    else
        std::cerr <<"Imposible d'ouvrir le fichier image.\n"<<std::endl;

}
template<typename T>
std::ostream& operator<< (std::ostream&  os ,const Image<T>& a )
{
    os<<a.mNumber<<std::endl
    <<a->width<<std::endl
    <<a.height<<std::endl
    <<a.data ;
    return os;
}

template<typename T>
std::istream& operator>> (std::istream&  is , Image<T>& a )
{
    is>> a.data;
    return is;
}


#endif /* Image_h */
