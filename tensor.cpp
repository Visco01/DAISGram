#include <iostream>
#include <string>
#include <random>
#include <math.h>
#include <fstream>

#include "dais_exc.h"
#include "tensor.h"

#define PI 3.141592654
#define FLT_MAX 3.402823466e+38F /* max value */
#define FLT_MIN 1.175494351e-38F /* min positive value */

using namespace std;

Tensor::Tensor(){
    data = nullptr;
    r = 0;
    c = 0;
    d = 0;
}

void Tensor::allocate_memory(){
    data = new float**[r];

    for(int i = 0; i < r; i++){
        data[i] = new float*[c];
        for(int j = 0; j < c; j++){
            data[i][j] = new float[d];
        }
    }
}

Tensor::Tensor(int r, int c, int d, float v) {
    this->r = r;
    this->c = c;
    this->d = d;

    allocate_memory();

    for(int i = 0; i < r; i++)
        for(int j = 0; j < c; j++)
            for(int k = 0; k < d; k++)
                data[i][j][k] = v;
}

Tensor::~Tensor(){
    for(int i = 0; i < r; i++)
        for(int j = 0; j < c; j++)
            delete data[i][j];
}

void Tensor::init_progressive(){

}

float Tensor::operator()(int i, int j, int k) const{
    if((i < 0 || i >= r) || (j < 0 || j >= c) || (k < 0 || k >= d)) throw(index_out_of_bound());
    else return data[i][j][k];
}

float& Tensor::operator()(int i, int j, int k){
    if((i < 0 || i >= r) || (j < 0 || j >= c) || (k < 0 || k >= d)) throw(index_out_of_bound());
    else return data[i][j][k];
}

Tensor::Tensor(const Tensor& that){
    this->r = that.r;
    this->c = that.c;
    this->d = that.d;

    allocate_memory();

    for(int i = 0; i < r; i++)
        for(int j = 0; j < c; j++)
            for(int k = 0; k < d; k++)
                data[i][j][k] = that(i, j, k);
}

//OPERATORI

Tensor Tensor::operator-(const Tensor &rhs){
    if(r != rhs.r || c != rhs.c || d != rhs.d) throw(dimension_mismatch());
    Tensor res{r, c, d};

    for(int i = 0; i < r; i++)
        for(int j = 0; j < c; j++)
            for(int k = 0; k < d; k++)
                res(i, j, k) = data[i][j][k] - rhs(i, j, k);
    
    return res;
}

Tensor Tensor::operator+(const Tensor &rhs){
    if(r != rhs.r || c != rhs.c || d != rhs.d) throw(dimension_mismatch());
    Tensor res{r, c, d};

    for(int i = 0; i < r; i++)
        for(int j = 0; j < c; j++)
            for(int k = 0; k < d; k++)
                res(i, j, k) = data[i][j][k] + rhs(i, j, k);
    
    return res;
}

Tensor Tensor::operator*(const Tensor &rhs){
    if(r != rhs.r || c != rhs.c || d != rhs.d) throw(dimension_mismatch());
    Tensor res{r, c, d};

    for(int i = 0; i < r; i++)
        for(int j = 0; j < c; j++)
            for(int k = 0; k < d; k++)
                res(i, j, k) = data[i][j][k] * rhs(i, j, k);
    
    return res;
}

Tensor Tensor::operator/(const Tensor &rhs){
    if(r != rhs.r || c != rhs.c || d != rhs.d) throw(dimension_mismatch());
    Tensor res{r, c, d};

    for(int i = 0; i < r; i++)
        for(int j = 0; j < c; j++)
            for(int k = 0; k < d; k++)
                res(i, j, k) = data[i][j][k] / rhs(i, j, k);
    
    return res;
}

Tensor Tensor::operator-(const float &rhs){
    Tensor res{r, c, d};

    for(int i = 0; i < r; i++)
        for(int j = 0; j < c; j++)
            for(int k = 0; k < d; k++)
                res(i, j, k) = data[i][j][k] - rhs;
    
    return res;
}

Tensor Tensor::operator+(const float &rhs){
    Tensor res{r, c, d};

    for(int i = 0; i < r; i++)
        for(int j = 0; j < c; j++)
            for(int k = 0; k < d; k++)
                res(i, j, k) = data[i][j][k] + rhs;
    
    return res;
}

Tensor Tensor::operator*(const float &rhs){
    Tensor res{r, c, d};

    for(int i = 0; i < r; i++)
        for(int j = 0; j < c; j++)
            for(int k = 0; k < d; k++)
                res(i, j, k) = data[i][j][k] * rhs;
    
    return res;
}

Tensor Tensor::operator/(const float &rhs){
    Tensor res{r, c, d};

    for(int i = 0; i < r; i++)
        for(int j = 0; j < c; j++)
            for(int k = 0; k < d; k++)
                res(i, j, k) = data[i][j][k] / rhs;
    
    return res;
}

Tensor & Tensor::operator=(const Tensor &other){
    this->r = other.r;
    this->c = other.c;
    this->d = other.d;

    allocate_memory();

    for(int i = 0; i < r; i++)
        for(int j = 0; j < c; j++)
            for(int k = 0; k < d; k++)
                data[i][j][k] = other(i, j, k);
    
    return *this;
}

//OPERAZIONI

/**
 * Random Initialization
 *
 * Perform a random initialization of the tensor
 *
 * @param mean The mean
 * @param std  Standard deviation
 */
void Tensor::init_random(float mean, float std){
    if(data){

        std::default_random_engine generator;
        std::normal_distribution<float> distribution(mean,std);

        for(int i=0;i<r;i++){
            for(int j=0;j<c;j++){
                for(int k=0;k<d;k++){
                    this->operator()(i,j,k)= distribution(generator);
                }
            }
        }

    }else{
        throw(tensor_not_initialized());
    }
}

void Tensor::init(int r, int c, int d, float v){

}

void Tensor::clamp(float low, float high){

}

void Tensor::rescale(float new_max){

}

Tensor Tensor::padding(int pad_h, int pad_w){

}

Tensor Tensor::subset(unsigned int row_start, unsigned int row_end, unsigned int col_start, unsigned int col_end, unsigned int depth_start, unsigned int depth_end){

}

Tensor Tensor::concat(const Tensor &rhs, int axis){

}

Tensor Tensor::convolve(const Tensor &f){

}

/* UTILITY */

int Tensor::rows(){

}

int Tensor::cols(){
    
}

int Tensor::depth(){
    
}

float Tensor::getMin(int k){

}

float Tensor::getMax(int k){

}

void Tensor::showSize(){
    
}

/* IOSTREAM */

ostream& operator<<(ostream& stream, const Tensor & obj){

}

void Tensor::read_file(string filename){

}

void Tensor::write_file(string filename){

}
