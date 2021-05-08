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
    for(int i = 0; i < r; i++){
        for(int j = 0; j < c; j++)
            delete[] data[i][j];
        delete[] data[i];
    }

    delete[] data;
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
    if(data){
        for(int i = 0; i < r; i++)
            for(int j = 0; j < c; j++)
                for(int k = 0; k < d; k++)
                    data[i][j][k] = v;
    }else{
        throw(tensor_not_initialized());
    }
}

void Tensor::clamp(float low, float high){
    for(int i = 0; i < r; i++){
        for(int j = 0; j < c; j++){
            for(int k = 0; k < d; k++){
                if(data[i][j][k] < low) data[i][j][k] = low;
                if(data[i][j][k] > high) data[i][j][k] = high;
            }
        }
    }
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
    return r;
}

int Tensor::cols(){
    return c;
}

int Tensor::depth(){
    return d;
}

float Tensor::getMin(int k){
    float res = data[0][0][k];

    for(int i = 0; i < r; i++){
        for(int j = 1; j < c; j++){
            if(data[i][j][k] < res) res = data[i][j][k];
        }
    }

    return res;
}

float Tensor::getMax(int k){
    float res = data[0][0][k];

    for(int i = 0; i < r; i++){
        for(int j = 1; j < c; j++){
            if(data[i][j][k] > res) res = data[i][j][k];
        }
    }

    return res;
}

void Tensor::showSize(){
    cout << r << " x " << c << " x " << d << endl;
}

/* IOSTREAM */

ostream& operator<<(ostream& stream, const Tensor & obj){
    for(int k = 0; k < obj.d; k++){
        for(int i = 0; i < obj.r; i++){
            stream << "[ ";
            for(int j = 0; j < obj.c; j++){
                stream << obj.data[i][j][k] << " ";
            }
            stream << "]\n";
        }
        stream << "\n";
    }
    return stream;
}

void Tensor::read_file(string filename){
    ifstream f{filename};
    if(!f) throw(unable_to_read_file());

    string line;
    f >> r >> c >> d;
    
    allocate_memory();

    for(int k = 0; k < d; k++)
        for(int i = 0; i < r; i++)
            for(int j = 0; j < c; j++)
                f >> data[i][j][k];

    f.close();
}

void Tensor::write_file(string filename){
    ofstream f{filename};
    f << r << "\n" << c << "\n" << d << "\n";

    for(int k = 0; k < d; k++)
        for(int i = 0; i < r; i++)
            for(int j = 0; j < c; j++)
                f << data[i][j][k] << "\n";
    
    f.close();
}
