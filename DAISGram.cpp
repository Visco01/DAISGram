#include <iostream>
#include <string>

#include "dais_exc.h"
#include "tensor.h"
#include "libbmp.h"
#include "DAISGram.h"

using namespace std;

DAISGram::DAISGram(){

}

DAISGram::~DAISGram(){

}

/**
 * Load a bitmap from file
 *
 * @param filename String containing the path of the file
 */
void DAISGram::load_image(string filename){
    BmpImg img = BmpImg();

    img.read(filename.c_str());

    const int h = img.get_height();
    const int w = img.get_width();

    data = Tensor(h, w, 3, 0.0);

    for(int i=0;i<img.get_height();i++){
        for(int j=0;j<img.get_width();j++){
            data(i,j,0) = (float) img.red_at(j,i);
            data(i,j,1) = (float) img.green_at(j,i);
            data(i,j,2) = (float) img.blue_at(j,i);
        }
    }
}


/**
 * Save a DAISGram object to a bitmap file.
 *
 * Data is clamped to 0,255 before saving it.
 *
 * @param filename String containing the path where to store the image.
 */
void DAISGram::save_image(string filename){

    data.clamp(0,255);

    BmpImg img = BmpImg(getCols(), getRows());

    img.init(getCols(), getRows());

    for(int i=0;i<getRows();i++){
        for(int j=0;j<getCols();j++){
            img.set_pixel(j,i,(unsigned char) data(i,j,0),(unsigned char) data(i,j,1),(unsigned char) data(i,j,2));
        }
    }

    img.write(filename);

}

/**
 * Get rows
 *
 * @return returns the number of rows in the image
 */
int DAISGram::getRows(){
    return data.rows();
}

/**
 * Get columns
 *
 * @return returns the number of columns in the image
 */
int DAISGram::getCols(){
    return data.cols();
}

/**
 * Get depth
 *
 * @return returns the number of channels in the image
 */
int DAISGram::getDepth(){
    return data.depth();
}

/**
 * Brighten the image
 *
 * It sums the bright variable to all the values in the image.
 *
 * Before returning the image, the corresponding tensor should be clamped in [0,255]
 *
 * @param bright the amount of bright to add (if negative the image gets darker)
 * @return returns a new DAISGram containing the modified object
 */
DAISGram DAISGram::brighten(float bright){
    DAISGram res{};
    res.data = data + bright;
    res.data.clamp(0, 255);
    return res;
}

/**
 * Create a grayscale version of the object
 *
 * A grayscale image is produced by substituting each pixel with its average on all the channel
 *
 * @return returns a new DAISGram containing the modified object
 */
DAISGram DAISGram::grayscale(){
    DAISGram res = *this;
    int average = 0;
    for(int i = 0; i < data.rows(); i++){
        for(int j = 0; j < data.cols(); j++){
            for(int k = 0; k < data.depth(); k++){
                average += res.data(i, j, k);
            }
            for(int k = 0; k < data.depth(); k++){
                res.data(i, j, k) = average / data.depth();
            }
            average = 0;
        }
    }
        
    return res;
}

/**
 * Create a Warhol effect on the image
 *
 * This function returns a composition of 4 different images in which the:
 * - top left is the original image
 * - top right is the original image in which the Red and Green channel are swapped
 * - bottom left is the original image in which the Blue and Green channel are swapped
 * - bottom right is the original image in which the Red and Blue channel are swapped
 *
 * The output image is twice the dimensions of the original one.
 *
 * @return returns a new DAISGram containing the modified object
 */
DAISGram DAISGram::warhol(){
    DAISGram original = *this;
    DAISGram green = *this;
    DAISGram red = *this;
    DAISGram blue = *this;

    //GREEN <-> RED
    for(int i = 0; i < data.rows(); i++){
        for(int j = 0; j < data.cols(); j++){
            float supp;

            //GREEN <-> RED
            supp = green.data(i, j, 0);
            green.data(i, j, 0) = green.data(i, j, 1);
            green.data(i, j, 1) = supp;

            //GREEN <-> BLUE
            supp = red.data(i, j, 1);
            red.data(i, j, 1) = red.data(i, j, 2);
            red.data(i, j, 2) = supp;

            //RED <-> BLUE
            supp = blue.data(i, j, 0);
            blue.data(i, j, 0) = blue.data(i, j, 2);
            blue.data(i, j, 2) = supp;
        }
    }

    DAISGram original_green, red_blue, res;
    original_green.data = original.data.concat(green.data, 1);
    red_blue.data = red.data.concat(blue.data, 1);
    res.data =  original_green.data.concat(red_blue.data, 0);

    return res;
}

/**
 * Sharpen the image
 *
 * This function makes the image sharper by convolving it with a sharp filter
 *
 * filter[3][3]
 *    0  -1  0
 *    -1  5 -1
 *    0  -1  0
 *
 * Before returning the image, the corresponding tensor should be clamped in [0,255]
 *
 * @return returns a new DAISGram containing the modified object
 */
DAISGram DAISGram::sharpen(){
    Tensor filter{3, 3, 3};
    
    for(int i = 0; i < 3; i++){
        for(int j = 0; j < 3; j++){
            for(int k = 0; k < 3; k++){
                if(i == 1 && j == 1){
                    filter(i, j, k) = 5;
                }else if(i % 2 != 0 || j % 2 != 0){
                    filter(i, j, k) = -1;
                }else{
                    filter(i, j, k) = 0;
                }
            }
        }
    }

    DAISGram res{};
    res.data = data.convolve(filter);
    res.data.clamp(0, 255);
    return res;
}

/**
 * Emboss the image
 *
 * This function makes the image embossed (a light 3D effect) by convolving it with an
 * embossing filter
 *
 * filter[3][3]
 *    -2 -1  0
 *    -1  1  1
 *     0  1  2
 *
 * Before returning the image, the corresponding tensor should be clamped in [0,255]
 *
 * @return returns a new DAISGram containing the modified object
 */
DAISGram DAISGram::emboss(){
    Tensor filter{3, 3, 3};

    int z = -3;
    for(int i = 0; i < 3; i++){
        z++;
        int temp = z-1;
        for(int j = 0; j < 3; j++){
            temp++;
            for(int k = 0; k < 3; k++){
                if(i == 1 && j == 1){
                    filter(i, j, k) = 1;
                }else{
                    filter(i, j, k) = temp;
                }
            }
        }
    }

    DAISGram temp = *this;
    DAISGram res;
    res.data = temp.data.convolve(filter);
    res.data.clamp(0, 255);
    return res;
}

/**
 * Smooth the image
 *
 * This function remove the noise in an image using convolution and an average filter
 * of size h*h:
 *
 * c = 1/(h*h)
 *
 * filter[3][3]
 *    c c c
 *    c c c
 *    c c c
 *
 * @param h the size of the filter
 * @return returns a new DAISGram containing the modified object
 */
DAISGram DAISGram::smooth(int h){
    float c = 1.0f / (h * h);
    Tensor filter{h, h, 3, c};

    DAISGram temp = *this;
    DAISGram res;
    res.data = temp.data.convolve(filter);
    return res;
}

/**
 * Edges of an image
 *
 * This function extract the edges of an image by using the convolution
 * operator and the following filter
 *
 *
 * filter[3][3]
 * -1  -1  -1
 * -1   8  -1
 * -1  -1  -1
 *
 * Remeber to convert the image to grayscale before running the convolution.
 *
 * Before returning the image, the corresponding tensor should be clamped in [0,255]
 *
 * @return returns a new DAISGram containing the modified object
 */
DAISGram DAISGram::edge(){
    Tensor filter{3, 3, 3};
    
    for(int i = 0; i < 3; i++){
        for(int j = 0; j < 3; j++){
            for(int k = 0; k < 3; k++){
                if(i == 1 && j == 1){
                    filter(i, j, k) = 8;
                }else{
                    filter(i, j, k) = -1;
                }
            }
        }
    }
    
    DAISGram temp = *this;
    temp = temp.grayscale();
    DAISGram res;
    res.data = temp.data.convolve(filter);
    return res;
}

/**
 * Blend with anoter image
 *
 * This function generate a new DAISGram which is the composition
 * of the object and another DAISGram object
 *
 * The composition follows this convex combination:
 * results = alpha*this + (1-alpha)*rhs
 *
 * rhs and this obejct MUST have the same dimensions.
 *
 * @param rhs The second image involved in the blending
 * @param alpha The parameter of the convex combination
 * @return returns a new DAISGram containing the blending of the two images.
 */
DAISGram DAISGram::blend(const DAISGram & rhs, float alpha){
    DAISGram check = rhs;
    if(getRows() != check.getRows() || getCols() != check.getCols() || getDepth() != check.getDepth()) throw(dimension_mismatch());
    DAISGram res;
    res.data = data * alpha + check.data * (1 - alpha);
    return res;
}


/**
 * Generate Random Image
 *
 * Generate a random image from nois
 *
 * @param h height of the image
 * @param w width of the image
 * @param d number of channels
 * @return returns a new DAISGram containing the generated image.
 */
void DAISGram::generate_random(int h, int w, int d){
    data = Tensor(h,w,d,0.0);
    data.init_random(128,50);
    data.rescale(255);
}

//PARTI OPZIONALI

/**
 * Green Screen
 *
 * This function substitutes a pixel with the corresponding one in a background image
 * if its colors are in the surrounding (+- threshold) of a given color (rgb).
 *
 * (rgb - threshold) <= pixel <= (rgb + threshold)
 *
 *
 * @param bkg The second image used as background
 * @param rgb[] The color to substitute (rgb[0] = RED, rgb[1]=GREEN, rgb[2]=BLUE)
 * @param threshold[] The threshold to add/remove for each color (threshold[0] = RED, threshold[1]=GREEN, threshold[2]=BLUE)
 * @return returns a new DAISGram containing the result.
 */
DAISGram DAISGram::greenscreen(DAISGram & bkg, int rgb[], float threshold[]){
    if(data.rows() != bkg.data.rows() || data.cols() != bkg.data.cols() ||data.depth() != bkg.data.depth()) throw(dimension_mismatch());
    DAISGram res = *this;
    int counter = 0;
    
    for(int i = 0; i < res.getRows(); i++){
        for(int j = 0; j < res.getCols(); j++){
            for(int k = 0; k < res.getDepth(); k++){
                if(res.data(i, j, k) >= (rgb[k] - threshold[k]) && res.data(i, j, k) <= (rgb[k] + threshold[k])){
                    counter++;
                }
            }

            if(counter == 3){
                res.data(i, j, 0) = bkg.data(i, j, 0);
                res.data(i, j, 1) = bkg.data(i, j, 1);
                res.data(i, j, 2) = bkg.data(i, j, 2);
            }
            counter = 0;
        }
    }

    return res;
}

/**
 * Equalize
 *
 * Stretch the distribution of colors of the image in order to use the full range of intesities.
 *
 * See https://it.wikipedia.org/wiki/Equalizzazione_dell%27istogramma
 *
 * @return returns a new DAISGram containing the equalized image.
 */
DAISGram DAISGram::equalize(){
    DAISGram res = *this;
    Tensor occurrencies{16, 16, 3};
    int rows = res.getRows();
    int cols = res.getCols();

    int a = 0, b = 0;

    for(int z = res.data.getMin(0); z < 256; z++){
        occurrencies(a, b, 0) = z;
        for(int i = 0; i < rows; i++){
            for(int j = 0; j < cols; j++){
                if(res.data(i, j, 0) == z){
                    occurrencies(a, b, 1)++;
                }
            }
        }
        b++;
        if(b == 16){
            a++;
            b = 0;
        }
    }

    rows = occurrencies.rows();
    cols = occurrencies.cols();

    float cdf = 0.0f;
    for(int i = 0; i < rows; i++){
        for(int j = 0; j < cols; j++){
            if(occurrencies(i, j, 1) != 0){
                cdf += occurrencies(i, j, 1);
                occurrencies(i, j, 2) = cdf;
            }
        }
    }

    float cdf_min = occurrencies(0, 0, 2);
    int dimension = res.getRows() * res.getCols();
    int L = 256;

    rows = res.getRows();
    cols = res.getCols();
    int depth = res.getDepth();

    a = 0;
    b = 0;

    for(int i = 0; i < rows; i++){
        for(int j = 0; j < cols; j++){
            for(int k = 0; k < depth; k++){
                if(res.data(i, j, k) != 0){
                    cdf = occurrencies.getCDF(res.data(i, j, k));

                    if(cdf != 0){
                        res.data(i, j, k) = (int) ( ((cdf - cdf_min) / (dimension - cdf_min)) * (L - 1) );
                    }
                }
            }
        }
    }

    return res;
}
