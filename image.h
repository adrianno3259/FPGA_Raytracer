#ifndef IMAGE_H
#define IMAGE_H

#include <vector>
#include <iostream>
#include "color.h"

class Image
{
public:
    std::vector< std::vector< Color> > pixels;
    int w, h;

    Image();// imagem preta, 200 x 200

    Image(const int s);
    Image(const int wi, const int he);

    Image( std::vector< std::vector<Color> > pxs);

    Color getPixel(const int x, const int y);

    void setPixel(const int x, const int y, const Color& c);

    void save(const char* filepath);

};


#endif // IMAGE_H

