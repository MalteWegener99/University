#include <tuple>
#include <string>
#include <algorithm>
#include <functional>
#include <math.h>
#include <iostream>

class Color{
public:
    Color();
    Color(int r, int g, int b){
        R = std::clamp<int>(r, 0, 255);
        G = std::clamp<int>(g, 0, 255);
        B = std::clamp<int>(b, 0, 255);
    };
    Color operator+(Color& RHS){
        return Color(R+RHS.R, G+RHS.G, B+RHS.B);
    };
    Color operator-(Color& RHS){
        return Color(R-RHS.R, G-RHS.G, B-RHS.B);
    };
    Color operator/(double RHS){
        return Color(
            std::round(R/RHS),
            std::round(G/RHS),
            std::round(B/RHS)
        );
    };
    Color interpolate(Color& RHS, double val){
        int Rd = RHS.R - R;
        int Gd = RHS.G - G;
        int Bd = RHS.B - B;
        return Color(
            R + std::round(val*Rd),
            G + std::round(val*Gd),
            B + std::round(val*Bd)
        );
    };
    unsigned int R, G, B;
};

class Colormap{
public:
    Colormap();
    Colormap(std::vector<Color> x){
        keyframes = x;
        h = (double)1/(x.size()-1);
    };
    Color get_val(double val){
        val = std::clamp<double>(val, 0, 1);
        int i = (int)(val/h);
        double x = (val-i*h)/h;
        return keyframes[i].interpolate(keyframes[i+1], x);
    };

    std::vector<Color> keyframes;
    double h;
};
