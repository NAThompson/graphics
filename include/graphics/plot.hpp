#ifndef GRAPHICS_PLOT_HPP
#define GRAPHICS_PLOT_HPP
#include "lodepng.h"
#include <stdexcept>
#include <vector>
#include <array>
#include <algorithm>
#include <string>
#include <cmath>
#include <iostream>
#include <utility>

namespace graphics {
template<typename Real>
class plot {
public:

    plot(std::pair<Real, Real> const & domain, unsigned width, unsigned height = 0, std::pair<Real, Real> range = {std::numeric_limits<Real>::quiet_NaN(), std::numeric_limits<Real>::quiet_NaN()})
    {
        a_ = domain.first;
        b_ = domain.second;
        if (b_ <= a_)
        {
            throw std::domain_error("b > a is required.");
        }
        if (width < 10)
        {
            throw std::domain_error("Must be at least 10 pixels wide!");
        }
        width_ = width;

        if (height == 0)
        {
            height_ = static_cast<unsigned>(std::floor(double(width_)/1.6180339));
        }
        else
        {
            height_ = height;
        }
        clip_max_ = range.second;
        clip_min_ = range.first;
        if (!std::isnan(clip_max_))
        {
            if (!std::isnan(clip_min_))
            {
                if (clip_min_ >= clip_max_)
                {
                    throw std::domain_error("Big trouble: clip_min_ is >= clip_max_!");
                }
            }
        }
        // turn it black:
        img_.resize(4*width_*height_, 0);
        for (size_t i = 1; 4*i < img_.size(); ++i)
        {
            img_[4*i-1] = 255;
        }
    }

    template<class F>
    plot& add_fn(F const & f, std::array<unsigned char, 4> rgba, unsigned thickness = 0)
    {
        if (thickness == 0)
        {
            thickness = width_/100;
        }
        std::vector<Real> ys(width_);
        //std::vector<Real> dydx(width_);
        for (unsigned ix = 0; ix < width_; ++ix)
        {
            // We query the function at the side of the pixels, not the center.
            // This allows us to fill in every pixel if the derivative is really large and create a continuous looking image:
            Real x = a_ + (b_-a_)*ix/Real(width_ - 1);
            ys[ix] = f(x);
            //dydx[ix] = boost::math::differentiation::finite_difference_derivative(f, x);
        }
        auto [it1, it2] = std::minmax_element(ys.begin(), ys.end());
        Real ymin = *it1;
        Real ymax = *it2;
        if (!std::isnan(clip_max_))
        {
            ymax = clip_max_;
        }
        if (!std::isnan(clip_min_))
        {
            ymin = clip_min_;
        }

        if (std::isnan(clip_max_) && std::isnan(clip_min_)  && (ymin == ymax) )
        {
            throw std::domain_error("You can graph constants, but if so, you have to set clips.");
        }
        
        for (unsigned ix = 0; ix < width_ - 1; ++ix)
        {
            Real yi = ys[ix];
            Real yip1 = ys[ix+1];
            if (yi > ymax || yi < ymin)
            {
                continue;
            }
             // iy = 0 is the top of the image, so use ymax - y, not y - ymin:
            Real yipx = (ymax-yi)*Real(height_-1)/(ymax - ymin);
            Real yip1px = (ymax-yip1)*Real(height_-1)/(ymax - ymin);
            unsigned iymin;
            unsigned iymax;
            // ys[ix] is associated with the left side, ys[ix+1] is associated with the right side.
            // Set every pixel in between.
            if (yipx <= yip1px)
            {
                iymin = static_cast<unsigned>(std::floor(yipx));
                iymax = static_cast<unsigned>(std::ceil(yip1px));
            }
            else
            { 
                iymin = static_cast<unsigned>(std::floor(yip1px));
                iymax = static_cast<unsigned>(std::ceil(yipx));
            }
            // If the function is (say) constant, then the graph width can get *very* thin.
            bool toggle = true;
            //unsigned t = thickness*std::min(std::max(Real(1), std::abs(dydx[ix])), Real(8));
            while (iymax - iymin < thickness)
            {
                if (toggle)
                {
                    ++iymax;
                    toggle = false;
                }
                else {
                    if (iymin > 0)
                    {
                        --iymin;
                        toggle = true;
                    }
                    else
                    {
                        break;
                    }
                }
            }
            
            for (unsigned iy = iymin; iy < iymax && iy < height_; ++iy)
            {
                img_[4 * width_ * iy + 4 * ix + 0] = rgba[0];
                img_[4 * width_ * iy + 4 * ix + 1] = rgba[1];
                img_[4 * width_ * iy + 4 * ix + 2] = rgba[2];
                img_[4 * width_ * iy + 4 * ix + 3] = rgba[3];
            }
        }
        // Now we set a clip, since we can't go back and re-render this function:
        if (std::isnan(clip_max_))
        {
            clip_max_ = ymax;
        }
        if (std::isnan(clip_min_))
        {
            clip_min_ = ymin;
        }
        return *this;
    }

    plot& write_gridlines(std::array<unsigned char, 4> rgba = {255, 255, 255, 255})
    {
        if (clip_min_ <= 0 && clip_max_ >= 0)
        {
            Real y = 0;
            Real ypx = (clip_max_-y)*Real(height_-1)/(clip_max_ - clip_min_);
            unsigned iy = static_cast<unsigned>(std::floor(ypx));

            for (unsigned ix = 0; ix < width_; ++ix)
            {
                img_[4 * width_ * iy + 4 * ix + 0] = rgba[0];
                img_[4 * width_ * iy + 4 * ix + 1] = rgba[1];
                img_[4 * width_ * iy + 4 * ix + 2] = rgba[2];
                img_[4 * width_ * iy + 4 * ix + 3] = rgba[3];
            }
        }

        if (a_ <= 0 && b_ > 0)
        {
            Real xpx = (0 - a_)*Real(width_-1)/(b_ - a_);
            unsigned ix = static_cast<unsigned>(std::round(xpx));
            for (unsigned iy = 0; iy < height_; ++iy)
            {
                img_[4 * width_ * iy + 4 * ix + 0] = rgba[0];
                img_[4 * width_ * iy + 4 * ix + 1] = rgba[1];
                img_[4 * width_ * iy + 4 * ix + 2] = rgba[2];
                img_[4 * width_ * iy + 4 * ix + 3] = rgba[3];
            }
        }

        return *this;
    }

    unsigned write(std::string filename)
    {
        unsigned error = lodepng::encode(filename, img_, width_, height_, LodePNGColorType::LCT_RGBA, 8);
        if(error)
        {
            std::cerr << "Error encoding png: " << lodepng_error_text(error) << "\n";
        }
        return error;
    }

private:
    Real a_;
    Real b_;
    unsigned width_;
    unsigned height_;
    std::vector<unsigned char> img_;
    Real clip_max_;
    Real clip_min_;
};

}
#endif