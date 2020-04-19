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
    plot& add_fn(F const & f, std::array<unsigned char, 3> rgb, unsigned thickness = 0)
    {
        if (thickness == 0)
        {
            thickness = width_/100;
        }
        std::vector<Real> ys(width_);
        for (unsigned ix = 0; ix < width_; ++ix)
        {
            // Query at the center of the pixel.
            // a is identified with the center of the leftmost pixel. (ix = 0).
            // b is identified with the center of the rightmost pixel. (ix = width - 1).
            // If there are 2 pixels, then the distance between their centers is (b-a) in data space, and 1 in pixel space.            
            Real dx = (b_ - a_)/Real(width_ - 1);
            Real x = a_ + ix*dx;
            ys[ix] = f(x);
        }
        Real ymin;
        Real ymax;
        if (std::isnan(clip_max_))
        {
            if (std::isnan(clip_min_))
            {
                auto [it1, it2] = std::minmax_element(ys.begin(), ys.end());
                ymin = *it1;
                ymax = *it2;
            }
            else
            {
                ymin = clip_min_;
                ymax = *std::max_element(ys.begin(), ys.end());
            }
        }
        else
        {
            ymax = clip_max_;
            ymin = clip_min_;
        }

        if (std::isnan(clip_max_) && std::isnan(clip_min_)  && (ymin == ymax) )
        {
            throw std::domain_error("You can graph constants, but if so, you have to set the range.");
        }

        // Set the first pixel:
        Real y = ys[0];
        if (y <= ymax && y >= ymin)
        {
            Real ypx = (ymax-y)*Real(height_-1)/(ymax - ymin);
            unsigned iy = static_cast<unsigned>(std::round(ypx));
            img_[4 * width_ * iy + 0] = rgb[0];
            img_[4 * width_ * iy + 1] = rgb[1];
            img_[4 * width_ * iy + 2] = rgb[2];
            img_[4 * width_ * iy + 3] = 255;
        }

        // Set the last pixel:
        y = ys.back();
        if (y <= ymax && y >= ymin)
        {
            Real ypx = (ymax-y)*Real(height_-1)/(ymax - ymin);
            unsigned iy = static_cast<unsigned>(std::round(ypx));
            unsigned ix = width_ - 1;
            img_[4 * width_ * iy + 4 * ix + 0] = rgb[0];
            img_[4 * width_ * iy + 4 * ix + 1] = rgb[1];
            img_[4 * width_ * iy + 4 * ix + 2] = rgb[2];
            img_[4 * width_ * iy + 4 * ix + 3] = 255;
        }

        for (unsigned ix = 1; ix < width_ - 1; ++ix)
        {
            Real yleft = ys[ix-1];
            Real y = ys[ix];
            Real yright = ys[ix+1];
            // If the function goes beyond the bounds, we want it to be graphed right up to the edge.
            // So there's a bit of if logic to figure out how far to extrapolate.
            Real local_ymin = std::min(std::min(yleft, y), yright);
            if (local_ymin > ymax)
            {
                continue;
            }
            Real local_ymax = std::max(std::max(yleft, y), yright);
            if (local_ymax < ymin)
            {
                continue;
            }
            // This happens at a pole:
            if (yleft < ymin && yright > ymax)
            {
                continue;
            }
            // This happens at a pole:
            if (yleft > ymax && yright < ymin)
            {
                continue;
            }
            if (yleft > ymax)
            {
                yleft = ymax;
            }
            if (yright > ymax)
            {
                yright = ymax;
            }

            // iy = 0 is the top of the image, so use ymax - y, not y - ymin.
            // ymax is associated with the *center* of the pixel at iy = 0. (Imagine marking an axis at ymax to see this.)
            // ymin is associated with the *center* of the pixel at iy = height - 1.
            // So if there are 2 pixels, dy = 1 in pixel space.
            //Real ypx = (ymax-y)*Real(height_-1)/(ymax - ymin);
            Real ypxleft = (ymax-yleft)*Real(height_-1)/(ymax - ymin);
            Real ypxright = (ymax-yright)*Real(height_-1)/(ymax - ymin);
            unsigned iymin;
            unsigned iymax;
            // Make sure the graph is continuous:
            if (ypxleft <= ypxright)
            {
                iymin = static_cast<unsigned>(std::floor(ypxleft));
                iymax = static_cast<unsigned>(std::ceil(ypxright));
            }
            else
            { 
                iymin = static_cast<unsigned>(std::floor(ypxright));
                iymax = static_cast<unsigned>(std::ceil(ypxleft));
            }
            
            for (unsigned iy = iymin; iy < iymax && iy < height_; ++iy)
            {
                img_[4 * width_ * iy + 4 * ix + 0] = rgb[0];
                img_[4 * width_ * iy + 4 * ix + 1] = rgb[1];
                img_[4 * width_ * iy + 4 * ix + 2] = rgb[2];
                img_[4 * width_ * iy + 4 * ix + 3] = 255;
            }
            // Now we need to mark the pixels along the normal to the tangent, whose equation is
            // y - f(x_i) = (-1/f'(x_i))*(x-x_i).
            // We can get a good approximation to this by using the symmetric finite difference.
            // I think this is closely related to Xiaolin Wu's Line algorithm:
            // https://en.wikipedia.org/wiki/Xiaolin_Wu's_line_algorithm
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

    plot& write_gridlines(std::array<unsigned char, 3> rgb = {255, 255, 255})
    {
        if (clip_min_ <= 0 && clip_max_ >= 0)
        {
            Real y = 0;
            Real ypx = (clip_max_-y)*Real(height_-1)/(clip_max_ - clip_min_);
            unsigned iy = static_cast<unsigned>(std::floor(ypx));
            for (unsigned ix = 0; ix < width_; ++ix)
            {
                img_[4 * width_ * iy + 4 * ix + 0] = rgb[0];
                img_[4 * width_ * iy + 4 * ix + 1] = rgb[1];
                img_[4 * width_ * iy + 4 * ix + 2] = rgb[2];
                img_[4 * width_ * iy + 4 * ix + 3] = 255;
            }
        }

        if (a_ <= 0 && b_ > 0)
        {
            Real xpx = (0 - a_)*Real(width_-1)/(b_ - a_);
            unsigned ix = static_cast<unsigned>(std::round(xpx));
            for (unsigned iy = 0; iy < height_; ++iy)
            {
                img_[4 * width_ * iy + 4 * ix + 0] = rgb[0];
                img_[4 * width_ * iy + 4 * ix + 1] = rgb[1];
                img_[4 * width_ * iy + 4 * ix + 2] = rgb[2];
                img_[4 * width_ * iy + 4 * ix + 3] = 255;
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