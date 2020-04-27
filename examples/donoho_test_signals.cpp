#include <iostream>
#include <iomanip>
#include <cmath>
#include <list>
#include <vector>
#include <array>
#include <chrono>
#include <graphics/plot.hpp>
#include <boost/math/special_functions/daubechies_scaling.hpp>
#include <boost/math/special_functions/daubechies_wavelet.hpp>
#include <boost/math/quadrature/trapezoidal.hpp>
#include <boost/math/constants/constants.hpp>

using std::floor;
using std::pow;
using std::abs;
using std::sqrt;
using std::ceil;
using boost::math::quadrature::trapezoidal;
using boost::math::sign;
using std::sin;
using boost::math::constants::pi;


constexpr const int p = 3;
auto phi = boost::math::daubechies_scaling<double, p>();
auto psi = boost::math::daubechies_wavelet<double, p>();

// Definition of this function is given in Donoho, 1994, Table 1, Formula for test functions.
template<class Real>
Real bumps(Real x)
{
    if (x <= 0 || x >= 1)
    {
        return 0;
    }
    static constexpr const std::array<Real, 11> t{0.1, 0.13, 0.15, 0.23, 0.25, 0.40, 0.44, 0.65, 0.76, 0.78, 0.81};
    static constexpr const std::array<Real, 11> h{4.0, 5.0, 3.0, 4.0, 5.0, 4.2, 2.1, 4.3, 3.1, 5.1, 4.2};
    static constexpr const std::array<Real, 11> w{0.005, 0.005, 0.006, 0.01, 0.01, 0.03, 0.01, 0.01, 0.005, 0.008, 0.005};

    Real f = 0;
    for (size_t i = 0; i < 11; ++i) {
        Real z = abs((x-t[i])/w[i]);
        f += h[i]*pow(1+z, -4);
    }
    return f;
}

template<class Real>
Real heavisine(Real x)
{
    return 4*sin(4*pi<Real>()*x) - sign(x - 0.3)  - sign(0.72 - x);
}

template<typename Real>
Real bump_function(Real x)
{
    if (x <= -1 || x >= 1)
    {
        return Real(0);
    }
    Real t = 1/(1-x*x);
    return std::exp(-t);
}

template<typename Real, int p>
class wavelet_series {
public:
    wavelet_series(std::function<Real(Real)> const f, 
                   std::pair<Real, Real> const & support,
                   std::pair<int, int> jrange,
                   boost::math::daubechies_wavelet<Real, p> psi,
                   boost::math::daubechies_scaling<Real, p> phi)
    {
        using std::floor;
        using std::ceil;
        using std::pow;
        using std::sqrt;
        using boost::math::quadrature::trapezoidal;

        x0_ = support.first;
        xf_ = support.second;
        jmin_ = jrange.first;
        jmax_ = jrange.second;
        // The support of phi_jmax,k is [2^jmax(k), 2^jmax(2p-1+k)], the support of f is [x0, xf].
        // If 2^jmax*k >= xf, then the coefficient is zero, so the last non-zero is 2^jmax*kmax < xf.
        int kmax = std::floor(xf_/pow(2, jmax_));
        // If 2^jmax(2p-1+k) <= x0, then the coefficient is zero.
        // So the first nonzero occurs when k > x0/2^jmax + 1 - 2p
        int kmin = std::ceil(x0_/pow(2, jmax_) + 1 - 2*p);
        for (int k = kmin; k <= kmax; ++k)
        {
            auto integrand = [&](Real x)->Real {
                Real t = x/pow(2.0, jmax_) - k;
                if (t <= 0 || t >= 2*p - 1)
                {
                    return 0;
                }
                return f(x)*phi(t)/sqrt(pow(2.0, jmax_));
            };
            Real tol = 10*std::numeric_limits<Real>::epsilon();
            Real a = std::max(pow(2, jmax_)*k, x0_);
            Real b = std::min(pow(2, jmax_)*(2*p-1+k), xf_);
            if (b < a) {
                std::cout << "b < a wft???\n";
                std::cout << "[a, b] = [" << a << ", " << b << "]\n";
                std::cout << "Occured at k = " << k << "\n";
            }
            Real ajk = trapezoidal(integrand, a, b, tol);
            a_.push_back(ajk);
        }

        // Now for the detail spaces:
        for (int j = jmin_; j <= jmax_; ++j)
        {
            // The support of psi_jk is [2^j(-p+1+k), 2^j(p+k)].
            // If 2^j(-p+1+k) >= xf, the coefficient is zero.
            // So the last nonzero coefficient is  k < 2^{-j}xf - 1 + p.
            int kmax = std::floor(pow(2.0, -j)*xf_ + p - 1);
            // If 2^j(p+k) <= x0, then the coefficient is zero.
            // So the first nonzero coefficient is
            // 2^j(p+k) > x0, so kmin > x0/2^j - p
            int kmin = std::ceil(pow(2, -j)*x0_ - p);
            std::vector<Real> d;
            for (int k = kmin; k <= kmax; ++k)
            {
                auto integrand = [&](Real x)->Real {
                    Real t = x/pow(2.0, j) - k;
                    if (t <= -p + 1 || t >= p) {
                        return 0;
                    }
                    return f(x)*psi(t)/sqrt(pow(2.0, j));
                };
                Real tol = 10*std::numeric_limits<Real>::epsilon();
                Real a = std::max(std::pow(2, j)*(-p+1+k), x0_);
                Real b = std::min(std::pow(2, j)*(p+k), xf_);
                if (b < a) {
                    std::cout << "b < a in detail coeffs; wft???\n";
                }
                Real djk = trapezoidal(integrand, a, b, tol);
                d.push_back(djk);
            }
            l_.push_back(d);
        }
    }

    Real operator()(Real x) const {
        if (x <= x0_ || x >= xf_)
        {
            return Real(0);
        }
        Real sum = 0;
        int kmin = std::ceil(x0_/pow(2, jmax_) + 1 - 2*p);
        // The analysis terms:
        Real tmin = x/pow(2, jmax_) - kmin;
        for (int kidx = 0; kidx < static_cast<int>(a_.size()); ++kidx)
        {
            Real t = tmin - kidx;
            if (t <= 0) {
                break;
            }
            if (a_[kidx] != 0)
            {
                sum += a_[kidx]*phi(t)/sqrt(pow(2.0, jmax_));
            }
        }
        // The detail terms:
        int j = jmin_;
        for (auto const & d : l_)
        {
            kmin = std::ceil(pow(2, -j)*x0_ - p);
            tmin = x/pow(2, j) - kmin;
            for (int kidx = 0; kidx < static_cast<int>(d.size()); ++kidx)
            {
                Real t = tmin - kidx;
                if (t <= -p + 1) {
                    break;
                }
                if (d[kidx] != 0) {
                    sum += d[kidx]*psi(t)/sqrt(pow(2,j));
                }
            }
            ++j;
        }
        return sum;
    }

    Real hoyer_sparsity() const {
        Real l1 = 0;
        Real l2_sq = 0;
        for (auto & t : a_) {
            l1 += abs(t);
            l2_sq += t*t;
        }
        Real terms = a_.size();

        for (auto const & d : l_) {
            terms += d.size();
            for (size_t k = 0; k < d.size(); ++k)
            {
                l1 += abs(d[k]);
                l2_sq += d[k]*d[k];
            }
        }
        return (sqrt(terms) - l1/sqrt(l2_sq))/(sqrt(terms)-1);
    }

    size_t bytes() {
        size_t b = a_.size();
        for (auto const & d : l_) {
            b += d.size();
        }
        return b*sizeof(Real) + sizeof(this);
    }

    void hard_threshold(Real threshold)
    {
        for (auto & a : a_) {
            if (abs(a) < threshold) {
                a = 0;
            }
        }
        for (auto & detail_level : l_) {
            for (auto & d : detail_level) {
                if (abs(d) < threshold)
                {
                    d = 0;
                }
            }
        }
    }

    size_t thresholded_bytes() {
        size_t nonzero_terms = 0;
        for (auto const & a : a_) {
            if (a != 0) {
                ++nonzero_terms;
            }
        }
        for (auto const & detail_level : l_) {
            for (auto & d : detail_level) {
                if (d != 0)
                {
                    ++nonzero_terms;
                }
            }
        }
        return nonzero_terms*sizeof(Real) + sizeof(this);
    }

private:
    Real x0_;
    Real xf_;
    int jmin_;
    int jmax_;
    std::vector<Real> a_;
    std::list<std::vector<Real>> l_;
};


int main()
{
    //std::pair<double, double> support{0, 2*p-1};
    //std::pair<double, double> support{-p+1, p};
    std::pair<double, double> support{0, 1};

    int width = 2800;
    auto plot = graphics::plot<double>(support, width, 1800);

    plot.add_fn(bumps<double>, {255,0,0});
    plot.write_gridlines();
    std::cout << std::setprecision(std::numeric_limits<double>::max_digits10);

    std::pair<int, int> jrange{-12, 0};
    auto ws = wavelet_series<double, p>(bumps<double>, support, jrange, psi, phi);
    std::cout << "Hoyer sparsity of wavelet series at p = " << p << " is " << ws.hoyer_sparsity() << "\n";
    size_t prebytes = ws.bytes();
    std::cout << "Occupies " << prebytes << " bytes before thresholding.\n";
    //ws.hard_threshold(1.5e-3);
    ws.hard_threshold(1e-4);
    size_t postbytes = ws.thresholded_bytes();
    std::cout << "Occupies " << postbytes << " bytes after thresholding.\n";
    std::cout << "compression ratio = " << prebytes/postbytes << "\n";

    auto start = std::chrono::system_clock::now();
    plot.add_fn(ws, {0, 255, 0});

    auto end = std::chrono::system_clock::now();

    auto elapsed =
    std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    std::cout << "Wrote " << width << " pixels in " << elapsed.count() << " ms.\n";

    plot.write_gridlines();
    plot.write("bumps.png");
}