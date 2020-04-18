#include <graphics/plot.hpp>
#include <iostream>
#include <cmath>

template<typename Real>
void test_plot()
{
    unsigned width = 4000;
    // when height = 0, it gets reset to the golden ratio:s
    unsigned height = 0;
    std::pair<Real, Real> domain(-5, 5);
    std::pair<Real, Real> range(-4, 4);
    auto graph = graphics::plot<Real>(domain, width, height, range);
    std::array<unsigned char, 4> green{0, 255, 0, 255};
    std::array<unsigned char, 4> red{255, 0, 0, 255};
    std::array<unsigned char, 4> blue{0, 0, 255, 255};
    graph.add_fn([](Real x) { return sin(x);}, green);
    graph.add_fn([](Real x) { if (x==0) {return Real(0); } return std::sin(1/x);} , red);
    graph.add_fn([](Real x) { return 1/x;}, blue);
    graph.write_gridlines().write("plot.png");
}

int main()
{
    test_plot<float>();
}
