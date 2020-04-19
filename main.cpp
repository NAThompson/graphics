#include <graphics/plot.hpp>
#include <iostream>
#include <cmath>

template<typename Real>
void test_plot()
{
    unsigned width = 4000;
    // when height = 0, it gets reset to the golden ratio
    unsigned height = width;
    std::pair<Real, Real> domain(-5, 5);
    std::pair<Real, Real> range(-5, 5);
    auto graph = graphics::plot<Real>(domain, width, height, range);
    graph.write_gridlines();
    std::array<unsigned char, 3> green{0, 255, 0};
    std::array<unsigned char, 3> red{255, 0, 0};
    std::array<unsigned char, 3> blue{0, 0, 255};
    auto f = [](Real x) { return std::sin(x);};
    graph.add_fn(f, green);

    auto g = [](Real x)->Real {
        if (x == 0)
        {
            return std::numeric_limits<Real>::quiet_NaN();
        }
        return 1/x;
    };
    graph.add_fn(g, red);
    auto g2 = [&](Real x)->Real {
        return -g(x);
    };
    std::array<unsigned char, 3> yellow{255, 255, 0};
    graph.add_fn(g2, yellow);

    auto h = [](Real x)->Real {
        if (x == 0)
        {
            return Real(0);
        }
        return sin(1/x);
    };
    graph.add_fn(h, blue);

    // Nothing gets displayed from this; that's good right?
    auto g3 = [](Real)->Real {
        return std::numeric_limits<Real>::quiet_NaN();
    };
    graph.add_fn(g3, blue);
    graph.write("plot.png");
}

int main()
{
    test_plot<float>();
}
