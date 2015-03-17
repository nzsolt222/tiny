# Tiny libraries

| Filename                      | Description                                       | Example                       |
|-------------------------------|---------------------------------------------------|-------------------------------|
| [bezier_curve.hpp](bezier_curve/bezier_curve.hpp) | Calculate the points of an n deegre bezier curve. | [Bezier curve example](#bezier-curve-example) |
| [hermite_spline.hpp](hermite_spline/hermite_spline.hpp) | Calculate the points of a Hermite spline. | [Hermite spline example](#hermite-spline-example) |

## Bezier curve example
```cpp
#include <iostream>
#include "bezier_curve.hpp"

struct Point
{
    Point(int x, int y) : x(x), y(y){};
    int x;
    int y;
};

Point operator*(double scalar, const Point& value)
{
    return Point(value.x * scalar, value.y * scalar);
}

Point operator+(const Point& lhs, const Point& rhs)
{
    return Point(lhs.x + rhs.x, lhs.y + rhs.y);
}

int main()
{
    // set the controll points, and precision
    nzs::BezierCurve<Point> curve({{0, 0}, {10, 35}, {15, 20}}, 10);

    // curve.set_precision(100);
    // curve.set_control_points({{9, 2}, {3, 2}, {15, 20}});
    // curve.update(); // if you make any changes, you have to call the update
    // function

    // print the points of the curve
    for (const auto& point : curve)
    {
        std::cout << point.x << " " << point.y << std::endl;
    }

    return 0;
}
```

## Hermite spline example
```cpp
#include "hermite_spline.hpp"

#include <iostream>
#include <vector>

struct Point
{
    Point(int x, int y) : x(x), y(y){};
    int x;
    int y;
};

int main()
{
    // create a hermite spline
    std::vector<nzs::PointParam> point_param = {
        // x, y, t
        nzs::PointParam(200, 200, -1.5), nzs::PointParam(300, 600, +0.5),
        nzs::PointParam(600, 300, +2),
    };

    std::vector<nzs::TangentParam> tangent_param = {
        // x, y
        nzs::TangentParam(50, 50), nzs::TangentParam(50, 50),
    };

    int precision = 100;
    nzs::HermiteSpline<Point> spline(std::move(point_param),
                                     std::move(tangent_param), precision);

    // if you change something you have to call the HermiteSpline::update()
    // function.

    for (const auto& point : spline)
    {
        std::cout << point.x << " " << point.y << std::endl;
    }

    return 0;
}
```
