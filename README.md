# Tiny libraries

| Filename                      | Description                                       | Example                       |
|-------------------------------|---------------------------------------------------|-------------------------------|
| [bezier_curve.hpp](bezier_curve/bezier_curve.hpp) | Calculate the points of an n deegre bezier curve. | [Bezier curve example](#bezier-curve-example) |

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
