// The MIT License (MIT)
//
// Copyright (c) 2015 Nagy Zsolt
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

#ifndef NZS_BEZIER_CURVE_HPP
#define NZS_BEZIER_CURVE_HPP

#include <vector>
#include <cmath>
#include <cassert>

namespace nzs
{
template <class Point>
class BezierCurve
{
public:
    using PointContainer_t = std::vector<Point>;
    using const_iterator = typename PointContainer_t::const_iterator;
    using iterator = typename PointContainer_t::iterator;

    BezierCurve(PointContainer_t control_points, std::size_t precision)
        : control_points_(std::move(control_points)),
          precision_(precision),
          have_to_update_(true)
    {
        assert(control_points_.size() >= 2 &&
               "The number of the control points are less than two.");
        update();
    }

    // apply the chnages
    void update()
    {
        if (have_to_update_)
        {
            calculate();
            have_to_update_ = false;
        }
    }

    inline void set_precision(std::size_t precision)
    {
        precision_ = precision;
        have_to_update_ = true;
    }

    inline std::size_t get_precision() const noexcept { return precision_; }

    inline const PointContainer_t &get_control_points() const noexcept
    {
        return control_points_;
    }

    inline const Point &get_control_point(int index) const
    {
        return control_points_.at(index);
    }

    inline void set_control_points(PointContainer_t control_points)
    {
        assert(control_points_.size() >= 2 &&
               "The number of the control points are less than two.");
        control_points_ = std::move(control_points);
        have_to_update_ = true;
    }

    inline void set_control_points(int index, Point new_point)
    {
        control_points_.at(index) = std::move(new_point);
        have_to_update_ = true;
    }

    inline const_iterator cbegin() const noexcept
    {
        return curve_points_.cbegin();
    }

    inline const_iterator cend() const noexcept { return curve_points_.cend(); }

    inline const_iterator begin() const noexcept { return cbegin(); }

    inline const_iterator end() const noexcept { return cend(); }

    inline iterator begin() noexcept { return curve_points_.begin(); }

    inline iterator end() noexcept { return curve_points_.end(); }

private:
    PointContainer_t control_points_;
    PointContainer_t curve_points_;
    std::size_t precision_;
    bool have_to_update_;

    inline int binomial_coefficient(int n, int k) const noexcept
    {
        int res = 1;

        if (k > n - k)
        {
            k = n - k;
        }

        for (int i = 0; i < k; ++i)
        {
            res *= (n - i);
            res /= (i + 1);
        }

        return res;
    }

    inline void calculate()
    {
        curve_points_.clear();

        for (double t = 0; t <= 1; t += (1. / precision_))
        {
            curve_points_.push_back(get_point(t));
        }
        // make sure t [0, 1]
        curve_points_.push_back(get_point(1));
    }

    // n degree Bézier curve general definition
    inline Point get_point(double t)
    {
        int n = control_points_.size() - 1;
        // prevent usage of default constructor
        Point final_point = binomial_coefficient(n, 0) *
                            std::pow(1 - t, n - 0) * std::pow(t, 0) *
                            control_points_[0];
        // i = 1 because of the above row
        for (int i = 1; i <= n; ++i)
        {
            final_point = final_point +
                          binomial_coefficient(n, i) * std::pow(1 - t, n - i) *
                              std::pow(t, i) * control_points_[i];
        }
        return final_point;
    }
};

}  // nzs

#endif  // NZS_BEZIER_CURVE_HPP
