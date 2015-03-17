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

#ifndef NZS_HERMITE_SPLINE_HPP
#define NZS_HERMITE_SPLINE_HPP

#include <vector>
#include <cassert>
#include <initializer_list>
#include <type_traits>
#include <cmath>

namespace nzs
{
namespace details
{
template <class T, class U>
auto divide(T t, U u) -> decltype(t / u);

template <class T, class U>
struct result_of_divide
{
    using type = typename std::result_of<decltype(divide<T, U>)&(T, U)>::type;
};

template <class T, class U>
auto mult(T t, U u) -> decltype(t *u);

template <class T, class U>
struct result_of_multiply
{
    using type = typename std::result_of<decltype(mult<T, U>)&(T, U)>::type;
};

template <typename value_type>
class DynamicMatrix
{
public:
    using Row_t = std::vector<value_type>;
    using MatrixData_t = std::vector<Row_t>;

    DynamicMatrix(std::size_t row, std::size_t column,
                  std::initializer_list<value_type> list)
        : DynamicMatrix(row, column)
    {
        assert(list.size() == row * column);

        int row_count = 0;
        int column_count = 0;
        for (const auto &value : list)
        {
            matrix_data_[row_count][column_count] = value;
            ++column_count;
            if (column_count == column_)
            {
                column_count = 0;
                ++row_count;
            }
        }
    }

    DynamicMatrix(std::size_t row, std::size_t column)
        : matrix_data_(row, Row_t(column)), row_(row), column_(column)
    {
    }

    inline DynamicMatrix(DynamicMatrix &) = default;
    inline DynamicMatrix &operator=(DynamicMatrix &) = default;

    DynamicMatrix(DynamicMatrix &&matrix)
        : matrix_data_(std::move(matrix.matrix_data_)),
          row_(matrix.row_),
          column_(matrix.column_)
    {
    }

    DynamicMatrix &operator=(DynamicMatrix &&matrix)
    {
        std::swap(matrix_data_, matrix.matrix_data_);
        std::swap(row_, matrix.row_);
        std::swap(column_, matrix.column_);

        return *this;
    }

    inline Row_t &operator[](std::size_t index) { return matrix_data_[index]; }

    inline const Row_t &operator[](std::size_t index) const
    {
        return matrix_data_[index];
    }

    inline std::size_t get_row_size() const noexcept { return row_; }

    inline std::size_t get_column_size() const noexcept { return column_; }

    inline auto get_inverse() const
        -> DynamicMatrix<typename result_of_divide<double, value_type>::type>
    {
        assert(row_ == column_ &&
               "number of columns of the left matrix is not the same as the "
               "number of rows of the right matrix");
        DynamicMatrix<typename result_of_divide<double, value_type>::type>
            inverse_matrix(row_, column_);

        auto determinant = 1.0 / get_determinant();

        for (std::size_t column = 0; column < column_; ++column)
        {
            for (std::size_t row = 0; row < row_; ++row)
            {
                auto minor = get_minor(column, row);
                inverse_matrix[row][column] =
                    determinant * minor.get_determinant();
                if ((row + column) % 2 == 1)
                {
                    inverse_matrix[row][column] *= -1;
                }
            }
        }

        return inverse_matrix;
    }

    inline DynamicMatrix get_minor(int minor_row, int minor_column) const
    {
        assert(row_ == column_ &&
               "number of columns is not the same as the number of rows");
        assert(row_ > 2 &&
               "the minimum number of the row or the column is two");

        DynamicMatrix return_value(row_ - 1, column_ - 1);
        int column_count = 0;
        int row_count = 0;

        for (int row = 0; row < row_; ++row)
        {
            for (int column = 0; column < column_; ++column)
            {
                if (row != minor_row && column != minor_column)
                {
                    return_value[column_count][row_count] =
                        matrix_data_[row][column];
                    ++row_count;
                    if (row_count == row_ - 1)
                    {
                        row_count = 0;
                        ++column_count;
                    }
                }
            }
        }

        return return_value;
    }

    inline value_type get_determinant() const
    {
        assert(row_ == column_ &&
               "number of columns is not the same as the number of rows");

        if (row_ == 2 && column_ == 2)
        {
            return (matrix_data_[0][0] * matrix_data_[1][1]) -
                   (matrix_data_[0][1] * matrix_data_[1][0]);
        };

        if (row_ == 1 && column_ == 1)
        {
            return matrix_data_[0][0];
        }

        double determinant_value = 0;
        for (std::size_t column = 0; column < column_; ++column)
        {
            auto minor = get_minor(0, column);
            // the recusion is here
            determinant_value += (column % 2 == 1 ? -1.0 : 1.0) *
                                 matrix_data_[0][column] *
                                 minor.get_determinant();
        }
        return determinant_value;
    }

private:
    MatrixData_t matrix_data_;
    std::size_t row_;
    std::size_t column_;
};

template <typename lhs_value_type, typename rhs_value_type>
inline DynamicMatrix<
    typename result_of_multiply<lhs_value_type, rhs_value_type>::type>
operator*(const DynamicMatrix<lhs_value_type> &lhs,
          const DynamicMatrix<rhs_value_type> &rhs)
{
    assert(lhs.get_column_size() == rhs.get_row_size() &&
           "number of columns of the left matrix is not the same as the number "
           "of rows of the right matrix");
    DynamicMatrix<
        typename result_of_multiply<lhs_value_type, rhs_value_type>::type>
        return_maxtrix(lhs.get_row_size(), rhs.get_column_size());

    for (long unsigned int row1 = 0; row1 < lhs.get_row_size(); ++row1)
    {
        for (long unsigned int column2 = 0; column2 < rhs.get_column_size();
             ++column2)
        {
            for (long unsigned int column1 = 0; column1 < lhs.get_column_size();
                 ++column1)
            {
                return_maxtrix[row1][column2] +=
                    lhs[row1][column1] * rhs[column1][column2];
            }
        }
    }
    return return_maxtrix;
}

}  // details

class BaseParam
{
public:
    inline void set_x(double x) noexcept { x_ = x; }

    inline void inc_x(double value) noexcept { x_ += value; }

    inline double get_x() const noexcept { return x_; }

    inline void set_y(double y) noexcept { y_ = y; }

    inline void inc_y(double value) noexcept { y_ += value; }

    inline double get_y() const noexcept { return y_; }

protected:
    double x_;
    double y_;

    BaseParam(int x, int y) : x_(x), y_(y) {}
};

class PointParam : public BaseParam
{
public:
    PointParam(double x, double y, double t) : BaseParam(x, y), t_(t) {}

    inline void set_t(double t) noexcept { t_ = t; }

    inline double get_t() const noexcept { return t_; }

private:
    double t_;
};

class TangentParam : public BaseParam
{
public:
    TangentParam(double x, double y) : BaseParam(x, y) {}
};

template <class Point>
class HermiteSpline
{
public:
    using const_iterator = typename std::vector<Point>::const_iterator;

    HermiteSpline(std::vector<PointParam> point_param,
                  std::vector<TangentParam> tangent_param,
                  std::size_t precision)
        : point_param_(std::move(point_param)),
          tangent_param_(std::move(tangent_param)),
          n_(point_param_.size() + tangent_param_.size()),
          G_(2, n_),
          M_(n_, n_),
          M_inverse_(n_, n_),
          precision_(precision),
          M_matrix_changed_(true),
          spline_changed_(true)

    {
        assert(point_param_.size() >= 2 &&
               "the minimum number of PointParam: 2");
        assert(tangent_param_.size() <= point_param_.size() &&
               "the maximum number of TangentParam is number of PointParam");
        assert(point_param_.size() >= tangent_param_.size() &&
               "the minimum number of PointParam is number of TangentParam");

        calculate_G_matrix();
        calculate_M_matrix();
        update();
    }

    // recalculate the spline points if something changed
    inline void update()
    {
        if (!spline_changed_)
        {
            return;
        }

        if (M_matrix_changed_)
        {
            M_inverse_ = M_.get_inverse();
            M_matrix_changed_ = false;
        }

        // remove previously calculated points
        spline_points_.clear();

        const auto from_t = point_param_.back().get_t();
        const auto to_t = point_param_.front().get_t();
        const auto inc = (from_t - to_t) / precision_;
        const auto C = G_ * M_inverse_;

        for (double t = to_t; t <= from_t; t += inc)
        {
            spline_points_.push_back(get_spline_point(C, t));
        }
        // make sure [to_t, from_t]
        spline_points_.push_back(get_spline_point(C, from_t));

        spline_changed_ = false;
    }

    inline std::size_t get_precision() const noexcept { return precision_; }

    inline void set_precision(std::size_t precision) noexcept
    {
        precision_ = precision;
        spline_changed_ = true;
    }

    inline int get_point_param_size() const noexcept
    {
        return point_param_.size();
    }

    inline const PointParam &get_point_param(int index) const
    {
        return point_param_.at(index);
    }

    inline void set_point_param(int index, const PointParam &p)
    {
        auto old = point_param_.at(index);
        point_param_.at(index) = p;

        // we have to recalculate the G matrix if the position changed
        if (old.get_x() != p.get_x() || old.get_y() != p.get_y())
        {
            G_[0][index] = p.get_x();
            G_[1][index] = p.get_y();
            spline_changed_ = true;
        }

        // we have to recalculate the M matrix if the t parameter changed
        if (old.get_t() != p.get_t())
        {
            M_matrix_changed_ = true;
            spline_changed_ = true;
            for (int row = 0; row < n_; ++row)
            {
                M_[row][index] = get_M_matrix_element(row, index);
                if (index + point_param_.size() < n_)
                {
                    M_[row][index + point_param_.size()] =
                        get_M_matrix_element(row, index + point_param_.size());
                }
            }
        }
    }

    inline int get_tangent_param_size() const noexcept
    {
        return tangent_param_.size();
    }

    inline const TangentParam &get_tangent_param(int index) const
    {
        return tangent_param_.at(index);
    }

    inline void set_tangent_param(int index, const TangentParam &t)
    {
        auto old = tangent_param_.at(index);
        tangent_param_.at(index) = t;

        // if the G matrix changed we have to recalculate the hermite spline
        if (old.get_x() != t.get_x() || old.get_y() != t.get_y())
        {
            G_[0][point_param_.size() + index] = t.get_x();
            G_[1][point_param_.size() + index] = t.get_y();
            spline_changed_ = true;
        }
    }

    inline const_iterator begin() const noexcept { return cbegin(); }

    inline const_iterator end() const noexcept { return cend(); }

    inline const_iterator cbegin() const noexcept
    {
        return spline_points_.cbegin();
    }

    inline const_iterator cend() const noexcept
    {
        return spline_points_.cend();
    }

private:
    using MatrixElementType = double;

    std::vector<PointParam> point_param_;
    std::vector<TangentParam> tangent_param_;
    const int n_;
    details::DynamicMatrix<MatrixElementType> G_;
    details::DynamicMatrix<MatrixElementType> M_;
    details::DynamicMatrix<MatrixElementType> M_inverse_;
    std::size_t precision_;
    std::vector<Point> spline_points_;
    bool M_matrix_changed_;
    bool spline_changed_;

    /*
        {
            point x, point x, point x, ..., tangent x, tangent x, tangent x, ...
            point y, point y, point y, ..., tangent y, tangent y, tangent y, ...
        }
     */
    inline void calculate_G_matrix()
    {
        for (std::size_t i = 0; i < point_param_.size(); ++i)
        {
            G_[0][i] = point_param_[i].get_x();
            G_[1][i] = point_param_[i].get_y();
        }

        for (std::size_t i = 0; i < tangent_param_.size(); ++i)
        {
            G_[0][i + point_param_.size()] = tangent_param_[i].get_x();
            G_[1][i + point_param_.size()] = tangent_param_[i].get_y();
        }
    }

    /*
    {
        point.t ^ n - 1, ...,  n-1 * tagent.t ^ n - 2
        point.t ^ n - 2, ...,  n-2 * tagent.t ^ n - 3
        point.t ^ n - 3, ...,  n-3 * tagent.t ^ n - 4
        ...,             ...,  ...
        point.t ^ 1,     ...,  n-(n-1) * tagent.t ^ n - (n-1) = 1
        point.t ^ 0,     ...,  n-n * tagent.t ^ n - n = 0
    }
    */
    inline void calculate_M_matrix()
    {
        for (std::size_t column = 0; column < n_; ++column)
        {
            for (int row = 0; row < n_; ++row)
            {
                M_[row][column] = get_M_matrix_element(row, column);
            }
        }
    }

    inline MatrixElementType get_M_matrix_element(int row, int column) const
    {
        if (column < point_param_.size())
        {
            return pow(point_param_[column].get_t(), n_ - row - 1);
        }
        else
        {
            return (n_ - row - 1) *
                   pow(point_param_[column - point_param_.size()].get_t(),
                       n_ - row - 2);
        }
    }

    /*
    {
        t ^ n-1
        t ^ n-2
        t ^ n-3
        ...
        t ^ n-(n-1) = t
        t ^ n-n = 1
    }
    */
    inline details::DynamicMatrix<MatrixElementType> get_T_matrix(
        double t) const
    {
        auto t_mat = details::DynamicMatrix<MatrixElementType>(n_, 1);

        for (int row = 0; row < n_; ++row)
        {
            t_mat[row][0] = pow(t, n_ - row - 1);
        }
        return t_mat;
    }

    inline Point get_spline_point(
        const details::DynamicMatrix<MatrixElementType> &C, double t) const
    {
        auto t_mat = get_T_matrix(t);
        auto final_point = C * t_mat;
        return Point(final_point[0][0], final_point[1][0]);
    }
};

}  // nzs

#endif  // NZS_HERMITE_SPLINE_HPP
