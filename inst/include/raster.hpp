#ifndef POPS_RASTER_HPP
#define POPS_RASTER_HPP

/*
 * PoPS model - native raster manipulation
 *
 * Copyright (C) 2015-2019 by the authors.
 *
 * Authors: Vaclav Petras <wenzeslaus gmail com>
 *          Completely rewritten by Vaclav Petras based on
 *          version by Zexi Chen <zchen22 ncsu edu>.
 *
 * The code contained herein is licensed under the GNU General Public
 * License. You may obtain a copy of the GNU General Public License
 * Version 2 or later at the following locations:
 *
 * http://www.opensource.org/licenses/gpl-license.html
 * http://www.gnu.org/copyleft/gpl.html
 */

#include <iostream>
#include <ostream>
#include <string>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <initializer_list>
#include <stdlib.h>

using std::string;
using std::cerr;
using std::endl;

namespace pops {

/*! Iterate over two ranges and apply a binary function which modifies
 *  the first parameter.
 */
template<class InputIt1, class InputIt2, class BinaryOperation>
BinaryOperation for_each_zip(InputIt1 first1, InputIt1 last1, InputIt2 first2, BinaryOperation f) {
    for (; first1 != last1; ++first1, ++first2) {
        f(*first1, *first2);
    }
    return f;
}

/*! Representation of a raster image.
 *
 * The object support raster algebra operations:
 *
 * ```
 * Raster<int> a = {{1, 2}, {3, 4}};
 * auto b = 2 * (a + 1);
 * ```
 *
 * The raster algebra operations sometimes overlap with matrix
 * operations, e.g. for plus operator or multiplication by scalar.
 * However, in some cases, the behavior is different, e.g.,
 * multiplication of the rasters results in a new raster with cell
 * values which are result of multiplying cell values in the relevant
 * positions of the two raster.
 *
 * ```
 * Raster<int> a = {{1, 2}, {3, 4}};
 * auto b = 2 * (a + 1);
 * ```
 *
 * The template parameter Number is the numerical type of the raster,
 * typically int, float, or double.
 *
 * The internal storage is directly accessible which comes with a great
 * responsibility. Although for the computations themselves, direct
 * access is not needed, it gives a lot of advantages when
 * initializing the values as well as when putting them to some storage
 * when the computation is done. The values need to be stored in
 * one liner array with row-major oder, i.e. individual value can be
 * accessed using the following:
 *
 * ```
 * row * total_number_of_columns + column
 * ```
 */
template<typename Number>
class Raster
{
protected:
    unsigned cols_;
    unsigned rows_;
    Number *data_;
public:
    Raster()
    {
        cols_ = 0;
        rows_ = 0;
        data_ = NULL;
    }

    Raster(const Raster& other)
    {
        cols_ = other.cols_;
        rows_ = other.rows_;
        data_ = new Number[cols_ * rows_];
        std::copy(other.data_, other.data_ + (cols_ * rows_), data_);
    }

    /*! Initialize size using another raster, but use given value
     *
     * The values in the other raster are not used.
     */
    Raster(const Raster& other, Number value)
    {
        cols_ = other.cols_;
        rows_ = other.rows_;
        data_ = new Number[cols_ * rows_]{value};
    }

    Raster(Raster&& other)
    {
        cols_ = other.cols_;
        rows_ = other.rows_;
        data_ = other.data_;
        other.data_ = nullptr;
    }

    Raster(int rows, int cols)
    {
        this->cols_ = cols;
        this->rows_ = rows;
        this->data_ = new Number[cols_ * rows_];
    }

    // TODO: size is unsigned?
    Raster(int rows, int cols, Number value)
    {
        this->cols_ = cols;
        this->rows_ = rows;
        this->data_ = new Number[cols_ * rows_]{value};
    }

    // maybe remove from the class, or make it optional together with
    // a reference
    Raster(std::initializer_list<std::initializer_list<Number>> l)
        : Raster(l.size(), l.begin()->size())
    {
         unsigned i = 0;
         unsigned j = 0;
         for (const auto& subl : l)
         {
            for (const auto& value : subl)
            {
               data_[cols_ * i + j] = value;
               ++j;
            }
            j = 0;
            ++i;
         }
    }

    ~Raster()
    {
        if (data_) {
            delete[] data_;
        }
    }

    unsigned cols() const
    {
        return cols_;
    }

    unsigned rows() const
    {
        return rows_;
    }

    /*! Returns pointer for direct access the underlying array.
     *
     * The values are stored in row-major order.
     * See the class description for details.
     */
    Number* data() noexcept
    {
        return data_;
    }

    /*! Returns pointer for direct access the underlying array.
     *
     * Same as the non-const version but used when the object is const.
     */
    const Number* data() const noexcept
    {
        return data_;
    }

    void fill(Number value)
    {
        std::fill(data_, data_ + (cols_ * rows_), value);
    }

    void zero()
    {
        std::fill(data_, data_ + (cols_ * rows_), 0);
    }

    template<class UnaryOperation>
    void for_each(UnaryOperation op)
    {
        std::for_each(data_, data_ + (cols_ * rows_), op);
    }

    const Number& operator()(unsigned row, unsigned col) const
    {
        return data_[row * cols_ + col];
    }

    Number& operator()(unsigned row, unsigned col)
    {
        return data_[row * cols_ + col];
    }

    Raster& operator=(const Raster& other)
    {
        if (this != &other)
        {
            if (data_)
                delete[] data_;
            cols_ = other.cols_;
            rows_ = other.rows_;
            data_ = new Number[cols_ * rows_];
            std::copy(other.data_, other.data_ + (cols_ * rows_), data_);
        }
        return *this;
    }

    Raster& operator=(Raster&& other)
    {
        if (this != &other)
        {
            if (data_)
                delete[] data_;
            cols_ = other.cols_;
            rows_ = other.rows_;
            data_ = other.data_;
            other.data_ = nullptr;
        }
        return *this;
    }

    Raster operator+(const Raster& image) const
    {
        if (this->cols_ != image.cols() || this->rows_ != image.rows()) {
            cerr << "The height or width of one image do not match with that of the other one!" << endl;
            return Raster();
        }
        else {
            auto re_width = this->cols_;
            auto re_height = this->rows_;
            auto out = Raster(re_height, re_width);

            for (unsigned i = 0; i < re_height; i++) {
                for (unsigned j = 0; j < re_width; j++) {
                    out.data_[i * cols_ + j] = this->data_[i * cols_ + j] + image.data_[i * cols_ + j];
                }
            }
            return out;
        }
    }

    Raster operator-(const Raster& image) const
    {
        if (this->cols_ != image.cols() || this->rows_ != image.rows()) {
            cerr << "The height or width of one image do not match with that of the other one!" << endl;
            return Raster();
        }
        else {
            auto re_width = this->cols_;
            auto re_height = this->rows_;
            auto out = Raster(re_height, re_width);

            for (unsigned i = 0; i < re_height; i++) {
                for (unsigned j = 0; j < re_width; j++) {
                    out.data_[i * cols_ + j] = this->data_[i * cols_ + j] - image.data_[i * cols_ + j];
                }
            }
            return out;
        }
    }

    Raster operator*(const Raster& image) const
    {
        if (cols_ != image.cols() || rows_ != image.rows()) {
            throw std::runtime_error("The height or width of one image do"
                                     " not match with that of the other one.");
        }
        auto out = Raster(rows_, cols_);

        std::transform(data_, data_ + (cols_ * rows_), image.data_, out.data_,
                       [](const Number& a, const Number& b) { return a * b; });
        return out;
    }

    Raster operator/(const Raster& image) const
    {
        if (cols_ != image.cols() || rows_ != image.rows()) {
            throw std::runtime_error("The height or width of one image do"
                                     " not match with that of the other one.");
        }
        auto out = Raster(rows_, cols_);

        std::transform(data_, data_ + (cols_ * rows_), image.data_, out.data_,
                       [](const Number& a, const Number& b) { return a / b; });
        return out;
    }

    Raster operator*(double value) const
    {
        auto out = Raster(rows_, cols_);

        std::transform(data_, data_ + (cols_ * rows_), out.data_,
                       [&value](const Number& a) { return a * value; });
        return out;
    }

    Raster operator/(double value) const
    {
        auto out = Raster(rows_, cols_);

        std::transform(data_, data_ + (cols_ * rows_), out.data_,
                       [&value](const Number& a) { return a / value; });
        return out;
    }

    Raster& operator+=(Number value)
    {
        std::for_each(data_, data_ + (cols_ * rows_),
                      [&value](Number& a) { a += value; });
        return *this;
    }

    Raster& operator-=(Number value)
    {
        std::for_each(data_, data_ + (cols_ * rows_),
                      [&value](Number& a) { a -= value; });
        return *this;
    }

    Raster& operator*=(double value)
    {
        std::for_each(data_, data_ + (cols_ * rows_),
                      [&value](Number& a) { a *= value; });
        return *this;
    }

    Raster& operator/=(double value)
    {
        std::for_each(data_, data_ + (cols_ * rows_),
                      [&value](Number& a) { a /= value; });
        return *this;
    }

    Raster& operator+=(const Raster& image)
    {
        for_each_zip(data_, data_ + (cols_ * rows_), image.data_,
                     [](Number& a, Number& b) { a += b; });
        return *this;
    }

    Raster& operator-=(const Raster& image)
    {
        for_each_zip(data_, data_ + (cols_ * rows_), image.data_,
                     [](Number& a, Number& b) { a -= b; });
        return *this;
    }

    Raster& operator*=(const Raster& image)
    {
        for_each_zip(data_, data_ + (cols_ * rows_), image.data_,
                     [](Number& a, Number& b) { a *= b; });
        return *this;
    }

    Raster& operator/=(const Raster& image)
    {
        for_each_zip(data_, data_ + (cols_ * rows_), image.data_,
                     [](Number& a, Number& b) { a /= b; });
        return *this;
    }

    bool operator==(const Raster& other) const
    {
        // TODO: assumes same sizes
        for (unsigned i = 0; i < cols_; i++) {
            for (unsigned j = 0; j < cols_; j++) {
                if (this->data_[i * cols_ + j] != other.data_[i * cols_ + j])
                    return false;
            }
        }
        return true;
    }

    bool operator!=(const Raster& other) const
    {
        // TODO: assumes same sizes
        for (unsigned i = 0; i < cols_; i++) {
            for (unsigned j = 0; j < cols_; j++) {
                if (this->data_[i * cols_ + j] != other.data_[i * cols_ + j])
                    return true;
            }
        }
        return false;
    }

    friend inline Raster operator*(double factor, const Raster& image)
    {
        return image * factor;
    }

    friend inline Raster pow(Raster image, double value) {
        image.for_each([value](Number& a){a = std::pow(a, value);});
        return image;
    }
    friend inline Raster sqrt(Raster image) {
        image.for_each([](Number& a){a = std::sqrt(a);});
        return image;
    }

    friend inline std::ostream& operator<<(std::ostream& stream, const Raster& image) {
        stream << "[[";
        for (unsigned i = 0; i < image.rows_; i++) {
            if (i != 0)
                stream << "],\n [";
            for (unsigned j = 0; j < image.cols_; j++) {
                if (j != 0)
                    stream << ", ";
                stream << image.data_[i * image.cols_ + j];
            }
        }
        stream << "]]\n";
        return stream;
    }
};

} // namespace pops

#endif // POPS_RASTER_HPP
