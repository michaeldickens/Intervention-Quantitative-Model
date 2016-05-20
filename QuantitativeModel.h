/* 
 * 
 * QuantitativeModel.h
 * -------------------
 * 
 * Author: Michael Dickens <mdickens93@gmail.com>
 * Created: 2016-05-18
 * 
 */

#define _USE_MATH_DEFINES // for M_PI

#include <cmath>
#include <functional>
#include <iostream>
#include <map>
#include <memory>
#include <vector>

#define NUM_BUCKETS 1000
#define STEP 1.25
#define EXP_OFFSET 100
#define NORMAL_90TH_PERCENTILE 1.2815515655

class Distribution {
public:
    enum class Type { buckets, lognorm };

    std::function<double(double)> pdf;
    std::vector<double> buckets;

    /* using enum instead of subclasses because C++ is stupid */
    Type type;

    /* params for log-normal */
    double p_m;
    double p_s;

    Distribution();
    Distribution(double p_m, double p_s);
    Distribution(std::function<double(double)> pdf);
    virtual double operator[](int index) const;
    virtual double get(int index) const;
    Distribution operator+(const Distribution& other) const;
    Distribution operator*(const Distribution& other) const;
    double mean() const;
    double variance(double mean1) const;
    double integrand(const Distribution& measurement, int index, bool ev) const;
    double integral(const  Distribution& measurement, bool ev) const;
    double posterior(const Distribution& measurement) const;
};
