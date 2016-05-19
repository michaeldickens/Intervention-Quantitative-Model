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
#include <vector>

class Distribution {
public:
    std::function<double(double)> pdf;
    std::vector<double> buckets;

    /* May be one of "buckets", "lognorm" */
    std::string type;

    Distribution();
    Distribution(std::function<double(double)> pdf);
    virtual double& operator[](int index);
    virtual double operator[](int index) const;
    virtual double get(int index) const;
    Distribution *operator+(const Distribution *other) const;
    Distribution *operator*(const Distribution *other) const;
    double mean();
    double variance(double mean1);
    double integrand(Distribution *measurement, int index, bool ev);
    double integral(Distribution *measurement, bool ev);
    double posterior(Distribution *measurement);
};

class LognormDist : public Distribution {
public:
    
    /*
     * p_m: exponent of mean of underlying Normal dist
     * p_s: base-10 standard deviation of underlying Normal dist
     */
    double p_m;
    double p_s;

    LognormDist(double p_m, double p_s);
};
