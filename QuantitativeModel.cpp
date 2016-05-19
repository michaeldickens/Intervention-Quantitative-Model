/* 
 * 
 * QuantitativeModel.cpp
 * ---------------------
 * 
 * Author: Michael Dickens <mdickens93@gmail.com>
 * Created: 2016-05-18
 * 
 * Port of my quantitative model spreadsheet for people who don't have
 * Excel.
 * 
 */

#include "QuantitativeModel.h"

using namespace std;

#define NUM_BUCKETS 1000
#define STEP 1.25
#define EXP_OFFSET 100
#define NORMAL_90TH_PERCENTILE 1.2815515655

function<double(double)> lomax_pdf(double x_m, double alpha)
{
    return [x_m, alpha](double x){
        return alpha / x_m / pow(1 + x / x_m, alpha + 1);
    };
}

function<double(double)> lognorm_pdf(double p_m, double p_s)
{
    double mu = log(p_m);
    double s = log(10) * p_s;
    return [mu, s](double x){
        return 1 / (x * s * sqrt(2 * M_PI)) *
            exp(-0.5 * pow((log(x) - mu) / s, 2));
    };
}

/*
 * Returns the bucket containing the given x value.
 */
int bucket_index(double x)
{
    // return (int) (log(x) / log(STEP) - 0.5) + EXP_OFFSET; // I believe this is more mathematically accurate but Excel does in the other way
    return (int) (log(x) / log(STEP)) + EXP_OFFSET;
}

/*
 * Gives probability density in the (geometric) middle of the
 * bucket. Sort-of the opposite of `bucket_index`.
 */
double bucket_prob(int index)
{
    return pow(STEP, index - EXP_OFFSET + 0.5);
}

/*
 * Gives probability density at the bottom of the bucket.
 */
double bucket_min_prob(int index)
{
    return pow(STEP, index - EXP_OFFSET);
}

/*
 * Gets the difference in x-values from the beginning to the end of a bucket.
 */
double get_delta(int index)
{
    return pow(STEP, index - EXP_OFFSET + 1) - pow(STEP, index - EXP_OFFSET);
}

/*
 * Combine an array of distributions using some binary operator (folds
 * left).
 */
Distribution *foldl_dists(
    function<Distribution *(const Distribution *, const Distribution *)> f,
    const Distribution *dists[], int length)
{
    Distribution *curr = NULL;
    Distribution *next = NULL;
    for (int i = 1; i < length; i++) {
        if (i == 1) {
            next = f(dists[0], dists[i]);
        } else {
            next = f(curr, dists[i]);
            delete curr;
        }
        curr = next;
    }
    return curr;
}

Distribution *sum_dists(const Distribution *dists[], int length)
{
    return foldl_dists([](const Distribution *x, const Distribution *y) {
            return *x + y;
        }, dists, length);
}

Distribution *product_dists(const Distribution *dists[], int length)
{
    return foldl_dists([](const Distribution *x, const Distribution *y) {
            return *x * y;
        }, dists, length);
}

Distribution::Distribution() : buckets(NUM_BUCKETS, 0)
{
    this->type = "buckets";
}

/*
 * Constructs buckets given a probability density function (PDF).
 */
Distribution::Distribution(function<double(double)> pdf) : Distribution()
{
    for (int i = 0; i < NUM_BUCKETS; i++) {
        buckets[i] = pdf(bucket_prob(i));
    }
}

/*
 * Should only be used for distributions that use buckets.
 */
double& Distribution::operator[](int index)
{
    return buckets[index];
}

double Distribution::operator[](int index) const
{
    return buckets[index];
}

double Distribution::get(int index) const
{
    if (type == "buckets") {
        return buckets[index];
    } else {
        return pdf(bucket_prob(index));
    }
}

/*
 * Calculates the sum of two probability distributions.
 */
Distribution *Distribution::operator+(const Distribution *other) const
{
    Distribution *res = new Distribution();
    for (int i = 0; i < NUM_BUCKETS; i++) {
        for (int j = 0; j < NUM_BUCKETS; j++) {
            int index = bucket_index(bucket_prob(i) + bucket_prob(j));
            double mass = get(i) * get_delta(i) * other->get(j) * get_delta(j);
            if (index >= NUM_BUCKETS) {
                index = NUM_BUCKETS - 1;
            }
            res->buckets[index] += mass / get_delta(index);
        }
    }
    return res;
}

/*
 * Calculates the product of two probability distributions.
 */
Distribution *Distribution::operator*(const Distribution *other) const
{
    if (this->type == "lognorm" && other->type == "lognorm") {
        LognormDist *this2 = (LognormDist *) this;
        LognormDist *other2 = (LognormDist *) other;
        double new_p_m = this2->p_m * other2->p_m;
        double new_p_s = sqrt(pow(this2->p_s, 2) + pow(other2->p_s, 2));
        return new LognormDist(new_p_m, new_p_s);
    }
    
    Distribution *res = new Distribution();
    for (int i = 0; i < NUM_BUCKETS; i++) {
        for (int j = 0; j < NUM_BUCKETS; j++) {
            int index = bucket_index(bucket_prob(i) * bucket_prob(j));
            double mass = get(i) * get_delta(i) * other->get(j) * get_delta(j);
            if (index >= NUM_BUCKETS) {
                index = NUM_BUCKETS - 1;
            } else if (index < 0) {
                index = 0;
            }
            res->buckets[index] += mass / get_delta(index);
        }
    }
    return res;
}

double Distribution::mean() const
{
    double mu = 0;
    for (int i = 0; i < NUM_BUCKETS; i++) {
        mu += bucket_prob(i) * get(i) * get_delta(i);
    }
    return mu;
}

double Distribution::variance(double mean1) const
{
    double sigma2 = 0;
    for (int i = 0; i < NUM_BUCKETS; i++) {
        sigma2 += pow(bucket_prob(i), 2) * get(i) * get_delta(i);
    }
    return sigma2 - pow(mean1, 2);
}

double Distribution::integrand(const Distribution *measurement, int index, bool ev) const
{
    double u = bucket_min_prob(index);
    double prior = this->get(index);
    double update = 0;
    if (measurement->type == "buckets") {
        /* approximate `measurement` with log-normal dist */
        // TODO: just do this directly numerically
        double mean1 = measurement->mean();
        double var = measurement->variance(mean1);
        double mu = log(mean1 / sqrt(1 + var / pow(mean1, 2)));
        double sigma = sqrt(log(1 + var / pow(mean1, 2)));
        update = lognorm_pdf(u, sigma / log(10))(exp(mu));
    } else {
        LognormDist *meas = (LognormDist *) measurement;
        update = lognorm_pdf(u, meas->p_s)(meas->p_m);
    }

    double res = prior * update;
    if (ev) {
        res *= u;
    }
    return res;
}

double Distribution::integral(const Distribution *measurement, bool ev) const
{
    double total = 0;
    double x_lo = pow(STEP, -EXP_OFFSET);
    double x_hi = x_lo * STEP;
    double y_lo = integrand(measurement, 0, ev);
    double y_hi;
    double avg, delta;
    for (int i = 1; i <= NUM_BUCKETS; i++) {
        y_hi = integrand(measurement, i, ev);
        avg = (y_lo + y_hi) / 2;
        delta = x_hi - x_lo;
        total += avg * delta;

        x_lo = x_hi;
        x_hi = x_hi * STEP;
        y_lo = y_hi;
    }
    return total;
}

double Distribution::posterior(const Distribution *measurement) const
{
    double c = this->integral(measurement, false);
    return this->integral(measurement, true) / c;
}

/*
 * Takes p_m as exp(mu) and p_s and base-10 standard deviation.
 */
LognormDist::LognormDist(double p_m, double p_s)
{
    type = "lognorm";
    this->p_m = p_m;
    this->p_s = p_s;
    this->pdf = lognorm_pdf(p_m, p_s);
}

/*
 * Converts an 80% credence interval into a log-normal distribution.
 */
LognormDist CI(double lo, double hi)
{
    double p_m = sqrt(lo * hi);
    double p_s = sqrt(log(hi / lo) / log(10) / 2 / NORMAL_90TH_PERCENTILE);
    return LognormDist(p_m, p_s);
}

void test_sum_dists()
{
    LognormDist one(5.464e33, sqrt(2.73));
    LognormDist two(5.12e36, sqrt(3.12));
    LognormDist three(1.619e40, sqrt(3.51));
    LognormDist four(3.84e38, sqrt(4.68));
    LognormDist five(8.095e43, sqrt(3.51));
    LognormDist six(1.265e46, sqrt(3.51));
    LognormDist *dists[] = {
        &one, &two, &three, &four, &five, &six
    };

    LognormDist *sum = (LognormDist *) sum_dists((const Distribution **) dists, 6);
    double mean1 = sum->mean();
    double var = sum->variance(mean1);
    double mu = log(mean1 / sqrt(1 + var / pow(mean1, 2)));
    double sigma = sqrt(log(1 + var / pow(mean1, 2)));
    // mean1, var should be 1.42e50, 2.40e108
    cout << mean1 << ", " << var << endl;
    // exp(mu), (sigma / log(10))^2 should be 1.3e46, 3.51
    cout << exp(mu) << ", " << pow(sigma / log(10), 2) << endl;
    delete sum;
}

double thl_posterior_direct(Distribution prior)
{
    LognormDist years_factory_farming_prevented_per_1000_dollars = CI(700, 13000);
    LognormDist chicken_wellbeing = CI(6, 6);
    LognormDist chicken_sentience_adjustment = CI(0.3, 0.3);
    LognormDist *arr[] = {
        &years_factory_farming_prevented_per_1000_dollars,
        &chicken_wellbeing,
        &chicken_sentience_adjustment,
    };
    LognormDist *utility_estimate = (LognormDist *) sum_dists((const Distribution **) arr, 3);
    double post = prior.posterior(utility_estimate);
    delete utility_estimate;
    return post;
}

int main(int argc, char *argv[])
{
    // TODO: put user input vars into hash tables and read from file
    
    printf("%f\n", thl_posterior_direct(LognormDist(1, 0.75)));
    return 0;
}
