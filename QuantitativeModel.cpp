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

#include <iostream>
#include <vector>

using namespace std;

#define NUM_BUCKETS 1000
#define STEP 1.25
#define EXP_OFFSET 200

/*
 * Returns the bucket containing the given x value.
 */
int bucket_index(int x) {
    return (int) (log(x) / log(STEP)) + EXP_OFFSET;
}

/*
 * Gives probability density in the (geometric) middle of the
 * bucket. Sort-of the opposite of `bucket_index`.
 */
double bucket_prob(int index) {
    return pow(STEP, index - EXP_OFFSET + 0.5);
}

/*
 * Gets the difference in x-values from the beginning to the end of a bucket.
 */
double get_delta(int index) {
    return pow(STEP, index - EXP_OFFSET + 1) - pow(STEP, index - EXP_OFFSET)
}

class Distribution {
    vector<double> buckets(NUM_BUCKETS, 0);

public:
    Distribution() {}

    /*
     * Constructs buckets given a probability density function (PDF).
     */
    Distribution(double (*pdf)(double)) {
        for (int i = 0; i < NUM_BUCKETS; i++) {
            this[i] = pdf(bucket_prob(i));
        }
    }

    double& operator[](size_t index) {
        return this.buckets[index];
    }

    /*
     * Calculates the sum of two probability distributions.
     */
    Distribution operator+(const Distribution *other) {
        Distribution res;
        for (int i = 0; i < NUM_BUCKETS; i++) {
            for (int j = 0; j < NUM_BUCKETS; j++) {
                index = bucket_index(bucket_prob(i) + bucket_prob(j));
                mass = this[i] * get_delta(i) * other[j] * get_delta[j];
                if (index >= num_buckets) {
                    index = num_buckets - 1;
                }
                res[index] += mass / get_delta(index);
            }
        }
        return res;
    }

    /*
     * Calculates the product of two probability distributions.
     */
    Distribution operator*(const Distribution *other) {
        Distribution res;
        for (int i = 0; i < NUM_BUCKETS; i++) {
            for (int j = 0; j < NUM_BUCKETS; j++) {
                int index = bucket_index(bucket_prob(i), bucket_prob(j));
                int mass = this[i] * get_delta(i) * other[j] * get_delta(j);
                if (index >= NUM_BUCKETS) {
                    index = num_buckets - 1;
                } else if (index < 0) {
                    index = 0;
                }
                res[index] += mass / get_delta(index);
            }
        }
        return res;
    }
}
