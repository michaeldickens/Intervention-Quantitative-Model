/* 
 * 
 * QuantitativeModel.cpp
 * ---------------------
 * 
 * Author: Michael Dickens <mdickens93@gmail.com>
 * Created: 2016-05-19
 * 
 * Port of my quantitative model spreadsheet for people who don't have
 * Excel.
 * 
 */

#include "QuantitativeModel.h"

using namespace std;

typedef map<string, Distribution> Table;

double Distribution::posterior(const Distribution& measurement) const
{
    double c = integral(measurement, false);
    return integral(measurement, true) / c;
}

/*
 * Converts an 80% credence interval into a log-normal distribution.
 */
Distribution CI(double lo, double hi)
{
    double p_m = sqrt(lo * hi);
    double p_s = sqrt(log(hi / lo) / log(10) / 2 / NORMAL_90TH_PERCENTILE);
    Distribution res(p_m, p_s);
    return res;
}

Distribution CI(double val)
{
    return CI(val, val);
}

Table read_input(string filename)
{
    Table table;
    ifstream file(filename);
    string key, comments, low_CI, high_CI;
    while (file.good()) {
        getline(file, key, ',');
        getline(file, low_CI, ',');
        getline(file, high_CI, ',');
        getline(file, comments);
        table[key] = CI(stof(low_CI), stof(high_CI));
    }
    return table;
}

void globals(Table& table)
{
    table["factory-farmed animal wellbeing"] = CI(6);
    table["factory-farmed animal sentience adjustment"] = CI(0.3);
    table["utility per factory-farmed animal"] =
        table["factory-farmed animal wellbeing"]
        * table["factory-farmed animal sentience adjustment"];
    table["utility per cage removed"] = CI(0.3);

    table["P(stay on earth)"] = CI(0.2);
    table["P(we reduce WAS on balance)"] = CI(0.7);
    table["P(fill universe with biology)"] = CI(0.4);
    table["P(society doesn't care about animals)"] = CI(0.8);
    table["P(we have factory farming)"] = CI(0.2);
    table["P(we spread WAS)"] = CI(0.4);
    table["P(we make suffering simulations)"] = CI(0.3);
    table["P(fill universe with computers)"] = CI(0.4);
    table["P(hedonium)"] = CI(0.05);
    table["P(ems)"] = CI(0.3);
    table["P(paperclip)"] = CI(0.649);
    table["P(dolorium)"] = CI(0.001);
}

// TODO: not currently called
void ev_far_future(Table& table)
{
    table["P(humans exist)"] = table["P(fill universe with biology)"];
    table["P(hedonium exists)"] =
        table["P(fill universe with computers)"] * table["P(hedonium)"];
    table["p(ems exist)"] =
        table["p(fill universe with computers)"] * table["p(ems)"];

    
    table["p(paperclips exist)"] =
        table["p(fill universe with computers)"] * table["p(paperclip)"];
    table["p(dolorium exists)"] =
        table["p(fill universe with computers)"] * table["p(dolorium)"];
}

double thl_posterior_direct(Table& table, const Distribution& prior)
{
    table["THL years factory farming prevented per $1000"] = CI(700, 13000);

    Distribution utility_estimate =
        table["THL years factory farming prevented per $1000"]
        * table["utility per factory-farmed animal"];
    printf("%f %f\n", utility_estimate.p_m, pow(utility_estimate.p_s, 2));
    return prior.posterior(utility_estimate);
}

double cage_free_posterior_direct(Table& table, const Distribution& prior)
{
    table["cage-free total expenditures ($M)"] = CI(2, 3);
    table["years until cage-free would have happened anyway"] = CI(5, 10);
    table["millions of cages prevented"] = CI(100, 150);
    table["proportion of change attributable to campaigns"] = CI(0.7, 1);
    table["cage-free years per cage prevented"] = CI(1, 1);
    Distribution utility_estimate =
        table["cage-free total expenditures ($M)"].reciprocal()
        * table["years until cage-free would have happened anyway"]
        * table["millions of cages prevented"]
        * table["proportion of change attributable to campaigns"]
        * table["cage-free years per cage prevented"]
        * table["utility per cage removed"]
        * 1000;
    Distribution r = table["cage-free total expenditures ($M)"].reciprocal();
    printf("%f %f\n", utility_estimate.p_m, pow(utility_estimate.p_s, 2));
    return prior.posterior(utility_estimate);
}

int main(int argc, char *argv[])
{
    Table table;
    globals(table);
    printf("%f\n", thl_posterior_direct(table, Distribution(1, 0.75))); // 194.8
    printf("%f\n", cage_free_posterior_direct(table, Distribution(1, 0.75))); // 2531
    return 0;
}
