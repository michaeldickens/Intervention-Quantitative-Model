"""

EV.py
----------------

Author: Michael Dickens <mdickens93@gmail.com>
Created: 2016-03-07

Probability distribution of intervention cost-effectiveness.

""" 

from integrate import *
from mpmath import mp
import numpy as np
import scipy.integrate as integrate
import scipy.special as special
import scipy.stats as stats

mp.dps = 50

"""
Converts an 80% CI into a lognormal distribution with mean and
sigma that can be input into `combine_normals()`.
"""
def CI(lo, hi=None, confidence=0.8):
    if hi == None:
        hi = lo
    mean = np.sqrt(lo * hi)
    sigma = np.log10(hi / lo) / 2 / stats.norm.ppf(1 - (1 - confidence)/2)
    return (mean, sigma)

def normal_to_CI(mean, sigma, confidence=0.8):
    [(mu, s)] = to_normal([(mean, sigma)])
    lo = stats.norm.ppf((1 - confidence)/2, mu, s)
    hi = stats.norm.ppf(1 - (1 - confidence)/2, mu, s)
    return (10**lo, 10**hi)

def log10s(xs):
    return [ int(round(np.log10(x))) for x in xs ]

"""
Takes a list of pairs (mean, stdev as orders of magnitude) and
transforms it to describe a normal distribution.
"""
def to_normal(pairs):
    return [ (np.log10(m), s) for (m, s) in pairs ]

"""
Converts a value from base-10 orders of magnitude to base-e.
"""
def base_e(s):
    return s * np.log(10)

"""
Takes `mu` parameters as medians for a log-normal distribution. Takes
the log of all the `mu`s and then returns the exponent of the combined
`mu`.

Takes standard deviations in terms of orders of magnitude
(i.e. base 10).
"""
def product_normals(*pairs):
    pairs = to_normal(pairs)
    
    mean = sum([ m / s**2 for (m, s) in pairs ]) / \
           sum([ 1.0 / s**2 for (m, s) in pairs ])
    sigma = 1 / np.sqrt(sum([ 1 / s**2 for (m, s) in pairs ]))
    return (10**mean, sigma)

"""
Takes means and sigmas in the same format as `product_normals`.
"""
def sum_normals(*pairs):
    pairs = to_normal(pairs)

    mean = sum([ m for (m, s) in pairs ])
    sigma = np.sqrt(sum([ s**2 for (m, s) in pairs ]))
    return (10**mean, sigma)

def log_posterior_EV(log_mu_p, s_p, log_m, s_m):
    return product_normals([(log_mu_p, s_p), (log_m, s_m)])

def posterior_EV(mu_p, s_p, m, s_m):
    return np.exp(log_posterior_EV(np.log(mu_p), s_p, np.log(m), s_m))

"""
Far future posterior: Takes estimate in terms of QALYs per
$1. Converts to QALYs per $1000 (i.e. number of times better than
GiveDirectly) before estimating the posterior.

"""
def ff_posterior(estimate):
    estimate = sum_normals(estimate, CI(1000))
    print normal_to_CI(*estimate)
    alpha = 1.5
    return pareto_posterior(alpha, estimate)

def net_EV(pos, neg):
    return ff_posterior(pos) - ff_posterior(neg)

"""
Prior is set of plausibly good charitable interventions. Assumes
that half are better than GiveDirectly.
"""
def direct():
    prior = (1, 1)
    print "GiveDirectly", product_normals(prior, (1, 0.25))[0]
    print "AMF", product_normals(prior, (11, 0.5))[0]
    print "cage-free campaigns", product_normals(prior, (22e3, 1))[0] 
    print "THL + far future", product_normals(prior, (13e3, 1.3))[0] \
        + product_normals(prior, (1e8, 2))[0]
    print "AI safety", product_normals(prior, (1e8, 2))[0]
    # Note: If you increase sigma to 3, the point at which AI safety looks better than corporate campaigns is 1e22 utility.

def career():
    print "E2G", product_normals((50, 0.5), (80, 0.25))[0]
    print "GiveWell", product_normals((50, 0.5), (110, 0.78))[0]
    # (0.78 means 80% confidence that the true value is +/- one order of magnitude)

def fundratios():
    # prior fundratio is 2:1
    print "REG", product_normals((3, 0.75), (10, 0.5))[0]
    print "GWWC", product_normals((3, 0.75), (600, 1.5))[0]
    print "GiveWell", product_normals((3, 0.75), (20, 0.8))[0]

# TODO: Add distributions for negative effects and subtract, like for
# GiveDirectly.
def far_future_estimates():
    # far future outcome probabilities
    p_stay_on_earth = CI(0.2)
    p_biological = CI(0.7)
    p_computronium = CI(0.1)
    p_not_care_about_animals = CI(0.8)
    p_WAS_exists = CI(0.4)
    p_suffering_simulations = CI(0.3)
    p_hedonium = CI(0.2)
    p_paperclip = CI(0.75)
    p_dolorium = CI(0.05)
    
    # computronium
    computronium_accessible_stars = CI(1e13, 1e14)
    biological_accessible_stars = CI(1e12, 1e14)
    usable_wattage = CI(1e20, 1e25)
    hedonium_brains_per_watt = CI(1e2)
    paperclip_brains_per_watt = CI(1e-2, 1e1)
    dolorium_brains_per_watt = hedonium_brains_per_watt
    years_of_future = CI(1e11, 1e12)
    humans_per_star = CI(1e10, 1e12)

    # AI safety
    # p_AI_related_extinction = (0.1, 0.5) # mean/sigma, not CI
    p_AI_related_extinction = CI(0.03, 0.3)
    size_of_FAI_community = CI(1.0/10000, 1.0/200)
    ai_safety_multiplicative_effect = CI(0.5, 3)
    ai_bad_scenarios_averted = CI(0.1, 0.7)
    ai_hours_research = CI(1e-9 * 2000, 1e-6 * 2000)
    cost_per_AI_researcher = CI(1.0/1.5e5, 1.0/7e4)

    # Effect of alleviating factory farming...
    factory_farming_today = CI(1e-10)
    cost_per_animal = CI(1.0/50, 1.0/5)

    # ...on factory farming (biological condition only)
    p_factory_farming_exists = CI(0.2)
    p_people_care_about_factory_farming = CI(0.5)
    factory_farming_per_star = CI(1e10, 1e13)
    qalys_per_farm_year = CI(0.2, 4)
    
    # ...on WAS (biological condition only)
    p_people_care_about_WAS = CI(0.1)
    qalys_per_WAS_year = CI(1)
    wild_animals_per_star = CI(1e14, 1e18) # adjusted by sentience

    # ...on subroutines (biological condition)
    p_biological_create_suffering_subroutines = CI(0.1)
    p_people_care_about_subroutines = CI(0.1)
    qalys_per_subroutine_year = qalys_per_WAS_year
    biological_subroutines_per_wild_animal = CI(0.001, 1)

    # AI-targeted values spreading
    p_build_FAI = CI(0.1)
    p_FAI_bad_for_animals = CI(0.2)
    p_animals_exist_given_FAI = CI(0.1)
    ai_researcher_values_propagation = CI(1, 3)
    cost_to_convince_AI_researcher = CI(1e-6, 1e-3)

    # GiveDirectly
    dollars_to_end_poverty = CI(1e-12, 1e-11)
    p_extinction = CI(0.3)
    end_poverty_x_risk_reduction = CI(0.01, 0.2)
    end_poverty_x_risk_increase = CI(0.008, 0.16)
    
    hedonium_value = sum_normals(
        computronium_accessible_stars,
        usable_wattage,
        years_of_future,
        hedonium_brains_per_watt,
    )

    paperclip_neg_value = sum_normals(
        computronium_accessible_stars,
        usable_wattage,
        years_of_future,
        paperclip_brains_per_watt,
    )

    dolorium_neg_value = sum_normals(
        computronium_accessible_stars,
        usable_wattage,
        years_of_future,
        dolorium_brains_per_watt,
    )

    biological_human_value = sum_normals(
        biological_accessible_stars,
        humans_per_star,
        years_of_future,
    )

    factory_farming_in_future = sum_normals(
        biological_accessible_stars,
        factory_farming_per_star,
        years_of_future,
    )

    wild_animals_in_future = sum_normals(
        biological_accessible_stars,
        wild_animals_per_star,
        years_of_future,
    )

    wild_animals_neg_value = sum_normals(
        wild_animals_in_future,
        qalys_per_WAS_year,
    )

    biological_subroutines_in_future = sum_normals(
        biological_accessible_stars,
        wild_animals_per_star,
        biological_subroutines_per_wild_animal,
        years_of_future,
    )

    biological_subroutines_neg_value = sum_normals(
        biological_subroutines_in_future,
        qalys_per_subroutine_year,
    )
    
    # Multiply by value of far future, cost per researcher to get
    # final result
    ai_model1 = sum_normals(
        p_AI_related_extinction,
        size_of_FAI_community,
        ai_safety_multiplicative_effect,
        ai_bad_scenarios_averted,
    )

    ai_model2 = sum_normals(
        p_AI_related_extinction,
        ai_hours_research,
    )

    ai_average = product_normals(ai_model1, ai_model2)
    ai_average = (ai_average[0], ai_average[1] * 1.1) # adjust for models being non-independent
    ai_safety = sum_normals(
        ai_average,
        hedonium_value,
        cost_per_AI_researcher,
    )

    ai_safety_1 = sum_normals(
        ai_model1,
        hedonium_value,
        cost_per_AI_researcher,
    )

    factory_farming = sum_normals(
        p_factory_farming_exists,
        p_people_care_about_factory_farming,
        factory_farming_in_future,
        qalys_per_farm_year,
        factory_farming_today,
        cost_per_animal,
    )

    was = sum_normals(
        p_WAS_exists,
        p_people_care_about_WAS,
        wild_animals_neg_value,
        factory_farming_today,
        cost_per_animal,
    )

    # TODO: redo this
    biological_suffering_subroutines = sum_normals(
        p_biological_create_suffering_subroutines,
        p_people_care_about_subroutines,
        biological_subroutines_neg_value,
        factory_farming_today,
        cost_per_animal,
    )

    (hm, hs) = hedonium_suffering_subroutines
    (bm, bs) = biological_suffering_subroutines
    suffering_subroutines_overall = ((hm + bm) / 2, np.sqrt((hs**2 + bs**2)/2))

    ai_targeted = sum_normals(
        p_build_FAI,

        p_FAI_bad_for_animals,
        p_animals_exist_given_FAI,
          # number of suffering subroutines
          hedonium_value,
          p_hedonium_create_suffering,
          proportion_suffering_minds,
        size_of_FAI_community,
        ai_researcher_values_propagation,
        cost_to_convince_AI_researcher,
    )

    givedirectly = sum_normals(
        hedonium_value,
        dollars_to_end_poverty,
        p_extinction,
        end_poverty_x_risk_reduction,
    )

    givedirectly_neg = sum_normals(
        hedonium_value,
        dollars_to_end_poverty,
        p_extinction,
        end_poverty_x_risk_increase,
    )

    # Expected values after accounting for probability of event
    earth_EV = sum_normals(
        p_stay_on_earth,
        CI(0.4),
        wild_animals_per_star,
        qalys_per_WAS_year,
    )
    bio_humans_EV = sum_normals(p_biological, biological_human_value)
    was_neg_EV = sum_normals(
        p_biological,
        p_not_care_about_animals,
        p_WAS_exists,
        wild_animals_neg_value,
    )
    bio_simulations_neg_EV = sum_normals(
        p_biological,
        p_not_care_about_animals,
        p_suffering_simulations,
        biological_subroutines_neg_value,
    )
    hedonium_EV = sum_normals(p_computronium, p_hedonium, hedonium_value)
    paperclip_neg_EV = sum_normals(p_computronium, p_paperclip, paperclip_neg_value)
    dolorium_neg_EV = sum_normals(p_computronium, p_dolorium, dolorium_neg_value)

    for ev in ['earth_EV', 'bio_humans_EV', 'was_neg_EV', 'bio_simulations_neg_EV',
               'hedonium_EV', 'paperclip_neg_EV', 'dolorium_neg_EV']:
        lo, hi = normal_to_CI(*eval(ev))
        print "| %s | %.1e | %.1e |" % (ev, lo, hi)

    # print "AI safety", ff_posterior(ai_safety)
    # print "factory farming", ff_posterior(factory_farming)
    # print "WAS", ff_posterior(was)
    # print "subroutines", ff_posterior(suffering_subroutines_overall)
    # print "AI targeted", ff_posterior(ai_targeted)
    # print "GiveDirectly", net_EV(givedirectly, givedirectly_neg)
    

"""
Helper function for computing combined Lomax/Lognormal distribution.
"""
def core_pdf(measurement, u, alpha):
    m, s = measurement
    return stats.lomax.pdf(u, alpha, scale = 1 / (2**(1.0/alpha) - 1)) * \
        stats.lognorm.pdf(m, base_e(s), scale=u)

def measurement_int(measurement, alpha, hi=np.inf):
    return integral(lambda u: core_pdf(measurement, u, alpha),
                    0.01, hi)[0]

"""
Computes posterior expected value.
"""
def pareto_posterior(alpha, measurement):
    c = measurement_int(measurement, alpha)
    return (integral(lambda u: u * core_pdf(measurement, u, alpha),
                     0.01, np.inf)[0]) / c

def lognorm_posterior(prior, measurement):
    return product_normals(prior, measurement)[0]

"""
Computes P(U < hi | M = m).
"""
def posterior_cdf(measurement, hi):
    return measurement_int(measurement, hi) / measurement_int(measurement)

"""
Prints a table of posteriors for a bunch of different parameters.
"""
def print_table(f, param):
    ms = [1, 11, 2e4, 1e39]
    ss = [0.25, 0.5, 1, 2]
    print "| s | GiveDirectly | AMF | animals | x-risk |"
    print "|---|--------------|-----|---------|--------|"
    for s in ss:
        print "| %.2f |" % s,
        for m in ms:
            p = f(param, (m, s))
            if p > 10000:
                print "%.2e |" % p,
            else:
                print "%.2f |" % p,
        print ""

def print_tables():
    print "### Pareto prior\n"
    for alpha in [1.1, 1.5, 2]:
        print "&alpha; = %.1f:\n" % alpha
        print_table(pareto_posterior, alpha)
        print "<br />\n"
        
    print "### Log-normal prior\n"
    for prior in [(1, 0.5), (1, 1), (1, 1.5)]:
        print "&sigma; = %.1f:\n" % prior[1]
        print_table(lognorm_posterior, prior)
        print "<br />\n"

