"""

EV.py
----------------

Author: Michael Dickens <mdickens93@gmail.com>
Created: 2016-03-07

Probability distribution of intervention cost-effectiveness.

""" 

import numpy as np
import scipy.integrate as integrate
import scipy.special as special
import scipy.stats as stats

robust_sigma = 0.5
mid_sigma =    1
weak_sigma =   2

def lognorm_pdf(x):
    global sigma
    return stats.lognorm.pdf(x, sigma)

def lognorm_cdf(x):
    global sigma
    return stats.lognorm.cdf(x, sigma)

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
Takes means as means for a log-normal distribution. Takes the log
of all the means and then returns the exponent of the combined mean.

Takes standard deviations in terms of orders of magnitude
(i.e. base 10). Internally converts them to base e.
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

def far_future_estimates():
    prior = (1, 1)

    # Value of far future
    accessible_stars = CI(1e13, 1e14)
    usable_wattage = CI(1e22, 1e25)
    brains_per_watt = CI(1e2)
    years_of_future = CI(1e11, 1e12)
    p_brains_are_happy = CI(0.1)
    humans_per_star = CI(1e10, 1e12)

    # AI safety
    p_AI_related_extinction = CI(0.1, 0.5)
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
    p_people_care_about_factory_farming = CI(1)
    factory_farming_in_future = CI(1e33, 1e37)
    qalys_per_farm_year = CI(0.2, 4)
    
    # ...on WAS (biological condition only)
    p_WAS_exists = CI(0.4)
    p_people_care_about_WAS = CI(0.1)
    qalys_per_WAS_year = CI(1)
    wild_animals_per_star = CI(1e16, 1e20)

    # ...on subroutines (hedonium condition)
    p_hedonium_create_suffering = CI(0.01)
    proportion_suffering_minds = CI(0.01, 0.1)
    p_people_care_about_subroutines = CI(0.1)
    qalys_per_subroutine_year = CI(1)

    # ...on subroutines (biological condition)
    p_biological_create_suffering = CI(0.5)
    biological_subroutines_per_wild_animal = CI(0.01, 1)

    # AI-targeted values spreading
    p_build_FAI = CI(0.1)
    p_FAI_bad_for_animals = CI(0.2)
    p_animals_exist_given_FAI = CI(0.1)
    ai_researcher_values_propagation = CI(1, 3)
    cost_to_convince_AI_researcher = CI(1e-6, 1e-3)
    
    hedonium_value = sum_normals(
        accessible_stars,
        usable_wattage,
        brains_per_watt,
        years_of_future,
        p_brains_are_happy,
    )

    biological_value = sum_normals(
        accessible_stars,
        humans_per_star,
        years_of_future,
    )

    wild_animals_in_future = sum_normals(
        accessible_stars,
        wild_animals_per_star,
        years_of_future,
    )

    biological_subroutines_in_future = sum_normals(
        accessible_stars,
        wild_animals_per_star,
        biological_subroutines_per_wild_animal,
        years_of_future,
    )
    
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
        cost_per_AI_researcher
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
        qalys_per_WAS_year,
        wild_animals_in_future,
        factory_farming_today,
        cost_per_animal,
    )

    hedonium_suffering_subroutines = sum_normals(
        hedonium_value,
        p_hedonium_create_suffering,
        proportion_suffering_minds,
        p_people_care_about_subroutines,
        qalys_per_subroutine_year,
        factory_farming_today,
        cost_per_animal,
    )

    biological_suffering_subroutines = sum_normals(
        p_biological_create_suffering,
        p_people_care_about_subroutines,
        biological_subroutines_in_future,
        qalys_per_subroutine_year,
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

    # print "factory farming", normal_to_CI(*factory_farming)
    # print "WAS", normal_to_CI(*was)

    print "AI safety", log10s(normal_to_CI(*product_normals(prior, ai_safety)))
    print "factory farming", log10s(normal_to_CI(*product_normals(prior, factory_farming)))
    print "WAS", log10s(normal_to_CI(*product_normals(prior, was)))
    print "subroutines", log10s(normal_to_CI(*product_normals(prior, suffering_subroutines_overall)))
    print "AI targeted", log10s(normal_to_CI(*product_normals(prior, ai_targeted)))
    

far_future_estimates()
