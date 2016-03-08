"""

charity.py
----------

Author: Michael Dickens <mdickens93@gmail.com>
Created: 2016-02-15

"""

from math import log
import scipy.integrate as integrate

class Model():

    """
    Gives total utility of donating summed over all years.
    """
    def utility(self, t_f = None):
        return int(self.diff(t_f)[0] / 1000)
    
class DonationModel(Model):

    """
    a: charity income growth rate
    c: your donation size (may be a constant or a function of time)
    d: rate of diminishing marginal utility of money
    g: charity operations growth rate
    M: maximum goal units achievable per period
    
    k: starting capital
    m: ratio of max goal units to current goal units
    U: utility per goal unit
    """
    def __init__(self, k, m, U, a=1.1, c=1e5, d=2.0/3, g=2.0):
        self.a = a
        self.c = c
        self.d = d
        self.g = g

        self.k = k = float(k)
        self.M = M = k**d * m
        self.U = U

        self.goal_units = dict() # dynamic programming bitches

    """
    The problem does harm every year. Every year when the charity
    operates, harm done by the problem reduces. This continues until
    the problem reaches 0. Disutility is the sum of harm done over all
    years.
    """
    def disutility_with_donation(self, c, t_f = None):
        a = self.a
        d = self.d
        g = self.g
        k = self.k
        M = self.M
        U = self.U
        if t_f == None:
            t_f = 999
        
        accum_progress = 0
        total_disutility = 0
        for t in range(t_f):
            if hasattr(c, '__call__'):
                c1 = c(t)
            else:
                c1 = c
            progress = min(k*a**t + c1, k*g**t)**d
            accum_progress = min(M, accum_progress + progress)
            total_disutility += U * (M - accum_progress)
            if accum_progress >= M:
                break
        
        return (total_disutility, t)

    def diff(self, t_f = None):
        w = self.disutility_with_donation(self.c, t_f)
        wo = self.disutility_with_donation(0, t_f)
        return (wo[0] - w[0], w[1], wo[1])

class FundraisingModel(Model):

    def __init__(self, donation_model, moved_per_employee, starting_employees,
                 counterfactual, starting_capital,
                 a = 1.1, k = 1e6):
        self.donation_model = donation_model
        self.moved_per_employee = moved_per_employee
        self.starting_employees = starting_employees
        self.counterfactual = counterfactual
        self.starting_capital = starting_capital
        self.a = a
        self.k = k
        
        self.capital_spending_rate = 1.0 / 20
        self.employee_cost = 5e4

        self.time_maxed = log(self.starting_capital * self.capital_spending_rate \
                              / (self.employee_cost * self.starting_employees))/log(a)
        self.U = donation_model.U

    def diff(self, t_f = None):
        w = self.donation_model.disutility_with_donation(self.moved_with_hire, t_f)
        wo = self.donation_model.disutility_with_donation(self.moved_without_hire, t_f)
        return (wo[0] - w[0], w[1], wo[1])

    def num_employees(self, t):
        # current growth rate is about 1.3, can probably continue
        # similar rate (?)
        return self.starting_employees * 1.2**t

    """
    Maximum number of employees, as constrained by the budget.
    """
    def max_employees(self, t):
        return self.starting_capital * self.capital_spending_rate \
            / self.employee_cost

        # assumes growing budget
        # a = self.a
        # k = self.k
        # return (self.starting_capital * self.capital_spending_rate \
            # + k * a**t) / self.employee_cost
            
    def moved_with_hire(self, t):
        return min(self.num_employees(t) + 1, self.max_employees(t)) \
            * self.moved_per_employee

    def moved_without_hire(self, t):
        return min(self.num_employees(t), self.max_employees(t)) \
            * self.moved_per_employee

career_length = 20

veg = DonationModel(1e7, 100, 5)
print "veg:", veg.utility(t_f=career_length)

miri = DonationModel(1e7, 10000, 20)
print "MIRI:", miri.utility(t_f=career_length)

# malaria = DonationModel(1e8, 10, 1)
# print "malaria:", utility(malaria)

malaria = DonationModel(1e8, 100, 1)
print "malaria:", malaria.utility(t_f=career_length)

counterfactual = 0.0
gw = FundraisingModel(malaria, moved_per_employee = 5e6,
                      starting_employees = 20,
                      counterfactual = counterfactual,
                      starting_capital = 2e8 # 1/10 of GV money
)
print "GiveWell (%s ctf):" % (counterfactual), gw.utility(t_f=career_length)
