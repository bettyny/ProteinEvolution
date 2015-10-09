# -*- coding: utf-8 -*-
"""
Created on Mon Sep 21 14:07:22 2015

@author: sshank
"""

import numpy as np
import scipy.stats as spst
from multiprocessing import Pool

class SinglePolyModel(object):
    # models with a single mutant allele segregating against a population
    def __init__(self, population_size, selective_coefficient):
        self.population_size = population_size
        self.selective_coefficient = selective_coefficient
        self.fitness = 1.+selective_coefficient
    
    def simulation(self):
        # return gene frequency from a single simulation
        self.population = np.ones(self.population_size)
        f = self.fitness
        self.population[0] = f
        current_freq = 1./self.population_size
        freq = [current_freq]
        while current_freq > 0 and current_freq < 1:
            self.update_population(current_freq)
            current_freq = np.sum(self.population == f)/float(self.population_size)
            freq.append(current_freq)
        return np.array(freq)
            

def seeded_simulation(input_tuple):
    np.random.seed(input_tuple[0])
    return input_tuple[1].simulation()

def multiple_simulations(model, number_of_simulations, initial_seed=1):
    p = Pool(12)
    arg_pairs = zip(initial_seed*number_of_simulations+np.array(range(number_of_simulations)),
                    number_of_simulations*[model])
    return p.map(seeded_simulation, arg_pairs)
    
def simulated_fixation_probability(all_sims):
    return sum([sim[-1] == 1 for sim in all_sims])/float(len(all_sims))


class MoranModel(SinglePolyModel):
    # Moran process
    def __init__(self, population_size, selective_coefficient):
        super(MoranModel, self).__init__(population_size, selective_coefficient)
    
    def update_population(self, freq):
        N = self.population_size
        f = self.fitness
        probability = N*freq*f/(N*freq*f+N*(1-freq))
        member_to_replace = np.random.randint(N)
        new_member_is_mutant = np.random.binomial(1, probability)
        new_member_fitness = 1.+self.selective_coefficient*new_member_is_mutant
        self.population[member_to_replace] = new_member_fitness
        
    def fixation_probability(self):
        N = self.population_size
        s = self.selective_coefficient
        return (1.-np.exp(-s))/(1.-np.exp(-N*s))


class WrightFisherModel(SinglePolyModel):
    def __init__(self, population_size, selective_coefficient):
        super(WrightFisherModel, self).__init__(population_size, selective_coefficient)
        
    def update_population(self, freq):
        N = self.population_size
        f = self.fitness
        probability = N*freq*f/(N*freq*f+N*(1-freq))
        new_member_is_mutant = np.random.binomial(1, probability, N)
        self.population = 1.+self.selective_coefficient*new_member_is_mutant
        
    def classic_fixation_probability(self):
        N = self.population_size
        s = self.selective_coefficient
        return (1-np.exp(-2.*s))/(1-np.exp(-2.*N*s))
    
    def new_fixation_probability(self):
        N = self.population_size
        f = self.fitness
        A = np.zeros((N-1,N-1))
        b = np.zeros(N-1)
        
        for i in np.arange(N-1):
            p = lambda x: N*f*x/(N*f*x+N*(1.-x))
            binomial = spst.binom(N, p((i+1.)/float(N)))
            b[i] = binomial.pmf(N)
            for j in np.arange(N-1):
                A[i,j] = binomial.pmf(j+1)
        x = np.linalg.solve(np.eye(N-1)-A, b)
        return x[0]

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    model = WrightFisherModel(200, 1.1)
#    model = MoranModel(200, 1.)
#    number_of_simulations = 1000
#    all_sims = multiple_simulations(model, number_of_simulations)
#    for i in range(number_of_simulations):
#        plt.plot(all_sims[i], alpha=.1, color='blue')
#    plt.plot(model.simulation())
#    plt.show()
    
    N = 200
    s = .5
    number_of_simulations = 2000
    title_str = ''
    for model, color in [(WrightFisherModel(N, s), 'blue'), (MoranModel(N, s), 'red')]:
        all_sims = multiple_simulations(model, number_of_simulations, initial_seed=1)
        if color == 'blue':
            all_sims = [np.kron(sim, np.ones(N)) for sim in all_sims]
        for i in range(number_of_simulations):
            plt.plot(all_sims[i], alpha=.05, color=color)
        ax = plt.gca()
        title_str += 'Simulated: %.2f, predicted:%.2f (%s)' % \
            (simulated_fixation_probability(all_sims), model.fixation_probability(), color)
    ax.set_title(title_str)
    plt.show()
