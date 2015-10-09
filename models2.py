# -*- coding: utf-8 -*-
"""
Created on Fri Oct  9 12:21:16 2015

@author: Betty
"""

import numpy as np
from multiprocesing import Pool
import matplotlib.pyplot as plt

class PolyModel(object):
    def __init__(self, population_size, selective_coefficients, frequency = None):
        self.population_size = population_size
        self.selective_coefficients = selective_coefficients
        self.frequency = frequency
        self.fitness = selective_coefficients + 1.
        
    def update_population(self, freq = 0):
        N = self.population_size
        f = self.fitness
        self.probability = N*freq*f/(N*freq*f+N*(1-freq))
        #return probability
        
        
    def update_frequency(self, current_freq = 0):
        self.update_population(current_freq)
        #self.update_populationWright(current_freq) #swap with update frequency
        #current_freq = np.sum(self.population == f)/float(self.population_size)
        #frequency.append(current_freq)
        
    
        
    def simulation(self):
        self.population = np.ones(self.population_size)
        f = self.fitness
        self.population[0] = f
        current_freq = 1./self.population_size
        #frequency = [current_freq]
        while current_freq > 0 and current_freq < 1:
            self.update_frequency(current_freq)
            
    
class Factory:
    def getModel(type, population, coefficients, freq = None):
        if type == "Moran":
            return Moran(population, coefficients, freq)
        elif type == "Wright":
            return Wright(population, coefficients, freq)
        assert 1, "Bad assertion " + type
#updates the population and the frequency
class Wright(PolyModel):
    
    def __init__(self, population_size, selective_coefficients, frequency = None):
        super().__init__(population_size, selective_coefficients, frequency = None)
    
    def update_population(self, freq):
        super().update_population(freq)
        N = self.population
        new_member_is_mutant = np.random.binomial(1, self.probability, N)
        self.population_size = 1.+self.selective_coefficient*new_member_is_mutant
        
    def update_frequency(self, current_freq = 0):
        #super().update_frequency(n)
        f = self.fitness
        self.update_population(current_freq) #swap with update frequency
        current_freq = np.sum(self.population == f)/float(self.population_size)
        self.frequency.append(current_freq)
        
        
    
    def fixation_probability(self):
        N = self.population_size
        s = self.selective_coefficient
        return (1-np.exp(-2.*s))/(1-np.exp(-4.*N*s))
        
class WrightPopulation(PolyModel): #does not save the updated population. Just the frequency
    def __init__(self, population_size, selective_coefficients, frequency = None):
        super().__init__(population_size, selective_coefficients, frequency = None)
    
    def update_population():
        pass
    
    def update_frequency(freq):
        super().update_frequency(freq)
        #to be finished
        
    def fixation_probability(self):
        N = self.population_size
        s = self.selective_coefficient
        return (1-np.exp(-2.*s))/(1-np.exp(-4.*N*s))
"""This class of Moran updates the population as it updates the frequency. Update frequecny calls
parent update frequency"""       
class Moran(PolyModel):
    def __init__(self,population_size, selective_coefficients, frequency = None):
        super().__init__(population_size, selective_coefficients, frequency = None)
    
    def update_population(self, freq):
        super().update_population(freq)
        N = self.population_size
        member_to_replace = np.random.randint(N)
        new_member_is_mutant = np.random.binomial(1, self.probability)
        new_member_fitness = 1.+self.selective_coefficient*new_member_is_mutant
        self.population[member_to_replace] = new_member_fitness
      
    #takes current_freq, but is not used. Will edit that
      #used now
    def update_frequency(self, current_freq = None):
        #super().update_frequency(current_freq)
        f = self.fitness
        self.update_population(self.current_freq) #swap with update frequency
        current_freq = np.sum(self.population == f)/float(self.population_size)
        self.frequency.append(current_freq)
    
    def fixation_probability(self):
        N = self.population_size
        s = self.selective_coefficient
        return (1.-np.exp(-s))/(1.-np.exp(-N*s))
        
#########################################################################
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

if __name__ == '__main__':
    model = Factory()
    model.getModel("Moran", 200, 1.1)
    N = 200
    s = .5
    number_of_simulations = 2000
    title_str = ''
    
