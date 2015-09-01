from SALib.sample import saltelli
from SALib.analyze import sobol
from collections import OrderedDict
from pysb.integrate import odesolve
from pysb.bng import generate_equations
from functools import partial
import multiprocessing
import SALib.analyze
import numpy
import copy_reg
import types

def _pickle_method(m):
    if m.im_self is None:
        return getattr, (m.im_class, m.im_func.func_name)
    else:
        return getattr, (m.im_self, m.im_func.func_name)

copy_reg.pickle(types.MethodType, _pickle_method)

class sensitivity_analysis:
    def __init__(self,model):
        self.model                  = model
        self.tspan                  = None
        self.problem                = None
        self.param_sets             = None
        self.sim                    = None
        self.output                 = None
        self.sa_result              = None
        
    def sensitivity(self,tspan=None, parameters_ref=None, sp_SA=None, N=1 , verbose=True):
      
        if tspan is not None:
            self.tspan = tspan
        elif self.tspan is None:
            raise Exception("'time t' must be defined.")
        
        if sp_SA is None:
            raise Exception("A species to do the sensitivity analysis on must be defined")
        if sp_SA not in [str(sp) for sp in self.model.species] and sp_SA  not in [str(obs.name) for obs in self.model.observables]:
            raise Exception("Species is not in model species")

        ref = odesolve(self.model, self.tspan, param_values=parameters_ref)

        if verbose: print "Getting parameters information"
        self.pars_info(parameters_ref, N=N)
        
        if verbose: print "Simulating with parameters from sampling"
        p=multiprocessing.Pool(2)
        func = partial(self.sim_model,sp_SA=sp_SA)
        self.sim = p.map(func, self.param_sets)
        self.output = [sum((r - ref[sp_SA])**2) for r in self.sim]
          
        if verbose: print "Sensitivity analysis"
        self.sa_result = sobol.analyze(self.problem, numpy.array(self.output), print_to_console=True)
        
    def pars_info(self, parameters_ref=None, upper_bound=10, lower_bound=0.1, N=1):
        
    #     sam_method = getattr(sample, sample_method)
        if parameters_ref is not None:
            # accept vector of parameter values as an argument
            if len(parameters_ref) != len(self.model.parameters):
                raise Exception("parameters_ref must be the same length as model.parameters")
            if not isinstance(parameters_ref, numpy.ndarray):
                parameters_ref = numpy.array(parameters_ref)
        else:
            # create parameter vector from the values in the model
            parameters_ref = numpy.array([p.value for p in self.model.parameters])
        new_pars = OrderedDict((p.name, parameters_ref[i]) for i, p in enumerate(self.model.parameters))
        #parameter information
        self.problem = {'num_vars':len(new_pars), 'names':new_pars.keys(), 'bounds' : [[i*lower_bound,i*upper_bound] for i in new_pars.values()]}
        self.param_sets = saltelli.sample(self.problem,N)
        return self.param_sets, self.problem

    
    def sim_model(self,par, sp_SA=None):   
        y = odesolve(self.model,self.tspan,par) 
        return y[sp_SA]  

def run_SA(model, tspan, parameters_ref=None, sp_SA = None, N=1):
    SA = sensitivity_analysis(model)
    SA.sensitivity(tspan, parameters_ref, sp_SA, N=N)
    return SA.sa_result


        
        
