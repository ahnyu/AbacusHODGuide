import os
import sys
import numpy as np
import argparse
import yaml
from scipy.stats import truncnorm
from nautilus import Sampler,Prior
from abacusnbody.hod.abacus_hod import AbacusHOD
from abacusnbody.hod.GRAND_HOD import *
from mpi4py.futures import MPIPoolExecutor

loaded_ball = None

def load_abacus_hod_once(sim_params, HOD_params, clustering_params):
    global loaded_ball
    if loaded_ball is None:
        loaded_ball = AbacusHOD(sim_params, HOD_params, clustering_params)
    return loaded_ball


class wp_Data(object):
    """
    Dummy object for calculating a likelihood
    """
    def __init__(self, data_params, HOD_params, nrpmin = 7, nrpmax = 21, npimax = 8):
        """
        Constructor of the power spectrum data
        """
        num_dens_mean = {}
        num_dens_std = {}
        for key in HOD_params['tracer_flags'].keys():
            if HOD_params['tracer_flags'][key]:
                num_dens_mean[key] = data_params['tracer_density_mean'][key]
                num_dens_std[key] = data_params['tracer_density_std'][key]
        self.num_dens_mean = num_dens_mean
        self.num_dens_std = num_dens_std

        rpfac = 2
        nbinrp = int(48/rpfac)
        # load the power spectrum for all tracer combinations
        clustering = {}
        wpstds = {}
        rs = {}
        for key in data_params['tracer_combos'].keys():
            clustering[key] = np.loadtxt(data_params['tracer_combos'][key]['path2wp'])[nrpmin:nrpmax, 1] # wp
            wpstds[key] = np.loadtxt(data_params['tracer_combos'][key]['path2wp'])[nrpmin:nrpmax, 2] # wp
            rs[key] = np.loadtxt(data_params['tracer_combos'][key]['path2wp'])[nrpmin:nrpmax, 0]
        self.clustering = clustering
        self.wpstds = wpstds
        self.rs = rs

        if 'path2cov' in data_params['tracer_combos'][key]:
            # load the covariance matrix for all tracer combinations
            cov = {}
            icov = {}
            for key in data_params['tracer_combos'].keys():

                newcov = np.load(data_params['tracer_combos'][key]['path2cov'])['cov']
                rescale = data_params['tracer_combos'][key].get('rescale', True)
                if rescale:
                    rescaledcov = np.zeros(newcov.shape)

                    for i in range(newcov.shape[0]):
                        for j in range(newcov.shape[1]):
                            rescaledcov[i, j] = newcov[i, j]/np.sqrt(newcov[i, i]*newcov[j, j])*wpstds[key][i]*wpstds[key][j]
                    cov[key] = rescaledcov
                else:
                    cov[key] = newcov
                
                icov[key] = np.linalg.inv(cov[key])
            self.icov = icov
            self.cov = cov


    def compute_likelihood(self, theory_clustering, theory_density, mockcov = False, ic_down = 1, jointcov_inv = None, fullscale = True):
        """
        Computes the likelihood using information from the context
        """
        skiprp = 0
        if not fullscale:
            skiprp = 4
            
        lnprob = 0.
        if jointcov_inv is None:
            for key in self.clustering.keys():
                delta = self.clustering[key] - theory_clustering[key][skiprp:]
                if mockcov:
                    lnprob += np.einsum('i,ij,j', delta, self.icov[key][skiprp:][skiprp:], delta)
                else:
                    lnprob += np.sum(delta**2/self.wpstds[key]**2)
        else:
            delta = np.concatenate([self.clustering[key] - theory_clustering[key] for key in np.sort(list(self.clustering.keys()))])
            lnprob += np.einsum('i,ij,j', delta, jointcov_inv, delta)
#        print(delta)
        lnprob *= -0.5

        # likelihood due to number density
        for etracer in self.num_dens_mean.keys():
            lnprob += -0.5*((self.num_dens_mean[etracer] - theory_density[etracer]*ic_down)/self.num_dens_std[etracer])**2

        print(" <><> Likelihood evaluated, lnprob = ",lnprob)
        return lnprob
    
    
def loglike(p, sim_params, HOD_params, clustering_params, param_mapping, mytracers, Data, nthread, bounds, mockcov, jointcov_inv, fullscale):

    if loaded_ball is None:
        print('loading data')
    Ball = load_abacus_hod_once(sim_params, HOD_params, clustering_params)
    
    for tracer_type in mytracers: 
        for key in param_mapping[tracer_type].keys():
            mapping_idx = param_mapping[tracer_type][key]
            if key == 'logsigma':
                Ball.tracers[tracer_type]['sigma'] = 10**p[mapping_idx].item()
            else:
                Ball.tracers[tracer_type][key] = p[mapping_idx].item()

    # impose Mmin cut (remove unphysical HOD, help chains to converge)
    for tracer_type in mytracers: 
        if tracer_type == 'LRG' and 10**Ball.tracers[tracer_type]['logM_cut']*Ball.tracers[tracer_type]['kappa'] < 1e12:
#            print("LRG Mmin < 1e12")
            return -np.inf
        elif tracer_type == 'ELG' and N_cen_ELG_v1(2e11, Ball.tracers[tracer_type]['p_max'], 
                                                   Ball.tracers[tracer_type]['Q'], 
                                                   Ball.tracers[tracer_type]['logM_cut'], 
                                                   Ball.tracers[tracer_type]['sigma'], 
                                                   Ball.tracers[tracer_type]['gamma']) > 0.01:
#            print("ELG N(2e11) > 0.01")
            return -np.inf
        elif tracer_type == 'QSO' and 10**Ball.tracers[tracer_type]['logM_cut']*Ball.tracers[tracer_type]['kappa'] < 2e11:
#            print("QSO Mmn < 2e11")
            return -np.inf

    # we need to determine the expected number density 
    for tracer_type in mytracers:
        Ball.tracers[tracer_type]['ic'] = 1

    ngal_dict, fsat_dict = Ball.compute_ngal(Nthread = nthread)

    for tracer_type in mytracers:
        if fsat_dict[tracer_type] > 0.6:
            return -np.inf
        if not tracer_type == 'ELG':
            N_tracer = ngal_dict[tracer_type]
            Ball.tracers[tracer_type]['ic'] = \
                min(1, Data.num_dens_mean[tracer_type]*Ball.params['Lbox']**3/N_tracer)
        if tracer_type == 'ELG':
            N_tracer = ngal_dict[tracer_type]
            if abs(Data.num_dens_mean[tracer_type]*Ball.params['Lbox']**3/N_tracer - 1) > 1:
                return -np.inf

    mock_dict = Ball.run_hod(Ball.tracers, Ball.want_rsd, Nthread = nthread, verbose = False)

    # put a satellite fraction cut
    theory_density = {}
    for tracer_type in mytracers:
        if mock_dict[tracer_type]['Ncent'] < 0.4*len(mock_dict[tracer_type]['x']):
            print(tracer_type, 'fsat > 0.6')
            return -np.inf
        theory_density[tracer_type] = len(mock_dict[tracer_type]['x'])/Ball.params['Lbox']**3

    clustering = Ball.compute_wp(mock_dict, Ball.rpbins, Ball.pimax, Ball.pi_bin_size, Nthread = nthread)

    lnP = Data.compute_likelihood(clustering, theory_density, mockcov = mockcov, jointcov_inv = jointcov_inv, fullscale = fullscale)

    return lnP  
            
def prepPrior(bounds):
    prior = Prior()
    for idx, bound_values in bounds:
        if len(bound_values) == 2:  # Flat prior: [lower, upper]
            lower, upper = bound_values
            prior.add_parameter(dist='flat', lower=lower, upper=upper)
        elif len(bound_values) == 4:  # Truncated normal: [lower, upper, mean, std]
            lower, upper, mean, std = bound_values
            a, b = (lower - mean) / std, (upper - mean) / std
            prior.add_parameter(dist=truncnorm(a, b, loc=mean, scale=std))
    return prior


def main(path2config):

    config = yaml.safe_load(open(path2config))
    sim_params = config['sim_params']
    HOD_params = config['HOD_params']
    clustering_params = config['clustering_params']   
    data_params = config['data_params']
    fit_config_params = config['fit_config_params']
    fit_params = config['fit_params']  

    mockcov = fit_config_params.get('mockcov', False)
    if mockcov:
        fit_config_params['chainsPrefix'] += '_mockcov'
    jointcov = fit_config_params.get('joint', False)
    jointcov_inv = None
    if jointcov:
#        fit_config_params['chainsPrefix'] += '_jointmockcov'
        jointcov_inv = np.linalg.inv(np.loadtxt(jointcov))
        
    fullscale = fit_config_params.get('fullscale', True)
    # read data parameters
    if fullscale:
        nrpmin = 3
    else:
        nrpmin = 7
    newData = wp_Data(data_params, HOD_params, nrpmin = nrpmin)
    
    prefix_check = fit_config_params['path2output']+fit_config_params['chainsPrefix']+'.hdf5'

    nparams = sum(len(v) for v in fit_params.values())
    param_mapping = {}
    bounds = []
    mytracers = []

    for tracertype, params in fit_params.items():
        mytracers.append(tracertype)
        param_mapping[tracertype] = {}
        for param, values in params.items():
            mapping_idx = values[0]
            param_mapping[tracertype][param] = mapping_idx
            # Store the entire value array, with the type inferred by the length
            bounds.append((mapping_idx, values[1:]))
    bounds.sort()
    prior=prepPrior(bounds)
    print('start loading AbacusHOD object')
    
    nthread=64
    
    sampler = Sampler(prior=prior, 
                      likelihood=loglike, 
                      n_live=2000,
                      filepath=prefix_check,
                      pass_dict=False,
                      pool=MPIPoolExecutor(),
                      resume=True,
                      likelihood_kwargs={'sim_params':sim_params,
                                         'HOD_params':HOD_params,
                                         'clustering_params':clustering_params,
                                         'param_mapping':param_mapping, 
                                         'mytracers':mytracers, 
                                         'Data':newData, 
                                         'nthread':nthread, 
                                         'bounds':bounds, 
                                         'mockcov':mockcov, 
                                         'jointcov_inv':jointcov_inv, 
                                         'fullscale':fullscale})
    print('run nested')
    sampler.run(verbose=True)
    
class ArgParseFormatter(argparse.RawDescriptionHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass
    
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=ArgParseFormatter)
    parser.add_argument('--path2config', dest='path2config', type=str, help='Path to config file.')
    args = vars(parser.parse_args())    
    main(**args)

    
    
    

    
    
    
    
    
    
