#!/usr/bin/env python
from Cef_tools import CEF_obj
import cPickle as pickle
import pystan
import numpy as np
from itertools import izip
from multiprocessing import Pool
import getopt
import sys, os

if __name__ == '__main__':
    n_processes = 2
    clustering_attribute = 'Cell_type'
    outfiles_path = ''

    optlist, args = getopt.gnu_getopt(sys.argv[1:], "hi:o:c:p:", ["help", "input=","output=","clus_attr=","processes="])

    if optlist== [] and args == []:
        print 'pystancef -i [INPUTFILE] -o [OUTPUTFOLDER] -c [CLUSTER_ATTRIBUTE] -p [THREADS]'
        sys.exit()
    for opt, a in optlist:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt in ('-i', '--input'):
            input_path = a
        elif opt in ("-o", "--output"):
            outfiles_path = a
        elif opt in ('-c','--clus_attr'):
            clustering_attribute = str(a)
        elif opt in ('-p', '--processes'):
            n_processes = int(a)
        else:
            assert False, "%s option is not supported" % opt

    cef = CEF_obj()

    try:
        print 'Loading CEF'
        cef.readCEF(input_path)
    except:
        print 'Error in loading CEF file'

    try:
        print 'Loading Model'
        sm = pickle.load(open(os.path.join(outfiles_path,'pystan_model_compiled.pickle'),'rb'))
    except IOError:  
        print 'Compiled model was not found. Defining and compiling model.'
        sten_code = '''
        # Bayesian Negative Binomial regression for single-cell RNA-seq

        data {
            int<lower=0> N;                 # number of outcomes
            int<lower=0> K;                 # number of predictors
            matrix<lower=0>[N,K] x;         # predictor matrix 
            int y[N];                       # outcomes
        }

        parameters {
            vector<lower=1>[K] beta;  # coefficients
            real<lower=0.001> r;  # overdispersion
        }

        model {	
            vector[N] mu;
            vector[N] rv;

            # priors
            r ~ cauchy(0, 1);
            beta ~ pareto(1, 1.5);

            # vectorize the scale
            for (n in 1:N) {
                rv[n] <- square(r + 1) - 1;
            }

            # regression
            mu <- x * (beta - 1) + 0.001;
            y ~ neg_binomial(mu ./ rv, 1 / rv[1]);
        }
        '''

        sm = pystan.StanModel(model_code=sten_code)

        print 'Saving model for future use.'
        pickle.dump(sm, open(os.path.join(outfiles_path,'pystan_model_compiled.pickle'),'wb'))




    print 'Formating the model input'
    for i,v in izip(cef.col_attr_names, cef.col_attr_values):
        if i == clustering_attribute:
            predictor_list = v
        if 'total' in i.lower():
            total_molecules = [float(j) for j in v ]
            
    for i,v in izip(cef.row_attr_names, cef.row_attr_values):
        if 'gene' in i.lower():
            gene_names = v

    total = sum(total_molecules)/len(total_molecules)
    total_molecules_norm = [j/total for j in total_molecules]
    predictors = ['Size']
    for i in predictor_list:
        if i not in predictors:
            predictors.append(i)

    predictor_matrix = []
    for i, c_p in enumerate( predictor_list ):
        predictor_matrix.append( [total_molecules_norm[i]] + [float(c_p==p) for p in predictors[1:]] )


    N = len( predictor_matrix )
    K = len( predictors )

    def one_gene_model(name_gene, gene_vector, predictor_matrix, N, K):
        path = os.path.join( outfiles_path, "beta_%s.npy" % name_gene )
        my_data = {'N': N, 'K': K, 'x': predictor_matrix, 'y': gene_vector}
        
        fit = sm.sampling(data = my_data, iter=3000, chains=1, seed='19900715',warmup=2000, n_jobs=1)

        # extract the traces
        traces = fit.extract(permuted=True)

        # save the traces with numpy
        
        np.save(path, traces['beta'])
        print path

    print 'Passing the inputs to mutliple threads'
    
    try:
        p = Pool(processes=n_processes)
        for name_gene, gene_vector in izip(gene_names, cef.matrix):
            p.apply_async(one_gene_model,(name_gene, gene_vector, predictor_matrix, N, K))
        p.close()
        p.join()
    except KeyboardInterrupt:
        p.terminate()
        

