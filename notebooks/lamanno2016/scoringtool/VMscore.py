#!/usr/bin/env python

# Copyright (c) 2015 Gioele La Manno and Sten Linnarsson
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import sys
import cPickle as pickle
import getopt
import numpy as np
import pandas as pd

optlist, args = getopt.gnu_getopt(sys.argv[1:], "hi:o:", ["help", "input=","output="])
if (optlist== [] and args == []):
    print '''
    usage:
    VMscore.py -i inputpath -o outputpath
    Input file containing single cell transcriptomes (UMI)
    formatted as a tab-delimited text file with 
    a row header (official gene symbol) and 
    a column header (unique_cell_id)
    '''
    sys.exit()
for opt, a in optlist:
    if opt in ("-h", "--help"):
        print '''
        usage:
        VMscore.py -i inputpath -o outputpath
        Input file containing single cell transcriptomes (UMI)
        formatted as a tab-delimited text file with 
        a row header (official gene symbol) and 
        a column header (unique_cell_id)
        '''
        sys.exit()
    elif opt in ('-i', '--input'):
        input_path = a
    elif opt in ("-o", "--output"):
        output_path = a

try:
	input_path == output_path
except NameError:
	 sys.exit()

dict_model = pickle.load(open('VMmodel.pickle'))
normalizer = dict_model['normalizer']
genenames = dict_model['genenames']
classesnames = dict_model['classesnames']
Logistic = dict_model['Logistic']

df_in = pd.read_csv(input_path, sep='\t', index_col=0)
df_in = np.log2( df_in.ix[genenames,:].fillna(0) +1)

prob = Logistic.predict_proba((df_in.values/normalizer).T)
df_out = pd.DataFrame(prob.T,index=classesnames,columns=df_in.columns)

df_out.to_csv(output_path,sep='\t')
