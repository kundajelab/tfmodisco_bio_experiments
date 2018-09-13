import numpy as np
import pandas as pd
import pysam
import pdb

import hyperas
from hyperas.distributions import uniform, choice
from hyperopt import Trials, STATUS_OK, tpe
from hyperas import optim
from hyperas.distributions import choice, uniform

import keras;
from keras.models import Sequential
from keras.layers.core import Dropout, Dense, Activation, Flatten
from keras.layers.convolutional import Conv1D, MaxPooling1D
from keras.optimizers import Adadelta, SGD, RMSprop, Adam;
import keras.losses;
from keras.layers.normalization import BatchNormalization
from keras import backend as K
import tabulate
import sys
import keras_genomics
import time

K.set_image_data_format('channels_last')
print(K.image_data_format())



#define training, validation, test data
def data():
    np.random.seed(91011)  # for reproducibility
    data=np.load('../data.npz')
    x_train=data['arr_0']
    y_train=data['arr_1']
    x_validate=data['arr_2']
    y_validate=data['arr_3']
    x_test=data['arr_4']
    y_test=data['arr_5']
    return x_train,y_train,x_validate,y_validate,x_test,y_test

def data_from_bed():

    #def load_bed(fname,ref_fasta="/srv/scratch/annashch/hg19.genome.fa"):
    def load_bed(fname,ref_fasta="/srv/scratch/ktian/kundajelab/ENCODE/genome/hg19.fa"):
        ltrdict = {'a':[1,0,0,0],
                   'c':[0,1,0,0],
                   'g':[0,0,1,0],
                   't':[0,0,0,1],
                   'n':[0,0,0,0],
                   'A':[1,0,0,0],
                   'C':[0,1,0,0],
                   'G':[0,0,1,0],
                   'T':[0,0,0,1],
                   'N':[0,0,0,0]}
        ref=pysam.FastaFile(ref_fasta)
        data=pd.read_csv(fname,header=None,sep='\t',index_col=[0,1,2], compression='gzip')
        seqs=[ref.fetch(i[0],i[1],i[2]) for i in data.index]
        seqs=np.array([[ltrdict.get(x,[0,0,0,0]) for x in seq] for seq in seqs])
        #x=np.expand_dims(seqs,1)
        x=seqs
        y=np.asarray(data)
        print("x.shape=", x.shape, "y.shape=", y.shape)
        return x,y

    x_train,y_train=load_bed("../hyper_splits/train_small.tsv.gz")
    x_validate,y_validate=load_bed("../hyper_splits/valid_small.tsv.gz")
    x_test,y_test=load_bed("../hyper_splits/test_small.tsv.gz")

    #x_train,y_train=load_bed("../hyper_splits/train.tsv.gz")
    #x_validate,y_validate=load_bed("../hyper_splits/valid.tsv.gz")
    #x_test,y_test=load_bed("../hyper_splits/test.tsv.gz")
    np.savez('../data.npz',x_train,y_train,x_validate,y_validate,x_test,y_test)
    return x_train,y_train,x_validate,y_validate,x_test,y_test

def create_model(x_train,y_train,x_validate,y_validate,x_test,y_test):
    import time
    import tabulate
    start = int(time.time())
    np.random.seed(91011)
    try:
        '''
        b) with and without a hidden fully-connected layer, 
        c) number of units in the hidden fc layer < 200, 
        c) different learning rates for adam (explore on a log scale - 0.001, 0.0001, etc), 
        d) maxpooling widths in the 10-60 range, 
        e) conv widths in the 10-40 range.
        '''
        model=Sequential()
        kernel_size1 = {{choice([13, 15, 19, 25, 39])}}
        kernel_size2 = kernel_size1 - 2
        model.add(Conv1D(filters={{choice([50,100,250,500,1000])}},kernel_size=(kernel_size1),input_shape=(1000,4)))
        model.add(BatchNormalization(axis=-1))
        model.add(Activation('relu'))

        ## a) number of layers between 1 and 4, 

        #decide on how many conv layers in model 
        n_conv = {{choice([0,1,2,3])}}

        filter_dim=[kernel_size1,kernel_size2,kernel_size2]

        for i in range(n_conv):
            model.add(Conv1D(filters={{choice([50,100,250,500,1000])}},kernel_size=(filter_dim[i])))
            model.add(BatchNormalization(axis=-1))
            model.add(Activation('relu'))

        model.add(MaxPooling1D(pool_size=({{choice([10,30,60])}})))

        model.add(Flatten())
        n_dense = {{choice([0,1,2,3])}}
        for i in range(n_dense):
            model.add(Dense({{choice([50,100,200])}}))
            model.add(BatchNormalization(axis=-1))
            model.add(Activation('relu'))
            model.add(Dropout({{choice([0.2,0.4,0.6])}}))

        model.add(Dense(4))
        model.add(Activation("sigmoid"))

        adam=keras.optimizers.Adam(lr={{choice([0.01, 0.001, 0.0001])}})

        model.compile(loss=keras_genomics.losses.ambig_binary_crossentropy, optimizer=adam, metrics=['accuracy'])
        print("compiled!")
        sys.stdout.flush()


        # added to collect optimization results
        if 'results' not in globals():
            global results
            results = []

        result = model.fit(x_train,y_train,
                           batch_size=200,
                           epochs=20,
                           verbose=2,
                           validation_data=(x_validate,y_validate))
        print("trained!")
        sys.stdout.flush()

        loss,acc = model.evaluate(x_validate,y_validate,verbose=2)
        print("Validation loss:",loss,"Validation acc:",acc)
        sys.stdout.flush()

        # added to collect results
        valLoss = result.history['val_loss'][-1]
        parameters = space
        parameters["loss"] = valLoss
        parameters["time"] = int(time.time() - start)
        results.append(parameters)
        print(parameters)
        if len(results) % 10 == 0 :
            tab = tabulate.tabulate(results, headers="keys", tablefmt="fancy_grid", floatfmt=".8f")
            print(tab.encode('utf-8'))
        else:
            tab = tabulate.tabulate(results[-1:], headers="keys", tablefmt="fancy_grid", floatfmt=".8f")
            print(tab.encode('utf-8'))
        print("model %d done ----------------" % len(results))
        sys.stdout.flush()

    except:
        loss=1000
        acc=0
        print("failed to run model")
        sys.stdout.flush()

        model=None

    return{'loss':loss,'status':STATUS_OK,'model':model}

def write_result(outf, best_eval, best_run, best_model): #,real_param_values):
    outf.write("Evalutation of best performing model:\n")
    outf.write(str(best_eval) + "\n")
    #outf.write("Best performing model chosen hyper-parameters:\n")
    #outf.write(str(real_param_values) + "\n")
    outf.write("Best performing model chosen hyper-parameters index:\n")
    outf.write(str(best_run) + "\n")
    outf.write("Summary of best model:\n")
    best_model.summary(print_fn=lambda x: outf.write(x + '\n'))
 
if __name__=="__main__":
    #x_train,y_train,x_validate,y_validate,x_test,y_test=data_from_bed()
    x_train,y_train,x_validate,y_validate,x_test,y_test=data()
    #model = create_model0(x_train,y_train,x_validate,y_validate,x_test,y_test)

    print("seed=91011")
    best_run,best_model=optim.minimize(model=create_model,
                                       data=data,
                                       algo=tpe.suggest,
                                       max_evals=50,
                                       trials=Trials(),
                                       rseed=91011)
    #from hyperas.utils import eval_hyperopt_space
    #space = get_space()
    #real_param_values = eval_hyperopt_space(space, best_run)

    best_eval = best_model.evaluate(x_test, y_test)

    outf=open('hyperas.results.txt','w')
    write_result(outf, best_eval, best_run, best_model)

    outf = sys.stdout
    write_result(outf, best_eval, best_run, best_model)

    import pdb
    pdb.set_trace()

    #dded to output results
    #print(tabulate.tabulate(results, headers="keys", tablefmt="fancy_grid", floatfmt=".8f"))


