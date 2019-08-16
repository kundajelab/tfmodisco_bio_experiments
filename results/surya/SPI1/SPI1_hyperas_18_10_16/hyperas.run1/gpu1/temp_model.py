#coding=utf-8

try:
    import numpy as np
except:
    pass

try:
    import pandas as pd
except:
    pass

try:
    import pysam
except:
    pass

try:
    import pdb
except:
    pass

try:
    import hyperas
except:
    pass

try:
    from hyperas.distributions import uniform, choice
except:
    pass

try:
    from hyperopt import Trials, STATUS_OK, STATUS_FAIL, tpe, mongoexp
except:
    pass

try:
    from hyperas import optim
except:
    pass

try:
    from hyperas.distributions import choice, uniform
except:
    pass

try:
    import keras
except:
    pass

try:
    from keras.models import Sequential
except:
    pass

try:
    from keras.layers.core import Dropout, Dense, Activation, Flatten
except:
    pass

try:
    from keras.layers.convolutional import Conv1D, MaxPooling1D
except:
    pass

try:
    from keras.optimizers import Adadelta, SGD, RMSprop, Adam
except:
    pass

try:
    import keras.losses
except:
    pass

try:
    from keras.layers.normalization import BatchNormalization
except:
    pass

try:
    from keras import backend as K
except:
    pass

try:
    import tabulate
except:
    pass

try:
    import sys
except:
    pass

try:
    import keras_genomics
except:
    pass

try:
    import time
except:
    pass

try:
    import time
except:
    pass

try:
    import tabulate
except:
    pass

try:
    import pdb
except:
    pass
from hyperopt import fmin, tpe, hp, STATUS_OK, Trials

np.random.seed(1234)  # for reproducibility
data=np.load('../data_big.npz')
x_train=data['arr_0']
y_train=data['arr_1']
x_validate=data['arr_2']
y_validate=data['arr_3']
x_test=data['arr_4']
y_test=data['arr_5']


def keras_fmin_fnct(space):

    start = int(time.time())
    np.random.seed(1234)
    parameters = space
    try:
        '''
        b) with and without a hidden fully-connected layer, 
        c) number of units in the hidden fc layer < 200, 
        c) different learning rates for adam (explore on a log scale - 0.001, 0.0001, etc), 
        d) maxpooling widths in the 10-60 range, 
        e) conv widths in the 10-40 range.
        '''
        model=Sequential()
        kernel_size1 = space['kernel_size1']
        kernel_size2 = kernel_size1 - 2
        model.add(Conv1D(filters=space['filters'],kernel_size=(kernel_size1),input_shape=(1000,4)))
        model.add(BatchNormalization(axis=-1))
        model.add(Activation('relu'))

        ## a) number of layers between 1 and 4, 

        #decide on how many conv layers in model 
        n_conv = space['n_conv']

        filter_dim=[kernel_size1,kernel_size2,kernel_size2]

        for i in range(n_conv):
            model.add(Conv1D(filters=space['filters_1'],kernel_size=(filter_dim[i])))
            model.add(BatchNormalization(axis=-1))
            model.add(Activation('relu'))

        model.add(MaxPooling1D(pool_size=(space['filters_2'])))

        model.add(Flatten())
        n_dense = space['n_conv_1']
        for i in range(n_dense):
            model.add(Dense(space['Dense']))
            model.add(BatchNormalization(axis=-1))
            model.add(Activation('relu'))
            model.add(Dropout(space['Dropout']))

        model.add(Dense(4))
        model.add(Activation("sigmoid"))

        adam=keras.optimizers.Adam(lr=space['lr'])

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
        status=STATUS_OK

    except:
        loss=1000
        status=STATUS_FAIL
        acc=0
        print("failed to run model")
        sys.stdout.flush()

        model=None

    return{'loss':loss,'status':status,'model_param':parameters}

def get_space():
    return {
        'kernel_size1': hp.choice('kernel_size1', [13, 15, 19, 25, 39]),
        'filters': hp.choice('filters', [50,100,250,500]),
        'n_conv': hp.choice('n_conv', [0,1,2,3]),
        'filters_1': hp.choice('filters_1', [50,100,250,500]),
        'filters_2': hp.choice('filters_2', [20,40,60]),
        'n_conv_1': hp.choice('n_conv_1', [0,1,2,3]),
        'Dense': hp.choice('Dense', [50,100,200]),
        'Dropout': hp.choice('Dropout', [0.2,0.4,0.6]),
        'lr': hp.choice('lr', [0.001, 0.0001]),
    }
