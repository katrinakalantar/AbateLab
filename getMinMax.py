__author__ = 'KATRINA'

import sys
import numpy as np


def getMinMax(file):
    min_ = []
    max_ = []
    f = open(file,'r').readlines()
    for line in f:
        if 'minimum' in line:
            min_.append(float(line.split(':')[1].strip()))
        elif 'maximum' in line:
            max_.append(float(line.split(':')[1].strip()))

    print(min_)
    print(max_)

    print('mean min: ' + str(np.mean(min_)))
    print('minimum min: ' + str(min(min_)))
    print('mean max: ' + str(np.mean(max_)))
    print('maximum max: ' + str(max(max_)))


file = sys.argv[1]
getMinMax(file)