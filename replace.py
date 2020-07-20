import numpy as np
import math
import pandas as pd

def replace(a):
    for i in range(len(a)):
        if a[i] == 0:
            a[i] = 1
        elif a[i] == 1:
            a[i] = 0
    return a