#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys, os
import numpy as np
import hood
from hoomd import azplugins
import gsd

from hoomd_util import *

if __name__=='__main__':

    time_start = time.time()

    input_file = sys.argv[1]
