#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys, os
import numpy as np
import hoomd
from hoomd import azplugins
import gsd

from hoomd_util import *

hoomd.context.initialize("")

# Read the "starting_config.gsd"
system = hoomd.init.read_gsd('starting_config.gsd')

snapshot = system.take_snapshot()

# Define the rigid bodies
## Types for the rigid bodies
type_rigid_1 = [aa_type[prot_id[i]] for i in range(284, 371)]
type_rigid_2 = [aa_type[prot_id[i]] for i in range(421, 453)]

rigid = hoomd.md.constrain.rigid()

rigid.set_param('R',
                types = type_rigid_1. positions=)
