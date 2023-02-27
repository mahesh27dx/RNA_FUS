import pathlib
import argparse
import numpy as np

valid_input_formats = ['.gsd']
valid_output_formats = ['.gsd']

parser = argparse.ArgumentParser(
    description=
    "Initialize the system from input file and perform MD runs")
# parser.add_argument(
#     'ifile',
#     type= str
#     help = f'Input file (allowed formats: {valid_input_formats})'
# )
# parser.add_argument(
#     'ofile',
#     type = str
#     help = f'Output file (allowed formats: {valid_output_formats})'
# )
parser.add_argument('-dt', '--dt', type=float, default=0.001, help='Time step size')
parser.add_argument('-time' ,'--time', type=int, default=10000, help='Simulation run time')
parser.add_argument('-temp' ,'--temp', type=int, default=300, help='temperature for the current simulation run')

args = parser.parse_args()

# ifile = pathlib.Path(args.ifile)
# ofile = pathlin.Path(args.ofile)

stat_file = 'input_files/stats_module.dat'


dt = args.dt
simulation_steps = args.time
T = args.temp
a = np.savetxt(T, "output_files/test" + str(args.temp) + '.txt')
print(dt)
print(simulation_steps)
print(T)
