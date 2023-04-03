# mr-seq
# c.chenes
#
# demo of the mr-seq abilities to generate an mr-seq using extended pypulseq for mr-sim
import os
from mr_seq.sequences import sequences
from mr_seq.sequence_generator import SequenceGenerator

# Init a sequence generator
seq_gen = SequenceGenerator()

print(f"Demo of the mr-seq module (v{SequenceGenerator.version()})")
print(f"Currently available sequences are: {seq_gen.get_sequences()}")

# 1. Generate a GRE seq with parameters from the gre_demo_configuration file

GRE_matrix_shape = (1, 7, 7)
title = "GRE_demo"
outdir = "demo"
physics_parameters = "demo/gre_demo_configuration.yml"
timestep = 1e-6
verbal = True

gre_res = seq_gen.generate(
    "gre", GRE_matrix_shape, title, outdir, physics_parameters, timestep, verbal
)

# 2. Generate a "border" seq with parameters from the border_demo_configuration file

matrix_shape = (1, 7, 7)
title = "border_demo"
outdir = "border_demo"
physics_parameters = "border_demo_configuration.yml"
timestep = 1e-6

border_res = seq_gen.generate(
    "border", matrix_shape, title, outdir, physics_parameters, timestep, verbal
)
