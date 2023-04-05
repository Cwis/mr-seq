# mr-seq
# c.chenes
#
# demo of the mr-seq abilities to generate an mr-seq using extended pypulseq for mr-sim

import os
from mr_seq.sequences import sequences
from mr_seq.sequence_generator import SequenceGenerator

# Init a sequence generator and retrieve available sequences
seq_gen = SequenceGenerator()
available_sequences = seq_gen.get_sequences()

print(f"Demo of the mr-seq module (v{SequenceGenerator.version()})")
print(f"Currently available sequences are: {available_sequences}")

# 1. Generate a GRE seq with parameters from the gre_demo_configuration file

sequence_type = "gre"
matrix_shape = (1, 7, 7)
title = "GRE_demo"
output_directory = f"out/demo/{title}"
physics_parameters = "demo/gre_demo_configuration.yml"
clock = 20e-6
verbal = True

gre_res = seq_gen.generate(
    sequence_type,
    matrix_shape,
    title,
    output_directory,
    physics_parameters,
    clock,
    verbal,
)

# 2. Generate a "border" seq with parameters from the border_demo_configuration file

sequence_type = "border"
matrix_shape = (1, 7, 7)
title = "border_demo"
output_directory = f"out/demo/{title}"
physics_parameters = "border_demo_configuration.yml"
clock = 1e-6
verbal = True

border_res = seq_gen.generate(
    sequence_type,
    matrix_shape,
    title,
    output_directory,
    physics_parameters,
    clock,
    verbal,
)
