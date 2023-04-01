"""Class to generate a sequence with extended pypulseq

sequence_type supported: gre

date: 2023-03-31
author: cc
"""
import os
from importlib.metadata import version

from mr_seq.sequences.sequences import Sequences


class SequenceGenerator:
    PRECISION = 5
    MS_CONVERT = 1e3

    available_sequences: Sequences

    def __init__(self):
        # read available sequences from configuration file
        self.available_sequences = Sequences()

    @staticmethod
    def version():
        v = version("mr-seq")
        return v

    def get_sequences(self):
        return self.available_sequences.SEQUENCE_TYPES #return seqs.SEQUENCE_TYPES

    @classmethod
    def generate(
        cls,
        sequence_type: str,
        matrix_shape: tuple,
        title: str,
        outdir: str,
        physics_configuration_filename: str,
        timestep: int = -1,
        verbal: bool = False,
    ):
        return True
        # return seqs.generate_sequence(
        #     sequence_type,
        #     matrix_shape,
        #     title,
        #     outdir,
        #     physics_configuration_filename,
        #     timestep,
        #     verbal,
        # )