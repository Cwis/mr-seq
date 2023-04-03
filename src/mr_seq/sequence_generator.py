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
        return self.available_sequences.SEQUENCE_TYPES

    def generate(
        self,
        sequence_type: str,
        matrix_shape: tuple,
        title: str,
        outdir: str,
        physics_configuration_filename: str,
        timestep: int = -1,
        verbal: bool = False,
    ):
        # Check seq type available
        if sequence_type not in self.available_sequences.SEQUENCE_TYPES:
            print(
                f"ERROR: asked sequence type ('{sequence_type}') is not available, available sequences are: {self.available_sequences.SEQUENCE_TYPES}"
            )
            return False

        # Get the seq type generator
        sequence_generator = self.available_sequences.SEQUENCE_GENERATORS[
            self.available_sequences.SEQUENCE_TYPES.index(sequence_type)
        ]
        print(sequence_generator)

        # Get the seq type class
        SequenceClass = self.available_sequences.SEQUENCE_CLASSES[
            self.available_sequences.SEQUENCE_TYPES.index(sequence_type)
        ]
        sc = SequenceClass()
        print(sc.check())

        # Generate the sequence using the parameters and the loaded SequenceClass
        sc.generate(
            matrix_shape,
            title,
            outdir,
            physics_configuration_filename,
            timestep,
            verbal,
        )

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
