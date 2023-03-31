from mr_seq.sequences import gre_sequence


class Sequences:
    SEQUENCE_TYPES = ["gre", "spiral", "border"]

    @classmethod
    def generate_sequence(
        self,
        sequence_type: str,
        matrix_shape: tuple,
        title: str,
        outdir: str,
        physics_configuration_filename: str,
        timestep: int = -1,
        verbal: bool = False,
    ):
        if not sequence_type in self.SEQUENCE_TYPES:
            raise Exception(
                f"Unknown sequence type '{sequence_type}',"
                + "sequence types allowed are: {SEQUENCE_TYPES}"
            )
        if sequence_type == "gre":
            return gre_sequence.generate_gre(
                matrix_shape,
                title,
                outdir,
                physics_configuration_filename,
                timestep,
                verbal,
            )
        elif sequence_type == "spiral":
            return False
            # return self.generate_spiral_sequence(
            # matrix_shape, title, outdir, physics_configuration_filename, timestep
            # )
        elif sequence_type == "border":
            return False
            # return self.generate_border_sequence(
            # matrix_shape, title, outdir, physics_configuration_filename, timestep
            # )
