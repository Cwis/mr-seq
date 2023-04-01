from pathlib import Path
import yaml

from mr_seq.configuration import Configuration


class Sequences(Configuration):
    __default_sequences_filename = "configurations/sequences.yml"

    SEQUENCE_TYPES = None
    SEQUENCE_GENERATORS = None

    def check_configuration(self) -> bool:
        print("Check seq config")
        return True

    def __init__(self, configuration_filename: str = ""):
        # read the configuration file containing the available sequences

        self.configuration_filename: str = configuration_filename
        if configuration_filename == "":
            self.configuration_filename = self.__default_sequences_filename
        self.read_configuration()

        '''
        seq_configuration_path = Path(configuration_filename)
        if seq_configuration_path.is_file():
            with open(seq_configuration_path, "r") as seq_configuration_file:
                try:
                    seq_configuration = yaml.safe_load(seq_configuration_file)
                    # TODO read sequences from configuration file
                except yaml.YAMLError as exc:
                    raise Exception(
                        "Error reading the sequences configuration file: "
                        + f"{seq_configuration_file}"
                    ) from exc
        else:
            raise Exception(f"Error reading the sequence configuration file: invalid path ({seq_configuration_path})")
        '''