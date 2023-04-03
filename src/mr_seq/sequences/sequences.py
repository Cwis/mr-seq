from pathlib import Path

from mr_seq.configuration import Configuration


def importName(modulename, name):
    """Import a named object from a module in the context of this function."""
    try:
        module = __import__(modulename, globals(), locals(), [name])
    except ImportError:
        return None
    return vars(module)[name]


class Sequences(Configuration):
    __default_sequences_filename = "configurations/sequences.yml"

    SEQUENCE_TYPES = []
    SEQUENCE_GENERATORS = []
    SEQUENCE_FILES = []
    SEQUENCE_CLASSES = []
    SEQUENCE_PREFIX = "generate_"
    SEQUENCE_SUFFIX = "_sequence"
    SEQUENCE_CLASS = "Sequence_"

    def check_configuration(self) -> bool:
        # Check all required keys are in the file
        if any([key not in self.configuration for key in self.configuration_keys]):
            print(
                f"Configuration check fail: sequences config must contain all keys: '{self.configuration_keys}'"
            )
            return False

        # Check sequences directory exists
        seq_dir = Path(self.configuration["sequences_directory"])
        if not (seq_dir.is_dir()):
            print(f"ERROR: directory of sequences does not exist: {seq_dir}")
            return False

        # Get the names of all sequences and check that the corresponding generator file exists
        for seq in self.configuration["sequences"]:
            seq_name = seq["name"]
            seq_generator = f"{self.SEQUENCE_PREFIX}{seq['name']}"
            seq_module = f"{seq['name']}{self.SEQUENCE_SUFFIX}"
            seq_filename = f"{seq_module}.py"
            seq_class = f"{self.SEQUENCE_CLASS}{seq_name}"
            seq_generator_path = seq_dir / seq_filename
            if not seq_generator_path.is_file():
                print(
                    f"Configuration check fail: sequence '{seq_name}' generator ('{seq_generator_path}') not found."
                )
                continue  # Skip this sequences => this allows for incomplete seq config

            # Add the sequence and its generator to the available sequences
            self.SEQUENCE_TYPES.append(seq_name)
            self.SEQUENCE_GENERATORS.append(seq_generator)
            self.SEQUENCE_FILES.append(seq_filename)

            # Import the generator function from the sequence file
            module_name = "mr_seq.sequences." + seq_filename
            # importName(f"mr_seq.sequences.{seq_module}", seq_generator)
            SequenceClass = importName(f"mr_seq.sequences.{seq_module}", seq_class)
            self.SEQUENCE_CLASSES.append(SequenceClass)

            # Check the import worked
            sc = SequenceClass()
            if not sc.check():
                print(f"ERROR: cannot check the SequenceClass: {seq_class}")
                return False

        return True

    def __init__(self, configuration_filename: str = ""):
        # init keys
        self.configuration_keys = [
            "title",
            "sequences_directory",
            "sequences",
        ]

        # read the configuration file containing the available sequences
        self.configuration_filename: str = configuration_filename
        if configuration_filename == "":
            self.configuration_filename = self.__default_sequences_filename
        self.read_configuration()
