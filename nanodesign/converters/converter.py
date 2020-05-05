# Copyright 2016 Autodesk Inc.
# Modifications Copyright (C) 2019 Dietzlab (TUM), Elija Feigl
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

""" This module is used to convert DNA design files into other file formats.

    An input caDNAno design file is conveted into a DnaStructure object
    containing information about the design (e.g. lattice type, virtual helix
    definitions) and information derived from that design (e.g. strands,
    domains). caDNAno design files may contain deleted/inserted bases. By
    default the DnaStructure is not created with deleted/inserted bases. The
    DnaStructure is created with deleted/inserted bases by specifying the
    --modify command-line argument.
"""
import os
import re

import logging
from .cadnano.reader import CadnanoReader
from .cadnano.writer import CadnanoWriter
from .cadnano.convert_design import CadnanoConvertDesign
from .cando.writer import CandoWriter

from ..data.parameters import DnaParameters
from .dna_sequence_data import dna_sequence_data
# TODO (JMS, 10/26/16): revisit where the sequence data is kept?


class ConverterFileFormats(object):
    """ File format names to convert to/from. """
    UNKNOWN = "unknown"
    CADNANO = "cadnano"
    CANDO = "cando"
    CIF = "cif"
    PDB = "pdb"
    SIMDNA = "simdna"
    STRUCTURE = "structure"
    TOPOLOGY = "topology"
    VIEWER = "viewer"
    names = [CADNANO, CANDO, CIF, PDB, SIMDNA, STRUCTURE, TOPOLOGY, VIEWER]


class Converter(object):
    """ This class stores objects for various models created when reading from
        a file.

        Attributes:
            cadnano_design (CadnanoDesign): The object storing the caDNAno
                design information.
            cadnano_convert_design (CadnanoConvertDesign): The object used to
                convert a caDNAno design into a DnaStructure.
            dna_parameters (DnaParameters): The DNA physical parameters used to
                generate the geometry of a DNA structure
            dna_structure (DnaStructure): The object storing connectivity and
                geometry of a DNA structure.
            infile (String): The file name to convert.
            informat (String): The format of the file to convert, taken from
                ConverterFileFormats.
            modify (bool): If true then DnaStructure is created with
                deleted/inserted bases.
            logg (bool): disable logging
            outfile (String): The name of the file for converter output.
    """
    def __init__(self, modify=False, logg=True):
        self.cadnano_design = None
        self.dna_structure = None
        self.cadnano_convert_design = None
        self.infile = None
        self.informat = None
        self.outfile = None
        self.modify = modify
        self.dna_parameters = DnaParameters()
        self.logger = logging.getLogger(__name__)

        if not logg:
            logging.disable(logging.INFO)

    def read_cadnano_file(self, file_name, seq_file_name, seq_name):
        """ Read in a caDNAno file.

            Arguments:
                file_name (String): The name of the caDNAno file to convert.
                seq_file_name (String): The name of the CSV file used to assign
                    a DNA base sequence to the DNA structure.
                seq_name (String): The name of a sequence used to assign a DNA
                    base sequence to the DNA structure.
        """
        cadnano_reader = CadnanoReader()
        self.cadnano_design = cadnano_reader.read_json(file_name)
        self.cadnano_convert_design = CadnanoConvertDesign(self.dna_parameters)
        self.dna_structure = self.cadnano_convert_design.create_structure(
            self.cadnano_design, self.modify
        )

        # Read in staple sequences from a CSV format file.
        if (seq_file_name):
            _, file_extension = os.path.splitext(seq_file_name)

            if (file_extension == ".csv"):
                modified_structure = False
                sequence = cadnano_reader.read_csv(seq_file_name)
                self.cadnano_convert_design.set_sequence(
                    self.dna_structure, modified_structure, sequence
                )

        # Assign a sequence using a name.
        if (seq_name):
            if (seq_name not in dna_sequence_data):
                self.logger.error(
                    "The sequence name %s is not recognized.", seq_name
                )
            modified_structure = False
            self.cadnano_convert_design.set_sequence_from_name(
                self.dna_structure, modified_structure, seq_name
            )

    def write_topology_file(self, file_name):
        """ Write a DNA topology file.

            Arguments:
                file_name (String): The name of the topology file to write.
        """
        self.dna_structure.write_topology(file_name, write_json_format=True)

    def write_structure_file(self, file_name):
        """ Write a DNA structure file.
            Arguments:
                file_name (String): The name of the structure file to write.
        """
        self.dna_structure.write(file_name, write_json_format=True)

    def write_cando_file(self, file_name):
        """ Write a CanDo .cndo file.
            Arguments:
                file_name (String): The name of the CanDo file to write.
        """
        cando_writer = CandoWriter(self.dna_structure)
        cando_writer.write(file_name)

    def write_cadnano_file(self, file_name):
        """ Write a caDNAno JSON file.

            Arguments:
                file_name (String): The name of the caDNAno file to write.
        """
        cadnano_writer = CadnanoWriter(self.dna_structure)
        cadnano_writer.write(file_name)

    def perform_staple_operations(self, staples_arg):
        """ Perform operations on staples.

            Arguments:
                staples_arg (String): The argument to the staples command-line
                    option.
        """
        tokens = staples_arg.split(",", 1)
        operation = tokens[0]
        retain_staples = []

        # Parse retained staples IDs.
        if len(tokens) == 2:
            pattern = re.compile('\W')
            retain_tokens = pattern.split(tokens[1])
            if retain_tokens[0] == "retain":
                retain_colors = [int(color)
                                 for color in retain_tokens[1:]
                                 if color != ''
                                 ]
            retain_staples = self.dna_structure.get_staples_by_color(
                retain_colors
            )

        # Remove all staple strands except those given in retain_staples[].
        if operation == "delete":
            self.dna_structure.remove_staples(retain_staples)

        # Generaqte the maximal staple strand set except those given in
        # retain_staples[].
        elif operation == "maximal_set":
            self.dna_structure.generate_maximal_staple_set(retain_staples)

    def set_module_loggers(self, names):
        module_names = names.split(",")
        for module in module_names:
            logger = logging.getLogger(module)
            logger.setLevel(logging.DEBUG)
