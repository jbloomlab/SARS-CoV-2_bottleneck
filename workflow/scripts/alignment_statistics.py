"""
This module contains utilities for assessing short-read 
alignments. Most utilities are wrappers around various 
`pysam` functions. 
"""
__author__ = "Will Hannon"
__copyright__ = "Copyright 2021 Will Hannon"
__email__ = "wwh22@uw.edu"
__license__ = "MIT"

import pysam
import pandas as pd #data frames
import numpy as np #arrays
import os #interacting with files
from Bio import SeqIO #reading fasta format
import re #regular expressions
