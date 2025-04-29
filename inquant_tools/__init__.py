from pyopenms import *
import pandas as pd
import numpy as np
import time
from Bio import SeqIO
import re
from tqdm.auto import tqdm
from math import floor, ceil, dist, sqrt
import random
from bs4 import BeautifulSoup
import requests

class Break(Exception):
    """
    Custom exception to signal an intentional break in script execution.

    This exception is used to control the flow of the script, allowing for 
    early termination or interruption within the context of a `with` statement.

    Attributes
    ----------
    message : str, optional
        A string that can be passed when the exception is raised to provide
        additional context about the break event.
    
    """
    pass

aa_mass = {'A':71.03711,'C':103.00919,'D':115.02694,'E':129.04259,'F':147.06841,'G':57.02146,'H':137.05891,'I':113.08406,'K':128.09496,'L':113.08406,'M':131.04049,'N':114.04293,'P':97.05276,'Q':128.05858,'R':156.10111,'S':87.03203,'T':101.04768,'V':99.06841,'W':186.07931,'Y':163.06333, 'H2O':18.01528} # Source:
alkylation_mass = 57.02146 # Alkylation with Iodoacetamide og 2-Chloroacetamide, sources: https://analyticalsciencejournals.onlinelibrary.wiley.com/doi/10.1002/(SICI)1096-9888(200004)35:4%3C572::AID-JMS971%3E3.0.CO;2-2 og https://pubs.acs.org/doi/10.1021/acs.jproteome.7b00022
aa_mass['C'] += alkylation_mass
proton_da = 1.007276466812 # the Da for a single proton
neutron_da = 1.00866491578 # the Da for a single neutron
        
from inquant_tools import INQuant
