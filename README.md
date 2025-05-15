# INQuant Quantification algorithm

`inquant_tools` is a Python module designed for quantification for InstaNovo predictions. 
Developed by Annekatrine Kirketerp-Møller and Ida Sofie Goldschmidt as part of a 5 ECTS special course and 20 ECTS Bachelors thesis at DTU.

## Documentation
[The documentation for the code, including the API documentation can be found here](https://uadklab.github.io/INQuant/)

## Installation

Ensure you have the required dependencies installed:

```bash
pip install -r requirements.txt
```

Install with git clone:

```bash
git clone https://github.com/UadKLab/INQuant.git
```

## Usage

Simple example of how to initialize the and run `INQuant` with minimal user input:

```python
from inquant_tools import INQuant  

with INQuant(predictions_file="predictions.csv",
                mzml_file_list = ["file_1.mzml", "file_2.mzml", "file_3.mzml"],
                proteome_file="example_proteome.fasta"
            ) as INQ:          
    INQ.run()
```
