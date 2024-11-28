# INQuant Quantification algorithm

`INQuant` is a Python class designed for quantification for InstaNovo predictions. 
Developed by Annekatrine Kirketerp-Møller and Ida Sofie Goldschmidt as part of a 5 ECTS special course at DTU.

## Installation

Ensure you have the required dependencies installed:

```bash
pip install pyopenms pandas numpy
```

## Usage

### Initialization

Initialize the `INQuant` class with a predictions file and a list of mzML file(s):

```python
from INQuant import INQuant

predictions_file = 'path/to/predictions.csv'
mzML_file_list = ['path/to/file1.mzML', 'path/to/file2.mzML']

with INQuant(predictions_file, mzML_file_list) as Exp:
    Exp.run()
```

### Methods

- **`__init__(self, predictions_file, mzML_file_list, script_purpose=None)`**: Initializes the class with the necessary files and optional script purpose.
- **`__enter__(self)`**: Starts the timer for the context manager.
- **`__exit__(self, *args)`**: Writes the status file upon completion.
- **`time_flag(self, name)`**: Sets a time flag to record execution time.
- **`run(self, write_individual_area_files=False, normalize_abundance=True, output_format='Wide', overwrite_quant_file=None)`**: Main method to execute the quantification process.
- **`load_preds(self)`**: Loads prediction data and matches it with mzML files. Sorts predictions for matching preds/targets
- **`load_mzML(self, mzML_id, make_area_files=False)`**: Loads mzML file data and calculates noise.
- **`get_areas(self, make_file=False)`**: Calculates the area under the curve for precursors in the sample.
- **`write_quant_file(self, file=None, output_format='Long', normalize=False, overwrite=None)`**: Writes the quantification results to a CSV file with desired settings.

### Example

```python
from INQuant import INQuant

predictions_file = 'path/to/predictions.csv'
mzML_file_list = ['path/to/file1.mzML', 'path/to/file2.mzML']

with INQuant(predictions_file, mzML_file_list, script_purpose='Quantification of peptides') as Exp:
    Exp.run(write_individual_area_files=True, normalize_abundance=True, output_format='Wide')
```
