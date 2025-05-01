==================
Installation
==================

We recommend downloading the package directly from the GitHub repository. You can do this by cloning the repository or downloading it as a ZIP file.
To clone the repository, use the following command:

Git repository
----------------

.. code-block:: bash

    git clone https://github.com/UadKLab/INQuant.git

Dependencies
----------------

This project requires several Python libraries. Below is a list of the main dependencies:

- **pyopenms**: A Python library for working with OpenMS data structures and algorithms.
- **pandas**: Used for data manipulation and analysis, especially for handling structured data like tables and time series.
- **numpy**: A fundamental package for scientific computing, useful for performing operations on arrays and matrices.
- **time**: A standard Python module for time-related functions.
- **Bio (BioPython)**: A suite of tools for biological computation, used here for reading sequence data.
- **re**: A standard library for regular expression operations, which allows searching and manipulating strings.
- **tqdm**: A library for creating progress bars in Python loops.
- **math**: A standard library for mathematical operations such as `floor`, `ceil`, and `dist`.
- **random**: A standard Python library for generating pseudo-random numbers.
- **BeautifulSoup (bs4)**: A library for parsing HTML and XML documents.
- **requests**: A simple HTTP library for making web requests.

For a full list of dependencies, you can also install them using the `requirements.txt` file, using 

.. code-block:: bash

    pip install -r requirements.txt

From here, the algorithms are ready to be imported into your Python scripts.

