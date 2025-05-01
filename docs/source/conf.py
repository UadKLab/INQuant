# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'INQuant'
copyright = '2025, Annekatrine Kirketerp-Møller, Ida Sofie Goldschmidt'
author = 'Annekatrine Kirketerp-Møller, Ida Sofie Goldschmidt'
release = '2.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',   
    'sphinx.ext.autosummary',
]

autosummary_generate = True  # Optional: auto-create .rst files for each module

templates_path = ['_templates']
exclude_patterns = []

autodoc_member_order = 'bysource' 



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "pydata_sphinx_theme"
html_static_path = ['_static']
