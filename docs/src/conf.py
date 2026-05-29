# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

from datetime import datetime

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'Analysis Multiverse Documentation'
copyright = '2026 Brain Health Lab, St. Francis Xavier University'
author = 'Kate Redfern'
release = '0.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

html_logo = 'icons/BrainHealthLabLogo.png'
html_favicon = 'icons/mag_brain.ico'

today = '%Y-%m-%d'

extensions = ['sphinx_rtd_theme', 'myst_parser']

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

source_suffix = {
    '.rst': 'restructuredtext',
    '.txt': 'restructuredtext',
    '.md': 'markdown',
    }

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

html_theme_options = {
    'logo_only': True,
    'prev_next_buttons_location': 'both',
    'style_nav_header_background': '#D25299',
    'titles_only': True,
    }

# -- Options for Markup ------------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-markup

rst_prolog = f"""
.. |update_month| replace:: {datetime.today().strftime('%Y-%m-%d')}
.. _`update_month`: https://www.iso.org/iso-8601-date-and-time-format.html
"""