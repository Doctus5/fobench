project = "FoBench"
copyright = "2026, Sergio Diaz-Meza, Jonas Pätzel"
author = "Sergio Diaz-Meza, Jonas Pätzel"
release = "0.1"

import os
import sys
sys.path.insert(0, os.path.abspath("../../"))

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ["sphinx_last_updated_by_git",
              "sphinx.ext.autodoc",
              "sphinx.ext.coverage",
              "sphinx.ext.napoleon",
              "sphinx.ext.viewcode",
              "sphinx.ext.intersphinx",
              "sphinx.ext.autosummary",
              "sphinx_copybutton",
              "sphinx_togglebutton",
              "sphinx_autodoc_typehints",
              "sphinx_design"]

templates_path = ["_templates"]
exclude_patterns = []
# napoleon_include_init_with_doc = True
autoclass_content = 'both'
autodoc_default_options = {"members": True,
                           "private-members": True}
autosummary_generate = True

typehints_defaults = "braces-after"
always_use_bars_union = True

autosummary_recursive = True
autosummary_generate_overwrite = True
autosummary_imported_members = False

intersphinx_mapping = {
    "scipy": ("https://docs.scipy.org/doc/scipy/", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "matplotlib": ('https://matplotlib.org/stable/', None),
    "python": ('https://docs.python.org/3', None),
    "obspy": ("https://docs.obspy.org/", None),
    "pandas": ("https://pandas.pydata.org/docs/", None),
    "zarr": ("https://zarr.readthedocs.io/en/stable/", None),
}

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
html_logo = "_static/logo.png"
html_favicon = "_static/logo.png"
html_title = "FoBench Docs"

html_theme_options = {
    "version_selector": True,
    "style_external_links": False,
}