# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'scATS'
copyright = '2025, LuChenLab'
author = 'LuChenLab'
release = '0.5.4'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [ 
    'myst_parser',
    "sphinxcontrib.mermaid",
    "sphinx_copybutton",
]

source_suffix = {
    '.rst': 'restructuredtext',
    '.md': 'markdown',
}


myst_enable_extensions = [
    "tasklist",
    "deflist",
    "dollarmath",
]

templates_path = ['_templates']
exclude_patterns = []

language = 'en'
highlight_langeuage="r,python"

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"
html_static_path = ['_static']
html_css_files = [
    'custom.css',  # 添加你的自定义CSS文件，代码块高亮等
]

html_theme_options = {
    'analytics_id': 'G-XXXXXXXXXX',  #  Provided by Google in your dashboard
    'analytics_anonymize_ip': False,
    'logo_only': False,
    'prev_next_buttons_location': 'bottom',
    'style_external_links': False,
    'vcs_pageview_mode': '',
    'flyout_display': 'hidden',
    'version_selector': True,
    'language_selector': True,
    # Toc options
    'collapse_navigation': True,
    'sticky_navigation': True,
    'navigation_depth': 4,
    'includehidden': True,
    'titles_only': False
}