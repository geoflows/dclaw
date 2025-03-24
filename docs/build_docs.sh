# Build the docs with source build_docs.sh

# Depends on the following python packages. 

# sphinx = ">=7.2.6"
# sphinx_design = ">=0.5"
# sphinx-autodoc2 = "*"
# myst-parser = "*"
# sphinx-copybutton = "*"
# furo = "*"



sphinx-build -b html . ../public
