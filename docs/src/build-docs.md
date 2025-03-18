# Build the documentation

The documentation may be built locally using the following procedure.

```bash
pip install sphinx sphinx_design sphinx-autodoc2 \
    myst-parser sphinx-copybutton furo
cd docs
sphinx-build -b html . ../public
```

This will generate a directory `public` in the top-level directory. Open `public/index.html` in a browser to access documentation.

