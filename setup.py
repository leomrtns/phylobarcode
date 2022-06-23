from pathlib    import Path
import setuptools 
import sys
min_version = (3, 8)

if sys.version_info < min_version:
    error = """
Python {0} or above is required.

Make sure you have an up-to-date pip installed.  
""".format('.'.join(str(n) for n in min_version)), sys.exit(error)

base_dir = Path(__file__).parent.resolve()
version_file = base_dir / "phylobarcode/__version__.py"
readme_file = base_dir / "README.md"

# Eval the version file to get __version__; avoids importing our own package
with version_file.open() as f:
    exec(f.read())

with readme_file.open(encoding = "utf-8") as f:
    long_description = f.read()

setuptools.setup(
    name = "phylobarcode",
    version = __version__,
    author = "Leonardo de Oliveira Martins",
    author_email = "Leonardo.de-Oliveira-Martins@quadram.ac.uk",
    description = "in-silico amplicons optimised for phylogenetic signal",
    long_description = long_description,
    long_description_content_type = "text/markdown",
    keywords = "amplicon,ribosomal_operon",
    url = "https://github.com/quadram-institute-bioscience/phylobarcode",
    project_urls = {
        "Source": "https://github.com/quadram-institute-bioscience/phylobarcode",
    },
    packages = setuptools.find_packages(),
    include_package_data=True,
    package_data = {'phylobarcode': ['data/*']},
    data_files = [("", ["LICENSE"])],
    python_requires = '>={}'.format('.'.join(str(n) for n in min_version)),
    license='GPLv3+',
    install_requires=[
           #'biopython>1.70',
           #'parasail>1.2.3',
           #'numpy>1.22.2',
           #'pandas>1.4.1',
           #'xxhash >= 0.8.0', # equiv to python-xxhash 2.0.0
           'biopython',
           'parasail',
           'numpy',
           'pandas',
           'xxhash', 
           'scikit-learn'
       ],
    classifiers = [
        "Development Status :: 2 - Pre-Alpha",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        # Python 3 only
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10", # conda does not have parasail-python for py3.10, only 3.9...
    ],
    # Installs a "phylobarcode" program which calls phylobarcode.__main__.main()
    #   https://setuptools.readthedocs.io/en/latest/setuptools.html#automatic-script-creation
    entry_points = {
        "console_scripts": [ "phylobarcode = phylobarcode.pb_main:main" ]
    }
)
