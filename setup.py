from setuptools import setup

setup(
    name='fasttree',
    version='1.0.0',
    py_modules=['cli', 'main', 'counter', 'distance_measures', 'plot_tree', 'topology', 'tree'],
    install_requires=[
        'Click',
        'newick',
        'numpy',
    ],
    entry_points={
        'console_scripts': [
            'fasttree = cli:run_fasttree',
        ],
    },
)