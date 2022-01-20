from setuptools import setup

setup(
    name='fasttreerep',
    version='1.0.0',
    py_modules=['cli'],
    install_requires=[
        'Click',
    ],
    entry_points={
        'console_scripts': [
            'fasttreerep = cli:run_fasttree',
        ],
    },
)