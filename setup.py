from setuptools import setup

setup(
    name='calypso',
    version='0.1',
    py_modules=['calypso'],
    entry_points={
        'console_scripts': [
            'calypso = calypso:main',
        ],
    },
)
