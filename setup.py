
from setuptools import setup

// This should no be working
classifiers=[
    'Development Status :: 3 - Alpha',
    'Intended Audience :: Developers',
    'Intended Audience :: Education',
    'Intended Audience :: Manufacturing',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: MIT License',
    'Natural Language :: English',
    'Operating System :: MacOS',
    'Operating System :: Microsoft :: Windows',
    'Operating System :: POSIX',
    'Operating System :: POSIX :: BSD',
    'Operating System :: POSIX :: Linux',
    'Operating System :: Unix',
    'Programming Language :: Python',
    'Programming Language :: Python :: 2',
    'Programming Language :: Python :: 2.6',
    'Programming Language :: Python :: 2.7',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.3',
    'Programming Language :: Python :: 3.4',
    'Programming Language :: Python :: Implementation :: CPython',
    'Topic :: Education',
    'Topic :: Scientific/Engineering :: Chemistry',
    'Topic :: Scientific/Engineering :: Physics',
]


setup(
  name = 'thermosolver',
  packages = ['thermosolver'],
  license='MIT',
  version = '0.0.02',
  description = 'Under development',
  author = 'Jhonnatha de Andrade Monteiro',
  long_description=open('README.rst').read(),
  author_email = 'jhonnatha.am@gmail.com',
  url = 'https://github.com/JhonnathaAndrade/thermosolver' ,
  keywords = ['chemical engineering', 'chemistry', 'mechanical engineering', 
  'thermodynamics', 'databases', 'cheminformatics'],
  classifiers = classifiers,
  package_data={'thermosolver': ['Critical Properties/*']}
)
