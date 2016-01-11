from setuptools import setup

setup(name='odradio',
      version='0.0.1',
      description='HPTCAD',
      url='http://github.com/manstetten/odradio',
      author='Paul Manstetten',
      author_email='manstetten@iue.tuwien.ac.at',
      license='MIT',
      packages=['odradio'],
      install_requires=[
          'pandas',
          'numpy',
          'matplotlib'
      ],      
      zip_safe=False)