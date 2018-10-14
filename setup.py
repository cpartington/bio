from setuptools import setup, find_packages

setup(name="bio_info",
      version="0.1",
      description="Basic algorithms for DNA manipulation.",
      url="http://github.com/cpartington/bio-info",
      author="Christie Partington",
      packages=find_packages(),
      install_requires=[
          "progressbar2",
          "numpy"
      ])
