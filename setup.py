from setuptools import setup

exec(open('intloc/version.py').read())

with open("README.md", "r", encoding="utf-8") as ld:
    long_description = ld.read()

setup(name='integration_locator',
      version=__version__,
      description='A tool for the detection of DNA integrations in reference genomes using long sequencing reads',
      long_description=long_description,
      long_description_content_type='text/markdown',
      url='https://github.com/simakro/Integration_locator/tree/2021-3',
      author='Simon Magin',
      author_email='simon.magin@uk-essen.de',
      license='to be MIT but not yet',
      packages=[
        'intloc',
        'intloc.il_aux',
        'intloc.il_aux.ideo_gen',
        'intloc.il_resources',
        'intloc.il_test',
        'intloc.il_test.Test_data',
        'intloc.il_test.Test_results',
        'intloc.il_test.read_dir_test',
        'intloc.il_test.integrator_reports',
        'intloc.il_resources.genome_assemblies'
        ],
       include_package_data=True,
#      package_data={
#        'intloc.il_resources': ["*"],
#        'intloc.il_test': ["*"],
#        'intloc.il_test.Test_data': ["*"],
#        'intloc.il_test.Test_results': ["*"],
#        'intloc.il_test.read_dir_test': ["*"],
#        'intloc.il_test.integrator_reports': ["*"]
#      },
      # scripts=['intloc/integration_locator'],
      entry_points={
      "gui_scripts":['intloc_gui = intloc.intloc_gui:gui_main'],
      "console_scripts":[
          'integration_locator = intloc.integration_locator:main',
          'intloc = intloc.integration_locator:main'
          ]
      },
      install_requires=[
          'matplotlib',
          'svglib',
          'reportlab>=3.5.54',
          'pillow',
          'gooey',
          'scipy>=1.5',
          'pandas',
      ],
      python_requires='>=3.6',
      classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
      ],
      )