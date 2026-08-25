Download and Installation
=========================

We will use Conda to create a virtual environment to safely manage the dependencies for FoBench. You can create and activate a minimal environment for FoBench with the following commands:

.. code-block:: bash

   conda create --name fobench_env python=3.11
   conda activate fobench_env


With the environment activated, install Fobench using `pip`:

.. code-block:: bash

   pip install git+https://git.gfz-potsdam.de/sergioad/fobench.git
   
.. admonition:: Installing in development mode

    To install in editable mode, it is best to clone the repository to your target destination and then install:

    .. code-block:: bash

       git clone https://git.gfz-potsdam.de/sergioad/fobench.git
       cd fobench
       pip install -e .


If you already have an installed version of FoBench make sure to update regularly, using the following command:

.. code-block:: bash

   pip install --upgrade fobench
   
Once installed, you're ready to start using the code in your Python environment!
