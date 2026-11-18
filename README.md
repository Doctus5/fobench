# FoBench: An open source fibre optic sensing toolbox

<div align="center">
<picture>
    <img src="./docs/source/_static/welcome_icon_light.png" alt="FoBench Logo" width="200px">
</picture>
</div>


Welcome to **FoBench**! This repository contains code for reading and processing Fibre Optic Sensing data. Your help in testing and improving this code is invaluable, and any suggestions regarding style, features, or bug fixes are highly appreciated!

Follow this guide or check out the [documentation](https://doctus5.github.io/fobench/) to get started with FoBench.



---

## Table of Contents
1. [Installation](#installation)
   - [Install Anaconda/Miniconda/Conda](#install-anacondaminicondaconda)
   - [Create a Conda Environment](#create-a-conda-environment)
2. [Installing FoBench](#installing-fobench)
3. [Contributing](#contributing)

---

## Installation

### Install Anaconda/Miniconda/Conda

We will use Conda to create a virtual environment to safely manage the dependencies for FoBench.

1. Download the latest version of Miniconda from [this link](https://docs.conda.io/en/latest/miniconda.html#linux-installers).

2. Run the following commands to install Miniconda and update Conda:

   ```bash
   bash Miniconda3-latest-Linux-x86_64.sh
   conda update conda
   ```

3. Change the main Conda channel to `conda-forge` (this avoids issues related to licensing):

   ```bash
   conda config --add channels conda-forge
   ```

4. Close and reopen the terminal to apply changes.

---

### Create a Conda Environment

You can create a minimal environment for FoBench with the following command:

```bash
conda create --name fobench_env python=3.11
```

Once created, activate the environment:

```bash
conda activate fobench_env
```

When you're done, deactivate the environment by running:

```bash
conda deactivate
```

---

## Installing FoBench

With the environment activated, install FoBench using `pip`:

```bash
pip install git+https://github.com/Doctus5/fobench.git
```

If you already have an installed version of FoBench, we recommend to regularly update it to receive the latest functionality and bug-fixes. To update, use the following command:

```bash
pip install --upgrade fobench
```

Once installed, you're ready to start using the code in your Python environment!

---

## Contributing

Your feedback is crucial to the development of this code. If you have suggestions, please feel free to reach out or submit an issue. You can contact:

- Sergio Diaz-Meza: [sergioad@gfz-potsdam.de](mailto:sergioad@gfz-potsdam.de)
- Jonas Pätzel: [jonas.patzel@ulb.be](mailto:jonas.patzel@ulb.be)

All contributions are welcome, including code improvements, bug reports, and feature suggestions!