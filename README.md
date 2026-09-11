# KOBEsim

**KOBEsim** is a Bayesian scheduler written in Python for radial velocity (RV) surveys. Once an emerging signal is identified in your data, KOBEsim finds the optimum next observing date to maximize the efficiency of confirming or discarding the periodicity as a Keplerian signal.

For a description of the method, see [Balsalobre-Ruza et al. (2023)]([https://ui.adsabs.harvard.edu/abs/2022arXiv221011207B/abstract](https://ui.adsabs.harvard.edu/abs/2023A%26A...669A..18B/abstract)).

<p align="center">
<img src="https://user-images.githubusercontent.com/47603865/188274702-4b41f705-4c27-4493-a853-eda4283b92cc.png" width="600" />
</p>

---

## Requirements

KOBEsim requires:

* Python 3.10.x
* The Python packages listed in [`requirements.txt`](requirements.txt)
* [Bayev](https://github.com/exord/bayev), which is maintained as a separate repository

> **Important:** `bayev` must be downloaded separately and placed alongside the KOBEsim repository, as described below.

---

## Installation

### 1. Clone KOBEsim

Clone this repository to your computer:

```bash
git clone https://github.com/olgabalsa/KOBEsim.git
```

Alternatively, you can download the repository as a ZIP file from GitHub.

### 2. Download Bayev

KOBEsim makes use of **Bayev** ([Díaz et al. 2016](https://ui.adsabs.harvard.edu/abs/2016A%26A...585A.134D/abstract)) to compute the Bayes Factor based on the Perrakis estimator.

Clone the Bayev repository separately, inside KOBEsim repository:

```bash
cd KOBEsim
git clone https://github.com/exord/bayev.git
```

The repositories must be organized as follows:

```text
KOBEsim/
└── bayev/
```

The directory containing the Bayev code must be named `bayev`.

### 3. Install uv

Use uv to create the Python environment and ensure that the required Python version is used.

On macOS and Linux, install uv with:
```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```
Then **restart your terminal** or **follow the instructions displayed** by the installer.

### 4. Create the Python environment

From the KOBEsim directory, create a virtual environment using Python 3.10:
```bash
uv venv --python 3.10 kobesim-env
```
uv will automatically download Python 3.10 if it is not already available on your computer.

Activate the environment:
```bash
source kobesim-env/bin/activate
```

### 5. Install the Python dependencies

Install the required packages:
```bash
uv pip install -r requirements.txt
```

### 6. Check the installation

From the KOBEsim directory, run:
```bash
python run_KOBEsim.py --help
```

If the installation is correctly configured, this should display the available command-line options.

## Usage

See a usage example in this [Jupyter notebook](https://github.com/olgabalsa/KOBEsim/blob/main/example/run_example.ipynb).
KOBEsim requires three main inputs:

* **Observatory:** either the name of a supported observatory or its coordinates.
* **Target:** the name of the target star (to be resolve by SIMBAD).
* **Previous RV data:** a file containing the RV measurements obtained so far.

### Using a predefined observatory

For observatories included in KOBEsim, specify the observatory name with `-obs_n`:

```bash
python run_KOBEsim.py -obs_n CAHA -star KOBE-1 -file example/data/mock_rv_20Me_60d_first15.csv
```

### Providing observatory coordinates

Alternatively, you can provide the observatory coordinates directly using:

```text
latitude longitude height
```

where the latitude and longitude are given in degrees and the height is given in meters.

For example:

```bash
python run_KOBEsim.py -obs 37.22 -2.55 2168 -star KOBE-1 -file example/data/mock_rv_20Me_60d_first15.csv
```

### Additional options

KOBEsim allows the user to customize the observing strategy through additional command-line arguments.

For a more complete description of the available parameters and their recommended use, see the Appendix of [Balsalobre-Ruza et al. (2023)](https://ui.adsabs.harvard.edu/abs/2022arXiv221011207B/abstract).

To have the most updated list, you can display the available command-line options with:

```bash
python run_KOBEsim.py --help
```

---

## Input data

KOBEsim accepts RV time series provided either as a **FITS file** or as a **Text/CSV file**.

### Text/CSV format

For an example of an input file, see the [example](https://github.com/olgabalsa/KOBEsim/blob/main/example/create_dataset_example.ipynb).
For text or CSV files, the file must contain the following three columns:

| Column | Description                        |
| ------ | ---------------------------------- |
| `jd`   | Observation time, in Julian Date   |
| `rv`   | Radial velocity                    |
| `erv`  | Uncertainty on the radial velocity |

The column names must be exactly `jd`, `rv`, and `erv`. All three columns are required.

The file is read using `pandas.read_csv`, so columns should be comma-separated. Lines beginning with `#` are treated as comments and ignored.

For example:

```text
jd,rv,erv
2459001.123,12.34,1.25
2459003.456,10.87,1.18
2459007.891,13.21,1.32
```

The units of `rv` and `erv` should be consistent throughout the input file (e.g. m/s). The time values must be given as Julian Dates.

### FITS format

KOBEsim also accepts FITS files containing the following columns in the first extension:

* `OBJ_DATE_BJD`: observation time
* `SPECTRO_CCF_RV`: RV
* `SPECTRO_CCF_ERV`: RV uncertainty

The code converts `OBJ_DATE_BJD` to the `jd` convention used internally by adding `2400000`.

---

## Output

KOBEsim produces two main outputs:

1. A **CSV file** containing the tested observing dates (orbital phases), ranked according to their preference. The most preferred observing date appears first in the list.

2. A **plot showing the expected increase in the detection metric (Bayes Factor)** for each tested orbital phase (observing date). The color coding indicates the preference assigned to each date, making it easy to identify the most promising observing windows.

Both outputs are automatically saved in an `outputs/` directory created in the current working directory:

```text
outputs/
├── [CSV output file]
└── [output plots]
```

The CSV file can be used to inspect and further analyze the ranking of the proposed observing dates, while the plot provides a visual representation of the expected gain in the Bayes Factor across the tested orbital phases.


---

## Citation

If you use KOBEsim in your research, please cite:

**Balsalobre-Ruza, O., Lillo-Box, J., Berihuete, A., et al. 2023, A&A, 669, A18.**

```bibtex
@ARTICLE{2023A&A...669A..18B,
       author = {{Balsalobre-Ruza}, O. and {Lillo-Box}, J. and {Berihuete}, A. and {Silva}, A.~M. and {Santos}, N.~C. and {Castro-Gonz{\'a}lez}, A. and {Faria}, J.~P. and {Hu{\'e}lamo}, N. and {Barrado}, D. and {Demangeon}, O.~D.~S. and {Marfil}, E. and {Aceituno}, J. and {Adibekyan}, V. and {Azzaro}, M. and {Barros}, S.~C.~C. and {Bergond}, G. and {Galad{\'\i}-Enr{\'\i}quez}, D. and {Pedraz}, S. and {Santerne}, A.},
        title = "{KOBEsim: A Bayesian observing strategy algorithm for planet detection in radial velocity blind-search surveys}",
      journal = {\aap},
     keywords = {planets and satellites: detection, methods: statistical, techniques: radial velocities, stars: solar-type, Astrophysics - Earth and Planetary Astrophysics, Astrophysics - Instrumentation and Methods for Astrophysics},
         year = 2023,
        month = jan,
       volume = {669},
          eid = {A18},
        pages = {A18},
          doi = {10.1051/0004-6361/202243938},
archivePrefix = {arXiv},
       eprint = {2210.11207},
 primaryClass = {astro-ph.EP},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2023A&A...669A..18B},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```

Please also cite the Bayev method when appropriate:

**Díaz et al. (2016), A&A, 585, A134.**

---

## Contact

If you have questions or feedback, please contact:

**Olga Balsalobre-Ruza**
[o.balsaruza@gmail.com](mailto:o.balsaruza@gmail.com)
