# KOBEsim

This is a Bayesian scheduler in Python for RV surveys. Once there is an emerging signal in your measurements, KOBEsim finds the optimum next observing date to maximize the efficiency confirming or discarding that periodicity (O. Balsalobre-Ruza et al., submitted).



<p align="center">
<img src="https://github.com/olgabalsa/KOBEsim-code/files/9234002/logo_KOBEsim.pdf" width="250" />


## Installation

You only have to download this folder or clone it to your computer.

## Usage

To run the code you have to provide:
- The observatory
- The target
- A file with the previous data

The way of calling it from the terminal is the following:
```python
python run_KOBEsim.py -obs_n CAHA -star hd147379 -file 'path/RVdata.ascii'
```

Or with the observatory coordinates (latitude in deg, longitude in deg, and height in m):

```python
python run_KOBEsim.py -obs_n 37.22 -2.55 2168 -star hd147379 -file 'path/RVdata.ascii'
```

You can customize the setting by giving additional inputs (see Appendix of Balsalobre-Ruza et al., submitted).
