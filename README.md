# KOBEsim

This is a Bayesian scheduler written in Python for RV surveys. Once there is an emerging signal in your data, KOBEsim finds the optimum next observing date to maximize the efficiency confirming or discarding that periodicity ([Balsalobre-Ruza et al., 2023](https://ui.adsabs.harvard.edu/abs/2022arXiv221011207B/abstract)).



<p align="center">
<img src="https://user-images.githubusercontent.com/47603865/188274702-4b41f705-4c27-4493-a853-eda4283b92cc.png" width="600" />

## Warning
This documentation is under construction. We are working on provide more information on the code and its usage. If you have any question or suggestion, please contact me: obalsalobre@cab.inta-csic.es.
  
  
## Installation

1) Download this folder or clone it to your computer.
2) Our code makes use of Bayev ([Díaz et al. 2016](https://ui.adsabs.harvard.edu/abs/2016A%26A...585A.134D/abstract)). In order to run KOBEsim, you need to download IN THE SAME DIRECTORY THAN STEP 1 the folder available in: https://github.com/exord/bayev. This new folder must be named: *bayev*.

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
python run_KOBEsim.py -obs 37.22 -2.55 2168 -star hd147379 -file 'path/RVdata.ascii'
```

You can customize the setting by giving additional inputs (see Appendix of [Balsalobre-Ruza et al., 2023](https://ui.adsabs.harvard.edu/abs/2022arXiv221011207B/abstract)).
  
## Citation
Please, if you use this code cite our work as: "O. Balsalobre-Ruza, J. Lillo-Box, A. Berihuete, et al., 2023, A&A, 669, A18".

```
  @ARTICLE{2022arXiv221011207B,
       author = {{Balsalobre-Ruza}, O. and {Lillo-Box}, J. and {Berihuete}, A. and {Silva}, A.~M. and {Santos}, N.~C. and {Castro-Gonz{\'a}lez}, A. and {Faria}, J.~P. and {Hu{\'e}lamo}, N. and {Barrado}, D. and {Demangeon}, O.~D.~S. and {Marfil}, E. and {Aceituno}, J. and {Adibekyan}, V. and {Azzaro}, M. and {Barros}, S.~C.~C. and {Bergond}, G. and {Galad{\'\i}-Enr{\'\i}quez}, D. and {Pedraz}, S. and {Santerne}, A.},
        title = "{$\texttt{KOBEsim}$: a Bayesian observing strategy algorithm for planet detection in radial velocity blind-search surveys}",
      journal = {arXiv e-prints},
     keywords = {Astrophysics - Earth and Planetary Astrophysics, Astrophysics - Instrumentation and Methods for Astrophysics},
         year = 2022,
        month = oct,
          eid = {arXiv:2210.11207},
        pages = {arXiv:2210.11207},
archivePrefix = {arXiv},
       eprint = {2210.11207},
 primaryClass = {astro-ph.EP},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2022arXiv221011207B},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```

