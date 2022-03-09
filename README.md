# SPAE
Stellar Parameters, Abundances, and Errors

SPAE is a Python code that uses a state-of-the-art Bayesian method to self-consistently propagate uncertainties from the stellar atmosphere solutions in calculating individual abundances. It employs the LTE plane-parallel spectral analysis code MOOG (Sneden 1973; https://www.as.utexas.edu/~chris/moog.html) and Kurucz model atmosphere grids (http://kurucz.harvard.edu/grids.html) to derive the abundances of elements included in a user-provided linelist. 


#Cloning and Installing SPAE
Clone SPAE from Github to a directory on your local machine:
git clone https://gihub.com/simon-schuler/SPAE

Install SPAE into a Python environment:
pip install ./
