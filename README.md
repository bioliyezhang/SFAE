# SFAE : A tool to  estimate relative RNA abundances of different subcellular fractionations

## About SFAE
A tool to  estimate relative RNA abundances of different subcellular fractionations. We used constrained minimization method “trust-constr” implemented in scipy.optimize.minimize function in Scipy package in Python to estimate the CR ratios(Virtanen et al. 2020). The cost function we used to minimize the differences is 

![image](https://github.com/bioliyezhang/SFAE/blob/main/formula.png)

Where represented Cytosolic RNA abundance Ratio (CR), and was corresponding to TPM value of gene  in Cytosol, Nucleus and Whole Cell (WC) fractions respectively and  was the gene number involved in estimation after filtering in preprocessing step.  
We’ve tested in multiple datasets and such cost function will normally yield a bell shape in the prediction error measurement. When dealing subcellular fractionation with more than two fractions, similar ratios can be estimated simultaneously with the same function.

![image](https://github.com/bioliyezhang/SFAE/blob/main/concept.png)

## Installation
The SFAE tool is a package written by python, you could just download `SFAE_CR_estimator.py` and directly use it.

## How to use
python SFAE_CR_estimator.py --input_file --num_component 

`--input_file <arg>` .csv file format. TPM value matrix of different fractions(columns represent fractions and rows represent genes, TPM of whole cell should be the first column).

`--num_component <arg>` the number of cell divided into fractions.



## System Requirements
Python 3.7.6

## Support
For help running the program or any questions/suggestions/bug reports, please contact daixm@shanghaitech.edu.cn and zhangly@shanghaitech.edu.cn


