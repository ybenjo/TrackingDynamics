# TrackingDynamics
## About
[Tracking Dynamics of Topic Trends Using a Finite Mixture Model(KDD 2004)](http://dl.acm.org/citation.cfm?id=1016919)
## TODO

- Initialize \mu in approprietly.  

> , and set \mu_i^{(0)} the first x_t s that are different each other.  
  
  
## Algorithm
Fitting Gaussian distribution in streaming documents. Parameters in step t are updated using parameters in step t-1.  
## Usage

    make ttfmm  
    ./ttfmm input_file alpha lambda k_max  
    

- alpha[NUM]: Parameter of \alpha( > 0)
- lambda[NUM]: Parameter of \labmda(0 < lambda < 1)  
- k_max[NUM]: # of maximum topic( > 0)  

## Format
Each line is one document.

    timestamp \t word_1 \t word_2 \t ...  
    timestamp \t word_1 \t word_2 \t ...  


## Memo
- Small \lambda causes inf in \lambda^{-(t\_{new} - t\_{old})}.  
- Dont calculate denominator in p(i|x\_t) because |\Sigma| sometimes zero.  
- Using logsumexp in p(i|x\_t) to avoid overflow.  
