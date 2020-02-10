# About CombinationsOfSparsities
CombinationsOfSparsities is a Julia package for both polynomial optimization and polynomial systems by combining [term](https://github.com/wangjie212/TSSOS) and correalative sparsities for the SDP hierarchies based on Putina's and [Putinar-Vasilescu's Positivstellensatz](https://arxiv.org/abs/1911.11428). 
## New advances
- Obtaining an approximate global optimizer of a general polynomial optimization 

     ![POP](https://github.com/maihoanganh/PV_TSSOS/blob/master/images/POP.gif)

with guarantee in theory.
- Obtaining generically a feasible solution of a basic semialgebraic set 

    ![POP](https://github.com/maihoanganh/PV_TSSOS/blob/master/images/S.gif)

with possibly uncountably many solutions.
## Installation
To use CombinationsOfSparsities in Julia, run:
```ruby
pkg> add https://github.com/maihoanganh/CombinationsOfSparsities
```
# Dependency
MOSEK (SDP solver)

# References
For more details, please refer to:
1. N. H. A. Mai, J.-B. Lasserre, and V. Magron. Positivity certificates and polynomial optimization on non-compact semialgebraic sets, 2019. Submitted.
https://arxiv.org/abs/1911.11428
2. J. Wang, V. Magron, and J.-B. Lasserre. TSSOS: a moment-SOS hierarchy that exploits term sparsity, 2019. Submitted. 
https://arxiv.org/abs/1912.08899
3. N. H. A. Mai, J. Wang, V. Magron, and J.-B. Lasserre. Combinations of correlative and term sparsities for the semidefinite hierarchies based on Putinar's and Putinar-Vasilescu's representation, 2020. Forthcoming.
