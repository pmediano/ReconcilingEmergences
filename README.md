Reconciling emergences: Identifying causal emergence in multivariate data
=========================================================================

This code provides basic functionality to test hypotheses about causal
emergence in time series data. This repository contains the companion software
for the paper:

Rosas\*, Mediano\*, et al. (2020). _Reconciling emergences: An
information-theoretic approach to identify causal emergence in multivariate
data_. PLoS Computational Biology, 16(12): e1008289. DOI:
[10.1371/journal.pcbi.1008289](https://doi.org/10.1371/journal.pcbi.1008289)

Please cite the paper (and give us a shout!) if you use this software. Please
contact Pedro Mediano for bug reports, pull requests, and feature requests.


Basic theory
------------

The paper provides a mathematical account of emergence in complex dynamical
systems. In a nutshell, a macroscopic feature `V` of a multivariate system `X`
is said to be _causally emergent_ if it contains some information about the
future of `X` that no microscopic element of `X` has on its own. Furthermore,
it is possible to distinguish between two kinds of emergence: downward
causation (when macroscopic features have an effect on microscopic elements),
and causal decoupling (when macroscopic features have an effect on other
macroscopic features). These intuitions are formulated rigorously using the
[Partial Information Decomposition](https://arxiv.org/abs/1004.2515) and the
[Integrated Information Decomposition](https://arxiv.org/abs/1909.02297)
frameworks.

Unfortunately, information decomposition is in general difficult and
computationally expensive (for now). Fortunately, the theory also provides a
set of practical criteria (Ψ, Δ, and Γ) to detect emergence which can be
used when a given macroscopic feature `V` is known. In particular:

* If Ψ > 0, then `V` is causally emergent.
* If Δ > 0, then `V` shows downward causation.
* If Ψ > 0 and Γ = 0, then `V` shows causal decoupling.

Note, however, that the converse does not hold: a negative Ψ does not imply
that `V` is _not_ causally emergent -- see section III.A of the paper for
details.

Finally, note that for the assumptions of the theory to hold, the candidate
emergent feature V has to be a _supervenient_ feature of X -- in other words,
`V(t)` has to be a (possibly stochastic) function of `X(t,:)`, and nothing
else.


Code examples
-------------

As a simple example, consider the bivariate binary system from Fig. 1 (right)
from the paper:

```octave
T = 1000;
X = zeros([T,2]);
for t=2:T
  X(t,1) = xor(X(t-1,1), X(t-1,2));
  X(t,2) = rand < 0.5;
end

V = xor(X(:,1), X(:,2));
psi   = EmergencePsi(X, V);
delta = EmergenceDelta(X, V);
```

You will see that, as expected, Ψ is close to zero but Δ is positive and
close to one.

Note, however, that the values are not exactly zero or one, as any estimation
on a finite dataset is bound to have some bias and variance. To estimate the
variance and correct the bias, it is common to perform surrogate data tests by
time shuffling. Assuming that `T = size(X,1)` and that `nb_surr` is a
large(-ish) integer, this can be done by:

```octave
surr_psi = arrayfun(@(j) EmergencePsi(X(randperm(T),:), V(randperm(T))), 1:nb_surr);
```

With this, the bias-corrected Ψ is `psi - mean(surr_psi)` and its standard
deviation is `std(surr_psi)`.


Download and installation
-------------------------

None, really. Just download this folder (with Github's zip download or with
`git clone`), add it to your Octave/Matlab path, and enjoy.

Note, however, that if you use Octave you will need to install the `statistics`
package. You can do this simply by running `pkg install -forge io statistics`
(and remember to run `pkg load statistics` every time you start a new session).

Tests are provided in the `tests/` subfolder. To run them in Matlab, run
`runtests('tests/')` from this repository's root folder.


Use in Python
-------------

All functions in this repository are Octave-friendly, which means they can be
easily called from Python through the wholesome
[oct2py](https://oct2py.readthedocs.io/) package. With a functional Octave and
Python installation, you can simply run (from the root folder of the repo):

```python
import numpy.random as rn
from oct2py import Oct2Py

oc = Oct2Py()
oc.EmergencePsi(rn.randn(100,2), rn.randn(100,1))
```


Licence
-------

This software is distributed under the modified 3-clause BSD Licence.


Further reading
---------------

* P. Mediano\*, F. Rosas\*, _et al._ (2019). Beyond integrated information: A
  taxonomy of information dynamics phenomena.
  [arXiv:1909.02297](https://arxiv.org/abs/1909.02297)

* P. Williams and R. Beer (2010). Nonnegative decomposition of multivariate
  information. [arXiv:1004.2515](https://arxiv.org/abs/1004.2515)


\(C\) Pedro Mediano and Fernando Rosas, 2020

