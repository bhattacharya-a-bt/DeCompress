# DeCompress

Welcome to `DeCompress` - a semi-reference-free method to
deconvolve targeted panels of mRNA expression
into tissue compartment. A tissue compartment is a group
of cells of similar type or biological function (i.e. immune
or stroma or tumor compartments).
Please cite 
[Bhattacharya et al 2020](https://www.biorxiv.org/content/10.1101/2020.08.14.250902v2) 
if you use our package.

## Installation

You can install `DeCompress` using 
`devtools::install_github('bhattacharya-a-bt/DeCompress')`.
Make sure to have dependencies installed as well!


## Using DeCompress

`DeCompress` requires two input data matrices:
- the *target* matrix ($k \times n$) with  $k$ genes and $n$ samples, and
- the *reference* matrix ($K \times n$) with $K$ genes and $N$ samples with
$K > k$.

`DeCompress` also requires *a priori* knowledge of the number of tissue
compartments. If you don't know the number of compartments you wish
to deconvolve to, you can use the `findNumberCells()` function to get
an estimate from an singular value decomposition of the target matrix.

To illustrate `DeCompress`, let's generate some fake data.
Here, `pure` is a matrix of compartment-specific expression profiles,
`reference_props` and `target_props` are different $200 \times 4$ matrices
of proportions, and `reference` and `target` gives us the input expression
matrices.

```r
pure = matrix(abs(rnorm(1e4*4)),ncol=4)
reference_props = apply(matrix(abs(rnorm(200*4)),ncol=4),
        1,function(x) x/sum(x))
target_props = apply(matrix(abs(rnorm(200*4)),ncol=4),
        1,function(x) x/sum(x))

reference = pure %*% reference_props
target = pure[sample(1:1e4,400),] %*% target_props
```

### Step 1: Feature selection

Step 1 of DeCompress is to feature select a set of $K'$ genes from
the reference. 
