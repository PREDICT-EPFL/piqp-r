# Convert a plain matrix or simple triplet form matrix to a [Matrix::dsCMatrix](https://rdrr.io/pkg/Matrix/man/dsCMatrix-class.html) (implicit, upper) form

Convert a plain matrix or simple triplet form matrix to a
[Matrix::dsCMatrix](https://rdrr.io/pkg/Matrix/man/dsCMatrix-class.html)
(implicit, upper) form

## Usage

``` r
make_csc_symm_matrix(m)
```

## Arguments

- m:

  a matrix or a simple triplet form matrix

## Value

a list of row pointer, column pointer, and values corresponding to a
[Matrix::dsCMatrix](https://rdrr.io/pkg/Matrix/man/dsCMatrix-class.html)
object
