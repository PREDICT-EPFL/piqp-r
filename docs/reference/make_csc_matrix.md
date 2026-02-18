# Convert a plain matrix or simple triplet form matrix to a [Matrix::dgCMatrix](https://rdrr.io/pkg/Matrix/man/dgCMatrix-class.html) (implicit) form

Convert a plain matrix or simple triplet form matrix to a
[Matrix::dgCMatrix](https://rdrr.io/pkg/Matrix/man/dgCMatrix-class.html)
(implicit) form

## Usage

``` r
make_csc_matrix(x)
```

## Arguments

- x:

  a matrix or a simple triplet form matrix

## Value

a list of row pointer, column pointer, and values corresponding to a
[Matrix::dgCMatrix](https://rdrr.io/pkg/Matrix/man/dgCMatrix-class.html)
object
