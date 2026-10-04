# Print an OncoArm object

Prints the calibrated parameters of a treatment group together with the
implied medians of PFS and OS (overall and by response) and the implied
Pearson correlations among PFS, OS and response.

## Usage

``` r
# S3 method for class 'OncoArm'
print(x, digits = 4, ...)
```

## Arguments

- x:

  An object of class `OncoArm`.

- digits:

  Number of significant digits.

- ...:

  Not used.

## Value

`x`, invisibly.

## Examples

``` r
OncoArm(pfs.median = 6, orr = 0.3, resp.cor = 0.4, death.prop = 0.15,
        os.median = 15)
#> OncoArm
#>   OS model           : idm
#>   Response timing    : none
#>   PFS hazard         : 0.1155
#>   Response rate      : 0.3
#>   Copula correlation : 0.5376
#>   Death proportion   : 0.15
#>   Post-progression hazard (non-responders, responders): 0.08767, 0.08767
#> Implied medians (all, responders, non-responders)
#>   PFS: 6, 11.39, 4.393
#>   OS : 15, 20.47, 12.89
#> Implied correlations
#>   Corr(PFS, R) = 0.4, Corr(OS, R) = 0.2435, Corr(PFS, OS) = 0.6089
#>   Partial correlation of OS and R given PFS = 3.818e-17
```
