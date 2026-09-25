# Third-party notices

## ClinGen SVI computational calibration code (Pejaver et al., 2022)

The following files are an R translation (derivative work) of the MATLAB code in
https://github.com/vpejaver/clingen-svi-comp_calibration, commit
`c56d32d4db45389bbd386c909cd49a625e3e2dc7`, and are distributed under the original MIT license:

- `calibration/R/pejaver_port.R` and its frozen copy `executed_pipeline/scripts/04_pejaver_port.R`
- `calibration/R/plot_posteriors.R` and its frozen copy `executed_pipeline/scripts/06_plot_posteriors.R`
  (reproduces the figure layout of `plotting/plot_both_posteriors_pub.m`)

```
MIT License

Copyright (c) 2022 Vikas Pejaver

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

The data associated with that study (Zenodo 10.5281/zenodo.8347415, CC-BY-4.0) are used for
validation but are not redistributed here.

Please also cite: Pejaver V, et al. Calibration of computational tools for missense variant
pathogenicity classification and ClinGen recommendations for PP3/BP4 criteria.
*Am J Hum Genet.* 2022;109(12):2163–2177.
