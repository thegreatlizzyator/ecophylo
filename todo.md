---
output:
  html_document: default
  pdf_document: default
---
TODO list :

## Article manuscrit

### intro

- [x] main ideas
- [ ] develop paragraphs

### M&M

- [x] Write
- [ ] Reread

### Result

- [ ] Vignette
- [x] Conceptual figure

### Discussion

- [ ] highlight novelty of model
- [ ] discuss limiting hypothesis (neutral)
- [ ] discuss limiting hypothesis (spatial struct)
- [ ] discuss limiting hypothesis (Coalescent with interaction)
- [ ] examples of uses

## Vignette peda and examples

- [ ] check for python vignette interest
- [ ] input refs in Rmd.

### intro

- [x] Write
- [ ] Reread

### Installation

- [ ] Explain workflow (python in R)

- [x] Write
- [ ] Install tuto more friendly

### Sim phylogeny

- [x] Constant size model
- [x] Demo fluctuation
- [ ] Vicariance migration
- [ ] make a example for timeframes usage

### Intensive sim with priors

- [ ] find case study example

### Conclusion

- [ ] Write
- [ ] Reread

## Codebase packaging and doc

- [x] Deploy on Pypi
- [ ] Write common README.md for Pypi and github.
- [ ] get error test out in special test file

### dosimulate.py

- [ ] Doc
- [ ] Parameters type
- [ ] Parameters doc
- [ ] Returns
- [ ] Test/examples
- [ ] Check loop trycatch

### islmodel.py

- [ ] population_configuration (11010)
- [ ] migration_matrix (11010)
- [ ] mass_migration (10010)

### pastdemo.py

- [x] timeframes (11111)
- [x] demographic_events (11111)
- [ ] check the case where T = 0.
- [ ] check the case where size are 0 at one time.

### toPhylo.py

- [x] toPhylo function (11111)
- [ ] Explain models in more detail
- [ ] give the options to simulate

- [x] ubranch_mutation function (11111)
- [ ] protracted speciation


Achievement in function is coded by 5 values 0 or 1 if the following tests are done : 

Doc ; Parameters type ; Parameters doc ; Returns ; Test/examples