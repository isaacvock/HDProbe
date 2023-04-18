# HDProbe
R package for performing differential mutation analysis. Implements [bakR's](https://github.com/simonlabcode/bakR) hierarchical model, minus the mixture model.

## Installing HDProbe with devtools

The easiest way to install HDProbe is with the `devtools` R package, as follows:

```
install.packages("devtools") # if not already installed
devtools::install_github("isaacvock/HDProbe")
```

You can then load HDProbe in an R session like any other R package with `library(HDProbe)`.

## Using HDProbe for development purposes

If you want to be able to edit HDProbe and easily implement the edited, development version of HDProbe, then follow these instructions:

1. Clone repository
2. In R script, add a line where you navigate to HDProbe repository directory with setwd("/path/to/HDProbe")
3. In the next line, call devtools::load_all() (make sure you have the devtools package installed)
4. You can now call any function in HDProbe like you would with functions from an R package that you loaded with library("package-name")

In other words, the top of any script that you want to use HDProbe in should look like:
```
setwd("/path/to/HDProbe/")
devtools::load_all()
```

## 

