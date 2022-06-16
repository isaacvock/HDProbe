# HDProbe
R package for performing differential mutation analysis

## Using HDprobe
1. Clone repository
2. In R script, add a line where you navigate to HDProbe repository directory with setwd("/path/to/HDProbe")
3. In the next line, call devtools::load_all() (make sure you have the devtools package installed)
4. You can now call any function in HDProbe like you would with functions from an R package that you loaded with library("package-name")

In other words, the top of any script that you want to use HDProbe in should look like:
```
setwd("/path/to/HDProbe/")
devtools::load_all()
```
