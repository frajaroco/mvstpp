# mvstpp: Mark variograms for spatio-temporal point processes 

Basic toolkit for the exploration and analysis of the spatio-temporal point patterns through spatial and temporal mark variograms. This repository is based on `spatstat`, `splancs`, `stpp`, `geoR` and `KernSmooth` packages.

## Installation guide

The easiest way to install the development version of mvstpp from github is using the devtools package which can be installed run the next command:
```
install.packages('devtools', dependencies=TRUE)
```
and thereafter run the commands:
```
require(devtools)
install_github('frajaroco/mvstpp')
```
**warning:** If you have problems to load the repository and work under Microsoft Windows. It is recommended load the package `geoR` before to load the `msftpp` repository.

## References
- [Stoyan, D., Rodríguez-Cortés, F. J., Mateu, J., and Gille, W. (2017). Mark variograms for spatio-temporal point processes. *Spatial Statistics*, **20**:125-147.](http://www.sciencedirect.com/science/article/pii/S2211675317300696)
- [González, J. A., Rodríguez-Cortés, F. J., Cronie, O. and Mateu, J. (2016). Spatio-temporal point process statistics: a review. *Spatial Statiscts*, **18**:505-544.](http://www.sciencedirect.com/science/article/pii/S2211675316301130)

## CiteBibtex
```
@misc{r15,
	author = {Francisco J. Rodr\'{i}guez-Cort\'{e}s},
	title = {mvstpp: {M}ark {V}ariogram for {S}patio-{T}emporal {P}oint {P}rocesses},
	year = {2015},
	note = {GitHub repository},
	url = {\url{https://github.com/frajaroco/mvstpp}}
}
```
