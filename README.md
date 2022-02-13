# Three-component ground motion selection
This repo contains supplemental materials for the [ICOLD short course](https://www.eko.polimi.it/index.php/icold-bw2019/2019) on ground motion selection and modification (GMSM) that I developed in collaboration with [Prof. Anil Chopra](https://ce.berkeley.edu/people/faculty/chopra) and [Dr. Arnkjell Løkke](https://www.researchgate.net/profile/Arnkjell-Lokke). The short course was held in Milan, Italy on Sept. 12, 2019 and I gratefully acknowledge [Guido Mazzà](https://www.researchgate.net/profile/Guido-Mazza) and [Antonella Frigerio](https://www.researchgate.net/profile/Antonella-Frigerio) for the opportunity to participate. In this repo, [OpenSHA](https://opensha.org/) is used to perform probabilistic seismic hazard analysis whereas Matlab is used to construct [CMS-UHS Composite Spectra](https://www.11ncee.org/images/program/papers/11NCEE-000569.pdf) (at a user-specified return period) for selecting [multicomponent ground motion time series](https://peer.berkeley.edu/peer-strong-ground-motion-databases) as inputs to nonlinear RHAs of dams. To help orient users who are familiar with GMSM, elements of [Prof. Jack Baker's CS software](https://github.com/bakerjw/CS_Selection) were used and extended when applicable. Finally, step-by-step instructions for using the files are provided (both as separate PDF file and as set of screenshots towards end of presentation slides).


## References
Kwong, N.S., and A.K. Chopra. (2020). "Selecting, scaling, and orienting three components of ground motions for intensity-based assessments at far-field sites." Earthquake Spectra, 36(3), 1013-1037.

Baker, J.W., and C. Lee. (2018). "An improved algorithm for selecting ground motions to match a Conditional Spectrum." Journal of Earthquake Engineering, 22(4), 708-723.

Bozorgnia, Y., and K.W. Campbell. (2016). "Vertical ground motion model for PGA, PGV, and linear response spectra using the NGA-West2 database." Earthquake Spectra, 32(2), 979-1004.

Bozorgnia, Y., et al. (2014). "NGA-West2 research project." Earthquake Spectra, 30(3), 973-987.

Gulerce, Z., and N.A. Abrahamson. (2011). "Site-specific design spectra for vertical ground motion." Earthquake Spectra, 27(4), 1023-1047.

Field, E.H., T.H. Jordan, and C.A. Cornell. (2003). "OpenSHA: A developing community-modeling environment for seismic hazard analysis." Seismological Research Letters, 74(4), 406-419.