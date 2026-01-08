# AstroSlurp

Application to grab observable objects for a given location/time. Also uses many data sources for exploration of given targets. 

## Libraries:
### Directly Relevant
- astroplan
- astropy
- astroquery
- numpy
- dustmaps
- timezonefinder
- flask
- pytz
- datetime
- matplotlib
- io

### Possibly relevant modules:
- Skyfield (computes positions for the stars, planets, and satellites in orbit around the Earth)
- SunPy (tools for sun astronomy)
- photoutils (background estimation, star finding, source detection and extraction, aperture photometry, PSF photometry, image segmentation, centroids, radial profiles, and elliptical isophote fitting)
- regions (draw regions on images)
- specutils (representing, loading, manipulating, and analyzing astronomical spectroscopic data)
- plotly (data visualization)
- Starplot (creates star/planet/DSO charts, trajectory plotting of comets)

## Features
### Sky Query
- Returns all observable targets given latitude, longitude, elevation (in meters), date/time observation period begins, magnitude limit, and observable horizon (ie how high an object has to be in the sky to be seen due to obstacles on the observer's horizon).
- Process of Sky Query:
1. Divide sky above the observable horizon at ***any point in time*** between the beginning observation time and next nautical twilight using a grid search adjusted for latitude (using trigonometry).
2. Search these grids for ***all*** optically visible objects (galaxies, nebulas, novas, supernovas, galactic clusters, stars of all kinds, stellar clusters, etc) unless the user restricts the queried type. There are checkboxes for types but defaults to all.
3. Data sources: SIMBAD, Vizier, Gaia (all available through astroquery).
4. Displays a table with the following information: name, RA/Dec, type, magnitude, time during which object is observable. The name is a link to object details (using deep_dive_search.html).


### Data Explorer
#### Functionality
- View and analyze imagery in optical, near infrared, infrared, x-ray, radio, and gamma-ray wavelengths.
- Search for data and export in PDF, CSV, XLS, SQLite formats
#### Layout
- Left side
  - Occupies 1/4 of the horizontal space
  - Has several CSS containers the appearance of which can be changed in "static/style.css"
    - First container at top of column:  
    - Second container from top of column: user chooses to obtain data from a given source. The link is titled "Slurp Raw Data" and goes to a page called "slurp_raw_data.html" 
### Spectra
- Observe spectra for a given target or region (user definable but defaults to 1 arcminute)
- Returns images as well as a key for correlating bands with elements

### Gas, dust, and extinction
-Observe gas and dust densities around target or region (definable)
-Overlay dots representing stars hidden by gas or dust
-Can choose among images to overlay or switch between by pressing a button or use a slider for transparency

### Star cube
-View a 3d plot of objects around a point in space
-If there is an asteroid or comet, trace it's path

## Catalogues
### Via dustmaps module
-Bayestar: Queries the Bayestar 3D dust maps (Green, Schlafly, Finkbeiner et al. 2015, 2018). The maps cover the Pan-STARRS 1 footprint (dec > -30 deg) amounting to three-quarters of the sky.
-Chen2014: The 3D dust map of Chen et al. (2014), based on stellar photometry from the Xuyi Schmidt Telescope Photometric Survey of the Galactic Anticentre. The map covers 140 deg < l < 240 deg, -60 deg < b < 40 deg.
-csfd: Queries the Corrected SFD dust map of Chiang (2023). This map is based on SFD, but contains a correction to remove contamination from large-scale structure (i.e., external galaxies).
-gaia_tge: Queries the Gaia Total Galactic Extinction (Delchambre 2022) dust map, which contains estimates of monochromatic extinction, A0, in mags.
-iphas: The 3D dust map of Sale et al. (2014), based on IPHAS imaging in the Galactic plane. The map covers 30 deg < l < 115 deg, -5 deg < b < 5 deg.
-lenz2017: Queries the Lenz, Hensley & Doré (2017) dust map
-marshall: Galactic-plane 3D dust map of Marshall et al. (2006), based on 2MASS photometry.
-planck: Queries the Planck Collaboration (2016) GNILC dust map.

### Via astroquery module
catalogs
-IRSA: For all image queries, the radius may be optionally specified. If missing the radius defaults to 5 degrees. Note that radius may be specified in any appropriate unit, however it must fall in the range of 2 to 37.5 degrees.
-NED: spectra and imges
-OGLE: Galactic bulge photometry, variable stars, reddening maps, microlensing
-HST: NIR, UV, spectroscopy, astrometry, photometry
-SDSS: optical and IR spectra, Images
-UKIDSS: ultra deep NIR galaxies images optical/IR w redshift

archives
-CADC: Access to a ton of catalogues
---APASS (All Sky Automated Survey): optical and IR photometry data
---BLAST: Photometry (redshifts, FIR luminosities)
---BRITE-Constellation: photometry, especially for bright stars.
---CFHT: shit ton of wavelengths (IR, UV, EUV, X-ray, gamma ray)
---GPS: radio, millimeter
---Subaru: optical
---TESS: optical
---VGPS, WALLABY, VLASS: radio
---XMM: x-ray, UV, optical
-ESA HST
-ESA JWST: spectroscopy, NIR/IR images
-Fermi: gamma rays
-Gaia TAP+: astrometry, photometry, spectroscopy, deduced positions, parallaxes, proper motions, radial velocities, and brightnesses

Skyview gives access to everythign image wise except HST and JWST
