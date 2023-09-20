# SkyMat
The Kujawinski Lab uses a unique method for processing samples prepared by the [Widner et al. (2021)](https://doi-org.libproxy.mit.edu/10.1021/acs.analchem.0c03769) method. 
The method has gone through some upgrades (for example, double stable isotope labels) in the intervening years, but for anyone using the chemical method and peak-picking in Skyline, this is meant to be a set of codes that will get you from peak areas to concentrations, provided a standard curve. 

## How to use this repository.
Below, you will find a diagram of how to use these scripts in different applications of the method.

<img src="images/flow.jpg" width="800">

The main interface here will be `riSkyline.m`, and you should only use the version specific to the isotope label you're working with. If you have used both the 13C and the D5 internal standards, you'll use both versions here, which will generate two nearly-identical datasets. Ultimately, you will have _four_ standard curves this way (one per isotope and per ion mode), although only one ionization mode or isotope may be present for a given metabolite. 