The Fortran codes for the Random-Fatigue Limit (RFL) are in the file stdrfl.f.  Other codes needed are listed in Dependencies.  
The file stdrfl-Rwin32bit.zip contains R codes for fitting the model to data.  
The codes have been tested ONLY on Fortran 77 and R 32-bit.  

Writing the Fortran codes began in the mid-1990s as part of my Statistics PhD dissertation at Iowa State University.  I had a computer programming background when I began the project, but I had no prior experience with programming in Fortran.  My adviser Professor Bill Meeker shared with me some of his Fortran codes related to my research project.  I had to learn Fortran on the fly by adapting these codes for my research.  

The calculation of the RFL pdf and cdf requires numerical integration.  The main problem with calculating them with one integral was that the integration algorithm would often miss the intervals over which the pdf integrand is nonzero.  As a result, the pdf or cdf value returned is smaller than what it should be.  So, to help the algorithm find positive integrand values, the calculation is broken over several intervals (of the integrator value) so that interval endpoints are where the integrand has positive mass.  

The codes have gone through a significant amount of revisions over 10 years.  The original codes (not posted here) were based on the 5-parameter model as described in the paper:
Pascual, F. G. and Meeker, W. Q. (1999). Estimating Fatigue Curves with the Random Fatigue-Limit 
  Model. Technometrics, volume 41, pages 277-290.

The current version of the codes is based on a "standardized" RFL with 3 parameters as described in this paper:
Pascual, F. G. (2003). A Standardized Form of the Random Fatigue-Limit Model. Communications in Statistics - Simulation and Computation, 32, 1209-1228.
The codes were created mainly to efficiently perform numerical integration outside of R.  The Fortran subroutines 
  srfl0pdf, srfl0cdf, srfl0quan
are called from within R to compute the pdf, cdf, and quantile values.  

I would like to thank Professor Bill Meeker for guidance with this project.
