A Vlasov-Fokker-Planck code for the study of electron transport in a collisional, magnetised plasma.  The algorithm is based on a spherical harmonic expansion of the distribution function to deal with the angular part.  This allows for spatially three dimensional simulations to be performed on relatively modest computational resources.

The full algorithm and benchmarks are outlined in the paper:

"A code to solve the Vlasov–Fokker–Planck equation applied to particle transport in magnetic turbulence"

W A Hornsby et al 2010 Plasma Phys. Control. Fusion 52 075011

The Fast Fourier Transform library FFTW (Version 2.5) is currently a requirement to compile the code, modify the included makefile if the libraries arent in  /usr/local.

You are free to modify the code as you please.  All I ask is that you cite the above paper.


Any questions/requests can be made directly to the author of the code at dr.william.hornsby@gmail.com.

Many thanks.