-----------------------------------------------------------------------------------------------------
- By Carlos Rios, August 2015
- Revised by Yi-Siou Huang, April 2022

-----------------------------------------------------------------------------------------------------
- Instructions on how to run are found commented at the beginning of each code. 

- All the programs require transfer.matrix function to run properly. ALso the "library" of refractive indexes in import_refr_in.

- Any material that we want to use should be included in the file import_refr_in.m.

-----------------------------------------------------------------------------------------------------------
-PROGRAMS:

** multilayers: Calculates A, R and T for any stack of layers of any thickness. 

** multilayers_bothPhasesART: Calculates A, R and T for any stack including phase change materials (whichever is in the import.refr.in "library"). This program will calculate and plot together the A, R and T of both states. REMEMBER: change the refractive index to the appropriate material. 

** multilayers_ColorContrast_SweepThickness_PCM_ITO: Calculates color contrast in the sweep thickness of PCM and ITO.
   The scale strats from zero, no matter which thickness starting from min (Same as _Alpha, reserved as an original code).

** multilayers_ColorContrast_SweepThickness_PCM_ITO_Alpha: Calculates color contrast in the sweep thickness of PCM and ITO.
   The scale strats from zero, no matter which thickness starting from min.

** multilayers_ColorContrast_SweepThickness_PCM_ITO_Beta: Calculates color contrast in the sweep thickness of PCM and ITO.
   The scale strats from the min thickness.

** multilayers_PhaseDifference: Calculates all R,A,T, phi - the phase difference in reflection between Am and Cry

** multilayers_PhaseDifference_ColorGamut: Calculates color gamut for both PCM phases and phase difference of  a stack of lyaers.

** multilayers_PhaseDifference_ColorGamut_ColorSwatches: Calculates color gamut, sawtches for both PCM phases and phase difference of  a stack of lyaers.

** multilayers_PhaseDifference_trans: Calculates all R,A,T, phi, and psi, where phi and psi are the reflected and trasnmitted phase, respectively. IT uses the extended transfer matrix code called transfer_matrix_trans. 

** multilayers_SweepThicknessART: Calculates A, R and T FOR sweep thickness for PCM of a stack of layers

-----------------------------------------------------------------------------------------------------------------------
ANY PROBLEM? Contact
Carlos: carlos.riosocampo@gmail.com
Yi-Siou: yisiouh@umd.edu

