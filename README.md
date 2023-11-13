# transition-models

## RAS transition models
### kOmegaSSTGam
#### Description
    A One-Equation Local Correlation-Based Transition Model based on the k-omega-SST RAS model.
#### References
    Menter, F.R., Smirnov, P.E., Liu, T. et al. (2015). A One-Equation Local Correlation-Based Transition Model. Flow Turbulence Combust 95, 583–619.
### kOmegaSSTLMK1
#### Description
    Langtry-Menter 4-equation transitional SST model based on the k-omega-SST RAS model with Krause's correlations.
#### References
    Langtry, R. B., & Menter, F. R. (2009). Correlation-based transition modeling for unstructured parallelized computational fluid dynamics codes. AIAA journal, 47(12), 2894-2906.
    Menter, F. R., Langtry, R., & Volker, S. (2006). Transition modelling for general purpose CFD codes. Flow, turbulence and combustion, 77(1-4), 277-303.
    Langtry, R. B. (2006). A correlation-based transition model using local variables for unstructured parallelized CFD codes. Phd. Thesis, Universität Stuttgart.
    Krause M, Behr M, Ballmann J. (2008). Modeling of Transition Effects in Hypersonic Intake Flows Using a Correlation-Based Intermittency Model. 15th AIAA International Space Planes and Hypersonic
	Systems and Technologies Conference.
### kOmegaSSTLMK2
#### Description
    Langtry-Menter 4-equation transitional SST model based on the k-omega-SST RAS model with Krause's correlations for hypersonic flow.
#### References
    Langtry, R. B., & Menter, F. R. (2009). Correlation-based transition modeling for unstructured parallelized computational fluid dynamics codes. AIAA journal, 47(12), 2894-2906.
    Menter, F. R., Langtry, R., & Volker, S. (2006). Transition modelling for general purpose CFD codes. Flow, turbulence and combustion, 77(1-4), 277-303.
    Langtry, R. B. (2006). A correlation-based transition model using local variables for unstructured parallelized CFD codes. Phd. Thesis, Universität Stuttgart.
    Frauholz S, Reinartz B U, Müller S, et al. (2014). Transition Prediction for Scramjets Using γ-Reθt Model Coupled to Two Turbulence Models. Journal of Propulsion & Power, 1-19.
    Krause M. (2010). Numerical analysis of transition effects for scramjet intake flows. Phd. Thesis, RWTH Aachen University.
### kOmegaSSTLMK2Y
#### Description
    Langtry-Menter 4-equation transitional SST model based on the k-omega-SST RAS model combined Krause's correlations for hypersonic flows with You's modifications for separated transitional flows.
#### References
    Langtry, R. B., & Menter, F. R. (2009). Correlation-based transition modeling for unstructured parallelized computational fluid dynamics codes. AIAA journal, 47(12), 2894-2906.
    Menter, F. R., Langtry, R., & Volker, S. (2006). Transition modelling for general purpose CFD codes. Flow, turbulence and combustion, 77(1-4), 277-303.
    Langtry, R. B. (2006). A correlation-based transition model using local variables for unstructured parallelized CFD codes. Phd. Thesis, Universität Stuttgart.
    Frauholz S, Reinartz B U, Müller S, et al. (2014). Transition Prediction for Scramjets Using γ-Reθt Model Coupled to Two Turbulence Models. Journal of Propulsion & Power, 1-19.
    Krause M. (2010). Numerical analysis of transition effects for scramjet intake flows. Phd. Thesis, RWTH Aachen University.
    You Y, Luedeke H, Eggers T, et al. (2013). Application of the γ-Reθt Transition Model in High Speed Flows. Aiaa/3af International Space Planes & Hypersonic Systems & Technologies Conference. 

## DES transition models
### kOmegaSSTDES
#### Description
    k-omega-SST DES turbulence model modified according to Menter & Kuntz (2004).
#### Usage
    FSST  0;
    //- Zonal filter choice
    //
    // - 0: no filtering
    // - 1: (1 - F1)
    // - 2: (1 - F2)
#### References
    Strelets, M. (2001). Detached Eddy Simulation of Massively Separated Flows, 39th AIAA Aerospace Sciences Meeting and Exhibit, Reno, NV
    Menter, F.R., Kuntz, M. (2004). Adaptation of Eddy-Viscosity Turbulence	Models to Unsteady Separated Flow Behind Vehicles. In: McCallen, R., Browand, F., Ross, J. (eds) The Aerodynamics of Heavy Vehicles: Trucks, Buses, and Trains. Lecture Notes in Applied and Computational Mechanics, vol 19. Springer, Berlin, Heidelberg.
### kOmegaSSTGamDES
#### Description
    A One-Equation Local Correlation-Based Transition Model based on the k-omega-SST DES model.
#### References
    Menter, F.R., Smirnov, P.E., Liu, T. et al. (2015). A One-Equation Local Correlation-Based Transition Model. Flow Turbulence Combust 95, 583–619.
### kOmegaSSTLMDES
#### Description
    Langtry-Menter 4-equation transitional SST model based on the k-omega-SST DES model.
#### References
    Langtry, R. B., & Menter, F. R. (2009). Correlation-based transition modeling for unstructured parallelized computational fluid dynamics codes. AIAA journal, 47(12), 2894-2906.
    Menter, F. R., Langtry, R., & Volker, S. (2006). Transition modelling for general purpose CFD codes. Flow, turbulence and combustion, 77(1-4), 277-303.
    Langtry, R. B. (2006). A correlation-based transition model using local variables for unstructured parallelized CFD codes. Phd. Thesis, Universität Stuttgart.
    Nakhostin, S.M. (2019). Investigation of transitional turbulence models to predict drag crisis for flow over spheres and cylinder. Master's Thesis, Universitetet i Stavanger.
    
## tutorials
