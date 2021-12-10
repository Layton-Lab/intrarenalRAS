# intrarenalRAS

This repository includes a public-facing version of the Layton Lab’s intrarenal renin-angiotensin system model. The code is the companion to “The intrarenal renin-angiotensin system in hypertension: Insights from mathematical modelling”.

This file contains information about how to simulate this model and run the associated code.

To create Figs 3–8 from the manuscript, run the script Smith2021.m 

To run your own simulations, call run_model.m with the appropriate inputs.

Smith2021.m
Running this will create all figures from the manuscript by calling get_Figk.m (k = 3,…,8).
Optional inputs can be specified by the user to
(1)	Generate only select figures
(2)	Alter the dose, type, and length of Ang II infusion experiment simulated in Figs 6 and 8
(3)	Alter the percentage used in the sensitivity analysis 
(4)	Save desired figures as PDFs

get_Figk.m (k = 3,…,8)
Running this will generate Fig k from the manuscript
User inputs:
If k = 6 or 8, you need to specify the Ang II infusion parameters
(1)	Dose: (0, infinity)
(2)	Units: ‘ng/min’ or ‘ng/min/kg’ (scaled to a 284g rat)
(3)	Type: ‘SC’ or ‘IV’
(4)	Days: (0, infinity)
If k = 7 or 8, you will need to specify the fold change in each parameter, delta ([-1,1]) to be used in the sensitivity analysis
