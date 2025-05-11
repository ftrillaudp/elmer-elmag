# Single-phase transformer

300 W single phase transformer (shell type with windings on the center core leg).


Primary winding (HV):
- $N_p$ = 840
- Filling factor: 0.7033
- $V_p^r$ = 127 V
- $I_p^r$ = 2.4 A

Secondary winding (LV):
- $N_s$ = 7
- $V_s^r$ = 1 V
- $I_s^r$ = 300 A

Electrical conductivity of Cu wires: 5.96e7 S/m

Core:
- Manufacturer: Tempel (26N174)
- EI type (catalogue: FR-250 R6W/0251)
- Material: M19 or M55
- Number of sheets: 238
- Stacking factor: 0.94-0.99
- Permeability relative: 1500
- Electric conductivity: 1000


# Run the codes

1. ```python run.py```: run the short-circuit and load transient case
2. ```python ocrun.py```: run the open-circuit transient case
3. ```python frequency.py```: run the frequency load model

Issue with current density: no sinusoidal form as expected in open-circuit test and current in same direction in frequency domain results.
