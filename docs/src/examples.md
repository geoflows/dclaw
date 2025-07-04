# Examples

There are two idealized examples to illustrate the capabilites of D-Claw: gully and radial slide. In addition, we provide a list of field-scale application maintained in separate repositories.


## Idealized examples

### Gully

The gully example is located within `dclaw/examples/gully` and is an idealized example of a debris flow moving down a confined channel and then out onto an unconfined plane. To run this example, first [install D-Claw](#install-dclaw) and then navigate to the gully application directory and execute: 

```bash 
make .exe
make input
make .plots
```

### Radial Slide

The radial slide example is located within `dclaw/examples/radial_slide` and is an ideallized example of a landslide-tsunami from a radially symmetric island. To run this example, first [install D-Claw](#install-dclaw) and then navigate to the radial slide application directory and execute: 

```bash 
make .exe
make input
make .plots
```

## Field-scale examples

:::{admonition} Note
This is an incomplete list and will be updated as more examples are disseminated.
:::


### Simulations of Barry Arm landslide complex

Barnhart and Collins (2025) provide an example of setting up and running multiple simulations of the Barry Arm landslide complex in support of Barnhart and others (2025).

### Simulation of the 2024 Surprise Inlet landslides

Barnhart (2025) provides an example of setting up and running a single simulation of the 2024 Surprise Inlet landslides in support of Karasözen and others (2025).

### References

- Barnhart, K.R., 2025, Simulation of Surprise Inlet landslides, version 1.0.0: U.S. Geological Survey software release, <https://doi.org/10.5066/P13BEXYP>.
- Barnhart, K.R., and Collins, A.L., 2025, Barry Arm landslide complex tsunami modeling scenarios: U.S. Geological Survey software release, <https://doi.org/10.5066/P175RQZR>.
- Barnhart, K. R., George, D. L., Collins, A. L., Schaefer, L. N., and Staley, D.M., 2025, Uncertainty reduction for subaerial landslide-tsunami hazards. Journal of Geophysical Research: Earth Surface, 130, e2024JF007906. <https://doi.org/10.1029/2024JF007906>.
- Karasözen, E., West, M.E., Barnhart, K.R., Lyons, J.J., Nichols, T., Schaefer, L.N., Bahng, B., Ohlendorf, S., Staley, D.M. and Wolken, G.J., 2025, 2024 Surprise Inlet landslides---Insights from a prototype landslide‐triggered tsunami monitoring system in Prince William Sound, Alaska: Geophysical Research Letters, v.52, e2025GL115911, <https://doi.org/10.1029/2025GL115911>.