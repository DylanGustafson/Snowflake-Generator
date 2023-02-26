This is a personal tinkering project to generate snowflakes using cellular automata using my own customized implementation. Crystal growth is dictated by simple neighbor-counting rules that are analogous to ambient temperature, pressure, and humidity effects for real snowflakes. This program is not meant to be an accurate physics simulation, but rather a demonstration of emergence.

![Screenshot](https://raw.githubusercontent.com/DylanGustafson/Snowflake-Generator/main/Screenshot.png)

## Background
Snowflake crystals grow when supercooled water droplets collide with seed crystals and fuse solid. This process releases heat due to the enthalpy of fusion, warming up the immediate area and temporarily inhibiting further growth. When humidity is low, the crystal growth is slow enough that this warming effect is insignificant. Thus, due to the hexagonal crystal structure of ice, the crystal grows outward uniformly into a hexagon. When humidity is high and the warming effect is not negligible, crystal growth will be reduced in areas that are already more "filled in." This results in the natural formation of spikes and dendrites with hexagonal symmetry.

This program contains a 2D lattice of cells in a hexagonal packing arrangement. All cells start "off" except for a single cell at the center, representing the seed crystal. When computing a generation, the program counts the adjacent neighbors (distance 1) and "proximate" neighbors (distance sqrt 3) of each "off" cell to check if it qualifies to be turned on in the next generation. Adjacent and proximate neighbor totals in both the current and previous generation are considered. These qualifying limits on the neighbor totals indirectly represent the effects of ambient temperature, pressure, and humidity at the time of generation. These limits can be modified between generations, representing how a falling snowflake can encounter areas of different ambient conditions that change the way it continues to grow.

A classic example is the simple rule that cells can only activate if they have one adjacent "on" neighbor. Running this rule produces shapes alternating between hexagonal plate and stellar plate. Another example is one that requires one adjacent neighbor and zero proximate neighbors. This stricter rule produces straight spikes radiating from their origin. When paired with other rules it can form branches and dendrites.

There are endless combinations of rules to experiment with, however so far my favorite results are created by combining the "radiating spikes" rule (as a sufficient condition) with neighbor restriction rules (as necessary conditions). The radiating spikes rule encourages branching, hence the "Encourage branching" option, while the neighbor restrictions dictate how these branches grow. See the Snowflake.generate() function for details. I encourage you to experiment with other combos.

I have included my favorites in a preset menu on the GUI, although the best shapes come from changing the rules between generations. This is analogous to a snowflake falling through changing atmospheric conditions.

## Installation
This program runs on Python 3.x, and uses libraries `NumPy`, `Matplotlib`, and `Tkinter`.
Run the program normally with:

```
python3 snowflake.py
```

The program uses Matplotlib's `FigureCanvasTkAgg` to composite the plot into the GUI, and `RangeSlider` to visually edit the neighbor restrictions. Both of these can be a bit finicky, especially on MacOS. Updating Matplotlib may fix things, but if there are still issues, try running the program with any of the follow options, which can be combined.

To display the snowflake in a separate window:
```
python3 snowflake.py -s
```

To use simple entry widgets rather than range sliders:
```
python3 snowflake.py -e
```

To run the GUI window without any Matplotlib elements:
```
python3 snowflake.py -t
```

## Controls

The upper right controls dictate how the snowflake grows, while the controls in the lower right are display options. 

Adjacent and Proximate neighbor limits are controlled using the range sliders in the upper right. Directly below them is a menu of presets I have cataloged that generate some nice shapes. Next to that is a check button that can toggle the rule that encourages branching.

Use the "Generate" button to compute the next generation, the "Generate all" button to compute generations until the canvas fram is filled, and "clear" to restart at the first frame.

At the bottom right are play buttons that can cycle forward or backward through previously generated frames, and can be held down to create an animation. Above the play buttons are display options such as zoom and resolution.

## Presets
Neighbor restriction presets are loaded from presets.txt, which must be located in the same directory as snowflake.py to work properly. Presets can be added to this file by copying the format of the included lines.
