# Helper debug scripts

Collection of scripts designed to visualize the shape functions involved with
generating the simulated topography. It uses `numpy` and `matplotlib` to
generate interactive plots of such functions.

- `plot_ellipse.py` plots the portion of the half ellipse to be created by the
  vibration turning process (v0.1.0)
- `plot_ellipse_full.py` draws the same ellipse along with the entire
  half-ellipse (v0.1.0)
- `plot_texture.py` draws the tool shape as used by the `texture.*` files (that
  is, considering `f` and `ap` as independent parameters) (v0.1.0, but still
  useful)
- `plot_texture_shallow.py` draws the tool shape as used by the
  `texture_shallow.*` files (that is, considering `f` and `ap` as dependent
parameters) (v0.1.0, but still useful)
- `plot_tool_path_param.py` plots a parametrized version of the tool path
  functions
- `plot_tool_path_sinusoid.py` plots a position-based version of the tool path 
  functions, based on sinusoidal functions. Not used for the implementation,
  but displays a possible approach to the modelling problem.
