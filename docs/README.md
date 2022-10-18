# Helper debug scripts

Collection of scripts designed to visualize the shape functions involved with
generating the simulated topography. It uses `numpy` and `matplotlib` to
generate interactive plots of such functions.

- `plot_ellipse.py` plots the portion of the half ellipse to be created by the
  vibration turning process
- `plot_ellipse_full.py` draws the same ellipse along with the entire
  half-ellipse
- `plot_texture.py` draws the tool shape as used by the `texture.*` files (that
  is, considering `f` and `ap` as independent parameters)
- `plot_texture_shallow.py` draws the tool shape as used by the
  `texture_shallow.*` files (that is, considering `f` and `ap` as dependent
parameters)
