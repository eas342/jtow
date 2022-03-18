==========
Parameters
==========

:code:`jtow takes` parameters from a parameter file as a dictionary.

Here are some descriptions of what the parameters do.


:code:`custBias`
~~~~~~~~~~~~~~~~~

:code:`custBias` controls the bias subtraction.

* :code:`None` in Python or :code:`null` in the YAML parameter file to use the default bias from the pipeline.
* Path If you give it a path to a custom fits file (e.g. :code:`bias/jwst_nircam_superbias_0027.fits`), it will use that superbias.
* :code:`selfBias` : it will use the first frame available


:code:`saveROEBAdiagnostics`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:code:`saveROEBAdiagnostics` saves diagnostics from row-by-row, odd/even by amplifier (ROEBA) correction.

* :code:`True` saves the step immediately following ROEBA correction (full ramp with all integrations). If growing the mask with a kernel, it also saves the growth kernel and before/after growth.
* :code:`False` does not save these images.

:code:`jumpRejectionThreshold`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
:code:`jumpRejectionThreshold` sets the sigma rejection threshold for the JWST jump step to detect cosmic rays

:code:`ROEBAmaskGrowthSize`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
:code:`ROEBAmaskGrowthSize` sets the size (in pixels) of how large a smoothing kernel should be used to grow the ROEBA mask.
This allows the mask to go into faint wings of the image
