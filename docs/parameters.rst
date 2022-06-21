.. _parameter-descriptions:

==========
Parameters
==========

:code:`jtow takes` parameters from a parameter file as a dictionary.
The default parameters can be found at :ref:`the default parameters <default-parameters>` 

Here are some descriptions of what the parameters do.


:code:`custBias`
~~~~~~~~~~~~~~~~~

:code:`custBias` controls the bias subtraction.

* :code:`None` in Python or :code:`null` in the YAML parameter file to use the default bias from the pipeline.
* Path If you give it a path to a custom fits file (e.g. :code:`bias/jwst_nircam_superbias_0027.fits`), it will use that superbias.
* :code:`selfBias` : it will use the first frame available
* :code:`cycleBias`: it will cycle through biases in a pattern defined by biasCycle
* :code:`lineIntercept`: it will fit a line to the integration and use the intercept as the bias frame

:code:`biasCycle`
~~~~~~~~~~~~~~~~~
Controls how the bias is subtracted in the case that :code:`custBias` is equal to :code:`cycleBias`
For example, :code:`['A','B','B','B']` will do bias A, B, B, B, A, B, B, B.

:code:`biasCycleSearch`
~~~~~~~~~~~~~~~~~~~~~~~
Controls where to find the bias cycle files. For example :code:`data_path/superbias_nrca3_?.fits` will search for :code:`superbias_nrca3_A.fits` and :code:`superbias_nrca3_B.fits` if the biasCycle contains A and B.

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

:code:`saveBiasStep`
~~~~~~~~~~~~~~~~~~~~
Save the result of the bias subtraction step of the pipeline?

:code:`saveJumpStep`
~~~~~~~~~~~~~~~~~~~~
Save the jump step result before ramp fitting?

:code:`doLincor`
~~~~~~~~~~~~~~~~~~~~
Do the linearity correction? It should be True for correct results, but sometimes can be helpful to turn off for troubleshooting
