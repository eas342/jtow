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

:code:`ROEBACorrection`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If :code:`ROEBACorrection` is True, it row-by-row, odd/even by amplifier (ROEBA) correction on each group after bias subtraction instead of the reference pixel correction. If it is False, no ROEBA is done. The reference pixel step is run, which only takes effect if reference pixels are available. If ROEBCorrection is :code:`GROEBA`, then it uses a Gaussian Process train-1/f noise model

:code:`ROEBAK`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If :code:`ROEBAK` is True, it adds a step to ROEBA that attempts to remove the kTC noise. Only has an effect if :code:`ROEBACorrection` is True. If False, no extra kTC correction will be done. This does a slow-read correction of all groups in a frame, calculates the median pixel for all groups for integration-specific kTC correction. This kTC correction is subtracted from all groups in that integration to find the underlying 1/f noise better. However, the kTC correction is undone after the 1/f correction to be sure it doesn't not introduce any new noise, especially on the source pixels. This is a similar idea to using the :code:`custBias="lineIntercept"` but with undoing the bias.

:code:`colByCol`
~~~~~~~~~~~~~~~~
If :code:`colByCol` is True, ROEBA will use the median of each column to estimate the slow-read correction. If False, it will subtract the median of all odd columns from all odd columns within an amplifier and the median of all even columns from all even columns within an amplifier.

:code:`saveROEBAdiagnostics`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:code:`saveROEBAdiagnostics` saves diagnostics from row-by-row, odd/even by amplifier (ROEBA) correction. If just using regular reference pixel correction, then the regular reference pixel correction is saved.

* :code:`True` saves the step immediately following ROEBA correction (full ramp with all integrations). If growing the mask with a kernel, it also saves the growth kernel and before/after growth. If doing regular reference pixel correction instead of ROEBA, save the output of the reference pixel step.
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

:code:`simpleSlopes`
~~~~~~~~~~~~~~~~~~~~
Do a simple line fit rather than the most-optimal (right now ordinary least squares) fit? 

* :code:`None` ('null' in .yaml file), no simple slopes are done. Regular jwst pipeline ramp fits.
* :code:`'Both'` Both the jwst fits and simple slopes are performed.
* :code:`'Only'` Only simple slopes are calculated and the (slow) most-optimal fit is skipped.

:code:`rampFitWeighting`
~~~~~~~~~~~~~~~~~~~~~~~~
Pass the weighting scheme along to the ramp fit step. If :code:`rampFitWeighting` is 'optimal' it will use optimal weights. If :code:`rampFitWeighting` is 'unweighted', not weighting of groups up the ramp is performed - this approaches last minus first for bright targets.
More info is available at the  `JWST pipelin ramp fitting description page <https://jwst-pipeline.readthedocs.io/en/latest/jwst/ramp_fitting/description.html>`_

:code:`side_smoothing_length`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Set the side smoothing length for reference pixels. This is passed to the :code:`jwst` reference pixel step. This does not affect ROEBA so if :code:`ROEBACorrection` is set to :code:`True`, this will not matter (in the current version of jtow at least). For ROEBA, use :code:`smoothSlowDir`

:code:`smoothSlowDir`
~~~~~~~~~~~~~~~~~~~~~~
If :code:`None` ('null' in .yaml file), no smoothing is done with ROEBA. If set to an int or float, a Savgol filter is applied along the slow-read direction to smooth the ROEBA model.

:code:`useGrismRefpx`
~~~~~~~~~~~~~~~~~~~~~
Use the reference pixels for ROEBA corrections? If True, the left refpix will be used for F322W2 data and right refpix will be used for F444W data. If False, only the background pixels will be used from the rightmost amplifier (F322W2 data) or leftmost amplifier (F444W data).

:code:`custGroupDQfile`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Set a custom data quality array for after the saturation step to manually set some pixels as saturated. Could be useful if you want to treat all the pixels in a column the same or enforce that the data quality flags are not as variable. Should be the path to a FITS file with the same number of groups as the data and the same dimensions as the data. Will be combined with all integrations' group data quality flags with a bitwise or. If :code:`None` ("Null" in YAML file), this is not done.

:code:`skipJumpDet`
~~~~~~~~~~~~~~~~~~~~~~
Skip the jump detection step? Passed to the jump detection skip parameter. If True, it is skipped.
