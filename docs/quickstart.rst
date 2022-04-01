==========
Quickstart
==========

To run ROEBA (row-by-row, odd/even, by amplifier) corrections to reduce 1/f noise do the following.

1. Install jtow as follows. It works if you already have the `jwst pipeline installed in an environment. <https://github.com/spacetelescope/jwst#installation>`_.

.. sourcecode:: python

   pip install jtow
   
2. Make sure 2 required jwst products are in your folder. For example

.. sourcecode:: bash

   jw01074100001_02108_00002_nrcblong_rate.fits
   jw01074100001_02108_00002_nrcblong_uncal.fits
   
(The :code:`_rate.fits` file is needed to make a source mask to find background pixels)

3. Run the command line utility

.. sourcecode:: bash

   roeba_run "*_uncal.fits"

Help on the command line utility.
