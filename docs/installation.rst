.. highlight:: shell

============
Installation
============


Stable release
--------------

To install jtow, run this command in your terminal:

.. code-block:: console

    $ pip install jtow

This is the preferred method to install jtow, as it will always install the most recent stable release.

If you don't have `pip`_ installed, this `Python installation guide`_ can guide
you through the process.

.. _pip: https://pip.pypa.io
.. _Python installation guide: http://docs.python-guide.org/en/latest/starting/installation/


Recommended Environmental Variables
-----------------------------------
It is recommended to set the environmental variables to that files are stored where expected.
For example, suppose /big_disk is a large storage location.

.. code-block:: bash

    export CRDS_PATH=/big_disk/crds_cache
    export CRDS_SERVER_URL=https://jwst-crds.stsci.edu
    export TSHIRT_DATA=/big_disk/tshirt_data
    export JWSTDOWNLOAD_OUTDIR=/big_disk/jwst_flight_data
    export MAST_API_TOKEN="blahblahblah"

:code:`CRDS_PATH` will store JWST reference data

:code:`TSHIRT_DAT` will store photometric/spectroscopic extractions and time series

:code:`JWSTDOWNLOAD_OUTDIR` will store downloaded JWST data

:code:`MAST_API_TOKEN` is the API token to access data on MAST.


From sources
------------

The sources for jtow can be downloaded from the `Github repo`_.

You can either clone the public repository:

.. code-block:: console

    $ git clone git://github.com/eas342/jtow

Or download the `tarball`_:

.. code-block:: console

    $ curl -OJL https://github.com/eas342/jtow/tarball/master

Once you have a copy of the source, you can install it with:

.. code-block:: console

    $ python setup.py install


.. _Github repo: https://github.com/eas342/jtow
.. _tarball: https://github.com/eas342/jtow/tarball/master
