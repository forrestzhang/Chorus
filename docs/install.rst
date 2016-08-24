Install
=======

Ubuntu 14.04
------------

Install dependent package

.. code-block:: bash

    $ apt-get update && apt-get install -y cython3  build-essential \
        zlib1g-dev \
        zlibc \
        git \
        libboost-dev \
        autoconf \
        libncursesw5-dev \
        libncurses5 \
        ncurses-dev \
        libboost-thread-dev \
        python3-pip \
        samtools \
        unzip \
        python \
        curl \
        wget \
        python3-pyqt5 \
        libfreetype6-dev \
        libxft-dev \
        python3-matplotlib

    $ apt-get remove -y python3-matplotlib


Install jellyfish

.. code-block:: bash

    $ mkdir /opt/software

    $ cd /opt/software

    $ wget https://github.com/gmarcais/Jellyfish/releases/download/v2.2.3/jellyfish-2.2.3.tar.gz

    $ tar zxvf jellyfish-2.2.3.tar.gz

    $ mv jellyfish-2.2.3  jellyfish

    $ cd jellyfish

    $ ./configure && make && make install


Install bwa

.. code-block:: bash

    $ cd /opt/software

    $ git clone https://github.com/lh3/bwa.git

    $ cd bwa

    $ make


Install primer3-py

.. code-block:: bash

    $ cd /opt/software

    $ wget https://github.com/forrestzhang/primer3-py/archive/unicode.zip

    $ unzip unicode.zip

    $ cd primer3-py-unicode

    $ python3 setup.py install


Install Python dependent package

.. code-block:: bash

    $ pip3 install numpy pyfasta matplotlib

    $ pip3 install pandas==0.16.2


