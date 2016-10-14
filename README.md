LD31-analysis
=============

Data analysis scripts &amp; base code for XPCS experiments.

Includes code to
- count photons on a per-pixel basis
- determine the size distribution of speckles via autocorrelation
- identification of individual speckles
- fit a negative binomial distrbution to speckle photon analysis
- a combination of the above to computate image contrast

--------------
Installation on an LCLS machine

1) Create a python virtual environment you can install private libraries into
```
$ virtualenv xcsm9816
```

2) download the speckle repo (be sure you are on pslogin so you have internet)
```
$ git clone https://github.com/tjlane/speckle.git
```

3) install speckle
```
ssh psana
cd speckle
python setup.py install
```

4) try running some of the examples, e.g.
```
cd examples
python acf.py
```

the tests have some good example code to steal.

IMPORTANT NOTE: you should use the latest ana build 
```
$ sit_setup ana-0.19.15
```

