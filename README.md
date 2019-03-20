# Athena
#### A toolkit for exploring structural variation mutation rates and dosage sensitivity

Copyright (c) 2019, [Ryan L. Collins](mailto:rlcollins@g.harvard.edu) and the Talkowski Laboratory.  
Distributed under terms of the [MIT License](/LICENSE) (see `LICENSE`).  

---  

## Table of Contents

#### Getting started
  * [Run Athena with Docker](https://github.com/talkowski-lab/athena#run-from-docker)
  * [Manual Installation](https://github.com/talkowski-lab/athena#manual-installation)

#### The Athena workflow

#### Other
  * [About the name](https://github.com/talkowski-lab/athena#about-the-name)


---  

## Run from Docker

The recommended way to run Athena is from its dedicated Docker container. This will handle all dependencies and installation for you, and ensure you are running the latest version.

```
$ docker pull talkowski/athena
$ docker run --rm -it talkowski/athena
```

## Manual installation

If you would prefer to install Athena on your own system, you can do so with `pip`.

```
$ git clone https://github.com/talkowski-lab/athena.git
$ cd athena
$ pip install -e .
```

## The Athena workflow

_Coming soon!_

### About the name
This package is named after [Athena](https://en.wikipedia.org/wiki/Athena), the Greek goddess of wisdom, strategy, tactics, and mathematics. She was selected as the namesake for this package given that it relies on understanding the features that influence structural variation mutation rates (_wisdom_), incorporating those features into a statistical model (_mathematics_), and using these models to infer which components of the genome are vulnerable to changes in copy number (a kind of genomic _tactics_/_strategy_).