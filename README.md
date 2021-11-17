# Lie Group code

This repository gathers different pieces of code related to our paper "Lie symmetries applied to interval integration" 


## C++

### Compilation requirements

To compile most of these you will need the following:
- The [IBEX](http://www.ibex-lib.org/) library
- The [CODAC](http://codac.io/) library
- [OPTIONAL] The [CAPD](http://capd.ii.uj.edu.pl/) library to compile the fastest version of the 
  algorithm in (src/ipegenerator/article_figures_ipe)
- [OPTIONAL] [ipegenerator](https://github.com/JulienDamers/ipe_generator) library if you choose 
  visualisation option 2

### Visualization
- [OPTION 1] (default)  [VIBes](https://enstabretagnerobotics.github.io/VIBES/) for quicker 
  visualization. Open the VIBes viewer **and then** run the example to see the figure appearing

- [OPTION 2] the [ipegenerator](https://github.com/JulienDamers/ipe_generator) library if you 
  want to compile the codes used to create the figures of the article  (placed in 
  src/ipegenerator_version/article_figures). Use the cmake option __-DIPE=ON__ if you want to use 
  this 
  version. As mentioned above you need to install [ipegenerator](https://github.com/JulienDamers/ipe_generator) 
  and [CAPD](http://capd.ii.uj.edu.pl/) libraries for this  option. When run, it will create a 
  .pdf and .ipe (XML) files. You can open the .ipe file with [IPE](https://ipe.otfried.org/) to 
  edit it.


## Python (in progress)

### Repl.it

If you only want a "quick" live computation of the different test cases presented in the article 
you can access these different repl. Just click on "run" on the webpage

- [Test case 1 continuous](https://replit.com/@JulienDamers/Lie-symmetries-test-case-1-continuous)
- [Test case 1 discrete](https://replit.com/@JulienDamers/Lie-symmetries-test-case-1-discrete)
- [Test case 2 continuous](https://replit.com/@JulienDamers/Lie-symmetries-test-case-2-continuous)
- [Test case 2 discrete](https://replit.com/@JulienDamers/Lie-symmetries-test-case-2-discrete)
- [Test case 3 continuous](https://replit.com/@JulienDamers/Lie-symmetries-test-case-3-continuous)
- [Test case 3 discrete](https://replit.com/@JulienDamers/Lie-symmetries-test-case-3-discrete)
- [Test case 4 continuous](https://replit.com/@JulienDamers/Lie-symmetries-test-case-4-continuous)
- [Test case 4 discrete](https://replit.com/@JulienDamers/Lie-symmetries-test-case-4-discrete) 


### Local run

If you prefer to run it (or do some tests) on your local machine, you will need python3 
(version>=3.6)


#### Required modules

- pyIbex (version >= 1.9.2)
- codac (version >=0.1.7)
- vibes (version >=0.2.2)

On linux, you can install the three modules using the following command
```sh
python3 -m pip install codac
```

#### Visualisation

- There is only [VIBes](https://enstabretagnerobotics.github.io/VIBES/) available. As in the C++ 
  version, visualisation option 1, open the VIBes viewer **and then** run your python program to see
  your figure appearing
