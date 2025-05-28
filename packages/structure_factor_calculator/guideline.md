# General guideline
- Create a new branch. 
- Only work in the following folders/files:
  1. `packages/structure_factor_calculator/` folder: this is the core backend calculator. Your code
     should go in here.
  2. `tests/structure_factor_calculator/` folder: this is the folder for tests.
  3. `requirements.txt` file: this is the file for dependencies.
  4. `documentations/` folder (Optional): this is the folder for documentation. 
- Try not to change other files. 
- request for review when the code is ready. I will merge the code to the main branch. 


## Object oriented
Wrap the code in a class. 
- Implement standard methods: `__init__()`, `initialize()`
`is_initialized()`. 
- Implement methods for calculating the structure factor. 
- lengthy codes could be put in a separate file. 
- the package `structure_factor_calculator` is only used for structure factor calculation. The
  visualization functionality should be implemented in the the package `packages/visualizer/`. But
  this is not necessary for now.


## Testing
- Write tests for the code. Use pytest and put tests file in the
  `tests/structure_factor_calculator/` folder. 
- Ideally, each function should be tested with multiple test cases. 


# Bonus
After implementing the structure factor calculator, you can continue to work on the visualization
functionality in the the package `packages/visualizer/`. Don't forget to define a class and wrap all
the necessary functions into the class.
