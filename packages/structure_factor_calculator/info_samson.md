# Prerequisites
0. setup SSH key:
   - **Set Up SSH Key (only if not already done)**
     - Generate SSH key:
       ```bash
       ssh-keygen -t ed25519 -C "your_email@example.com"
       ```
     - Add SSH key to ssh-agent:
       ```bash
       eval "$(ssh-agent -s)"
       ssh-add ~/.ssh/id_ed25519
       ```
     - Add public key to GitHub:
       - Go to GitHub → Settings → SSH and GPG Keys → New SSH Key
       - Paste content of `~/.ssh/id_ed25519.pub`
   - **Clone the Repository**
     - Using SSH (recommended if SSH is set up):
       ```bash
       git clone git@github.com:HongXunyang/rixs_pre_toolbox.git
       ```
   - **Create a Branch**
     ```bash
     git switch -c <your_branch_name>
     ```
   - **Make a Test Changes and Commit**
     ```bash
     git add .
     git commit -m "Describe your changes"
     ```
   - **Push Changes to GitHub**
     ```bash
     git push origin <your_branch_name>
     ```

1. Install the dependencies:
   ```bash
   pip install -r requirements.txt
   ```


# General guideline
- Work only in your own branch. 
- Work only in the following folders/files:
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
- Push to github when all tests pass.


# Bonus
After implementing the structure factor calculator, you can continue to work on the visualization
functionality in the the package `packages/visualizer/`. Don't forget to define a class and wrap all
the necessary functions into the class.
